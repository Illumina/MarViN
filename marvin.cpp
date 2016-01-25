#define EIGEN_DONT_PARALLELIZE
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <iostream>
#include <queue>
#include <random>
#include <iomanip>
#include "marvin.hpp"
#include "EM.hpp"
#include <time.h>  
#include <limits>      
#include <Eigen/Eigenvalues> 

#include "htslib/hts.h"

using namespace std;
using namespace Eigen;


extern "C" {
#include "htslib/synced_bcf_reader.h"
}

void usage(){
  cerr << "Imputation from GLs" << endl;
  cerr << "Usage:" << endl;
  cerr << "./marvin -f input.vcf" << endl;
  cerr << "Expects input.vcf to contain GLs/PLs." << endl;
  cerr << "\n\t  Input params...\n" << endl;
  cerr << "\t -f  : 		input vcf" << endl;
  cerr << "\t -o  : 		output vcf" << endl;
  cerr << "\t -O  : 		output vcf type" << endl;
  cerr << "\t -Ogp:			output genotype probabilities" << endl;
  cerr << "\t -v:			verbose output" << endl;
  cerr << "\t -num_threads: 	threads in a single window (default 1)" << endl;
  cerr << "\n\t  Window params...\n" << endl;
  cerr << "\t -r: 			chromosome region (mandatory with sites file)" << endl;
  cerr << "\t -b: 			block size (mandatory with sites file)" << endl;
  cerr << "\t -ov: 			overlap (mandatory with sites file)" << endl;
  cerr << "\n\t  Panel params...\n" << endl;
  cerr << "\t -site:	vcf with sites in panel label" << endl;
  cerr << "\t -c:		collapse snps|indels|both|all|some|none" << endl;
  cerr << "\n\t  Run params...\n" << endl;
  cerr << "\t -max_its: 	number of iterations of algorithm (default 5)" << endl;
  cerr << "\t -inner_its:	number of iterations of inner loop (updating means, default 1)" << endl;
  cerr << "\t -EMits:		number of iterations for EM initialization" << endl;
  cerr << "\t -maxlr: 		max likelihood ratio" << endl;
  cerr << "\t -bias:		dosage < bias gt = 0, dosage > 2-bias gt = 2; (default 0.5)" << endl;
  cerr << "\n\t  Regularization params...\n" << endl;
  cerr << "\t -sigma_reg: 	use sigma regularization with params lambda, lambda2, pct" << endl;
  cerr << "\t -lambda: 		regularization parameter (0.06)" << endl;
  cerr << "\t -lambda2:		threshold steepness (4)" << endl;
  cerr << "\t -pct:			sig mid point (0.2)" << endl;

  exit(1);
	
}

int main(int argc, char* argv[])
{	
	
  //input vcf
  if( !cmdOptionExists(argv, argv+argc, "-f" ) ){ usage(); }
  string filename = getCmdOption(argv, argv+argc, "-f" );
	
  //Read in fixed reference panel data
  bool use_panel = false;
  string pfilename;
  if( cmdOptionExists(argv, argv+argc, "-site" ) ){ 
	use_panel = true;
    pfilename = getCmdOption(argv, argv+argc, "-site" ); 
  }

  int block = 0;  
  if( cmdOptionExists(argv, argv+argc, "-b" ) ){ block = atoi(getCmdOption(argv, argv+argc, "-b" )); }
  else{ usage(); }
  
  int overlap = 0;
  if( cmdOptionExists(argv, argv+argc, "-ov" ) ){ overlap = atoi(getCmdOption(argv, argv+argc, "-ov" )); }
  else{ usage(); }
	
  bool zero_missing = true; //cmdOptionExists(argv, argv+argc, "-zm" );
  bool verbose = cmdOptionExists(argv, argv+argc, "-v" );
  bool isec_only = cmdOptionExists(argv, argv+argc, "-isec" );
  
  //output vcf
  string out_filename = "out.vcf";
  string out_type = "w";
  if( cmdOptionExists(argv, argv+argc, "-o" ) ){ out_filename = getCmdOption(argv, argv+argc, "-o" ); } 
  if( cmdOptionExists(argv, argv+argc, "-O" ) ){ out_type += (string)( getCmdOption(argv, argv+argc, "-O" ) );} 
  bool glout = cmdOptionExists(argv, argv+argc, "-Ogp" );
	
  //regularization
  float lambda = 0.06;
  if( cmdOptionExists(argv, argv+argc, "-lambda" ) ){ lambda = atof(getCmdOption(argv, argv+argc, "-lambda" )); }
  float lambda2 = 4;
  if( cmdOptionExists(argv, argv+argc, "-lambda2" ) ){ lambda2 = atof(getCmdOption(argv, argv+argc, "-lambda2" )); }
  float pct = 0.2;
  if( cmdOptionExists(argv, argv+argc, "-pct" ) ){ pct = atof(getCmdOption(argv, argv+argc, "-pct" )); }
  bool sigma_reg = cmdOptionExists(argv, argv+argc, "-sigma_reg" );

  //number of iterations
  int max_its = 5;
  if( cmdOptionExists(argv, argv+argc, "-max_its" ) ){
    max_its = atoi( getCmdOption(argv, argv+argc, "-max_its" ) );
  }

  //maximum liklihood ratio
  float maxlr = -100;
  if( cmdOptionExists(argv, argv+argc, "-maxlr" ) ){
    maxlr = atof( getCmdOption(argv, argv+argc, "-maxlr" ) );
  }

  //How to interpret evidence from probabilities
  float bias = 0.5;
  if( cmdOptionExists(argv, argv+argc, "-bias" ) ){
    bias = atof( getCmdOption(argv, argv+argc, "-bias" ) );
  }
				
  //How many inner iterations of algorithm
  int inner_its = 1;
  if( use_panel ){ inner_its = 5; }
  if( cmdOptionExists(argv, argv+argc, "-inner_its" ) ){
    inner_its = atoi( getCmdOption(argv, argv+argc, "-inner_its" ) );
  }

  //How many steps of EM to start, too many can lose some information
  int EMits = 1;
  if( cmdOptionExists(argv, argv+argc, "-EMits" ) ){
    EMits = atoi( getCmdOption(argv, argv+argc, "-EMits" ) );
  }	

  //number of threads
  int num_threads = 1;
  if( cmdOptionExists(argv, argv+argc, "-num_threads" ) ){
    num_threads = atoi( getCmdOption(argv, argv+argc, "-num_threads" ) );
  }	

  omp_set_num_threads(num_threads);
  #pragma omp parallel
  {
    if(omp_get_thread_num() == 0){
      if( omp_get_num_threads() != 1){
	cerr << "there are " << omp_get_num_threads() << " threads" << endl;
	cerr << "parallelism gives only small speed up here, use extra cores for other windows" << endl;
      }
    }
  }
  
	int collapse = COLLAPSE_NONE;
	if( cmdOptionExists(argv, argv+argc, "-c" ) ){ 
		string ctype = (string)( getCmdOption(argv, argv+argc, "-c" ) );
		if( ctype == "none"){
			collapse = COLLAPSE_NONE; cerr << "collapse none" << endl;
		} else if (ctype == "snps") {
			collapse = COLLAPSE_SNPS; cerr << "collapse snps" << endl;
		} else if (ctype == "indels") {
			collapse = COLLAPSE_INDELS; cerr << "collapse indels" << endl;
		} else if (ctype == "any") {
			collapse = COLLAPSE_ANY; cerr << "collapse any" << endl;
		} else if (ctype == "some") {
			collapse = COLLAPSE_SOME; cerr << "collapse some" << endl;
		} else if (ctype == "both") {
			collapse = COLLAPSE_BOTH; cerr << "collapse both" << endl;
		} else {
			collapse = COLLAPSE_NONE; cerr << "collapse none" << endl;
		}
	} 

    //subsetting
	string regions;
	string use_regions;
	bool get_regions = cmdOptionExists(argv, argv+argc, "-r" );
	if(get_regions){ regions = getCmdOption(argv, argv+argc, "-r" ); }
	else{ usage(); }
	int minpos, maxpos;

	stringstream ss(regions);
	string chr, item;
	getline(ss, chr, ':'); getline(ss, item, ':');

	if(item==""){
		pair<int32_t, int32_t> fl = get_range(filename);
		minpos=fl.first;
		maxpos=fl.second;
	} else {
		stringstream ss2(item);
		string lt, rt;
		getline(ss2, lt, '-'); getline(ss2, rt, '-');			
		minpos = stoi(lt);
		maxpos = stoi(rt);
	}
	cerr << "Reading chromosome " << chr << " range " << minpos << "-" << maxpos << endl;
		
	htsFile *out_fh  = hts_open( out_filename.c_str(), out_type.c_str());
	bcf_hdr_t *new_hdr;
	
  	for(int pos=minpos; pos<maxpos; pos+=block){

		int pos_left = max(0,pos-overlap);
		int pos_right = pos+block+overlap;
		use_regions = chr + ":" + to_string( pos_left ) + "-" + to_string( pos_right );
		//cout << use_regions << endl;
		
		//input reader set up
		bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
		sr->collapse = collapse;
		
		sr->require_index = 1;
		if ( bcf_sr_set_regions(sr, use_regions.c_str(), false)<0 ){
		  cerr << "Failed to read the regions: " <<  use_regions << endl; exit(1);
		}
	
		if(!(bcf_sr_add_reader (sr, filename.c_str() ))){ 
		cerr << "Problem opening " << filename << endl; 
		bcf_sr_destroy(sr);	
		return 0;
		}
		if(use_panel){if(!(bcf_sr_add_reader (sr, pfilename.c_str() ))){ 
		  cerr << "Problem opening " << pfilename << endl; 
		  bcf_sr_destroy(sr);	
		  return 0;
		}}
	
		int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples
		float nm = 1.0/N; 
		float snm = sqrt(nm);
		int ngl = N*3;
		int ngl_arr = N*3;

		bcf1_t *line, *line2;///bcf/vcf line structure.
		vector< float* > gl;
		int M=0;
	
		vector<int> sites_to_skip;	//sites in panel but not sample
		int panel_uniq=0;
		int sample_number=0;
		int sample_uniq=0;
		int sample_ma=0;
		int common_var=0;
		int overlapping_indel=0;
		
  		int min_idx = 0;
		int max_idx = 0;
		
		//this goes before the while loop starts
		int32_t *PL; bool readPL = false;//stores the phred scaled likelihoods (if needed)
		int32_t *GT = (int32_t *)malloc(N*2*sizeof(int32_t));  //new int[N*2];
		
		//float GL is built on the fly per line.
		if ( bcf_hdr_id2int(sr->readers[0].header, BCF_DT_ID, "PL")!=-1 ){
			PL = (int32_t *)malloc(N*3*sizeof(int32_t)); //new int32_t[N*3]; 
			readPL = true;
		} else if ( bcf_hdr_id2int(sr->readers[0].header, BCF_DT_ID, "GL")<0 ) {
			cerr << "ERROR: neither FORMAT/GL nor FORMAT/PL were defined!"<<endl;
			exit(1);
		}    
	  


		while(bcf_sr_next_line (sr)) { 
			bool count = true;
			
			if( use_panel ){ count = (bcf_sr_has_line(sr,0) && bcf_sr_has_line(sr,1)); }
			if( count ){

			  //if using panel, in panel and sample
			  line =  bcf_sr_get_line(sr, 0);
			  bool read = ( line->n_allele == 2 );
			  if( use_panel && line->n_allele != 2 ){	//multiallelic in sample but not panel
				line2 = bcf_sr_get_line(sr, 1);
				read = ( line2->n_allele == 2 );
				++sample_ma;
			  }
			  if( use_panel && read ){
				 read = (line->pos+1 >= pos_left && line->pos+1 < pos_right);
				 if( !read ) ++overlapping_indel;
			  } 
			  if( read ){
				
				if(use_panel){
					if(line->pos+1 <= pos){ 	 ++min_idx; }
					if(line->pos+1 < pos+block){ ++max_idx; }
				}
				
				///reads GLs (or PLs if there are no GLs (or GTs if there are no PLs) )
				float *gl_farr= new float[N*3];//new set of GLs that will get pushed into vector.
				ngl = bcf_get_format_float(sr->readers[0].header, line, "GL", &gl_farr, &ngl_arr); 
				if(ngl<0) {//-1 or -3 FORMAT/GL is missing.
		
					if(ngl==-2) {
						cerr<<"FORMAT/GL had bad type at "<< line->pos+1<<" "<<ngl<<endl;
						exit(1);
					}
					ngl = bcf_get_format_int32(sr->readers[0].header, line, "PL", &PL, &ngl_arr);        
					if(ngl==-2) {
						cerr<<"FORMAT/PL had bad type at "<< line->pos+1<<" "<<ngl<<endl;
						exit(1);
					} else if(ngl<0) {
						if(bcf_get_genotypes(sr->readers[0].header, line, &GT, &ngl_arr)==2*N) {
							cerr << "WARNING: no genotype likelihoods available at "<< line->pos+1<<" (deriving from GT) "<<endl;

							for(int i=0; i<N; ++i) {
								for(int j=0;j<3;j++){ gl_farr[i*3+j] = -1000; }
								if(GT[i*2]!=bcf_gt_missing&&GT[i*2+1]!=bcf_gt_missing) {
									int g = bcf_gt_allele(GT[i*2])+bcf_gt_allele(GT[i*2+1]);
									gl_farr[i*3 + g] = 0.0;
								}
							}       
						} else {//flat likelihood.
							cerr << "WARNING: no genotype likelihoods or GT available at "<< line->pos+1<<" (flat GL assigned) "<<endl;
							for(int i=0; i<3*N; ++i) gl_farr[i] = 0.0;
						}
					}
					else {
						for(int i=0; i<3*N; ++i) gl_farr[i] = (float)(-0.1 * PL[i]);
					}
				}
				
				gl.push_back( gl_farr );

				++common_var;
				if( !use_panel ) ++M;
			} 		
			} 
			if( use_panel ){ 
					
				if( bcf_sr_has_line(sr,1) ){	//in panel
					
					if( !bcf_sr_has_line(sr,0) ){	//not in sample or missing
						sites_to_skip.push_back(M); //skip these?
						float *gl_farr= new float[N*3];
						for(int i=0; i<3*N; ++i){ gl_farr[i] = (float)1; }
						gl.push_back( gl_farr ); 
						++panel_uniq;
					} 
					
					line =  bcf_sr_get_line(sr, 1);
					if( (line->pos+1 >= pos_left && line->pos+1 < pos_right) ){
						if(line->pos+1 <= pos){      ++min_idx; }
						if(line->pos+1 < pos+block){ ++max_idx; }
						++M; 
					}
				} 
				if( bcf_sr_has_line(sr,0) ){ //in sample
					++sample_number; 
					if( !bcf_sr_has_line(sr,1) ){ //not in panel
						++sample_uniq;
					}
				} 
			}
	  }
		
	  if( verbose ){
		  if( pos == minpos ){ cerr << N << " samples in " << filename << endl; }
		  cerr << use_regions << endl;
		  cerr << "processing " << gl.size() << " variants on " << N << " individuals" << endl; 
		  if( use_panel ){ 
			cerr << M << "\tin panel" << endl;
			cerr << sample_number 	<< "\tin sample" << endl;
			cerr << common_var 		<< "\tcommon to both"  << endl;
			cerr << sample_ma 		<< "\tmulti-allelic sample sites that are biallelic in panel" << endl;
			cerr << overlapping_indel << "\tdeletion overlaps boundary" << endl;
			cerr << panel_uniq 		<< "\tunique to panel" << endl;
			cerr << sample_uniq 	<< "\tunique to sample" << endl;
			if( zero_missing ){ cerr << "\t Zeroing " << sites_to_skip.size() << " sites" << endl; }
			if( sample_ma ){ cerr << "Taking AF for first variant of multi-allelics!" << endl; }
			if( panel_uniq != 0 || sample_uniq != 0){ cerr << "Impute at intersection of sample and panel" << endl; }
		  } 
	  }
	  if(M != gl.size()){ cerr << "count error" << endl; exit(1); }
	  free( GT );
	  if( readPL ){ free(PL); }
	  bcf_sr_destroy(sr);	
  
		//GLs to GPs
		float norm;
		int sidx = 0;
		for(int i=0; i<M; ++i){
			for(int j=0; j<N; ++j){
				
				int j0 = 3*j;
				int j1 = j0+1;
				int j2 = j1+1;

				if( bcf_float_is_missing(gl[i][j0]) || bcf_float_is_missing(gl[i][j1]) || bcf_float_is_missing(gl[i][j2]) ){
					gl[i][j0] = 0; gl[i][j1] = 0; gl[i][j2] = 0;
				}

				gl[i][j0] = pow(10,gl[i][j0]);
				gl[i][j1] = pow(10,gl[i][j1]);
				gl[i][j2] = pow(10,gl[i][j2]);

				norm = gl[i][j0]+gl[i][j1]+gl[i][j2];
				if( norm != 0 ){
					gl[i][j0] /= norm;
					gl[i][j1] /= norm;
					gl[i][j2] /= norm;
				} 
							
				if( maxlr > 0 ){
					float maxl = gl[i][j0]; 
					if(gl[i][j1] > maxl){ maxl = gl[i][j1]; }
					if(gl[i][j2] > maxl){ maxl = gl[i][j2]; }
						
					gl[i][j0] = (maxl > maxlr*gl[i][j0]) ? maxl/maxlr : gl[i][j0]; 
					gl[i][j1] = (maxl > maxlr*gl[i][j1]) ? maxl/maxlr : gl[i][j1];
					gl[i][j2] = (maxl > maxlr*gl[i][j2]) ? maxl/maxlr : gl[i][j2];
					
					norm = gl[i][j0]+gl[i][j1]+gl[i][j2];
				
					gl[i][j0] /= norm;
					gl[i][j1] /= norm;
					gl[i][j2] /= norm;
				} 

				if(zero_missing &&  sidx < sites_to_skip.size() && i==sites_to_skip[sidx] ){ 
					gl[i][j0] = 0;
					gl[i][j1] = 0;
					gl[i][j2] = 0;
					++sidx;
				}
			}
		} 
		
		MatrixXf D(M, N);
		
		if(use_panel){ max_its = 1; } else { EMexpectation(gl, D, EMits); }

		VectorXf mu(M);
		VectorXf vcfmu(M);
		MatrixXf Sigma(M, M);
		MatrixXf lv(M,M);
		VectorXf sigs(M);
		MatrixXf Probs;
		if( use_panel ){
				
			//cerr << "Reading and resizing panel params" << endl;

			//input reader set up
			bcf_srs_t *sr2 =  bcf_sr_init() ; ///htslib synced reader.
			sr2->collapse = collapse;

			if(get_regions){
			sr2->require_index = 1;
				if ( bcf_sr_set_regions(sr2, use_regions.c_str(), false)<0 ){
				  cerr << "Failed to read the regions: " <<  use_regions << endl; exit(1);
				}
			}

			sr2->require_index = 1;
			if(!(bcf_sr_add_reader (sr2, pfilename.c_str() ))){ 
				cerr << "Problem opening " << pfilename << endl; 
				bcf_sr_destroy(sr2);	
				return 0;
			}

			float *one_ptr=NULL; //(float *)malloc(1*sizeof(float)); 
			float *left_cov=NULL; //(float *)malloc( M*sizeof(float));
			int32_t *left_cov_pos=NULL; //(int32_t *)malloc( M*sizeof(int32_t));
			float *right_cov=NULL; //(float *)malloc( M*sizeof(float));
			int32_t *right_cov_pos=NULL; //(int32_t *)malloc( M*sizeof(int32_t));
			
			/*int nval = 0;
			int left = 0;
			int right = 0;*/
				
			int idx = 0;
			lv = MatrixXf::Zero(M,M);
			while(bcf_sr_next_line (sr2)) { 
				if( bcf_sr_has_line(sr2,0) ){

					line =  bcf_sr_get_line(sr2, 0);
					if( line->n_allele == 2 && line->pos+1 >= pos_left && line->pos+1 < pos_right ){
						string ov_tag = "";
						if( pos+block < maxpos && line->pos+1 >= pos+block-overlap ){
							ov_tag="2";
						}
						int nval = 0;
						int ret = bcf_get_info_float(sr2->readers[0].header, line, "MU", &one_ptr, &nval);
						if( ret > 0 ) mu(idx) = one_ptr[0];
						nval = 0;
						string tag = "SIG" + ov_tag;
						ret = bcf_get_info_float(sr2->readers[0].header, line, tag.c_str(), &one_ptr, &nval);
						if( ret > 0 ) sigs(idx) = one_ptr[0];
						
						int left = 0;
						tag = "LCOR" + ov_tag;
						ret = bcf_get_info_float(sr2->readers[0].header, line, tag.c_str(), &left_cov, &left);
						left = 0;
						tag = "LCORP" + ov_tag;
						ret = bcf_get_info_int32(sr2->readers[0].header, line, tag.c_str(), &left_cov_pos, &left);
						for(int m=0; m<ret; ++m){ lv(idx, left_cov_pos[m] ) = left_cov[m]; }
						
						int right = 0;
						tag = "RCOR" + ov_tag;
						ret = bcf_get_info_float(sr2->readers[0].header, line, tag.c_str(), &right_cov, &right);
						right = 0;
						tag = "RCORP" + ov_tag;
						ret = bcf_get_info_int32(sr2->readers[0].header, line, tag.c_str(), &right_cov_pos, &right);
						//cout << M << " " << line->pos+1 << " " << ret << " " << right << " " << idx << " " << " "
						//<< right_cov_pos[ret-1] << " "
						//<< lv.rows() << " " << lv.cols() << endl; 
						for(int m=0; m<ret; ++m){ lv(idx, right_cov_pos[m] ) = right_cov[m]; }
						++idx;
					}
					
				}
			}
			bcf_sr_destroy(sr2);

			free( one_ptr );
			free(left_cov);
			free(left_cov_pos);
			free(right_cov);
			free(right_cov_pos);
						
			//bad guess
			prob_from_af(gl, mu, D);

			//zero the panel at missing sites
			if( zero_missing ){
			for(int i=0; i<sites_to_skip.size(); ++i){ 
				D.row( sites_to_skip[i] ) = VectorXf::Zero(N); 
				mu(  sites_to_skip[i] ) = 0;
				sigs(  sites_to_skip[i] ) = 0;
				lv.row( sites_to_skip[i] ) = VectorXf::Zero(N); 
				lv.col( sites_to_skip[i] ) = VectorXf::Zero(N); 
			}} 
			
		}
		
		for(int its = 0; its<max_its; ++its){

			 if(!use_panel) {
				
				for(int i=0; i<M; ++i){ mu(i) = D.row(i).sum()*nm; }
				//cerr << "calculating Sigma for iteration " << its << endl;				
				Sigma = D*D.transpose()*nm - mu*mu.transpose();

				for(int i=0; i<M;++i){ 
					if( sigma_reg ){
						Sigma(i,i) += lambda / (1.0 + exp( lambda2 * ( Sigma(i,i) - pct ) ) );
					} else {
						Sigma(i,i) += lambda;
					}
				}
			}
			
			MVNiterate(Sigma, mu, lv, sigs, !use_panel, inner_its, gl, D);

		}//end its 


	  //cerr <<"Printing genotypes to " << out_filename << endl; 
	  bcf_srs_t *reader = bcf_sr_init();
	  reader->collapse = collapse;

	  if(get_regions){
		reader->require_index = 1;
		if ( bcf_sr_set_regions(reader, use_regions.c_str(), false)<0 ){
		  cerr << "Failed to read the regions: " << use_regions << endl; exit(1);
		}
	  }

	  reader->require_index = 1;
	  if(!(bcf_sr_add_reader (reader, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }
	  if(use_panel){
		if(!(bcf_sr_add_reader (reader, pfilename.c_str() ))){ cerr << "Problem opening " << pfilename << endl; exit(1); }
	  }
	
	  vector< float* > gp;

	  if(pos == minpos ){
		  	new_hdr = bcf_hdr_dup(reader->readers[0].header);

			bcf_hdr_append(new_hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
			bcf_hdr_append(new_hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">");
			bcf_hdr_append(new_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
			if(glout) bcf_hdr_append(new_hdr, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">");
			
			bcf_hdr_write(out_fh, new_hdr);
	  }
	  if(glout){
		for(int i=0; i<M; ++i){ float *gp_arr = new float[N*3]; gp.push_back(gp_arr); }
		write_triple(Sigma, mu, lv, sigs, gl, D, gp);
	  }
	  //MVNiterate_probs(Sigma, mu, lv, sigs, !use_panel, inner_its, gl, D, gp);
	  //float lhood_sum = log_likelihood(0,D,gp,bias);

	  int nline = 0;
	  int nwrote = 0;
	  int32_t *gt_arr= new int32_t[N*2];
	  float *ds_arr= new float[N];
	  float *gp_arr;
	  if(glout) gp_arr = new float[N*3];
		
	  int ns = 0;
	  int np = 0;
	  int32_t *gq_arr = new int32_t[N];//FORMAT/GQ this is phred scale probablity the most likely genotype is wrong (good for filtering)
	  while(bcf_sr_next_line (reader)) { 
			
		bool count = true;
		if( use_panel ){ count = (bcf_sr_has_line(reader,0) && bcf_sr_has_line(reader,1)); }
		if( count ){	
		  line =  bcf_sr_get_line(reader, 0); 

		  if(line->n_allele == 2  && line->pos+1 >= pos_left && line->pos+1 < pos_right ){
					
		  for(int i=0; i<N; ++i){ 
								
			  if(D(nline,i) <= bias){	
				gt_arr[2*i] = bcf_gt_unphased(0); gt_arr[2*i+1] = bcf_gt_unphased(0);
			  } 
			  if(D(nline,i) > bias && D(nline,i) <= 2.0 - bias){
				gt_arr[2*i] = bcf_gt_unphased(0); gt_arr[2*i+1] = bcf_gt_unphased(1);
			  }
			  if(D(nline,i) > 2.0 - bias){
				gt_arr[2*i] = bcf_gt_unphased(1); gt_arr[2*i+1] = bcf_gt_unphased(1);
			  }							
							
			  ds_arr[i] = D(nline,i);
			  if( glout ){
				gp_arr[3*i] = gp[nline][3*i];
				gp_arr[3*i+1] = gp[nline][3*i+1];
				gp_arr[3*i+2] = gp[nline][3*i+2];
				float tmp_prob = gp[nline][3*i + bcf_gt_allele(gt_arr[2*i]) + bcf_gt_allele(gt_arr[2*i+1])];
				if(tmp_prob<0){ //avoids huge GQ values and -0 issues
					gq_arr[i] = (int32_t)(-10*log10(1 - pow(10,tmp_prob)));
				} else { 	      
					gq_arr[i] = 100;
				}
			  }	
		  } //end i
					 

		bcf_update_genotypes(new_hdr, line, (void*)gt_arr, 2*N);
		bcf_update_format_float(new_hdr,line,"DS",ds_arr, N);
		if( glout ) {
		  bcf_update_format_float(new_hdr,line,"GP",(void*)gp_arr, 3*N);
		  bcf_update_format_int32(new_hdr,line,"GQ",(void*)gq_arr, N);
		}
		bcf_unpack(line, BCF_UN_ALL);
		if( line->pos+1 >= pos && line->pos+1 < pos+block ){
			++nwrote; bcf_write(out_fh, new_hdr, line); 
		}
		++nline;
      } //n_alleles
    } //in both
		
    if( use_panel ){ 
      if( bcf_sr_has_line(reader,1) && !bcf_sr_has_line(reader,0) ){ 
		  line =  bcf_sr_get_line(reader, 0); 
		  if(line->n_allele == 2  && line->pos+1 >= pos_left && line->pos+1 < pos_right ){++nline;}
	  }
      if( !isec_only && bcf_sr_has_line(reader,0) && !bcf_sr_has_line(reader,1) ){ 
		  line =  bcf_sr_get_line(reader, 0); 
		  ++nwrote; bcf_write(out_fh, new_hdr, line);  
	  }
    }	
  }	
  //cerr << "wrote " << nwrote << " records" << endl;


  while(!gl.empty()) delete [] gl.back(), gl.pop_back();
  if(glout){
    while(!gp.empty()) delete [] gp.back(), gp.pop_back();	
    delete [] gp_arr;
  }
  bcf_sr_destroy(reader);
  
  delete [] gq_arr;
  delete [] gt_arr;
  delete [] ds_arr;
}//windows 

  hts_close(out_fh);
  bcf_hdr_destroy(new_hdr);	
	
  return 0;
}


