#define __STDC_LIMIT_MACROS
#include <iostream>
#include "marvin.hpp"
#include <queue>
#include <iomanip>
#include <stdint.h>
#include "htslib/hts.h"
#include <stdio.h>

using namespace std;
using namespace Eigen;

extern "C" {
#include "htslib/synced_bcf_reader.h"
}

void usage(){
	cerr << "Set up means and covariances for marvin_ref" << endl;	
	cerr << "Usage:" << endl;
	cerr << "./marvin_prep -f input_filename.vcf" << endl;
	cerr << "Expects input_filename.vcf to contain hard genotypes only with no missing" << endl;
	cerr << "\t -f : 		input vcf (mandatory)" << endl;
	cerr << "\t -o : 		site only vcf (sites.vcf)" << endl;
	cerr << "\t -O  : 		vcf type (w)" << endl;
	cerr << "\t -r: 		chromosome region (mandatory)" << endl;
	cerr << "\t -sigma_reg: use sigma regularization with params: lambda, lambda2, pct (false)" << endl;
	cerr << "\t -lambda: 	regularization parameter (0.06)" << endl;
	cerr << "\t -lambda2:	threshold steepness (4)" << endl;
	cerr << "\t -pct:		sig mid point (0.2)" << endl;
	cerr << "\t -b: 		block size (mandatory)" << endl;
	cerr << "\t -ov: 		overlap (mandatory)" << endl;
	cerr << "\t -max_ratio: only keep elements of cor s.t. |Cij/max(|Cij|)| > r (0.01)" << endl;
	exit(1);
}

//bcftools record is basically this but I don't want to have to keep track of pointers.
struct corr_container{
	
	int32_t rid;
	int32_t pos;
	int32_t qual;
	int32_t rlen;
	string id;
	string alleles;	
					
	float mu;
	float sig;
	
	vector<float> left_cov;
	vector<int32_t> left_cov_pos;
	vector<float> right_cov;
	vector<int32_t> right_cov_pos;
};

int main(int argc, char* argv[])
{
	
	//Input file - should be a multisample vcf
	if( !cmdOptionExists(argv, argv+argc, "-f" ) ){ usage(); }
	string filename = getCmdOption(argv, argv+argc, "-f" );
	//Site file - will be a site only vcf containing input file positions
	string sitefile = "sites.vcf";
	bool make_site_file = false;
	if( cmdOptionExists(argv, argv+argc, "-o" ) ){ make_site_file = true; sitefile = getCmdOption(argv, argv+argc, "-o" ); }
	string out_type = "w";
	if( cmdOptionExists(argv, argv+argc, "-O" ) ){ out_type += (string)( getCmdOption(argv, argv+argc, "-O" ) ); } 
	
	//window parameters
	int block = 0;
	if( cmdOptionExists(argv, argv+argc, "-b" ) ){ block = atoi(getCmdOption(argv, argv+argc, "-b" )); }
	else{ usage(); }
	int overlap = 0;
	if( cmdOptionExists(argv, argv+argc, "-ov" ) ){ overlap = atoi(getCmdOption(argv, argv+argc, "-ov" )); }
	else{ usage(); }
	float rat = 0.01;
	if( cmdOptionExists(argv, argv+argc, "-max_ratio" ) ){ rat = atof(getCmdOption(argv, argv+argc, "-max_ratio" )); }
		
	//regularization parameters
	bool sigma_reg = cmdOptionExists(argv, argv+argc, "-sigma_reg" );
	float lambda = 0.06;
	if( cmdOptionExists(argv, argv+argc, "-lambda" ) ){ lambda = atof(getCmdOption(argv, argv+argc, "-lambda" )); }
	float lambda2 = 4;
	if( cmdOptionExists(argv, argv+argc, "-lambda2" ) ){ lambda2 = atof(getCmdOption(argv, argv+argc, "-lambda2" )); }
	float pct = 0.2;
	if( cmdOptionExists(argv, argv+argc, "-pct" ) ){ pct = atof(getCmdOption(argv, argv+argc, "-pct" )); }

	//subsetting
	string regions;
	bool get_regions = cmdOptionExists(argv, argv+argc, "-r" );
	if(get_regions){ regions = getCmdOption(argv, argv+argc, "-r" ); }
	else{ usage(); }
	
	string use_regions;
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
	
	htsFile *out_fh  = hts_open( sitefile.c_str(), out_type.c_str());
	bcf_hdr_t *new_hdr;
		
	queue<corr_container> left_ov_buffer;
	queue<corr_container> right_ov_buffer;
	
	for(int pos=minpos; pos<maxpos; pos+=block){
		
		int pos_left = max(0,pos-overlap);
		int pos_right = pos+block+overlap;
		use_regions = chr + ":" + to_string( pos_left ) + "-" + to_string( pos_right );
		//Setup htslib reader
		bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
		sr->collapse = COLLAPSE_ANY;

		if(get_regions){
			sr->require_index = 1;
			if ( bcf_sr_set_regions(sr, use_regions.c_str(), false)<0 ){
				cerr << "Failed to read the region: " <<  use_regions << endl; exit(1);
			}
		}

		sr->require_index = 1;
		int	ridx = 0;
		if(!(bcf_sr_add_reader (sr, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }

		int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples
		
		float nm = 1.0/N; 
		float snm = sqrt(nm);
		
		int ngt = N*2;
		int ngt_arr = N*2;

		bcf1_t *line; ///bcf/vcf line structure.
		vector< int* > gt;
		int M=0;

		//read the input
		while(bcf_sr_next_line (sr)) { 
			
			line =  bcf_sr_get_line(sr, 0);
			
			//exclude haploid sites & multi allelics
			if( line->n_allele == 2 && line->pos+1 >= pos_left && line->pos+1 < pos_right){	
				
				int *gt_arr= new int[2*N];
				ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);   
				gt.push_back( gt_arr );
				++M;
			}
				
		}
		bcf_sr_destroy(sr);	
	
		MatrixXf D(M,N);
		VectorXf mu = VectorXf::Zero(M); //these are the allele freqs
		for(int i=0; i<M; ++i){
			for(int j=0; j<N; ++j){
				if(gt[i][2*j] == bcf_gt_missing || gt[i][2*j+1] == bcf_gt_missing ){
					cerr << "Found a missing site." << endl; exit(1); 
				}
				
				D(i,j) = (float)(bcf_gt_allele(gt[i][2*j]) + bcf_gt_allele(gt[i][2*j+1]));
				mu(i) += D(i,j);
					
			} 
			mu(i) *= nm; 
		}	
		//cerr << "Calculating Sigma." << endl;
		MatrixXf Sigma = D*D.transpose()*nm - mu*mu.transpose();
		
		for(int i=0; i<M;++i){ 
			if( sigma_reg ){
				Sigma(i,i) += lambda / (1.0 + exp( lambda2 * ( Sigma(i,i) - pct ) ) );
			} else {
				Sigma(i,i) += lambda;
			}
		} 		
		while(!gt.empty()) delete [] gt.back(), gt.pop_back(); 
		  
		//mu correction for very small sample sizes
		float theta = 0; for(int i=1; i<N; ++i){ theta += 1.0/(float)(i); }
		theta = (1.0/theta) / (N + (1.0/theta) );
		mu *= (1.0 - theta);
		mu += VectorXf::Ones(M) * (theta/2.0);
		
		//Only keep
		//cerr << "Calculating Sigma Inverse" << endl;
		MatrixXf SI = MatrixXf::Zero(M,M);
		VectorXf sigs = VectorXf::Ones(M);
		calc_lv_mem_noblock(Sigma, SI, sigs);	//do the necessary linear algebra

		bcf_srs_t *reader =  bcf_sr_init() ; ///htslib synced reader.
		reader->collapse = COLLAPSE_ANY;

		reader->require_index = 1;
		if ( bcf_sr_set_regions(reader, use_regions.c_str(), false)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
		if(!(bcf_sr_add_reader (reader, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }

		if(pos == minpos){

			bcf_hdr_t *hdr = bcf_hdr_dup(reader->readers[0].header);
			new_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
			bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures

			bcf_hdr_append(new_hdr, "##INFO=<ID=MU,Number=A,Type=Float,Description=\"Allele frequency\">");
			bcf_hdr_append(new_hdr, "##INFO=<ID=SIG,Number=A,Type=Float,Description=\"Variance\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=SIG2,Number=A,Type=Float,Description=\"Variance\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=LCOR,Number=G,Type=Float,Description=\"Left Correlation\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=LCOR2,Number=G,Type=Float,Description=\"Left Correlation\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=LCORP,Number=G,Type=Integer,Description=\"Left Correlation Position\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=LCORP2,Number=G,Type=Integer,Description=\"Left Correlation Position\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=RCOR,Number=G,Type=Float,Description=\"Right Correlation\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=RCOR2,Number=G,Type=Float,Description=\"Right Correlation\">"); 
			bcf_hdr_append(new_hdr, "##INFO=<ID=RCORP,Number=G,Type=Integer,Description=\"Right Correlation Position\">");
			bcf_hdr_append(new_hdr, "##INFO=<ID=RCORP2,Number=G,Type=Integer,Description=\"Right Correlation Position\">");
			string argstr = "##marvin_prep_Args= -r" + regions + " -b " + to_string(block) + " -ov" + to_string(overlap); 
			bcf_hdr_append(new_hdr, argstr.c_str()); 
					
			bcf_hdr_write(out_fh, new_hdr);
			bcf_hdr_destroy(hdr); 		

		}
		bcf1_t *rec = bcf_init1() ;

		if( M > 0 ){
			
			int idx = 0;
			int nline = 0;
			int nwrote = 0;
			int lb = 0;
			int rb = 0;
			int *gt_arr=(int *)malloc(N*2*sizeof(int));

			float mincor = SI.minCoeff();
			float maxcor = SI.maxCoeff();
			float mx = max( abs(mincor), abs(maxcor) );
				
			while(bcf_sr_next_line (reader)) { 
				line =  bcf_sr_get_line(reader, 0); 
				if( line->n_allele == 2 && line->pos+1 >= pos_left && line->pos+1 < pos_right){	
					//exclude haploid sites & multi allelics and annoying large indels
					
					bool buffer_right = false;
					if( pos+block < maxpos && line->pos+1 >= pos+block-overlap){ buffer_right = true; }
					bool buffer_left = false;
					if( pos > minpos && line->pos+1 < pos+overlap){ buffer_left = true; }
					
					/*cout << use_regions << "  " << line->pos+1;
					if( buffer_left ){ cout << " left buffer"; ++lb; }
					if( buffer_right ){ cout << " right buffer"; ++rb; }
					cout << endl;*/
					
					//bcf1_t *rec = bcf_init1();
					rec->rid = line->rid;
					rec->pos = line->pos;
					rec->qual = line->qual;
					rec->rlen = line->rlen;
			
					bcf_update_id(new_hdr, rec, line->d.id);
					string ref = line->d.allele[0];
					string alt = line->d.allele[1];
					string alleles = ref + "," + alt;
					bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
						
					float taf = mu(idx);
					float tsig = sigs(idx);
				
					int left = idx;
					int right = M-1-idx;
					
					vector<float> left_cov;
					vector<int32_t> left_cov_pos;
					int ct = 0; int st = 0;
					while(ct < left){ 
						if( abs( SI(idx, st+ct) ) > rat*mx ){
							left_cov.push_back( SI(idx, st+ct) );
							left_cov_pos.push_back( st+ct );
						}
						++ct;
					}
					vector<float> right_cov;
					vector<int32_t> right_cov_pos;
					ct = 0; st = idx+1;
					while(ct < right){ 
						if( abs( SI(idx, st+ct) ) > rat*mx ){
							right_cov.push_back( SI(idx, st+ct) );
							right_cov_pos.push_back( st+ct );
						}
						++ct;
					}
					
								

					if( buffer_left || buffer_right ){ 
						corr_container cc;
						
						cc.rid = line->rid;
						cc.pos = line->pos;
						cc.qual = line->qual;
						cc.rlen = line->rlen;
						cc.id =  line->d.id;
						
						string ref = line->d.allele[0];
						string alt = line->d.allele[1];
						cc.alleles = ref + "," + alt;
					
									
						cc.mu = taf;
						cc.sig = tsig;
						
						cc.left_cov = left_cov;
						cc.left_cov_pos = left_cov_pos;
						cc.right_cov = right_cov;
						cc.right_cov_pos = right_cov_pos;
					
						if( buffer_left ){ left_ov_buffer.push(cc); }
						else{ right_ov_buffer.push(cc); }
					}
					
					if( left_ov_buffer.size() > 0 && left_ov_buffer.size() == right_ov_buffer.size() ){
						while( !left_ov_buffer.empty() ){
							bcf_clear1(rec);
							
							rec->rid = left_ov_buffer.front().rid;
							rec->pos = left_ov_buffer.front().pos;
							rec->qual = left_ov_buffer.front().qual;
							rec->rlen = left_ov_buffer.front().rlen;
					
							bcf_update_id(new_hdr, rec, left_ov_buffer.front().id.c_str());
							bcf_update_alleles_str(new_hdr, rec, left_ov_buffer.front().alleles.c_str());
					
							bcf_update_info_float(new_hdr,rec,"MU", &left_ov_buffer.front().mu, 1);	
							bcf_update_info_float(new_hdr,rec,"SIG",&left_ov_buffer.front().sig,1);
							bcf_update_info_float(new_hdr,rec,"SIG2",&right_ov_buffer.front().sig,1);
							
							bcf_update_info_float(new_hdr,rec,"LCOR",&left_ov_buffer.front().left_cov[0],left_ov_buffer.front().left_cov.size() );
							bcf_update_info_int32(new_hdr,rec,"LCORP",&left_ov_buffer.front().left_cov_pos[0],left_ov_buffer.front().left_cov_pos.size() );
							bcf_update_info_float(new_hdr,rec,"RCOR",&left_ov_buffer.front().right_cov[0],left_ov_buffer.front().right_cov.size() );
							bcf_update_info_int32(new_hdr,rec,"RCORP",&left_ov_buffer.front().right_cov_pos[0],left_ov_buffer.front().right_cov_pos.size() );
							
							bcf_update_info_float(new_hdr,rec,"LCOR2",&right_ov_buffer.front().left_cov[0],right_ov_buffer.front().left_cov.size() );
							bcf_update_info_int32(new_hdr,rec,"LCORP2",&right_ov_buffer.front().left_cov_pos[0],right_ov_buffer.front().left_cov_pos.size() );
							bcf_update_info_float(new_hdr,rec,"RCOR2",&right_ov_buffer.front().right_cov[0],right_ov_buffer.front().right_cov.size() );
							bcf_update_info_int32(new_hdr,rec,"RCORP2",&right_ov_buffer.front().right_cov_pos[0],right_ov_buffer.front().right_cov_pos.size() );
						
							//cout << "left  L=" << left_ov_buffer.front().left_cov.size() << " R=" << left_ov_buffer.front().right_cov.size() << endl;
							//cout << "right L=" << right_ov_buffer.front().left_cov.size() << " R=" << right_ov_buffer.front().right_cov.size() << endl;
							bcf_unpack(rec, BCF_UN_ALL);	
						
							bcf_write1(out_fh, new_hdr, rec); //++nwrote;
							bcf_clear1(rec);
						
							left_ov_buffer.pop();
							right_ov_buffer.pop();
						}
					}
					
					if( !buffer_left && !buffer_right ){
						
						bcf_update_info_float(new_hdr,rec,"MU", &taf, 1);	
						bcf_update_info_float(new_hdr,rec,"SIG",&tsig,1);
						bcf_update_info_float(new_hdr,rec,"LCOR",&left_cov[0],left_cov.size() );
						bcf_update_info_int32(new_hdr,rec,"LCORP",&left_cov_pos[0],left_cov_pos.size() );
						bcf_update_info_float(new_hdr,rec,"RCOR",&right_cov[0],right_cov.size() );
						bcf_update_info_int32(new_hdr,rec,"RCORP",&right_cov_pos[0],right_cov_pos.size() );
						bcf_unpack(rec, BCF_UN_ALL);	
						//cout << "L=" << left_cov.size() << " R=" << right_cov.size() << endl;
					
						bcf_write1(out_fh, new_hdr, rec); //++nwrote;
						bcf_clear1(rec);
					}
					++idx;
				}
				++nline;
			}
			free(gt_arr);
		} //M>0
		bcf_destroy1(rec);
		bcf_sr_destroy(reader);
		
		if( pos == minpos ){ cerr << N << " samples in " << filename << endl; }
		cerr << M << " variants in " << use_regions << endl;
		
		
	}//blocks
	hts_close(out_fh);
	bcf_hdr_destroy(new_hdr); 		
}

