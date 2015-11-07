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
	cerr << "Expects input_filename.vcf to contain hard genotypes only" << endl;
	cerr << "\t -f : 		input vcf" << endl;
	cerr << "\t -o : 		site only vcf" << endl;
	cerr << "\t -O  : 		vcf type" << endl;
	cerr << "\t -om: 		output allele frequencies" << endl;
	cerr << "\t -os: 		(optional) output covariance matrix" << endl;
	cerr << "\t -ow: 		output MVN update matrix" << endl;
	cerr << "\t -ov: 		output MVN variances" << endl;
	cerr << "\t -sigma_reg: use sigma regularization with params lambda, lambda2, pct" << endl;
	cerr << "\t -lambda: 	regularization parameter (0.06)" << endl;
	cerr << "\t -lambda2:	threshold steepness (4)" << endl;
	cerr << "\t -pct:		sig mid point (0.2)" << endl;
	cerr << "\t -r: 		chromosome region (default all)" << endl;
	exit(1);
	
}

int main(int argc, char* argv[])
{
	
	//Input file - should be a multisample vcf
	if( !cmdOptionExists(argv, argv+argc, "-f" ) ){ usage(); }
	char* filename = getCmdOption(argv, argv+argc, "-f" );
	
	//Site file - will be a site only vcf containing input file positions
	string sitefile = "sites.vcf";
	bool make_site_file = false;
	if( cmdOptionExists(argv, argv+argc, "-o" ) ){ make_site_file = true; sitefile = getCmdOption(argv, argv+argc, "-o" ); }
	string out_type = "w";
	if( cmdOptionExists(argv, argv+argc, "-O" ) ){ out_type += (string)( getCmdOption(argv, argv+argc, "-O" ) ); } 
	
	//Output matrices - mu=allele freq, sigma=covariance, omega=matrix required by marvin.cpp, var=list of variances
	string mfilename = "mu.dat";
	string sfilename = "sigma.dat"; bool print_sig = false;
	string wfilename = "omega.dat"; 
	string vfilename = "var.dat"; 
	
	if( cmdOptionExists(argv, argv+argc, "-om" ) ){ mfilename = getCmdOption(argv, argv+argc, "-om" ); }
	if( cmdOptionExists(argv, argv+argc, "-os" ) ){ print_sig = true; sfilename = getCmdOption(argv, argv+argc, "-os" ); }
	if( cmdOptionExists(argv, argv+argc, "-ow" ) ){ wfilename = getCmdOption(argv, argv+argc, "-ow" ); }
	if( cmdOptionExists(argv, argv+argc, "-ov" ) ){ vfilename = getCmdOption(argv, argv+argc, "-ov" ); }
		
	//regularization parameter
	float lambda = 0.06;
	if( cmdOptionExists(argv, argv+argc, "-lambda" ) ){ lambda = atof(getCmdOption(argv, argv+argc, "-lambda" )); }
	float lambda2 = 4;
	if( cmdOptionExists(argv, argv+argc, "-lambda2" ) ){ lambda2 = atof(getCmdOption(argv, argv+argc, "-lambda2" )); }
	float pct = 0.2;
	if( cmdOptionExists(argv, argv+argc, "-pct" ) ){ pct = atof(getCmdOption(argv, argv+argc, "-pct" )); }
  	bool sigma_reg = cmdOptionExists(argv, argv+argc, "-sigma_reg" );

	//subsetting
	char* regions;
	bool get_regions = cmdOptionExists(argv, argv+argc, "-r" );
	if(get_regions) regions = getCmdOption(argv, argv+argc, "-r" );
	
	
	//Setup htslib reader
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	sr->collapse = COLLAPSE_ANY;

	if(get_regions){
		sr->require_index = 1;
		if ( bcf_sr_set_regions(sr, regions, false)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
	}

	sr->require_index = 1;
	int	ridx = 0;
	if(!(bcf_sr_add_reader (sr, filename ))){ cerr << "Problem opening " << filename << endl; exit(1); }

	
	int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples
	
	float nm = 1.0/N; 
	float snm = sqrt(nm);
	
	int ngt = N*2;
	int ngt_arr = N*2;
	cerr << N << " samples in " << filename << endl;
	
	bcf1_t *line;///bcf/vcf line structure.
	vector< int* > gt;
	int M=0;

	//set up htslib for site only vcf output
	htsFile *out_fh;
	bcf_hdr_t *new_hdr;
	bcf1_t *rec;
	if( make_site_file ){ 
		cerr << "Printing coefficients to " << sitefile << endl; 
		///creates a new subsetted header (with 0 samples) from src_header
		new_hdr = bcf_hdr_subset(sr->readers[0].header,0,NULL,NULL); 
		bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures
		rec = bcf_init1() ;
		out_fh  = hts_open(sitefile.c_str(), out_type.c_str());
		bcf_hdr_write(out_fh, new_hdr);
	}

	//read the input and write the site file
	while(bcf_sr_next_line (sr)) { 
		
		line =  bcf_sr_get_line(sr, 0);
		
		if( line->n_allele == 2){	//exclude haploid sites & multi allelics
			
			if( make_site_file ){ 
				rec->rid = line->rid;
				rec->pos = line->pos;
				rec->qual = line->qual;

				bcf_update_id(new_hdr, rec, line->d.id);
				string ref = line->d.allele[0];
				string alt = line->d.allele[1];
				string alleles = ref + "," + alt;
				bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
				
				bcf_unpack(rec, BCF_UN_ALL);				
				bcf_write1(out_fh, new_hdr, rec);
				bcf_clear1(rec);
			}
			
			int *gt_arr= new int[2*N];
			ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);   
			gt.push_back( gt_arr );
			++M;
		}
			
	}
	
	if( make_site_file ){ bcf_destroy1(rec); hts_close(out_fh); bcf_hdr_destroy(new_hdr); }
	bcf_sr_destroy(sr);	

	cerr << M << " variants on " << N << " individuals." << endl;
	
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
	cerr << "Calculating Sigma." << endl;
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
	
	//write mu file
    ofstream f(mfilename.c_str(), ios::binary);
	f.write((char *)mu.data(), sizeof(typename VectorXf::Scalar)*M);
	f.close();

	if(print_sig){	//optional as marvin does not require sigma
		cerr << "Printing Sigma" << endl;	
	
		ofstream outs(sfilename.c_str(), ios::binary);
		outs.write((char *)Sigma.data(), sizeof(typename MatrixXf::Scalar)*M*M);
		outs.close();
	}
	
	cerr << "Calculating Sigma Inverse" << endl;
	MatrixXf SI = MatrixXf::Zero(M,M);
	VectorXf sigs = VectorXf::Ones(M);
	calc_lv_mem_noblock(Sigma, SI, sigs);	//do the necessary linear algebra

	//print variances
	ofstream outv(vfilename.c_str(), ios::binary);
	outv.write((char *)sigs.data(), sizeof(typename VectorXf::Scalar)*M);
	outv.close();

	//print 'inverses'
	ofstream outi(wfilename.c_str(), ios::binary);
	outi.write((char *)SI.data(), sizeof(typename MatrixXf::Scalar)*M*M);
	outi.close();
	

	
}

