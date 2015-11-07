#ifndef MVN_HPP
#define MVN_HPP

#define EIGEN_DONT_PARALLELIZE

#include <iostream>
#include <Eigen/Dense>
#include <queue>
#include <iomanip>
#include <time.h>       
#include <fstream>
#include <string>
#include <omp.h>

using namespace std;
using namespace Eigen;

//No dependency on boost etc. Very simple command line reading
char* getCmdOption(char ** begin, char ** end, const string & option)
{
    char ** itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return NULL;
}
bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return find(begin, end, option) != end;
}

//un-normalised Gaussian function
float Normal(float x, float mu, float sigma){
	return exp( -(x - mu) * (x - mu) / (2.0 * sigma) ); //normalise by hand anyway
}

//Compute all inverses and dot products in one step.
void calc_lv_mem_noblock(
Ref<MatrixXf> Sigma, 
Ref<MatrixXf> lv, 
Ref<VectorXf> sigs
){
	
	int M = Sigma.rows();
	//sub matrixes to do all the Schur compliments
	MatrixXf inv = Sigma.partialPivLu().inverse(); 
	MatrixXf Sigp = Sigma; Sigp.diagonal() = VectorXf::Zero(M);
	Sigp *= inv; 
	lv = Sigp - Sigp.diagonal().cwiseQuotient( inv.diagonal() ).asDiagonal() * inv;
	sigs = Sigma.diagonal() - (lv * Sigma).diagonal();
	
}

void MVNiterate(Ref<MatrixXf> Sigma, 
Ref<VectorXf> mu, 
Ref<MatrixXf> lv,
Ref<VectorXf> sigs, 
bool calc_invs, 
int tot,
const vector< float* >& gl,
Ref<MatrixXf> D
 ){
	
	int M = D.rows();
	int N = D.cols();
	
	if(calc_invs){
		cerr << "calculating inverses " << endl;	
		calc_lv_mem_noblock(Sigma, lv, sigs); 
	} 

	cerr << "Iterating over samples" << endl;
	#pragma omp parallel for	//might as well use extra cores for other windows but here anyway for some cases
	for(int ind=0; ind<N; ++ind){

		vector<float> phap(3);

		for(int it=0; it<tot; ++it){
				
			VectorXf mus = mu + lv * (D.col(ind) - mu); //calculate one column 

			for(int i=0; i<M; ++i){	

				float norm = 0;
				for(int k=0; k<3; ++k){
					phap[k] = gl[i][3*ind+k] * Normal(k, mus(i), sigs(i) ); 				
					norm += phap[k];
				}

				if( phap[2] > phap[1] && phap[0] > phap[1] ){
					(phap[2] > phap[0]) ? D(i,ind) = 2: D(i,ind) = 0;
				} else {
					(norm>0) ? D(i,ind) = (phap[1] + 2.0 * phap[2])/norm : D(i,ind) = 0;
				}
				//D(i,ind) = distance(phap.begin(), max_element(phap.begin(), phap.end()));	//Max prob instead of expectation

			} //end i

		} //end it 
		
	}//end ind
	
}

//Routine for re-calculating correct final genotype probabilities
void write_triple(Ref<MatrixXf> Sigma, 
Ref<VectorXf> mu, 
Ref<MatrixXf> lv,
Ref<VectorXf> sigs, 
const vector< float* >& gl,
Ref<MatrixXf> D,
vector< float* >& gp
 ){
	int M = D.rows();
	int N = D.cols();
	cerr << "Calculating final genotype probabilities" << endl;
	//No inv calculation, already have it from last iteration

	cerr << "Iterating over samples" << endl;
	//#pragma omp parallel for
	for(int ind=0; ind<N; ++ind){

		//no multiple iterations - just recomputing probabilities
		VectorXf mus = mu + lv * (D.col(ind) - mu);

		for(int i=0; i<M; ++i){	

			float n0 = Normal(0, mus(i), sigs(i) );
			float n1 = Normal(1, mus(i), sigs(i) );
			float n2 = Normal(2, mus(i), sigs(i) );
				
			vector<float> phap(3,0);
			phap[0] = gl[i][3*ind] * n0;
			phap[1] = gl[i][3*ind+1] * n1;
			phap[2] = gl[i][3*ind+2] * n2;

			float norm = 0; for(int k=0; k<3; ++k) norm += phap[k];
			//log likelihood, minimum is -5.
			gp[i][3*ind] =   max((float)-5.0, round( log10( phap[0]/norm )*(float)1.0e5 )/(float)1.0e5 );
			gp[i][3*ind+1] = max((float)-5.0, round( log10( phap[1]/norm )*(float)1.0e5 )/(float)1.0e5 );
			gp[i][3*ind+2] = max((float)-5.0, round( log10( phap[2]/norm )*(float)1.0e5 )/(float)1.0e5 );

		} //end i	

	}//end ind
	
}
 
#endif
