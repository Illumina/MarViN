#ifndef EM_HPP
#define EM_HPP

#include <iostream>
#include <Eigen/Dense>
#include <queue>
#include <random>
#include <iomanip>

using namespace std;
using namespace Eigen;

void calcDosageE(const vector< float* >& gl, Ref<MatrixXf> Ynew){
	for(int i=0; i<Ynew.rows(); ++i){
	for(int j=0; j<Ynew.cols(); ++j){
		int j0 = 3*j; int j1 = j0+1; int j2 = j1+1;
		Ynew(i,j) = ( gl[i][j1] + 2*gl[i][j2] )/( gl[i][j0] + gl[i][j1] + gl[i][j2] );
	}}
}

void EMexpectation(const vector< float* >& gl, Ref<MatrixXf> Ynew, int max_its=1) {

	int M = Ynew.rows();
	int N = Ynew.cols();

	calcDosageE(gl, Ynew);

	//EM algorithm Expectation
	for(int it=0; it<max_its; ++it){

		VectorXf AF = VectorXf::Zero(M); 
		for(int i=0; i<M; ++i){ 
			for(int j=0; j<N; ++j){
				AF(i) += Ynew(i,j);
			}
			AF(i) /= (2*N);
		}
		
		//Calculate P(g | reads) ~ P(reads | g) P(g)
		for(int i=0; i<M; ++i){
			
			float pg0 = (1.0 - AF(i))*(1.0 - AF(i)); 	
			float pg1 = 2*(1.0 - AF(i))*AF(i); 					  
			float pg2 = AF(i)*AF(i); 							
				
			for(int j=0; j<N; ++j){
				
				float tp0 = gl[i][3*j] *   pg0; 
				float tp1 = gl[i][3*j+1] * pg1;  
				float tp2 = gl[i][3*j+2] * pg2;  

				Ynew(i,j) = ( tp1 + 2*tp2 )/(tp0 + tp1 + tp2);
			}
		}
		
	}


}

void prob_from_af(const vector< float* >& gl, const Ref<VectorXf> AF, Ref<MatrixXf> D){
	
int M = D.rows();
int N = D.cols();

//Calculate P(g | reads) ~ P(reads | g) P(g)
for(int i=0; i<M; ++i){
	float p = AF(i)/2.0;
	float pg0 = (1.0 - p)*(1.0 - p); 	
	float pg1 = 2*p*(1.0 - p); 					  
	float pg2 = p*p; 
			
	for(int j=0; j<N; ++j){

		float tp0 = gl[i][3*j] *   pg0; 
		float tp1 = gl[i][3*j+1] * pg1;  
		float tp2 = gl[i][3*j+2] * pg2;  
		
		D(i,j) = ( tp1 + 2*tp2 )/(tp0 + tp1 + tp2);
		
	}
}

	
}

void random_phase(Ref<MatrixXf> D, Ref<MatrixXf> D2, float rand_sub=0.1){
	
int M = D.rows();
int N = D.cols();
	cout << rand_sub << endl;
	for(int i=0; i<M; ++i){
	for(int j=0; j<N; ++j){
		double min = abs(D(i,j)); int kmin = 0;
		for(int k=1; k<3; ++k){ if( abs(D(i,j) - k)  < min ){ kmin = k; min = abs(D(i,j) - k); } }

		if(kmin == 0){ D(i,j) = 0; D2(i,j) = 0; }
		else if(kmin == 2){ D(i,j) = 1; D2(i,j) = 1; }
		else if(kmin == 1){ 
			//D(i,j) = rand() % 2; D2(i,j) = 1 - D(i,j); 
			if( rand() % 2 ){
				D(i,j) = 0.5 - rand_sub; D2(i,j) = 0.5 + rand_sub;
			} else {
				D(i,j) = 0.5 + rand_sub; D2(i,j) = 0.5 - rand_sub;
			}
		} 
		else { cerr << "Dosage error" << endl; exit(1); }
	}}
	
}


#endif
