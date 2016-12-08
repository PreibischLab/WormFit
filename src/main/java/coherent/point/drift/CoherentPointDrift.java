package coherent.point.drift;

import org.la4j.matrix.dense.Basic2DMatrix;

public class CoherentPointDrift {

	int D; // dimensionality of the point set
	int N; // # of points in the first point set
	int M; // # of points in the second point set
	float [] X; // first point set N x D
	float [] Y; // second point set M x D
	float [] I; // (in theory) identity matrix

	float [] W; // matrix of coefficients M x D 
	float sigma2; // TODO: 2 for squared

	float [] G; 

	float [] P; // probability matrix 

	// parameters
	float w; // amount of noise [0, 1]
	float beta; // Briefly speaking, parameter   defines the model of the smoothness regularizer (width of smoothing Gaussian filter in (20)).
	float lambda; //Parameter   represents the trade-off between the goodness of maximum likelihood fit and regularization.
	
	// TODO: convert stuff to the HPC
	Basic2DMatrix mX; // N D
	Basic2DMatrix mY; // M D 
	Basic2DMatrix mI; // 
	Basic2DMatrix mW; // M D 
	Basic2DMatrix mG; // M M
	Basic2DMatrix mP; // M N
	
	// straight to the non-rigid registration

	// constructor
	public CoherentPointDrift(){
	}

	public void runNonRigidRegistration(){
		// init part  
		W = new float[M*D];

		// TODO: part to read X and Y
		X = new float[N*D];
		Y = new float[M*D];

		float res = 0;

		for (int i = 0; i < N; ++i){
			for(int j = 0; j < M; ++j){
				for (int d = 0; d < D; ++d){
					res += (X[d + i*D] - Y[d + j*D])*(X[d + i*D] - Y[d + j*D]);
				}
			}
		}
		sigma2 = 1/ (D*N*M) * res;

		w = 0.1f; 
		beta = 1;
		lambda = 1;

		G = new float[M*M];

		// TODO: it is possible to make it more efficient 
		// dure to the symmetry
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				float val = 0;
				for(int d = 0; d < D; ++d)
					val += (Y[d + i*D] - Y[d + j*D])*(Y[d + i*D] - Y[d + j*D]);
				G[j + M*i] = (float)Math.exp(-val/(2*beta*beta));
			}
		}


		int counter = 0;
		while(counter++ < 10){

			P = new float [M*N]; 

			// if this thing is correct then you are a freaking hacker! 
			for(int i = 0; i < M; i++){
				for(int j = 0; j < N; ++j){
					float den = 0;
					float nom = 0;	
					float[] tmp = new float [D];

					// lower part first
					for (int k = 0; k < M; ++k){
						multVecMatrix(G, W, tmp, k);				
						for (int d = 0; d < D; ++d){
							nom += (float)Math.exp(- (X[j*D + d] - (Y[k*D + d] + tmp[d]))*X[j*D + d] - (Y[k*D + d] + tmp[d])/ (2*sigma2));
						}
					}

					nom += w/(1 - w)*(2 * Math.PI *sigma2)*D/2*M/N;

					// upper part after that
					multVecMatrix(G, W, tmp, i);				
					for (int d = 0; d < D; ++d){
						den = (float)Math.exp(- (X[j*D + d] - (Y[i*D + d] + tmp[d]))*X[j*D + d] - (Y[i*D + d] + tmp[d])/ (2*sigma2));
					}

					P[i*M + j] = den/nom;
				}
			}
			
			
		}

	}


	// y = G(m, :) * W
	public void multVecMatrix(float[] G, float[] W, float [] y, int m){
		for(int d = 0; d < D; ++d){
			float res = 0;
			for (int i = 0; i < M; ++i){
				res += G[m*M + i]*W[i*D + d];
			}
			y[d] = res;
		}
	}

	public void testLAPACK(){
		
	}
	

	public static void main(String[] args){

	}
}
