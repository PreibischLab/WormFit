package coherent.point.drift;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.la4j.Matrix;
import org.la4j.matrix.dense.Basic2DMatrix;

import ij.ImagePlus;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;

public class ApacheCPD {
	
	Img<FloatType> img; 
	ImagePlus imp; 
	
	// TODO: read this parameters from file or from variable
	final int D; // dimensionality of the point set
	final int N; // # of points in the first point set
	final int M; // # of points in the second point set

	double sigma2; // TODO: 2 for squared

	// parameters
	final double w; 	// amount of noise [0, 1]
	final double beta;  // Briefly speaking, parameter defines the model of the smoothness regularizer (width of smoothing Gaussian filter).
	final double lambda; // trade-off between the goodness of maximum likelihood fit and regularization.

	// TODO: convert stuff to the HPC
	// how to make this final
	 RealMatrix mX; // N D
	 RealMatrix mY; // M D
	 RealMatrix mW; // M D
	 RealMatrix mG; // M M
	 RealMatrix mP; // M N
	 RealMatrix mT; // M D transformation matrix

	// these guys are necessary for proper visualization
	double[] translate;
	double[] scale;
	double[] sigma;
		
	final int maxIteration;

	// constructor
	public ApacheCPD(Img<FloatType> img, int N, int M, int D, double w, double beta, double lambda, int maxIteration) {
		this.img = img;

		this.w = w;
		this.beta = beta;
		this.lambda = lambda;
		
		this.maxIteration = maxIteration;

		this.D = D; // dimensionality of the point set
		this.N = N; // # of points in the first point set
		this.M = M; // # of points in the second point set
		
		mX.createMatrix(N, D);
		mY.createMatrix(M, D);
		mW.createMatrix(M, D);
		mG.createMatrix(M, M);
		mP.createMatrix(M, N);
		mT.createMatrix(M, D);
		
		// TODO: if nothing is working set matrix mW to 0 explicitly
		// TODO: make this things calculated automatically
		translate = new double[] { 250, 250 };
		scale = new double[] { 150, 150 };
		sigma = new double[] { 3, 3 };
	}
	
	
	protected void normalize(RealMatrix M){
		int numColumns = M.getColumnDimension();
		int numRows = M.getRowDimension();
		
		double [] sumOverColumns = new double [numColumns];
		
		
		for (int d = 0; d < numColumns; ++d){
			double[] column = M.getColumn(d);
			sumOverColumns[d]= 0;
			for(int j = 0; j < numRows; ++j){
				sumOverColumns[d] += column[j];
			}
			
			M.getColumnVector(d).mapDivide(sumOverColumns[d]);
			
		}
		
		
	}
	
	
	public static void main(String[] args){
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}