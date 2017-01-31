package coherent.point.drift;

import org.apache.commons.math3.linear.MatrixUtils;
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

	public ApacheCPD(){
		mX = null; // N D
		mY = null;  // M D
		mW = null; // M D
		mG = null; // M M
		mP = null; // M N
		mT = null; // M D transformation matrix

		D = 1;
		M = 1;
		N = 1; 

		w = 1;
		beta = 1;
		lambda = 1;

		maxIteration = 1;
	}


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

		System.out.println("rows: " + numRows + " cols: " + numColumns);


		double [] sumOverColumns = new double [numColumns];

		printMatrix(M, 3);

		// zero mean
		for (int d = 0; d < numColumns; ++d){
			sumOverColumns[d]= 0;
			for(int j = 0; j < numRows; ++j){
				
				sumOverColumns[d] += M.getEntry(j, d);
			}
			sumOverColumns[d] /= numRows;
		}
		
		for (int d = 0; d < numColumns; ++d ){
			M.setColumnVector(d, M.getColumnVector(d).mapSubtract(sumOverColumns[d]));		
		}
		
		// std = 1
		// TODO: std is calculated in a wroong way
		double scaleValue = 0; 
		for (int  i = 0; i < numRows; ++i){
			for (int d =0; d < numColumns; ++d){
				//  TODO: looks correct but might be not, have a look into this again  
				scaleValue += (M.getRowVector(i).ebeMultiply(M.getRowVector(i))).getL1Norm();
			}
		}

		scaleValue /= numRows;
		scaleValue = Math.sqrt(scaleValue);

		M.scalarMultiply(1./scaleValue);
		// printMatrix(M, numRows);

	}


	// debug 
	public void printMatrix(RealMatrix M, int RowDimension){

		for (int j = 0; j < RowDimension; j++ ){
			for (int i = 0; i < M.getColumnDimension(); i++ ){
				System.out.print(M.getEntry(j, i) + " ");
			}
			System.out.println();
		}

	}


	public void test(){
		RealMatrix matrix = MatrixUtils.createRealMatrix(10, 3);
		
		for (int j = 0; j < matrix.getColumnDimension(); ++j){
			for (int i = 0; i < matrix.getRowDimension(); ++i){
				matrix.setEntry(i, j, i+j);
			}
		}
		
		// M.setEntry(j, d, d + 1);
		normalize(matrix);
	}

	public static void main(String[] args) {
		new ApacheCPD().test();
	}

}

