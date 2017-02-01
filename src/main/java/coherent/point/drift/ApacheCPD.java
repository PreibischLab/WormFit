package coherent.point.drift;

import org.apache.commons.math3.geometry.euclidean.oned.Vector1D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.RealVector;
import org.la4j.Matrix;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.operation.VectorVectorOperation;

import ij.ImagePlus;
import mdbtools.dbengine.functions.Upper;
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
	final RealMatrix mX; // N D
	final RealMatrix mY; // M D
	final RealMatrix mW; // M D
	final RealMatrix mG; // M M
	final RealMatrix mP; // M N
	final RealMatrix mT; // M D transformation matrix

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

		mX = MatrixUtils.createRealMatrix(N, D);
		mY = MatrixUtils.createRealMatrix(M, D);
		mW = MatrixUtils.createRealMatrix(M, D);
		mG = MatrixUtils.createRealMatrix(M, M);
		mP = MatrixUtils.createRealMatrix(M, N);
		mT = MatrixUtils.createRealMatrix(M, D);

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
		// zero mean
		for (int d = 0; d < numColumns; ++d)
			sumOverColumns[d]= 0;

		ColumnSumVisitor columnSumVisitor = new ColumnSumVisitor(sumOverColumns);
		columnSumVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		M.walkInOptimizedOrder(columnSumVisitor);
		columnSumVisitor.end();

		for(int d = 0; d < numColumns; ++d)		
			sumOverColumns[d] /= numRows;

		ColumnSubtractValueVisitor columnSubtractValueVisitor = new ColumnSubtractValueVisitor(sumOverColumns);
		columnSubtractValueVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		M.walkInOptimizedOrder(columnSubtractValueVisitor);
		columnSubtractValueVisitor.end();

		// std = 1
		MatrixSumElementsVisitor matrixSumSquaredElementsVisitor = new MatrixSumElementsVisitor(1);
		matrixSumSquaredElementsVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		M.walkInOptimizedOrder(matrixSumSquaredElementsVisitor);
		double scaleValue = matrixSumSquaredElementsVisitor.end();

		scaleValue /= numRows;
		scaleValue = Math.sqrt(scaleValue);

		MatrixVisitor matrixVisitor = new MatrixVisitor(scaleValue);
		matrixVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		M.walkInOptimizedOrder(matrixVisitor);
		matrixVisitor.end();

		// System.out.println("If this factor is not 1 then normalization is wrong! " + res/numRows);
	}

	// calculates sigma squared
	protected double getSigma2(RealMatrix X, RealMatrix Y){
		int M = Y.getRowDimension();
		int N = X.getRowDimension();

		int D = X.getColumnDimension(); 


		double [] sumMX = new double[D];
		double [] sumMY = new double[D];

		for (int d =0; d < D; d++){
			sumMX[d] = 0;
			sumMY[d] = 0;
		}

		ColumnSumVisitor columnSumVisitorMX= new ColumnSumVisitor(sumMX);
		columnSumVisitorMX.start(N, D, 0, N - 1, 0, D - 1);
		X.walkInRowOrder(columnSumVisitorMX);
		columnSumVisitorMX.end();

		ColumnSumVisitor columnSumVisitorMY = new ColumnSumVisitor(sumMY);
		columnSumVisitorMY.start(M, D, 0, M - 1, 0, D - 1);
		Y.walkInRowOrder(columnSumVisitorMY);
		columnSumVisitorMY.end();

		double sum = 0;
		for (int d  = 0; d < D; d++)
			sum += sumMX[d]*sumMY[d];

		double res = (double) (M * X.transpose().multiply(X).getTrace() + N * Y.transpose().multiply(Y).getTrace() - 2*sum) / (M*N*D);

		return res; 

	}

	// calculates G matrix 
	protected void calculateG(RealMatrix Y, RealMatrix G, double beta){
		long M = Y.getRowDimension();

		// TODO: G is symmetric any better way to store it here 
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				double val = (double) Y.getRowVector(i).subtract(Y.getRowVector(j)).getNorm(); // TODO: might be simplified ?!
				val *= -val/(2*beta*beta);
				G.setEntry(i, j, Math.exp(val));
			}	
		}
	}

	// upadate sigma squared value 
	protected double updateSigma2(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix T){
		double res = 0; 

		double[] terms = new double [3];

		int M = Y.getRowDimension();
		int N = X.getRowDimension();
		int D = X.getColumnDimension();


		RealVector bigOneM = new ArrayRealVector(M, 1);
		RealVector bigOneN = new ArrayRealVector(N, 1);


		// TODO: Check this part for errors
		terms[0] = (X.transpose().multiply(new DiagonalMatrix(P.transpose().operate(bigOneM).toArray()))).multiply(X).getTrace();
		terms[1] = ((P.multiply(X)).transpose().multiply(T)).getTrace()*(-2);
		terms[2] = T.transpose().multiply(new DiagonalMatrix(P.operate(bigOneN).toArray()).multiply(T)).getTrace();

		for (int j = 0; j < terms.length; ++j)
			res += terms[j];

		double N_P = 0; 
		MatrixSumElementsVisitor matrixSumElementsVisitor = new MatrixSumElementsVisitor(0);
		matrixSumElementsVisitor.start(M, N, 0, M - 1, 0, N - 1);
		P.walkInOptimizedOrder(matrixSumElementsVisitor);
		N_P = matrixSumElementsVisitor.end();

		res /= (N_P*D); 

		return res;

	}


	protected double calculateP(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix W, RealMatrix G, double w_, double sigmaSq){
		long M = Y.getRowDimension();
		long N = X.getRowDimension();
		long D = X.getColumnDimension();		

		double error = 0;
		for (int n = 0; n < N; ++n) {
			for (int m = 0; m < M; ++m) {
				double val = X.getRowVector(n).subtract(Y.getRowVector(m).add(W.preMultiply(G.getRowVector(m)))).getNorm();
				val *= -val / (2 * sigmaSq);
				val = Math.exp(val);				
				double denom = 0;				
				for (int k = 0; k < M; ++k) {
					double tmp = X.getRowVector(n).subtract(Y.getRowVector(k).add(W.preMultiply(G.getRowVector(k)))).getNorm();
					tmp *= -tmp / (2 * sigmaSq);
					denom += Math.exp(tmp);
				}

				denom += w_ / (1 - w_) * Math.pow(2 * Math.PI * sigmaSq, D / 2.0) * M / N;
				P.setEntry(m, n, val/denom);
			}
		} 
		return error; 
	}

	//-- Move all visitors to the new class file --//

	// DEBUG: Print out the matrix 
	public void printMatrix(RealMatrix M, int ColumnDimension, int RowDimension){
		for (int j = 0; j < RowDimension; j++ ){
			for (int i = 0; i < ColumnDimension; i++ ){
				System.out.print(M.getEntry(j, i) + " ");
			}
			System.out.println();
		}

	}

	// this is a matrix visitor that scales matrix
	public class MatrixVisitor implements RealMatrixChangingVisitor{

		final double scale; 

		public MatrixVisitor(double scale){
			this.scale = scale;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){

		}

		public double visit(int row, int col, double val){
			return val/scale; // val/scaleValue;
		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}

	// this is a matrix visitor that returns the sum over columns
	public class ColumnSumVisitor implements RealMatrixPreservingVisitor{
		final double [] sumOverColumns; 

		public ColumnSumVisitor(double[] sumOverColumns){
			this.sumOverColumns = sumOverColumns;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){

		}

		public void visit(int row, int col, double val){
			sumOverColumns[col] += val;

		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}

	// this is a matrix visitor that subtracts array values from corresponding columns
	public class ColumnSubtractValueVisitor implements RealMatrixChangingVisitor{
		final double [] values; 

		public ColumnSubtractValueVisitor(double[] values){
			this.values = values;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){

		}

		public double visit(int row, int col, double val){
			return val - values[col];
		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}

	// this is a matrix visitor that returns the squared sum over columns
	public class MatrixSumElementsVisitor implements RealMatrixPreservingVisitor{
		double sum; 
		double flag; // specifies if one want to sum squares or values

		public MatrixSumElementsVisitor(int flag){
			this.sum = 0;
			this.flag = flag;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){
		}

		public void visit(int row, int col, double val){
			sum += val* (flag == 0 ? 1 : val);

		}

		public double end(){
			return sum; // 0 for everything is fine 
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

