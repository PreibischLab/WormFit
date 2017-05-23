package coherent.point.drift;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.la4j.Matrix;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import mpicbg.imglib.util.Util;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.opencsv.CSVReader;
import util.opencsv.CSVWriter;

public class ApacheCPD {

	Img<FloatType> img; 
	ImagePlus imp; 

	// TODO: read this parameters from file or from variable
	final int D; // dimensionality of the point set
	final int N; // # of points in the first point set
	final int M; // # of points in the second point set

	double sigma2; // 2 for squared

	// parameters
	final double w; 	 // amount of noise [0, 1]
	final double beta;   // Briefly speaking, parameter defines the model of the smoothness regularizer (width of smoothing Gaussian filter).
	final double lambda; // trade-off between the goodness of maximum likelihood fit and regularization.

	final RealMatrix mX; // NxD
	final RealMatrix mY; // MxD
	final RealMatrix mW; // MxD
	final RealMatrix mG; // MxM
	final RealMatrix mP; // MxN
	final RealMatrix mT; // MxD transformation matrix

	// these variables are necessary for proper visualization
	double[] translate;
	double[] scale;
	double[] sigma;

	final int maxIteration;


	// default constructor   
	public ApacheCPD(){
		D = 1;
		M = 1;
		N = 1; 
		
		mX = MatrixUtils.createRealMatrix(N, D); 
		mY = MatrixUtils.createRealMatrix(M, D); 
		mW = MatrixUtils.createRealMatrix(M, D); 
		mG = MatrixUtils.createRealMatrix(M, M); 
		mP = MatrixUtils.createRealMatrix(M, N); 
		mT = MatrixUtils.createRealMatrix(M, D); 

		w = 1;
		beta = 1;
		lambda = 1;

		maxIteration = 1;
	}

	/** 
	 * Constructor 
	 * @param img - displays the result of fitting
	 * @param N - total number of data points
	 * @param M - total number of GMM centroids
	 * @param D - dimensionality
	 * @param w - amount of noise [0, 1]
	 * @param beta - defines the model of the smoothness regularizer (width of smoothing Gaussian filter)
	 * @param lambda - trade-off between the goodness of maximum likelihood fit and regularization
	 * @param maxIteration - number of maximum iterations for the algorithm  
	 * */
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
	
	// 
	/**
	 * Constructor  
	 * @param img - displays the result of fitting
	 * @param mX - data points
	 * @param mY - GMM centroids
	 * @param w - amount of noise [0, 1]
	 * @param beta - defines the model of the smoothness regularizer (width of smoothing Gaussian filter)
	 * @param lambda - trade-off between the goodness of maximum likelihood fit and regularization
	 * @param maxIteration - number of maximum iterations for the algorithm  
	 * */
 	public ApacheCPD(Img<FloatType> img, RealMatrix mX, RealMatrix mY, double w, double beta, double lambda, int maxIteration) {
		this.img = img;

		this.w = w;
		this.beta = beta;
		this.lambda = lambda;

		this.maxIteration = maxIteration;
		
		this.D = mX.getColumnDimension(); // dimensionality of the point set
		this.N = mX.getRowDimension(); // # of points in the first point set
		this.M = mY.getRowDimension(); // # of points in the second point set

		this.mX = mX;
		this.mY = mY;		
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
	
	
	
 	/**
 	 * Normalization of the matrix so that the mean = 0 and standard deviation = 1
 	 * @param mA - matrix to normalize
 	 * */
 	protected void normalize(RealMatrix mA){
		int numColumns = mA.getColumnDimension();
		int numRows = mA.getRowDimension();

		double [] sumOverColumns = new double [numColumns];
		// zero mean
		for (int d = 0; d < numColumns; ++d)
			sumOverColumns[d]= 0;

		ColumnSumVisitor columnSumVisitor = new ColumnSumVisitor(sumOverColumns);
		columnSumVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnSumVisitor);
		columnSumVisitor.end();

		for(int d = 0; d < numColumns; ++d)		
			sumOverColumns[d] /= numRows;

		ColumnSubtractValueVisitor columnSubtractValueVisitor = new ColumnSubtractValueVisitor(sumOverColumns);
		columnSubtractValueVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnSubtractValueVisitor);
		columnSubtractValueVisitor.end();

		// std = 1
		MatrixSumElementsVisitor matrixSumSquaredElementsVisitor = new MatrixSumElementsVisitor(1);
		matrixSumSquaredElementsVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(matrixSumSquaredElementsVisitor);
		double scaleValue = matrixSumSquaredElementsVisitor.end();

		scaleValue /= numRows;
		scaleValue = Math.sqrt(scaleValue);

		mA.setSubMatrix(mA.scalarMultiply(1./scaleValue).getData(), 0, 0);
	}

	/**
	 * Calculate sigma<sup>2</sup>
	 * @param mA - data points
	 * @param mB - GMM centroids
	 * @return sigma squared
	 * */
	protected double getSigma2(RealMatrix mA, RealMatrix mB){
		int M = mB.getRowDimension();
		int N = mA.getRowDimension();

		int D = mA.getColumnDimension(); 

		double [] sumMX = new double[D];
		double [] sumMY = new double[D];

		for (int d =0; d < D; d++){
			sumMX[d] = 0;
			sumMY[d] = 0;
		}

		ColumnSumVisitor columnSumVisitorMX = new ColumnSumVisitor(sumMX);
		columnSumVisitorMX.start(N, D, 0, N - 1, 0, D - 1);
		mA.walkInRowOrder(columnSumVisitorMX);
		columnSumVisitorMX.end();

		ColumnSumVisitor columnSumVisitorMY = new ColumnSumVisitor(sumMY);
		columnSumVisitorMY.start(M, D, 0, M - 1, 0, D - 1);
		mB.walkInRowOrder(columnSumVisitorMY);
		columnSumVisitorMY.end();

		double sum = 0;
		for (int d  = 0; d < D; d++)
			sum += sumMX[d]*sumMY[d];

		double res = (double) (M * mA.transpose().multiply(mA).getTrace() + N * mB.transpose().multiply(mB).getTrace() - 2*sum) / (M*N*D);

		return res; 

	}

	/*
	 * Compute the G matrix
	 * */
	protected void calculateG(RealMatrix Y, RealMatrix G, double beta){
		long M = Y.getRowDimension();
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				double val = Y.getRowVector(i).getDistance(Y.getRowVector(j));
				val *= -val/(2*beta*beta);
				G.setEntry(i, j, Math.exp(val));
			}	
		}
	}

	/*
	 * Update sigma <sup>2</sup> 
	 * */
	protected double updateSigma2(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix T){
		double res = 0; 

		double[] terms = new double [3];

		int M = Y.getRowDimension();
		int N = X.getRowDimension();
		int D = X.getColumnDimension();
		
		RealVector bigOneM = new ArrayRealVector(M, 1);
		RealVector bigOneN = new ArrayRealVector(N, 1);

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

	/*
	 * Calculate P for the non-rigid transformation
	 * */
	protected double calculatePnonRigid(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix W, RealMatrix G, double w_, double sigmaSq){
		long M = Y.getRowDimension();
		long N = X.getRowDimension();
		long D = X.getColumnDimension();		

		double error = 0;
		for (int n = 0; n < N; ++n) {
			for (int m = 0; m < M; ++m) {
				double val = X.getRowVector(n).getDistance(Y.getRowVector(m).add(W.preMultiply(G.getRowVector(m))));
				val *= -val / (2 * sigmaSq);
				val = Math.exp(val);				
				double denom = 0;				
				for (int k = 0; k < M; ++k) {
					double tmp = X.getRowVector(n).getDistance(Y.getRowVector(k).add(W.preMultiply(G.getRowVector(m))));
					tmp *= -tmp / (2 * sigmaSq);
					denom += Math.exp(tmp);
				}

				denom += w_ / (1 - w_) * Math.pow(2 * Math.PI * sigmaSq, D / 2.0) * M / N;
				P.setEntry(m, n, val/denom);
			}
		} 
		return error; 
	}


	public void printLog(int idx, double sigma, double error){
		System.out.println("ITERATION #" + idx + ":");
		System.out.println("Sigma squared: " + String.format(java.util.Locale.US, "%.2e", sigma)); // TODO: was sigma2 inititally
		System.out.println("Error: " + String.format(java.util.Locale.US, "%.2e", error));
	}

	/*
	 * use this method if you want to read data X, Y from file 
	 * */
	public int runNonRigidRegistration(int flag, String from, String to){
		readData(mX, mY, from, to);
		return runNonRigidRegistration(flag);
	}
	
	/*
	 * use this method if you already have data X, Y
	 * */
	public int runNonRigidRegistration(int flag, RealMatrix mX, RealMatrix mY){
		this.mX.setSubMatrix(mX.getData(), 0, 0);
		this.mY.setSubMatrix(mY.getData(), 0, 0);
		return runNonRigidRegistration(flag);
	}
		
	protected int runNonRigidRegistration(int flag){	
		addPoints(img);
		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
		imp.show();

		sigma2 = getSigma2(mX, mY);
		calculateG(mY, mG, beta);
		// TODO: FIXME: this is dumb copying to keep everything final 
		// is mT = ... faster ?
		mT.setSubMatrix(mY.add(mG.multiply(mW)).getData(), 0, 0); 
		addOverlay(imp, mT);
		double error= 1; 
		double errorOld = 1;

		int iter = 0; 

		// TODO: add error estimator here 
		while (iter++ < maxIteration && sigma2 > 1e-10){
			printLog(iter, sigma2, Math.abs((error - errorOld)/error));

			errorOld = error;
			error = calculatePnonRigid(mX, mY, mP, mW, mG, w, sigma2);
			error += lambda/2*(mW.transpose().multiply(mG).multiply(mW)).getTrace();

			RealVector bigOneN = new ArrayRealVector(N, 1);

			RealVector invP = mP.operate(bigOneN);
			invP.mapToSelf(new ElementwiseInverse());

			System.out.println("goes to this part");

			RealMatrix A = mG.add((new DiagonalMatrix(invP.toArray())).scalarMultiply(sigma2 * lambda));
			RealMatrix b = ((new DiagonalMatrix(invP.toArray())).multiply(mP).multiply(mX)).subtract(mY);

			DecompositionSolver solver = new QRDecomposition(A).getSolver();
			mW.setSubMatrix(solver.solve(b).getData(), 0, 0) ;

			mT.setSubMatrix(mY.add(mG.multiply(mW)).getData(), 0, 0); 
			sigma2 = updateSigma2(mX, mY, mP, mT);

			addOverlay(imp, mT);
		}

		System.out.println("Done");
		return 0; // TODO: is 0 here for the normal execution 
	}

	// TODO: extend the functionality of this one
	// should be used by affine and rigid registrations
	// AR for affine and rigid 
	public double calculatePAR(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix B, RealVector t, double w_, double sigmaSq){
		long M = Y.getRowDimension();
		long N = X.getRowDimension();
		long D = X.getColumnDimension();		

		double error = 0;
		for (int n = 0; n < N; ++n) {
			for (int m = 0; m < M; ++m) {
				double val = (X.getRowVector(n).subtract(B.operate(Y.getRowVector(m)).add(t))).getNorm();
				val *= -val / (2 * sigmaSq);
				val = Math.exp(val);				
				double denom = 0;				
				for (int k = 0; k < M; ++k) {
					double tmp = (X.getRowVector(n).subtract(B.operate(Y.getRowVector(k)).add(t))).getNorm();
					tmp *= -tmp / (2 * sigmaSq);
					denom += Math.exp(tmp);
				}

				denom += w_ / (1 - w_) * Math.pow(2 * Math.PI * sigmaSq, D / 2.0) * M / N;
				P.setEntry(m, n, val/denom);
			}
		} 
		return error; 
	}

	/*
	 * use this method if you want to read data X, Y from file 
	 * */
	public int runAffineRegistration(int flag, String from, String to){
		readData(mX, mY, from, to);
		return runAffineRegistration(flag);
	}
	
	/*
	 * use this method if you already have data X, Y
	 * */
	public int runAffineRegistration(int flag, RealMatrix mX, RealMatrix mY){
		this.mX.setSubMatrix(mX.getData(), 0, 0);
		this.mY.setSubMatrix(mY.getData(), 0, 0);
		return runAffineRegistration(flag);
	}

	public int runAffineRegistration(int flag){		
		addPoints(img);
		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
		imp.show();
		
		RealMatrix mB = MatrixUtils.createRealIdentityMatrix(D);
		// double t = 0; // some other parameter 

		RealVector t = MatrixUtils.createRealVector(new double[D]);
		t.set(0);


		sigma2 = getSigma2(mX, mY); // TODO: this looks fine but who knows

		double error= 1; 
		double errorOld = 1;

		int iter = 0;
		while (iter++ < maxIteration && sigma2 > 1e-10){
			printLog(iter, sigma2, Math.abs((error - errorOld)/error));

			computePaffine(mX, mY, mP, mB, t, w, sigma2);

			double N_P = 0; 
			MatrixSumElementsVisitor matrixSumElementsVisitor = new MatrixSumElementsVisitor(0);
			matrixSumElementsVisitor.start(M, N, 0, M - 1, 0, N - 1);
			mP.walkInOptimizedOrder(matrixSumElementsVisitor);
			N_P = matrixSumElementsVisitor.end();

			RealVector muX = (mX.transpose().multiply(mP.transpose()).operate(new ArrayRealVector(M, 1))).mapDivide(N_P);
			RealVector muY = (mY.transpose().multiply(mP).operate(new ArrayRealVector(N, 1))).mapDivide(N_P);

			// zero mean
			RealMatrix meanX = mX.subtract(new ArrayRealVector(M, 1).outerProduct(muX));
			RealMatrix meanY = mY.subtract(new ArrayRealVector(N, 1).outerProduct(muY));

			// TODO: use set instead
			// check if this freaking line is correct at all 
			mB = (meanX.transpose().multiply(mP.transpose().multiply(meanY))).multiply(MatrixUtils.inverse(meanY.transpose().multiply((new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(M, 1)).toArray()).multiply(meanY)))));

			t.setSubVector(0, muX.subtract(mB.operate(muY)));

			// twoo terms for the sum
			double [] terms = new double [2];

			terms[0] = meanX.transpose().multiply(new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(N, 1)).toArray()).multiply(meanX)).getTrace();
			terms[1] = meanX.transpose().multiply(mP.transpose().multiply(meanY).multiply(mB.transpose())).getTrace();

			sigma2 = (terms[0] - terms[1])/(N_P*D);
			
			// TODO: FIXME: line below is only to see the result
			mT.setSubMatrix((mY.multiply(mB.transpose())).add(new ArrayRealVector(M, 1).outerProduct(t)).getData(), 0, 0);
			addOverlay(imp, mT);

		}

		mT.setSubMatrix((mY.multiply(mB.transpose())).add(new ArrayRealVector(M, 1).outerProduct(t)).getData(), 0, 0);

		return 0;
	}

	public void computePrigid(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix R, double s, RealVector t, double w_, double sigmaSq){
		calculatePAR(X, Y, P, R.scalarMultiply(s), t, w_, sigmaSq);
	}
	
	public void computePaffine(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix B, RealVector t, double w_, double sigmaSq){
		calculatePAR(X, Y, P, B, t, w_, sigmaSq);
	}
	
	/*
	 * use this method if you want to read data X, Y from file 
	 * */
	public int runRigidRegistration(int flag, String from, String to){
		readData(mX, mY, from, to);
		return runRigidRegistration(flag);
	}
	
	/*
	 * use this method if you already have data X, Y
	 * */
	public int runRigidRegistration(int flag, RealMatrix mX, RealMatrix mY){
		this.mX.setSubMatrix(mX.getData(), 0, 0);
		this.mY.setSubMatrix(mY.getData(), 0, 0);
		return runRigidRegistration(flag);
	}
	
	public int runRigidRegistration(int flag){		
		addPoints(img);
		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
		imp.show();
		
		RealMatrix mR = MatrixUtils.createRealIdentityMatrix(D);
	
		RealVector t = MatrixUtils.createRealVector(new double[D]);
		t.set(0);
		double s = 1;
		
		sigma2 = getSigma2(mX, mY); // TODO: this looks fine but who knows

		double error= 1; 
		double errorOld = 1;

		int iter = 0;
		while (iter++ < maxIteration && sigma2 > 1e-10){
			printLog(iter, sigma2, Math.abs((error - errorOld)/error));

			computePrigid(mX, mY, mP, mR, s, t, w, sigma2);

			double N_P = 0; 
			MatrixSumElementsVisitor matrixSumElementsVisitor = new MatrixSumElementsVisitor(0);
			matrixSumElementsVisitor.start(M, N, 0, M - 1, 0, N - 1);
			mP.walkInOptimizedOrder(matrixSumElementsVisitor);
			N_P = matrixSumElementsVisitor.end();

			RealVector muX = (mX.transpose().multiply(mP.transpose()).operate(new ArrayRealVector(M, 1))).mapDivide(N_P);
			RealVector muY = (mY.transpose().multiply(mP).operate(new ArrayRealVector(N, 1))).mapDivide(N_P);

			// zero mean
			RealMatrix meanX = mX.subtract(new ArrayRealVector(M, 1).outerProduct(muX));
			RealMatrix meanY = mY.subtract(new ArrayRealVector(N, 1).outerProduct(muY));

			
			RealMatrix mA = meanX.transpose().multiply(mP.transpose().multiply(meanY));
			SingularValueDecomposition SVDmA = new SingularValueDecomposition(mA);
			
			double[] identity = new double[D];
			for(int d = 0; d < D - 1; d++)
				identity[d] = 1;
			identity[D - 1] = new LUDecomposition(SVDmA.getU().multiply(SVDmA.getV())).getDeterminant();
					
			mR.setSubMatrix(SVDmA.getU().multiply(new DiagonalMatrix(identity, false).multiply(SVDmA.getVT())).getData() , 0, 0);
			
			s = (mA.transpose().multiply(mR).getTrace()) / (meanY.transpose().multiply(new DiagonalMatrix( mP.transpose().operate(new ArrayRealVector(M, 1)).toArray()).multiply(meanY)).getTrace());
			t.setSubVector(0, mR.operate(muY).mapMultiply(s));

			// two terms for the sum
			double [] terms = new double [2];

			terms[0] = meanX.transpose().multiply(new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(N, 1)).toArray()).multiply(meanX)).getTrace();
			terms[1] = s*(mA.transpose().multiply(mR)).getTrace(); 
					
			sigma2 = (terms[0] - terms[1])/(N_P*D);
			
			// TODO: FIXME: line below is only to see the result
			mT.setSubMatrix(mY.multiply(mR.transpose()).scalarMultiply(s).add(new ArrayRealVector(M, 1).outerProduct(t) ).getData(), 0, 0);
			addOverlay(imp, mT);

		}

		mT.setSubMatrix(mY.multiply(mR.transpose()).scalarMultiply(s).add(new ArrayRealVector(M, 1).outerProduct(t) ).getData(), 0, 0);
		return 0;
	}
	
	//-- reading part move to another class --// 
	public void readData(RealMatrix X, RealMatrix Y, String from, String to) {
		readCSV(from, to); // reading is fine
		normalize(X);
		normalize(Y);
	}

	// TODO: move intput and visual stuff to another class file
	public void readCSV(String from, String to) {
		CSVReader reader = null;
		String[] nextLine;

		try {
			reader = new CSVReader(new FileReader(from), '\t');
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		try {
			int i = 0;
			while ((nextLine = reader.readNext()) != null) {
				mX.setEntry(i, 0, Double.parseDouble(nextLine[0]));
				mX.setEntry(i, 1, Double.parseDouble(nextLine[1]));
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			reader = new CSVReader(new FileReader(to), '\t');
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			int i = 0;
			while ((nextLine = reader.readNext()) != null) {
				mY.setEntry(i, 0, Double.parseDouble(nextLine[0]));
				mY.setEntry(i, 1, Double.parseDouble(nextLine[1]));
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public void writeCSV() {
		CSVWriter writer = null;
		String[] nextLine;

		try {
			writer = new CSVWriter(new FileWriter("src/main/resources/xRes.csv"), '\t');
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {
			for (int i = 0; i < N; ++i) {
				nextLine = mX.getRow(i).toString().split(" ");
				writer.writeNext(nextLine);
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			writer = new CSVWriter(new FileWriter("src/main/resources/yRes.csv"), '\t');
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {
			for (int i = 0; i < N; ++i) {
				nextLine = mY.getRow(i).toString().split(" ");
				writer.writeNext(nextLine);
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/*
	 * Shows coordinates defined by A in the imp image
	 * 
	 */
	public void addOverlay(ImagePlus imp, RealMatrix A) {
		int numDimensions = A.getColumnDimension();
		int numPoints = A.getRowDimension();

		Overlay overlay = imp.getOverlay();
		if (overlay == null) {
			overlay = new Overlay();
		}
		overlay.clear();

		double[] location = new double[numDimensions];
		for (int i = 0; i < numPoints; ++i) {
			for (int d = 0; d < numDimensions; ++d) {
				location[d] = A.getEntry(i, d) * scale[d] + translate[d];
			}
			final OvalRoi or = new OvalRoi(location[0] - sigma[0], location[1] - sigma[0], Util.round(2 * sigma[0]),
					Util.round(2 * sigma[1]));
			or.setStrokeColor(Color.RED);
			overlay.add(or);
		}
		imp.setOverlay(overlay);
	}

	// add points that has to be detected
	public void addPoints(Img<FloatType> img) {
		double[] location = new double[D];
		for (int i = 0; i < N; ++i) {
			for (int d = 0; d < D; ++d) {
				location[d] = (mX.getEntry(i, d) * scale[d] + translate[d]);
			}
			// TODO: Use function from klim.utils
			addGaussian(img, location, sigma);
		}
	}
	
	//----------------------------------------//

	private static class ElementwiseInverse implements UnivariateFunction {
		public double value(double x) {
			return 1.0/x;
		}
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

	public void testNonRigid(){
		new ImageJ();

		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;

		final Img<FloatType> img = new ArrayImgFactory<FloatType>().create(new long[] { 500, 500 }, new FloatType());

		// read data first
		// TODO: reading here
		// TODO:  make the calculation of these parameters automatic
		int D = 2; // dimensionality of the point set
		int M = 112; // # of points in the first point set
		int N = 133; // # of points in the second point set

		String path = "/home/milkyklim/Documents/imglib2Dev/WormFit/src/main/resources/";
		String from = path + "worm-folded.csv";
		String to = path + "worm-straight.csv";
		
		
		// pass it as arguments
		new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration).runNonRigidRegistration(0, from, to) ;
		// new CoherentPointDrift().readCSV();
	}

	public void testAffine(){
		new ImageJ();

		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;

		final Img<FloatType> img = new ArrayImgFactory<FloatType>().create(new long[] { 500, 500 }, new FloatType());

		// read data first
		// TODO: reading here
		int D = 2; // dimensionality of the point set
		int M = 91; // # of points in the first point set
		int N = 91; // # of points in the second point set
		
		String path = "/home/milkyklim/Documents/imglib2Dev/WormFit/src/main/resources/";
		String from = path + "fishy-fish-x.csv";
		String to = path + "fishy-fish-y.csv";
		

		// pass it as arguments
		new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration).runAffineRegistration(0, from, to);
		// new CoherentPointDrift().readCSV();
	}
	
	public void testRigid(){
		new ImageJ();

		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;

		final Img<FloatType> img = new ArrayImgFactory<FloatType>().create(new long[] { 500, 500 }, new FloatType());

		// read data first
		// TODO: reading here
		int D = 2; // dimensionality of the point set
		int M = 91; // # of points in the first point set
		int N = 91; // # of points in the second point set
		
		String path = "/home/milkyklim/Documents/imglib2Dev/WormFit/src/main/resources/";
		String from = path + "fishy-fish-rigid-x.csv";
		String to = path + "fishy-fish-rigid-y.csv";
		

		// pass it as arguments
		new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration).runRigidRegistration(0, from, to);
		// new CoherentPointDrift().readCSV();
	}
	
	
	// TODO: move this one to some other class
	// taken frome radial symmetry for test purposes only
	final public static void addGaussian( final Img< FloatType > image, final double[] location, final double[] sigma )
	{
		final int numDimensions = image.numDimensions();
		final int[] size = new int[ numDimensions ];
		
		final long[] min = new long[ numDimensions ];
		final long[] max = new long[ numDimensions ];
		
		final double[] two_sq_sigma = new double[ numDimensions ];
		
		for ( int d = 0; d < numDimensions; ++d )
		{
			size[ d ] = Util.getSuggestedKernelDiameter( sigma[ d ] ) * 2;
			min[ d ] = (int)Math.round( location[ d ] ) - size[ d ]/2;
			max[ d ] = min[ d ] + size[ d ] - 1;
			two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];
		}

		final RandomAccessible< FloatType > infinite = Views.extendZero( image );
		final RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
		final IterableInterval< FloatType > iterable = Views.iterable( interval );
		final Cursor< FloatType > cursor = iterable.localizingCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			
			double value = 1;
			
			for ( int d = 0; d < numDimensions; ++d )
			{
				final double x = location[ d ] - cursor.getIntPosition( d );
				value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
			}
			
			cursor.get().set( cursor.get().get() + (float)value );
		}
	}
	

	public static void main(String[] args) {
		// new ApacheCPD().test();
		new ApacheCPD().testNonRigid();
	}

}

