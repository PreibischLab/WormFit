package coherent.point.drift.refactored;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import coherent.point.drift.refactored.CustomVisitorsFunctions.ColumnDivideVisitor;
import coherent.point.drift.refactored.CustomVisitorsFunctions.ColumnSquareSumVisitor;
import coherent.point.drift.refactored.CustomVisitorsFunctions.ColumnSubtractValueVisitor;
import coherent.point.drift.refactored.CustomVisitorsFunctions.ColumnSumVisitor;
import coherent.point.drift.refactored.CustomVisitorsFunctions.MatrixSumElementsVisitor;

public class ApacheCPD {
	// coherent point drift parameters
	private final double w; 	 // amount of noise [0, 1]
	private final double beta; // Briefly speaking, parameter defines the model of the smoothness regularizer (width of smoothing Gaussian filter).
	private final double lambda; // trade-off between the goodness of maximum likelihood fit and regularization.
	private final int maxIteration; // how many steps of coherent point drift do we want

	final RealMatrix mX; // NxD
	final RealMatrix mY; // MxD

	// TODO: mov these away; it can be local
	// final RealMatrix mW; // MxD
	// final RealMatrix mG; // MxM
	// final RealMatrix mP; // MxN
	// final RealMatrix mT; // MxD transformation matrix

	// this params should never be changed after reading the files or creating the matrices
	private final int numDimensions; // dimensionality of the point set
	private final int numSpotsSample; // # of points in the first point set
	private final int numSpotsLineage; // # of points in the second point set

	// TODO: do I need this one really? 
	// private double sigma2; // 2 for squared

	public ApacheCPD(RealMatrix mX, RealMatrix mY, double w, double beta, double lambda, int maxIteration) {
		this.mX = mX;
		this.mY = mY; 

		this.w = w;
		this.beta = beta;
		this.lambda = lambda; 
		this.maxIteration = maxIteration;

		// these fields are set but should never be passed as params
		this.numDimensions = mX.getColumnDimension(); // dimensionality of the point set
		this.numSpotsSample = mX.getRowDimension(); // # of points in the first point set
		this.numSpotsLineage = mY.getRowDimension(); // # of points in the second point set

		// matrices used for computations
		// mW = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);
		// mG = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsLineage);
		// mP = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsSample);
		// mT = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);
	}

	// TODO: make another constructor to read the files? 

	// FIXED
	/**
	 * Normalization of the matrix so that the mean = 0 and standard deviation = 1
	 * @param mA - matrix to normalize
	 * */
	public static RealMatrix normalize(RealMatrix mA){
		int numColumns = mA.getColumnDimension();
		int numRows = mA.getRowDimension();

		double [] sumOverColumns = new double [numColumns];
		// zero mean => per coordinate!
		for (int d = 0; d < numColumns; ++d)
			sumOverColumns[d]= 0;

		ColumnSumVisitor columnSumVisitor = new CustomVisitorsFunctions().new ColumnSumVisitor(sumOverColumns);
		columnSumVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnSumVisitor);
		columnSumVisitor.end();
		
		for(int d = 0; d < numColumns; ++d)
			sumOverColumns[d] /= numRows;

		ColumnSubtractValueVisitor columnSubtractValueVisitor = new CustomVisitorsFunctions().new ColumnSubtractValueVisitor(sumOverColumns);
		columnSubtractValueVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnSubtractValueVisitor);
		columnSubtractValueVisitor.end();
		
		// std = 1
		// this is the old version of the calculation of the std; I suppose it is wrong 
		// because we want to have the std = 1 over each coordinate not over the whole matrix!
//		MatrixSumElementsVisitor matrixSumSquaredElementsVisitor = new CustomVisitorsFunctions().new MatrixSumElementsVisitor(1);
//		matrixSumSquaredElementsVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
//		mA.walkInOptimizedOrder(matrixSumSquaredElementsVisitor);
//		double scaleValue = matrixSumSquaredElementsVisitor.end();
//		scaleValue /= numRows;
//		scaleValue = Math.sqrt(scaleValue);
//		mA.setSubMatrix(mA.scalarMultiply(1./scaleValue).getData(), 0, 0);
		
		// this seems to be the correct version
		for (int d = 0; d < numColumns; ++d)
			sumOverColumns[d] = 0;
		ColumnSquareSumVisitor columnSquareSumVisitor = new CustomVisitorsFunctions().new ColumnSquareSumVisitor(sumOverColumns);
		columnSquareSumVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnSquareSumVisitor);
		columnSquareSumVisitor.end();
		
		for(int d = 0; d < numColumns; ++d)
			sumOverColumns[d] = Math.sqrt(sumOverColumns[d]/numRows);
		
		ColumnDivideVisitor columnMultiplyVisitor = new CustomVisitorsFunctions().new ColumnDivideVisitor(sumOverColumns);
		columnMultiplyVisitor.start(numRows, numColumns, 0, numRows - 1, 0, numColumns - 1);
		mA.walkInOptimizedOrder(columnMultiplyVisitor);
		columnMultiplyVisitor.end();
		
		return mA;
	}

	// FIXED
	/**
	 * Calculate sigma<sup>2</sup>
	 * @param mA - data points
	 * @param mB - GMM centroids
	 * */
	public static double getSigma2(RealMatrix mA, RealMatrix mB){
		int M = mB.getRowDimension();
		int N = mA.getRowDimension();

		int D = mA.getColumnDimension(); 

		double [] sumMX = new double[D];
		double [] sumMY = new double[D];

		for (int d = 0; d < D; d++){
			sumMX[d] = 0;
			sumMY[d] = 0;
		}

		ColumnSumVisitor columnSumVisitorMX = new CustomVisitorsFunctions().new ColumnSumVisitor(sumMX);
		columnSumVisitorMX.start(N, D, 0, N - 1, 0, D - 1);
		mA.walkInRowOrder(columnSumVisitorMX);
		columnSumVisitorMX.end();

		ColumnSumVisitor columnSumVisitorMY = new CustomVisitorsFunctions().new ColumnSumVisitor(sumMY);
		columnSumVisitorMY.start(M, D, 0, M - 1, 0, D - 1);
		mB.walkInRowOrder(columnSumVisitorMY);
		columnSumVisitorMY.end();

		double sum = 0;
		for (int d  = 0; d < D; d++)
			sum += sumMX[d]*sumMY[d];

		double res = (double) (M * mA.transpose().multiply(mA).getTrace() + N * mB.transpose().multiply(mB).getTrace() - 2*sum) / (M*N*D);

		return res; 

	}

	// FIXED
	/**
	 * Compute the G matrix
	 * */
	protected static void calculateG(RealMatrix Y, RealMatrix G, double beta){
		long M = Y.getRowDimension();
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				double val = Y.getRowVector(i).getDistance(Y.getRowVector(j));
				val *= -val/(2*beta*beta);
				G.setEntry(i, j, Math.exp(val));
			}	
		}
	}

	// FIXED
	/**
	 * Update sigma <sup>2</sup> 
	 * */
	protected static double updateSigma2(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix T){
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
		MatrixSumElementsVisitor matrixSumElementsVisitor = new CustomVisitorsFunctions().new MatrixSumElementsVisitor(0);
		matrixSumElementsVisitor.start(M, N, 0, M - 1, 0, N - 1);
		P.walkInOptimizedOrder(matrixSumElementsVisitor);
		N_P = matrixSumElementsVisitor.end();

		res /= (N_P*D); 

		return res;
	}

	// FIXED
	/**
	 * Calculate P for the non-rigid transformation
	 * */
	protected static double calculatePnonRigid(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix W, RealMatrix G, double w_, double sigmaSq){
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

	// FIXED // FIXME: Seems to be useless
	/**
	 * use this method if you want to read data X, Y from file 
	 * */
	public int runNonRigidRegistation(String from, String to, char delimiter) {
		RealMatrix mA = IOUtils.readPositionsFromCSV(from, delimiter);
		RealMatrix mB = IOUtils.readPositionsFromCSV(to, delimiter);
		normalize(mA);
		normalize(mB);
		this.mX.setSubMatrix(mA.getData(), 0, 0);
		this.mY.setSubMatrix(mB.getData(), 0, 0);
		
		final RealMatrix mW = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);
		final RealMatrix mG = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsLineage);
		final RealMatrix mP = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsSample);
		final RealMatrix mT = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);

		return runNonRigidRegistration(mX, mY, mW, mG, mP, mT, w, beta, lambda, maxIteration);
	}

	/**
	 * use this method if you already have data X, Y
	 * */
	public int runNonRigidRegistration(RealMatrix mA, RealMatrix mB){
		normalize(mA);
		normalize(mB);

		this.mX.setSubMatrix(mA.getData(), 0, 0);
		this.mY.setSubMatrix(mB.getData(), 0, 0);
		
		final RealMatrix mW = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);
		final RealMatrix mG = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsLineage);
		final RealMatrix mP = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsSample);
		final RealMatrix mT = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);

		return runNonRigidRegistration(mX, mY, mW, mG, mP, mT, w, beta, lambda, maxIteration);
	}

	// TODO: run the same thing but with parameters passed 
	protected int runNonRigidRegistration(RealMatrix mX, RealMatrix mY, RealMatrix mW, RealMatrix mG, RealMatrix mP, RealMatrix mT, double w, double beta, double lambda, int maxIteration){
		// important mapping
		// M -> numSpotsLineage
		// N -> numSpotsSample
		// D -> numDimensions

		// this params should never be changed after reading the files or creating the matrices
		final int numDimensions = mX.getColumnDimension(); ; // dimensionality of the point set
		final int numSpotsSample = mX.getRowDimension(); // # of points in the first point set
		final int numSpotsLineage = mY.getRowDimension(); // # of points in the second point set

//		final RealMatrix mW = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);
//		final RealMatrix mG = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsLineage);
//		final RealMatrix mP = MatrixUtils.createRealMatrix(numSpotsLineage, numSpotsSample);
//		final RealMatrix mT = MatrixUtils.createRealMatrix(numSpotsLineage, numDimensions);

		// TODO: do I need this one really? 
		double sigma2; // 2 for squared

		// START HERE

		sigma2 = getSigma2(mX, mY);
		calculateG(mY, mG, beta);
		// TODO: FIXME: this is dumb copying to keep everything final 
		// is mT = ... faster ?
		mT.setSubMatrix(mY.add(mG.multiply(mW)).getData(), 0, 0); 
		// addOverlay(imp, mT, radius);
		double error= 1; 
		double errorOld = 1;

		int iter = 0; 

		// TODO: add error estimator here 
		while (iter++ < maxIteration && sigma2 > 1e-10){
			Utils.printLog(iter, sigma2, Math.abs((error - errorOld)/error));

			errorOld = error;
			error = calculatePnonRigid(mX, mY, mP, mW, mG, w, sigma2);
			error += lambda/2*(mW.transpose().multiply(mG).multiply(mW)).getTrace();

			RealVector bigOneN = new ArrayRealVector(numSpotsSample, 1);

			RealVector invP = mP.operate(bigOneN);
			invP.mapToSelf(new CustomVisitorsFunctions.ElementwiseInverse()); // TODO: is this one correct

			// System.out.println("goes to this part");

			RealMatrix A = mG.add((new DiagonalMatrix(invP.toArray())).scalarMultiply(sigma2 * lambda));
			RealMatrix b = ((new DiagonalMatrix(invP.toArray())).multiply(mP).multiply(mX)).subtract(mY);

			DecompositionSolver solver = new QRDecomposition(A).getSolver();
			mW.setSubMatrix(solver.solve(b).getData(), 0, 0);

			mT.setSubMatrix(mY.add(mG.multiply(mW)).getData(), 0, 0);
			sigma2 = updateSigma2(mX, mY, mP, mT);
		}

		System.out.println("Non-rigid CPD: Done");
		return 0; // TODO: is 0 here for the normal execution 
	}

	// TODO: extend the functionality of this one
	// should be used by affine and rigid registrations
	// AR for affine and rigid 
	protected double calculatePAR(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix B, RealVector t, double w_, double sigmaSq){
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


// FIXME: 
//	public int runAffineRegistration(int flag){		
//		addPoints(img, mX);
//		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
//		imp.show();
//
//		RealMatrix mB = MatrixUtils.createRealIdentityMatrix(numDimensions);
//		// double t = 0; // some other parameter 
//
//		RealVector t = MatrixUtils.createRealVector(new double[numDimensions]);
//		t.set(0);
//
//
//		sigma2 = getSigma2(mX, mY); // TODO: this looks fine but who knows
//
//		double error= 1; 
//		double errorOld = 1;
//
//		int iter = 0;
//		while (iter++ < maxIteration && sigma2 > 1e-10){
//			printLog(iter, sigma2, Math.abs((error - errorOld)/error));
//
//			computePaffine(mX, mY, mP, mB, t, w, sigma2);
//
//			double N_P = 0; 
//			MatrixSumElementsVisitor matrixSumElementsVisitor = new MatrixSumElementsVisitor(0);
//			matrixSumElementsVisitor.start(numSpotsLineage, numSpotsSample, 0, numSpotsLineage - 1, 0, numSpotsSample - 1);
//			mP.walkInOptimizedOrder(matrixSumElementsVisitor);
//			N_P = matrixSumElementsVisitor.end();
//
//			RealVector muX = (mX.transpose().multiply(mP.transpose()).operate(new ArrayRealVector(numSpotsLineage, 1))).mapDivide(N_P);
//			RealVector muY = (mY.transpose().multiply(mP).operate(new ArrayRealVector(numSpotsSample, 1))).mapDivide(N_P);
//
//			// zero mean
//			RealMatrix meanX = mX.subtract(new ArrayRealVector(numSpotsLineage, 1).outerProduct(muX));
//			RealMatrix meanY = mY.subtract(new ArrayRealVector(numSpotsSample, 1).outerProduct(muY));
//
//			// TODO: use set instead
//			// check if this freaking line is correct at all 
//			mB = (meanX.transpose().multiply(mP.transpose().multiply(meanY))).multiply(MatrixUtils.inverse(meanY.transpose().multiply((new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(numSpotsLineage, 1)).toArray()).multiply(meanY)))));
//
//			t.setSubVector(0, muX.subtract(mB.operate(muY)));
//
//			// twoo terms for the sum
//			double [] terms = new double [2];
//
//			terms[0] = meanX.transpose().multiply(new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(numSpotsSample, 1)).toArray()).multiply(meanX)).getTrace();
//			terms[1] = meanX.transpose().multiply(mP.transpose().multiply(meanY).multiply(mB.transpose())).getTrace();
//
//			sigma2 = (terms[0] - terms[1])/(N_P*numDimensions);
//
//			// TODO: FIXME: line below is only to see the result
//			mT.setSubMatrix((mY.multiply(mB.transpose())).add(new ArrayRealVector(numSpotsLineage, 1).outerProduct(t)).getData(), 0, 0);
//			addOverlay(imp, mT, radius);
//
//		}
//
//		mT.setSubMatrix((mY.multiply(mB.transpose())).add(new ArrayRealVector(numSpotsLineage, 1).outerProduct(t)).getData(), 0, 0);
//
//		return 0;
//	}

	protected void computePrigid(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix R, double s, RealVector t, double w_, double sigmaSq){
		calculatePAR(X, Y, P, R.scalarMultiply(s), t, w_, sigmaSq);
	}

	protected void computePaffine(RealMatrix X, RealMatrix Y, RealMatrix P, RealMatrix B, RealVector t, double w_, double sigmaSq){
		calculatePAR(X, Y, P, B, t, w_, sigmaSq);
	}

// FIXME:
//	public int runRigidRegistration(int flag){
//		addPoints(img, mX);
//		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
//		imp.show();
//
//		RealMatrix mR = MatrixUtils.createRealIdentityMatrix(numDimensions);
//
//		RealVector t = MatrixUtils.createRealVector(new double[numDimensions]);
//		t.set(0);
//		double s = 1;
//
//		sigma2 = getSigma2(mX, mY); // TODO: this looks fine but who knows
//
//		double error= 1; 
//		double errorOld = 1;
//
//		int iter = 0;
//		while (iter++ < maxIteration && sigma2 > 1e-10){
//			printLog(iter, sigma2, Math.abs((error - errorOld)/error));
//
//			computePrigid(mX, mY, mP, mR, s, t, w, sigma2);
//
//			double N_P = 0; 
//			MatrixSumElementsVisitor matrixSumElementsVisitor = new MatrixSumElementsVisitor(0);
//			matrixSumElementsVisitor.start(numSpotsLineage, numSpotsSample, 0, numSpotsLineage - 1, 0, numSpotsSample - 1);
//			mP.walkInOptimizedOrder(matrixSumElementsVisitor);
//			N_P = matrixSumElementsVisitor.end();
//
//			RealVector muX = (mX.transpose().multiply(mP.transpose()).operate(new ArrayRealVector(numSpotsLineage, 1))).mapDivide(N_P);
//			RealVector muY = (mY.transpose().multiply(mP).operate(new ArrayRealVector(numSpotsSample, 1))).mapDivide(N_P);
//
//			// zero mean
//			RealMatrix meanX = mX.subtract(new ArrayRealVector(numSpotsLineage, 1).outerProduct(muX));
//			RealMatrix meanY = mY.subtract(new ArrayRealVector(numSpotsSample, 1).outerProduct(muY));
//
//
//			RealMatrix mA = meanX.transpose().multiply(mP.transpose().multiply(meanY));
//			SingularValueDecomposition SVDmA = new SingularValueDecomposition(mA);
//
//			double[] identity = new double[numDimensions];
//			for(int d = 0; d < numDimensions - 1; d++)
//				identity[d] = 1;
//			identity[numDimensions - 1] = new LUDecomposition(SVDmA.getU().multiply(SVDmA.getV())).getDeterminant();
//
//			mR.setSubMatrix(SVDmA.getU().multiply(new DiagonalMatrix(identity, false).multiply(SVDmA.getVT())).getData() , 0, 0);
//
//			s = (mA.transpose().multiply(mR).getTrace()) / (meanY.transpose().multiply(new DiagonalMatrix( mP.transpose().operate(new ArrayRealVector(numSpotsLineage, 1)).toArray()).multiply(meanY)).getTrace());
//			t.setSubVector(0, mR.operate(muY).mapMultiply(s));
//
//			// two terms for the sum
//			double [] terms = new double [2];
//
//			terms[0] = meanX.transpose().multiply(new DiagonalMatrix(mP.transpose().operate(new ArrayRealVector(numSpotsSample, 1)).toArray()).multiply(meanX)).getTrace();
//			terms[1] = s*(mA.transpose().multiply(mR)).getTrace(); 
//
//			sigma2 = (terms[0] - terms[1])/(N_P*numDimensions);
//
//			// TODO: FIXME: line below is only to see the result
//			mT.setSubMatrix(mY.multiply(mR.transpose()).scalarMultiply(s).add(new ArrayRealVector(numSpotsLineage, 1).outerProduct(t) ).getData(), 0, 0);
//			addOverlay(imp, mT, radius);
//
//		}
//
//		mT.setSubMatrix(mY.multiply(mR.transpose()).scalarMultiply(s).add(new ArrayRealVector(numSpotsLineage, 1).outerProduct(t) ).getData(), 0, 0);
//		return 0;
//	}	




// FIXME:
	/**
	 * use this method if you want to read data X, Y from file 
	 * */
//	public int runAffineRegistration(int flag, String from, String to){
//		readData(mX, mY, from, to);
//		return runAffineRegistration(flag);
//	}

// FIXME:
	/**
	 * use this method if you already have data X, Y
	 * */
//	public int runAffineRegistration(int flag, RealMatrix mX, RealMatrix mY){
//		this.mX.setSubMatrix(mX.getData(), 0, 0);
//		this.mY.setSubMatrix(mY.getData(), 0, 0);
//		return runAffineRegistration(flag);
//	}

// FIXME:
	/**
	 * use this method if you want to read data X, Y from file 
	 * */
//	public int runRigidRegistration(int flag, String from, String to){
//		readData(mX, mY, from, to);
//		return runRigidRegistration(flag);
//	}

// FIXME:
	/**
	 * use this method if you already have data X, Y
	 * */
//	public int runRigidRegistration(int flag, RealMatrix mX, RealMatrix mY){
//		this.mX.setSubMatrix(mX.getData(), 0, 0);
//		this.mY.setSubMatrix(mY.getData(), 0, 0);
//		return runRigidRegistration(flag);
//	}


	public static void main(String[] args) {
		// TODO: you have to scale this data otherwise
		// you won't see the convergence process
		// new ApacheCPD().testNonRigid();
		new ApacheCPDTests().testNonRigid();
	}
}

