package coherent.point.drift;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.la4j.Matrices;
import org.la4j.Matrix;
import org.la4j.Vector;
import org.la4j.Vectors;
import org.la4j.linear.GaussianSolver;
import org.la4j.linear.LeastSquaresSolver;
import org.la4j.matrix.dense.Basic2DMatrix;
import org.la4j.operation.MatrixVectorOperation;
import org.la4j.vector.VectorFactory;
import org.la4j.vector.dense.BasicVector;
import org.la4j.vector.functor.VectorAccumulator;
import org.la4j.vector.functor.VectorFunction;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import imglib.ops.operator.binary.Multiply;
import mpicbg.imglib.util.Util;
import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import test.TestGauss3d;
import util.ImgLib2Util;
import util.opencsv.CSVReader;
import util.opencsv.CSVWriter;

public class CoherentPointDrift {

	Img<FloatType> img;

	int D = 2; // dimensionality of the point set
	int N = 133; // # of points in the first point set
	int M = 112; // # of points in the second point set
	double[] X; // first point set N x D
	double[] Y; // second point set M x D
	double[] I; // (in theory) identity matrix

	double[] W; // matrix of coefficients M x D
	double sigma2; // TODO: 2 for squared

	double[] G;
	double[] P; // probability matrix

	// parameters
	double w; // amount of noise [0, 1]
	double beta; // Briefly speaking, parameter defines the model of the
	// smoothness regularizer (width of smoothing Gaussian
	// filter in (20)).
	double lambda; // Parameter represents the trade-off between the goodness of
	// maximum likelihood fit and regularization.

	// TODO: convert stuff to the HPC
	Basic2DMatrix mX; // N D
	Basic2DMatrix mY; // M D
	Basic2DMatrix mI; //
	Basic2DMatrix mW; // M D
	Basic2DMatrix mG; // M M
	Basic2DMatrix mP; // M N

	Basic2DMatrix mT; // M D transformation matrix

	// straight to the non-rigid registration

	// constructor
	public CoherentPointDrift() {
	}

	/**
	 * Inverts each element of vector (1/value).
	 */
	public static final VectorFunction INV2_FUNCTION = new VectorFunction() {
		@Override
		public double evaluate(int i, double value) {
			double eps = 0.000000001;
			double res = 0;
			if (value < eps) {
				System.out.println(i + "'th value is too small to invert!");
			} else {
				res = 1 / value;
			}
			return res;
		}
	};

	public static void debugResult(Matrix m, int idx) {
		for (int j = 0; j < 20; ++j) {
			double val = m.get(j, idx);
			System.out.println(j + " : " + val);
		}
	}

	public static void debugResult(Vector m) {
		for (int j = 0; j < 20; ++j) {
			double val = m.get(j);
			System.out.println(j + " : " + val);
		}
	}

	// these guys are necessary for proper coordinates adjustment
	double[] translate = new double[] { 250, 250 };
	double[] scale = new double[] { 150, 150 };
	double[] sigma = new double[] { 3, 3 };

	public void addPoints(Img<FloatType> img) {
		double[] location = new double[D];
		// double [] sigma = new double [D];

		// for (int d = 0; d <D; ++d ){
		// sigma[d] = 3;
		// }

		for (int i = 0; i < N; ++i) {
			for (int d = 0; d < D; ++d) {
				location[d] = (mX.get(i, d) * scale[d] + translate[d]);
				// System.out.print(location[d]);
			}
			// System.out.println();
			TestGauss3d.addGaussian(img, location, sigma);
		}
	}

	// display output
	
	public void addOverlay(ImagePlus imp) {
		// new ImageJ();
		Overlay overlay = imp.getOverlay();
		if (overlay == null) {
			System.out.println("addOverlay: overlay is null!");
			overlay = new Overlay();
			imp.setOverlay(overlay);
		}
		overlay.clear();

		// here should be the T coordinate of the fitted points
		for (int i = 0; i < M; ++i) {
			double[] xy = new double[D];
			for (int d = 0; d < D; ++d) {
				xy[d] = mT.get(i, d)*scale[d] + translate[d];
				// System.out.print(xy[d] + " ");
			}
			// System.out.println();
			final OvalRoi or = new OvalRoi(xy[0] - sigma[0], xy[1] - sigma[0], Util.round(2 * sigma[0]),
					Util.round(2 * sigma[1]));
			or.setStrokeColor(Color.RED);
			overlay.add(or);
		}
		imp.setOverlay(overlay);
		// imp.updateAndDraw();
	}

	// normalize with zero mean and variance 1
	// normalization along column
	// TODO: this function is horrible cause it doesn't use the API of la4j
	// but it is working!
	public  void normalize(Matrix M){
		double[] sumM = M.foldColumns(Vectors.asSumAccumulator(0.0));
		// TODO: Rewrite using vectorization
		for (int i =0; i < D; ++i){
			sumM[i] /= M.getColumn(i).length();
		}
		
		// System.out.println(sumM[0] + " " + sumM[1]);
		
		// TODO: Rewrite using vectorization
		for (int i = 0; i < M.rows(); ++i){
			for (int d = 0 ; d < D; ++d){
				double val = M.get(i, d) - sumM[d];
				M.set(i, d, val);
			}	
		}
		
		// debugResult(M, 0);	
		double scaleVal = 0;
		
		for (int i = 0; i < M.rows(); ++i){
			for (int d = 0; d < D; ++d){
				scaleVal += M.get(i, d)*M.get(i, d);
			}
		}
		
		scaleVal /= (M.rows());
		scaleVal = Math.sqrt(scaleVal);
		System.out.println(scaleVal);
		
		for (int i = 0; i < M.rows(); ++i){
			for (int d = 0 ; d < D; ++d){
				double val = M.get(i, d)/scaleVal;
				M.set(i, d, val);
			}	
		}		
	}

	public int runNonRigidRegistration(int flag) {
		// initialization
		// reading the files

		mX = new Basic2DMatrix(N, D);
		mY = new Basic2DMatrix(M, D);
		mI = new Basic2DMatrix(N, D); // TODO: ??
		mW = new Basic2DMatrix(M, D);
		mG = new Basic2DMatrix(M, M);
		mP = new Basic2DMatrix(M, N);

		mT = new Basic2DMatrix(M, D);

		mW.setAll(0);
		// TODO: this assholes also do normalization for some reasons
		readCSV(); // reading is fine
		normalize(mX);
		normalize(mY);
	
		img = new ArrayImgFactory<FloatType>().create(new long[]{500, 500}, new FloatType());
		addPoints(img);
		// img = ImgLib2Util.openAs32Bit(new File("/Users/kkolyva/Desktop/cpd_ex/worm-folded.tif"));

		// need a container because not yet implemented in I2
		final ImagePlus imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
		imp.show();
		// ImageJFunctions.show(img);
		// column-wise summation
		double[] sumMX = mX.foldColumns(Vectors.asSumAccumulator(0.0));
		double[] sumMY = mY.foldColumns(Vectors.asSumAccumulator(0.0));

		double sum = 0;

		for (int d = 0; d < D; ++d)
			sum += sumMX[d] * sumMY[d];

		sigma2 = (double) (M * mX.transpose().multiplyByItsTranspose().trace()
				+ N * mY.transpose().multiplyByItsTranspose().trace() - 2 * sum) / (M * N * D);

		System.out.println("sigma2 : " + sigma2);
		

		
		
		w = 0.1f;
		beta = 2;
		lambda = 3;

		// there might be a better way for this part
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < M; ++j) {
				double val = (double) mY.getRow(i).subtract(mY.getRow(j)).norm();
				val *= -val / (2 * beta * beta);
				mG.set(i, j, Math.exp(val));
			}
		}
		// stuff above looks valid
		// G looks correct at the moment
		
		

		mT = (Basic2DMatrix) mT.insert(mY.add(mG.multiply(mW)));
		addOverlay(imp);
		
		int counter = 0;
		while (counter++ < 30) {

			System.out.println("ITERATION #" + counter);

			// TODO: OPTIMIZE
			for (int m = 0; m < M; ++m) {
				for (int n = 0; n < N; ++n) {
					double val = (double) mX.getRow(n).subtract(mY.getRow(m).add(mG.getRow(m).multiply(mW))).norm();
					val *= -val / (2 * sigma2);
					val = (double) Math.exp(val);

					// TODO: is it possible to move this part out of the loop
					double denom = 0;
					for (int k = 0; k < M; ++k) {
						double tmp = (double) mX.getRow(n).subtract(mY.getRow(k).add(mG.getRow(k).multiply(mW))).norm();
						tmp *= -tmp / (2 * sigma2);
						denom += Math.exp(tmp);
					}
					denom += w / (1 - w) * Math.pow(2 * Math.PI * sigma2, D / 2) * M / N;
					mP.set(m, n, val / denom);
				}
			}

			// TODO: I think that should be N instead of M here
			// TODO: FIXED: check if this thing is correct at all
			BasicVector ones = BasicVector.constant(N, 1);
			Vector invP = mP.multiply(ones);
			// Up to here output is correct
			invP.update(INV2_FUNCTION); // you will reuse this value !!
			// Up to here output is correct
			// debugResult(invP);
			// if (true) return 0;

			// TODO: now you know that the parameters are fine!
			// check out the values in the left handside and right hand side
			// after that adjust the solver

			// the left hand side is calculated correctly
			// check the right hand-side

			System.out.println(sigma2 + " : " + lambda);

			// Problem pops up after this statement

			for (int d = 0; d < D; ++d) {
				// mW.setColumn(d, new
				// GaussianSolver(mG.add(invP.toDiagonalMatrix().multiply(sigma2
				// * lambda))).solve(
				// (invP.toDiagonalMatrix().multiply(mP).multiply(mX)).subtract(mY).getColumn(d)));
				mW.setColumn(d, new LeastSquaresSolver(mG.add(invP.toDiagonalMatrix().multiply(sigma2 * lambda)))
						.solve((invP.toDiagonalMatrix().multiply(mP).multiply(mX)).subtract(mY).getColumn(d)));
			}

			double N_P = (double) mP.sum();
			mT = (Basic2DMatrix) mT.insert(mY.add(mG.multiply(mW)));

			// debugResult(mT, 0);
			// if (true) return 0;

			// System.out.println("Np : " + N_P);

			double[] terms = new double[3];

			Vector tmpVec = mP.transpose().multiply(BasicVector.constant(M, 1));
			tmpVec.update(INV2_FUNCTION);

			terms[0] = (double) (mX.transpose()
					.multiply(mP.transpose().multiply(BasicVector.constant(M, 1)).toDiagonalMatrix().multiply(mX))
					.trace());
			terms[1] = (double) (mP.multiply(mX).transpose().multiply(mT).trace()) * (-2);
			terms[2] = (double) mT.transpose()
					.multiply(mP.multiply(BasicVector.constant(N, 1)).toDiagonalMatrix().multiply(mT)).trace();

			sigma2 = 0;
			for (int i = 0; i < terms.length; ++i) {
				sigma2 += terms[i];
				// System.out.println(terms[i]);
			}

			sigma2 /= (N_P * D);
			System.out.println(sigma2);

			// if (true) return 0;
			addOverlay(imp);
		}

		// writeCSV();

		// T is the transformation you are looking for
		System.out.println("hello");
		return 1;
	}

	// add parameters
	// matrix file name
	// ugly function to read mx and mY from CSV
	public void readCSV() {

		// int N = 112;
		// int M = 133;
		// int D = 2;
		// mX = new Basic2DMatrix(N, D);
		// mY = new Basic2DMatrix(M, D);

		CSVReader reader = null;
		String[] nextLine;

		try {
			reader = new CSVReader(new FileReader("/Users/kkolyva/Desktop/cpd_ex/worm-folded.csv"), '\t');
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		try {
			int i = 0;
			while ((nextLine = reader.readNext()) != null) {
				mX.set(i, 0, Double.parseDouble(nextLine[0]));
				mX.set(i, 1, Double.parseDouble(nextLine[1]));
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		try {
			reader = new CSVReader(new FileReader("/Users/kkolyva/Desktop/cpd_ex/worm-straight.csv"), '\t');
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		try {
			int i = 0;
			while ((nextLine = reader.readNext()) != null) {
				mY.set(i, 0, Double.parseDouble(nextLine[0]));
				mY.set(i, 1, Double.parseDouble(nextLine[1]));
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

	public static void main(String[] args) {
		new ImageJ();
		new CoherentPointDrift().runNonRigidRegistration(0);
		// new CoherentPointDrift().readCSV();
	}
}
