package coherent.point.drift;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.la4j.LinearAlgebra;
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
	ImagePlus imp;

	// TODO: read this parameters from file or from variable
	final int D; // dimensionality of the point set
	final int N; // # of points in the first point set
	final int M; // # of points in the second point set

	double sigma2; // TODO: 2 for squared

	// parameters
	double w; // amount of noise [0, 1]
	double beta; // Briefly speaking, parameter defines the model of the
	// smoothness regularizer (width of smoothing Gaussian
	// filter in (20)).
	double lambda; // Parameter represents the trade-off between the goodness of
	// maximum likelihood fit and regularization.

	// TODO: convert stuff to the HPC
	final Matrix mX; // N D
	final Matrix mY; // M D
	final Matrix mI; //
	final Matrix mW; // M D
	final Matrix mG; // M M
	final Matrix mP; // M N

	final Matrix mT; // M D transformation matrix

	// straight to the non-rigid registration

	// these guys are necessary for proper coordinates adjustment
	double[] translate;
	double[] scale;
	double[] sigma;

	// constructor
	public CoherentPointDrift(Img<FloatType> img, int N, int M, int D, double w, double beta, double lambda) {

		this.img = img;

		this.w = w;
		this.beta = beta;
		this.lambda = lambda;

		this.D = D; // dimensionality of the point set
		this.N = N; // # of points in the first point set
		this.M = M; // # of points in the second point set

		mX = new Basic2DMatrix(N, D);
		mY = new Basic2DMatrix(M, D);
		mI = new Basic2DMatrix(N, D); // TODO: ??
		mW = new Basic2DMatrix(M, D);
		mG = new Basic2DMatrix(M, M);
		mP = new Basic2DMatrix(M, N);

		mT = new Basic2DMatrix(M, D);

		mW.setAll(0);
		// settings for ouput image
		// these guys are necessary for proper coordinates adjustment
		// TODO: make this things calculated automatically
		translate = new double[] { 250, 250 };
		scale = new double[] { 150, 150 };
		sigma = new double[] { 3, 3 };
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

	// add points that has to be detected
	public void addPoints(Img<FloatType> img) {
		double[] location = new double[D];
		for (int i = 0; i < N; ++i) {
			for (int d = 0; d < D; ++d) {
				location[d] = (mX.get(i, d) * scale[d] + translate[d]);
			}
			// TODO: Use function from klim.utils
			TestGauss3d.addGaussian(img, location, sigma);
		}
	}

	public void readData(Matrix X, Matrix Y) {
		readCSV(); // reading is fine
		normalize(mX);
		normalize(mY);
	}

	/**
	 * Shows coordinates defined by A in the imp image
	 * 
	 */
	public void addOverlay(ImagePlus imp, Matrix A) {
		int numDimensions = A.columns();
		int numPoints = A.rows();

		Overlay overlay = imp.getOverlay();
		if (overlay == null) {
			overlay = new Overlay();
		}
		overlay.clear();

		double[] location = new double[numDimensions];
		for (int i = 0; i < numPoints; ++i) {
			for (int d = 0; d < numDimensions; ++d) {
				location[d] = A.get(i, d) * scale[d] + translate[d];
			}
			final OvalRoi or = new OvalRoi(location[0] - sigma[0], location[1] - sigma[0], Util.round(2 * sigma[0]),
					Util.round(2 * sigma[1]));
			or.setStrokeColor(Color.RED);
			overlay.add(or);
		}
		imp.setOverlay(overlay);
	}

	// normalize with zero mean and variance 1
	// normalization along column
	// TODO: this function is horrible cause it doesn't use the API of la4j
	// but it is working!
	public void normalize(Matrix M) {
		double[] sumM = M.foldColumns(Vectors.asSumAccumulator(0.0));
		// TODO: Rewrite using vectorization
		for (int i = 0; i < D; ++i) {
			sumM[i] /= M.getColumn(i).length();
		}

		// System.out.println(sumM[0] + " " + sumM[1]);

		// TODO: Rewrite using vectorization
		for (int i = 0; i < M.rows(); ++i) {
			for (int d = 0; d < D; ++d) {
				double val = M.get(i, d) - sumM[d];
				M.set(i, d, val);
			}
		}

		// debugResult(M, 0);
		double scaleVal = 0;

		for (int i = 0; i < M.rows(); ++i) {
			for (int d = 0; d < D; ++d) {
				scaleVal += M.get(i, d) * M.get(i, d);
			}
		}

		scaleVal /= (M.rows());
		scaleVal = Math.sqrt(scaleVal);
		System.out.println(scaleVal);

		for (int i = 0; i < M.rows(); ++i) {
			for (int d = 0; d < D; ++d) {
				double val = M.get(i, d) / scaleVal;
				M.set(i, d, val);
			}
		}
	}

	// calculates sigma squared
	public double getSigma2(Matrix X, Matrix Y) {
		long M = Y.rows();
		long N = X.rows();
		long D = X.columns(); // Y.columns();

		double res = 0;
		// column-wise summation
		double[] sumMX = X.foldColumns(Vectors.asSumAccumulator(0.0));
		double[] sumMY = Y.foldColumns(Vectors.asSumAccumulator(0.0));

		double sum = 0;
		for (int d = 0; d < D; ++d)
			sum += sumMX[d] * sumMY[d];

		res = (double) (M * X.transpose().multiplyByItsTranspose().trace()
				+ N * Y.transpose().multiplyByItsTranspose().trace() - 2 * sum) / (M * N * D);

		return res;
	}

	// calculates G matrix
	// G might be symmetric
	public void calculateG(Matrix Y, Matrix G, double beta) {
		long M = Y.rows();

		// there might be a better way for this part
		// TODO: make an upper/lower triangle matrix
		// G is symmetric
		for (int i = 0; i < M; ++i) {
			for (int j = 0; j < M; ++j) {
				double val = (double) Y.getRow(i).subtract(Y.getRow(j)).norm();
				val *= -val / (2 * beta * beta);
				G.set(i, j, Math.exp(val));
			}
		}
	}

	// update sigma squared value
	public double updateSigma2(Matrix X, Matrix Y, Matrix P, Matrix T) {
		double res = 0;

		double[] terms = new double[3];

		int M = Y.rows();
		int N = X.rows();
		int D = X.columns();

		terms[0] = (X.transpose()
				.multiply(P.transpose().multiply(BasicVector.constant(M, 1)).toDiagonalMatrix().multiply(X)).trace());
		terms[1] = (P.multiply(X).transpose().multiply(T).trace()) * (-2);
		terms[2] = T.transpose().multiply(P.multiply(BasicVector.constant(N, 1)).toDiagonalMatrix().multiply(T))
				.trace();

		// sigma2 = 0; // ?
		for (int i = 0; i < terms.length; ++i) {
			res += terms[i];
		}

		double N_P = P.sum();
		res /= (N_P * D);

		return res;
	}

	// TODO:
	public double calculateError(Matrix W, Matrix G, double lambda, double LOld) {
		double res = 0;

		return res;
	}

	// TODO: rewrite all assignments -> in place!

	public int runNonRigidRegistration(int flag) {
		// TODO: update to use the functionality of la4j
		readData(mX, mY);
		// fill image with points that has to be detected
		addPoints(img);
		imp = ImageJFunctions.wrapFloat(img, "Put an overlay on top!");
		imp.show();

		sigma2 = getSigma2(mX, mY);
		calculateG(mY, mG, beta);

		mY.add(mG.multiply(mW)).apply(LinearAlgebra.IN_PLACE_COPY_MATRIX_TO_MATRIX, mT);
		addOverlay(imp, mT);

		int counter = 0;
		while (counter++ < 30) {
			// TODO: should be info for current iteration
			// # + error
			System.out.println("ITERATION #" + counter);

			// TODO: OPTIMIZE
			// Need paper on fast gauss transform
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

			Vector invP = mP.multiply(BasicVector.constant(N, 1));
			// TODO: Rewrite function below
			// use in build functions
			invP.update(INV2_FUNCTION); // you will reuse this value !!

			// using least square method here
			for (int d = 0; d < D; ++d) {
				mW.setColumn(d, new LeastSquaresSolver(mG.add(invP.toDiagonalMatrix().multiply(sigma2 * lambda)))
						.solve((invP.toDiagonalMatrix().multiply(mP).multiply(mX)).subtract(mY).getColumn(d)));
			}

			mY.add(mG.multiply(mW)).apply(LinearAlgebra.IN_PLACE_COPY_MATRIX_TO_MATRIX, mT);
			sigma2 = updateSigma2(mX, mY, mP, mT);

			addOverlay(imp, mT);
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

		double w = 0.1;
		double beta = 2;
		double lambda = 3;

		final Img<FloatType> img = new ArrayImgFactory<FloatType>().create(new long[] { 500, 500 }, new FloatType());

		// read data first
		// TODO: reading here
		int D = 2; // dimensionality of the point set
		int N = 133; // # of points in the first point set
		int M = 112; // # of points in the second point set

		// pass it as arguments
		new CoherentPointDrift(img, N, M, D, w, beta, lambda).runNonRigidRegistration(0);
		// new CoherentPointDrift().readCSV();
	}
}
