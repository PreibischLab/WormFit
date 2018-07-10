package coherent.point.drift.refactored;

import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.real.FloatType;

import org.apache.commons.math3.linear.RealMatrix;

import ij.ImageJ;

public class ApacheCPDTests {
	

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

		// String path = "/home/milkyklim/Documents/imglib2Dev/WormFit/src/main/resources/";
		String path = "/Users/kkolyva/workspace/WormFit/src/main/resources/";
		String from = path + "worm-folded.csv";
		String to = path + "worm-straight.csv";
		
		// vis parameters
		// int radius = 2; 
		// long [] shift = new long [D];
//		for (int d = 0; d < D; d++) {
//			// extra padding to fit the whole cells
//			shift[d] = (img.dimension(d) + 2*radius + 5)/2;
//		}
		
		// pass it as arguments
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, '\t');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, '\t');
		
		// new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration).runNonRigidRegistration(0, from, to) ;
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runNonRigidRegistration(mX, mY);
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
		
		// vis parameters
		int radius = 10; 
		long [] shift = new long [D];
		for (int d = 0; d < D; d++) {
			// extra padding to fit the whole cells
			shift[d] = (img.dimension(d) + 2*radius + 5)/2;
		}

		// pass it as arguments
		// new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration, radius, shift).runAffineRegistration(0, from, to);
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
		
		// vis parameters
		int radius = 10; 
		long [] shift = new long [D];
		for (int d = 0; d < D; d++) {
			// extra padding to fit the whole cells
			shift[d] = (img.dimension(d) + 2*radius + 5)/2;
		}
		
		
		// pass it as arguments
		// new ApacheCPD(img, N, M, D, w, beta, lambda, maxIteration, radius, shift).runRigidRegistration(0, from, to);
		// new CoherentPointDrift().readCSV();
	}
}
