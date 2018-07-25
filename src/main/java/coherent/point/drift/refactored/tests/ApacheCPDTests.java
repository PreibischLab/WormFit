package coherent.point.drift.refactored.tests;

import org.apache.commons.math3.linear.RealMatrix;

import coherent.point.drift.refactored.ApacheCPD;
import coherent.point.drift.refactored.utils.IOUtils;

public class ApacheCPDTests {
	
	// non rigid test with for coherent point drift
	public static void testNonRigid(){
		double w = 0.1; 
		double beta = 2;
		double lambda = 3;
		int maxIteration = 10;

		// read data first
		String path = "./src/main/resources/";
		String from = path + "worm-folded.csv";
		String to = path + "worm-straight.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, '\t');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, '\t');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runNonRigidRegistration();
	}
	
	// affine test with for coherent point drift
	public static void testAffine(){
		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 10;
		
		// read data first
		String path = "./src/main/resources/";
		String from = path + "x.csv";
		String to = path + "y.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, ',');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, ',');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runAffineRegistration();
	}
	
// rigid test with for coherent point drift
	public static void testRigid(){
		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 10;
		
		// read data first
		String path = "./src/main/resources/";
		String from = path + "x.csv";
		String to = path + "y.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, ',');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, ',');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runRigidRegistration();
	}
	
	public static void main(String[] args) {
		testNonRigid();
	}
}
