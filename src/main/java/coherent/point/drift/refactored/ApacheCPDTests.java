package coherent.point.drift.refactored;

import org.apache.commons.math3.linear.RealMatrix;

public class ApacheCPDTests {
	
	// non rigid test with for coherent point drift
	public void testNonRigid(){
		double w = 0.1; 
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;

		// read data first
		String path = "./src/main/resources/";
		String from = path + "worm-folded.csv";
		String to = path + "worm-straight.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, '\t');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, '\t');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runNonRigidRegistration(mX, mY);
	}
	
	// affine test with for coherent point drift
	public void testAffine(){
		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;
		
		// read data first
		String path = "./src/main/resources/";
		String from = path + "fishy-fish-x.csv";
		String to = path + "fishy-fish-y.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, '\t');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, '\t');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runAffineRegistration(mX, mY);
	}
	
// rigis test with for coherent point drift
	public void testRigid(){
		double w = 0.1;
		double beta = 2;
		double lambda = 3;
		int maxIteration = 30;
		
		// read data first
		String path = "./src/main/resources/";
		String from = path + "fishy-fish-rigid-x.csv";
		String to = path + "fishy-fish-rigid-y.csv";
		
		RealMatrix mX = IOUtils.readPositionsFromCSV(from, '\t');
		RealMatrix mY = IOUtils.readPositionsFromCSV(to, '\t');
		
		new ApacheCPD(mX, mY, w, beta, lambda, maxIteration).runRigidRegistration(mX, mY);
	}
}
