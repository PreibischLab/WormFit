package coherent.point.drift.refactored.tests;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import coherent.point.drift.refactored.ApacheCPD;
import coherent.point.drift.refactored.IOUtils;
import coherent.point.drift.refactored.Utils;

public class MatrixOperationsTests {
	
	private static double eps = 1e-8;

	protected static RealMatrix simpleZeroMeanUnitVariance(double [][] mA, int n, int numDimensions){
		double [] sumOverColumns = new double [n];
	
		for (int d = 0; d < numDimensions; ++d) {
			sumOverColumns[d] = 0;
			for (int j = 0; j < n; j++)
				sumOverColumns[d] += mA[j][d];
		}
		
		for (int d = 0; d < numDimensions; ++d) {
			sumOverColumns[d] /= n;
		}
		
		for (int d = 0; d < numDimensions; ++d)
			for (int j = 0; j < n; j++)
				mA[j][d] -= sumOverColumns[d];
		// at this point the mean is 0!
		
		for (int d = 0; d < numDimensions; ++d) {
			sumOverColumns[d] = 0;
			for (int j = 0; j < n; j++)
				sumOverColumns[d] += mA[j][d]*mA[j][d];
		}
		
		for (int d = 0; d < numDimensions; ++d)
			sumOverColumns[d] = Math.sqrt(sumOverColumns[d]/n);
		
		
		for (int d = 0; d < numDimensions; ++d) {
			for (int j = 0; j < n; j++)
				mA[j][d] /= sumOverColumns[d];
		}
		
		return MatrixUtils.createRealMatrix(mA);
	}
	
	//FIXME: TEST
	public static void testNormalize(){
		RealMatrix matrix = MatrixUtils.createRealMatrix(3, 3);
		
		for (int j = 0; j < matrix.getColumnDimension(); ++j)
			for (int i = 0; i < matrix.getRowDimension(); ++i)
				matrix.setEntry(i, j, (j + 2)*(i+1));
		
		RealMatrix foundMatrix = simpleZeroMeanUnitVariance(matrix.getData(), 3, 3);
		RealMatrix realMatrix = ApacheCPD.normalize(matrix);
		
		double res = realMatrix.subtract(foundMatrix).getNorm();
		
		if (res > eps)
			throw new RuntimeException(String.format("ApacheCPD.normalize is not working properly; Norm is too high: %f", res));
		
		System.out.println("ApacheCPD.normalize works fine");
	}
	

	public static void main(String[] args) {
		testNormalize();
		System.out.println("Done!");
	}

}
