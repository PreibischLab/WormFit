package coherent.point.drift.refactored.tests;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import coherent.point.drift.refactored.utils.IOUtils;

public class IOUtilsTests {
	
	private static double eps = 1e-8;
	
	public static void testGetNumLines() {
		String filepath = "/Users/kkolyva/workspace/WormFit/src/main/resources/coherent-point-drift-tests/coordinates-3D.csv";
		char delimiter = '\t';
		
		int realNumLines = 4;
		int foundNumLines = IOUtils.getNumLines(filepath, delimiter);
		
		if (realNumLines != foundNumLines)
			throw new RuntimeException(String.format("IOUtils.getNumLines is not working properly; Lines should be %d but found %d", realNumLines, foundNumLines));
		
		System.out.println("IOUtils.getNumLines works fine");
	}
	
	public static void testGetNumDimensions() {
		String filepath = "/Users/kkolyva/workspace/WormFit/src/main/resources/coherent-point-drift-tests/coordinates-3D.csv";
		char delimiter = '\t';
		
		int realNumDimensions= 3;
		int foundNumDimensions = IOUtils.getNumDimensions(filepath, delimiter);
		
		if (realNumDimensions != foundNumDimensions)
			throw new RuntimeException(String.format("IOUtils.getNumDimensions is not working properly; Dimensions should be %d but found %d", realNumDimensions, foundNumDimensions));
		
		System.out.println("IOUtils.getNumDimensions works fine");
	}

	public static void testReadPositionsFromCSV() {
		String filepath = "/Users/kkolyva/workspace/WormFit/src/main/resources/coherent-point-drift-tests/coordinates-3D.csv";
		char delimiter = '\t';
		
		RealMatrix realMatrix = MatrixUtils.createRealMatrix(new double [][] {
			{ 466.670822281167, 485.9442086648983, 49.94341290893015 },
			{ 588.9463278667272, 576.7424645102742, 55.12458676346665 },
			{ 481.41731797614153, 574.8794734677087, 54.50010283833813 },
			{ 345.9463278667272, 334.9442086648983, 37.50010283833813 }
		});
		
		RealMatrix foundMatrix = IOUtils.readPositionsFromCSV(filepath, delimiter);
		
		if (realMatrix.subtract(foundMatrix).getNorm() > eps)
			throw new RuntimeException(String.format("IOUtils.readPositionsFromCSV is not working properly; Norm is too high: %f", realMatrix.subtract(foundMatrix).getNorm()));
		
		System.out.println("IOUtils.readPositionsFromCSV works fine");
	}
	
	public static void testWritePositionsToCSV() {
		String filepath = "/Users/kkolyva/workspace/WormFit/src/main/resources/coherent-point-drift-tests/test-coordinates-3D.csv";
		char delimiter = '\t';
		
		RealMatrix realMatrix = MatrixUtils.createRealMatrix(new double [][] {
			{ 466.670822281167, 485.9442086648983, 49.94341290893015 },
			{ 588.9463278667272, 576.7424645102742, 55.12458676346665 },
			{ 481.41731797614153, 574.8794734677087, 54.50010283833813 },
			{ 345.9463278667272, 334.9442086648983, 37.50010283833813 }
		});
		
		IOUtils.writePositionsToCSV(realMatrix, filepath, delimiter);
		RealMatrix foundMatrix = IOUtils.readPositionsFromCSV(filepath, delimiter);
		
		if (realMatrix.subtract(foundMatrix).getNorm() > eps)
			throw new RuntimeException(String.format("IOUtils.writePositionsToCSV is not working properly; Norm is too high: %f", realMatrix.subtract(foundMatrix).getNorm()));
		
		System.out.println("IOUtils.writePositionsToCSV works fine");
	}

}
