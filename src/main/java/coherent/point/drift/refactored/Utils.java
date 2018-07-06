package coherent.point.drift.refactored;

import org.apache.commons.math3.linear.RealMatrix;

public class Utils {

	public static void printMatrix(RealMatrix M, int ColumnDimension, int RowDimension){
		for (int j = 0; j < RowDimension; j++ ){
			for (int i = 0; i < ColumnDimension; i++ ){
				System.out.print(M.getEntry(j, i) + " ");
			}
			System.out.println();
		}
	}	
	
	public static void printLog(int idx, double sigma, double error){
		System.out.println("ITERATION #" + idx + ":");
		System.out.println("Sigma squared: " + String.format(java.util.Locale.US, "%.2e", sigma*sigma));
		// System.out.println("Error: " + String.format(java.util.Locale.US, "%.2e", error));
	}
	
}
