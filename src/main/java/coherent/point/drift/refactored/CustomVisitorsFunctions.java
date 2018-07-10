package coherent.point.drift.refactored;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.RealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrixPreservingVisitor;

public class CustomVisitorsFunctions {
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
			return val/scale;
		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}
	
	// this is a matrix visitor that returns the scaled value over columns
	public class ColumnDivideVisitor implements RealMatrixChangingVisitor{
		final double [] columnScale; 

		public ColumnDivideVisitor(double[] columnScale){
			this.columnScale = columnScale;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){
			// empty
		}

		public double visit(int row, int col, double val){
			return val/columnScale[col];
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
			// empty
		}

		public void visit(int row, int col, double val){
			sumOverColumns[col] += val;

		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}
	
	// this is a matrix visitor that returns the squared sum over columns
	public class ColumnSquareSumVisitor implements RealMatrixPreservingVisitor{
		final double [] sqSumOverColumns; 

		public ColumnSquareSumVisitor(double[] sqSumOverColumns){
			this.sqSumOverColumns = sqSumOverColumns;
		}

		public void start (int rows, int columns,
				int startRow, int endRow, int startColumn, int endColumn){
			// empty
		}

		public void visit(int row, int col, double val){
			sqSumOverColumns[col] += val*val;
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
			// empty
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
	
	public static class ElementwiseInverse implements UnivariateFunction {
		public double value(double x) {
			return 1.0/x;
		}
	}
}
