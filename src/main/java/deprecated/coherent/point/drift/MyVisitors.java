package coherent.point.drift;

import org.apache.commons.math3.linear.RealMatrixChangingVisitor;

public class  MyVisitors {
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
			return val/scale; // val/scaleValue;
		}

		public double end(){
			return 0; // 0 for everything is fine 
		}
	}

	// this is a matrix visitor that returns 
}
