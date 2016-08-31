package klim;

import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

public class ObjectSegmentation {
	// returns labeling for each connected component 
	public static <T extends RealType<T>> void setLabeling(RandomAccessibleInterval<BitType> bitImg, ImgLabeling<Integer, IntType> labeling) {
		final Iterator<Integer> labelIterator = AllConnectedComponents.getIntegerNames(0);
		ConnectedComponents.labelAllConnectedComponents(bitImg, labeling, labelIterator, ConnectedComponents.StructuringElement.EIGHT_CONNECTED );		
	}
	
	// find the max label //TODO: maybe this one is alreay implemented
	public static int getMaxLabel(final ImgLabeling<Integer, IntType> labeling){
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		int maxLabel = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			if(cursor.get().get() > maxLabel)
				maxLabel = cursor.get().get();
		}
		return maxLabel;
	} 

	// searches for the center of beads
	public static <T extends RealType<T>> void findBeads(final RandomAccessibleInterval<T> img, final ImgLabeling<Integer, IntType> labeling, PointSampleList<T> beads){
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();

		int numBeads = 0;
		numBeads = getMaxLabel(labeling) + 1;
		
		System.out.println("Number of beads found: " + numBeads);
		
		long[][] averageBeadsCoordinates = new long[numBeads][randomAccess.numDimensions()];
		long[] pixelsPerBead = new long[numBeads];

		while(cursor.hasNext()){
			cursor.fwd();
			int idx = cursor.get().getInteger();
			for (int d = 0; d < img.numDimensions(); d++){
				averageBeadsCoordinates[idx][d] += cursor.getLongPosition(d);
			}	
			pixelsPerBead[idx]++;
		}

		for (int j = 0; j < numBeads; ++j){
			for (int d = 0; d < img.numDimensions(); d++)
				averageBeadsCoordinates[j][d] /= pixelsPerBead[j];
			if (j != 0) // skip the background
				beads.add(new Point(averageBeadsCoordinates[j]), randomAccess.get());
			
		}		
	}
}
