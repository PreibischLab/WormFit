package klim;

import java.util.ArrayList;
import java.util.Iterator;

import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.labeling.AllConnectedComponents;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.labeling.Labeling;
import net.imglib2.labeling.NativeImgLabeling;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.view.Views;

public class ObjectSegmentation {
	
	private final long[] min; 
	private final long[] max; 
	private final long[] size;
	private long index;
	
	public ObjectSegmentation(long[] min, long[] max, long[] size, long index){
		this.min = new long[min.length];
		this.max = new long[max.length];
		this.size = new long[size.length];
		this.index = index;
		
		for (int d = 0; d < min.length; d++){
			this.min[d] = min[d];
			this.max[d] = max[d];
			this.size[d] = size[d];
		}
	}
	
	// returns labeling for each connected component 
	public static <T extends RealType<T>> void setLabeling(RandomAccessibleInterval<BitType> bitImg, ImgLabeling<Integer, IntType> labeling) {
		final Iterator<Integer> labelIterator = AllConnectedComponents.getIntegerNames(0);
		ConnectedComponents.labelAllConnectedComponents(bitImg, labeling, labelIterator, ConnectedComponents.StructuringElement.EIGHT_CONNECTED );		
	}
	
	public static <T extends RealType<T>> void setParameters(final RandomAccessible<T> img, final ImgLabeling<Integer, IntType> labeling){
		
		ArrayList<ObjectSegmentation> objects = new ArrayList<ObjectSegmentation>();
		
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();
		
		// long[] dummy = new long[randomAccess.numDimensions()];
		
		long curNum = 0;
		while(cursor.hasNext()){
			
			long[] max = new long[randomAccess.numDimensions()];
			long[] min = new long[randomAccess.numDimensions()];
			
			for (int d = 0; d < img.numDimensions(); d++){
				max[d] = 0;
				min[d] = Long.MAX_VALUE;
			}
			
			cursor.fwd();
			long curElement = cursor.get().get();
			// to prevent crashes
			if (curElement >= curNum){
				while(curElement >= curNum){
					objects.add(new ObjectSegmentation(min, max, max, curElement));
					++curNum; 
				}
			}
			try{
				long [] curPosition = new long[randomAccess.numDimensions()];
				cursor.localize(curPosition);
				int index = cursor.get().get(); 
				
				for (int d = 0; d < img.numDimensions(); d++){
					if (objects.get(index).max[d] < curPosition[d]){
						objects.get(index).max[d] = curPosition[d];
					}
					
					if (objects.get(index).min[d] > curPosition[d]){
						objects.get(index).min[d] = curPosition[d];
					}
				}
//				
				
			}
			catch(Exception e){
				System.out.println("Lol! You failed!");
			}
		}
		
		for (ObjectSegmentation obj : objects){
			for (int d = 0; d < img.numDimensions(); d++){
				System.out.print("[" + obj.min[d] + " " + obj.max[d] + "]");
			}
			System.out.println();
		}
		
		
	}
	
	public static <T extends RealType<T>> void segmentImage(final ImgLabeling<Integer, IntType> labeling, final RandomAccessible<T> img){
		// this one bellow is the final list you'll need 
		ArrayList<PointSampleList<T>> objectsList = new ArrayList<>();

		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();

		// number of objects in the picture calculated dynamically
		// this part divides pixels into chunks
		int curMax = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			int curElement = cursor.get().get();
			// increase the size of the objects list
			if (curElement >= curMax){
				while(curElement >= curMax){
					objectsList.add(new PointSampleList<T>(img.numDimensions()));
					++curMax; 
				}
			}
			try{
				randomAccess.setPosition(cursor);
				objectsList.get(curElement).add(new Point(randomAccess), randomAccess.get()); // .copy here ?
			}
			catch(Exception e){
				System.out.println(objectsList.get(curElement).size());
			}
		}
		
		
		// search for worm index 
		// taking into account that worm is the largest object
//		int idx = 0; 
//		long maxSize = 0;
//		// skip background i == 0 
//		for (int i = 1; i < objectsList.size(); ++i){
//			if (objectsList.get(i).size() > maxSize){
//				maxSize = objectsList.get(i).size();
//				idx = i;
//			}			
//		}		

		// copy the output
//		Cursor<T> it = objectsList.get(idx).cursor(); // copy worm
//		while(it.hasNext()){
//			it.fwd();
//			worm.add(new Point(it), it.get());
//		}
	}
}
