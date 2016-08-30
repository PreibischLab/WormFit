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
	private int numDimensions;

	public ObjectSegmentation(long[] min, long[] max, long[] size, long index){
		this.numDimensions = min.length;
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

	public static <T extends RealType<T>> void findBeads(final RandomAccessibleInterval<T> img, final ImgLabeling<Integer, IntType> labeling, PointSampleList<T> beads){
		Cursor<IntType> cursor = Views.iterable(labeling.getIndexImg()).cursor();
		RandomAccess <T> randomAccess = img.randomAccess();

		int numBeads = 0;
		numBeads = getMaxLabel(labeling) + 1;
		
		System.out.println(numBeads);
		long[][] averageBeadsCoordinates = new long[numBeads][randomAccess.numDimensions()];
		long[] numPossibleBeads = new long[numBeads];

		while(cursor.hasNext()){
			cursor.fwd();
			int idx = cursor.get().getInteger();
			for (int d = 0; d < img.numDimensions(); d++){
				averageBeadsCoordinates[idx][d] += cursor.getLongPosition(d);
			}	
			numPossibleBeads[idx]++;
		}

//		for (int j = 0; j < numBeads; ++j){
//			System.out.print(numPossibleBeads[j] + " ");
//		}
		
		System.out.println();
		for (int j = 0; j < numBeads; ++j){
			for (int d = 0; d < img.numDimensions(); d++){
				averageBeadsCoordinates[j][d] /= numPossibleBeads[j];
				System.out.print(averageBeadsCoordinates[j][d] + " ");
			}
			beads.add(new Point(averageBeadsCoordinates[j]), randomAccess.get());
			System.out.println();
		}
		
	}

	public static <T extends RealType<T>> ArrayList<ObjectSegmentation> setParameters(final RandomAccessible<T> img, final ImgLabeling<Integer, IntType> labeling){

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
			// TODO: switch to if else
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
			}
			catch(Exception e){
				System.out.println("Lol! You failed!");
			}
		}

		// DEBUG:
		//		for (ObjectSegmentation obj : objects){
		//			for (int d = 0; d < img.numDimensions(); d++){
		//				System.out.print("[" + obj.min[d] + " " + obj.max[d] + "]");
		//			}
		//			System.out.println();
		//			// ImageJFunctions.show(Views.interval(img, obj.min, obj.max));
		//		}



		// TODO: make a function
		for (int j = 0; j < objects.size(); ++j){
			for (int d = 0; d < objects.get(j).numDimensions; ++d){
				objects.get(j).size[d] = objects.get(j).max[d] - objects.get(j).min[d];
				if (objects.get(j).size[d] < 0)
					System.out.println("Size value is negative!");
			}
		}

		System.out.println("Size: " + objects.size());

		// TODO: rewrite as a proper function
		objects.remove(0); // drop background
		return objects;

	}

	public static <T extends RealType<T>>void removeNoise(ArrayList<ObjectSegmentation> beads, RandomAccessibleInterval<T> img){
		System.out.println("Denoizing");

		int numDimensions = beads.get(0).numDimensions;
		for (int j = 0; j < beads.size(); ++j){
			long [] position = new long[numDimensions];
			boolean isNoise = false; 
			for (int d = 0; d < beads.get(j).numDimensions; ++d){
				if(beads.get(j).size[d] <= 4)
					isNoise = true;
			}
			if (isNoise){
				beads.remove(j);
				j--;
			}
			// System.out.println(j);
		}

		for (ObjectSegmentation obj : beads){
			for (int d = 0; d < numDimensions; d++){
				System.out.print("[" + obj.min[d] + " " + obj.max[d] + "]");
			}
			System.out.println();
			ImageJFunctions.show(Views.rotate(Views.interval(img, obj.min, obj.max), 0, 2));
		}


	}

}
