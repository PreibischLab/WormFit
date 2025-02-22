package klim;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.RectangleShape;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class CentralCurveDetection {

	private static final boolean showImage = true;
	private static final boolean inProgress = false;

	// detects the start and end points
	public static <T extends RealType<T>, U extends RealType<U>> void detectEnds(PointSampleList<T> endPoints, PointSampleList<U> wormOutline, Img <T> notDst, T intensity){
		// TODO: automate search
		// Use second derivative (steep gradient) to detect edges
		// use sliding triangular or quad window
		// start 42 846
		// end 985 485
		endPoints.add(new Point(new long[]{42, 846}), intensity);
		endPoints.add(new Point(new long[]{985, 485}), intensity);



		if (!inProgress){
			// real search 
			// TODO: adjust to 3D case 
			//for (int d = 0; d < endPoints.numDimensions() - 1; ++d){
			ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
			Img< FloatType > tmpImg = imgFactory.create(notDst, new FloatType() );
			RandomAccess<FloatType> ra = tmpImg.randomAccess();		

			Cursor<U> wormOulineCursor = wormOutline.cursor(); 

			// wormOutline.firstElement();

			long idx = 0;
			while(wormOulineCursor.hasNext()){
				wormOulineCursor.fwd();
				if (idx % 5 ==0){
					ra.setPosition(wormOulineCursor);
					ra.get().set(255);
				}
				idx++;
				for (int d = 0; d < wormOutline.numDimensions(); ++d){
					System.out.print(wormOulineCursor.getIntPosition(d) + " ");
				}
				System.out.println();

			}

			ImageJFunctions.show(tmpImg).setTitle("DEBUG");

			//}

		}
	}
	// detect the direction of movement
	public static void detectDirection(final long[] pPos, final long [] cPos, final long [] minmax, int n){

		int numMoves = 0;
		final long[] dir = new long[n];

		for (int d = 0; d < n; ++d){
			dir[d] = cPos[d] - pPos[d];
			if (dir[d] == 1 ||  dir[d] == -1) // better use != 0
				numMoves++;
		}

		// treat 1D- 2D- 3D- movement cases differently	
		if (numMoves == 1){
			// works correctly for 1D: 2D: 3D: cases
			for(int d = 0; d < n; ++d){
				minmax[d] = cPos[d] + dir[d] + (dir[d] == 0 ? -1 : 0);
				minmax[d + cPos.length] = cPos[d] + dir[d] + (dir[d] == 0 ? 1 : 0);
			}
		}

		if (numMoves == 2){
			for(int d = 0; d < n; ++d){
				// TODO: rewrite this loop in the if else fashion
				// TODO: simplify formulas! dir[d] - 1
				// works correctly for 1D: 2D: 3D: cases
				long addMin = 0;
				long addMax = 0;
				if (dir[d] > 0){
					addMin = dir[d] - 1;
					addMax = dir[d] + 0;
				}
				else{
					if(dir[d] == 0){
						addMin = dir[d] - 1;
						addMax = dir[d] + 1;
					}
					else{
						// dir[d] < 0
						addMin = dir[d] - 0;
						addMax = dir[d] + 1;
					}
				}

				minmax[d] = cPos[d] + addMin;   
				minmax[d + cPos.length] = cPos[d] + addMax; 

				// OLD VERSION WORKS ONLY FOR 2D
				// minmax[d] = cPos[d] + dir[d] + (dir[d] > 0 ? -1 : 0); // 1 instead of 0 ? 
				// minmax[d + cPos.length] = cPos[d] + dir[d] + (dir[d] < 0 ? 1 : 0); // -1 instead of 0?

			}
		}

		if (numMoves == 3){
			for(int d = 0; d < n; ++d){
				// TODO: rewrite this loop in the if else fashion
				// TODO: simplify formulas! dir[d] - 1
				// works correctly for 1D: 2D: 3D: cases
				long addMin = 0;
				long addMax = 0;
				if (dir[d] > 0){
					addMin = dir[d] - 1;
					addMax = dir[d] + 0;
				}
				else{
					if(dir[d] == 0){ 
						// if this is correct, this case is never the case :) 
						addMin = dir[d] - 1;
						addMax = dir[d] + 1;
						System.out.println("Congrats! You've found the bug");
					}
					else{
						// dir[d] < 0
						addMin = dir[d] - 0;
						addMax = dir[d] + 1;
					}
				}

				minmax[d] = cPos[d] + addMin;   
				minmax[d + cPos.length] = cPos[d] + addMax; 
			}
		}
	}


	// this one defines the neighbors we should consider
	// TODO: Deprecate! 
	public static void setMinMax(long[] minmax, long[] cPos){
		long neighbors = 1;

		for (int i = 0; i < cPos.length; ++i){
			minmax[i] = cPos[i] - neighbors;
			minmax[i + cPos.length] = cPos[i] + neighbors;
		}
	}

	public static boolean isFar(long[] px, long[] py, long threshold){
		boolean res = false;

		for(int d = 0; d < px.length; ++d){
			if(Math.abs(px[d] - py[d]) > threshold){
				res = true; 
				break;
			}
		}

		return res;

	}

	// detects central line 
	public static <T extends Comparable<T> & RealType<T>> void detectCentralCurve(RandomAccessibleInterval<T> img, PointSampleList<T> endPoints, RandomAccessible<T> out, PointSampleList<T> points, T minValue){
		RandomAccess<T> rImg = img.randomAccess();
		RandomAccess<T> rOut = out.randomAccess();

		Cursor<T> tmpCursor = endPoints.cursor();
		tmpCursor.fwd();
		// tmpCursor.fwd();
		rImg.setPosition(tmpCursor); 

		tmpCursor.fwd();
		long [] endPoint = new long [img.numDimensions()];
		tmpCursor.localize(endPoint);

		// distance to track the progress
		// also keeps us away from the infinite loop
		long dist = 20;
		long idx = 0;

		// previous and current positions of the randomAccess
		long[] pRandomAccessPosition = new long[img.numDimensions()];
		long[] cRandomAccessPosition = new long[img.numDimensions()];
		// the interval we consider for neighbors
		long[] minmax = new long[2*img.numDimensions()];

		rImg.localize(pRandomAccessPosition);
		rImg.localize(cRandomAccessPosition);

		//TODO: fix this dirty hack to enforce the initial iteration
		pRandomAccessPosition[0]--;
		pRandomAccessPosition[1]--;

		while(isFar(endPoint, cRandomAccessPosition, dist) && idx <= 8000){
			rImg.localize(cRandomAccessPosition);

			//setMinMax(minmax, cRandomAccessPosition);
			detectDirection(pRandomAccessPosition, cRandomAccessPosition, minmax, img.numDimensions());

			Interval interval = Intervals.createMinMax(minmax);
			Cursor <T> cursor = Views.interval(img, interval).cursor();

			// DEBUG: dimensions look correct
			//			System.out.print("interval: ");
			//			for (int i = 0; i < interval.numDimensions(); ++i){
			//				System.out.print(interval.dimension(i) + " ");
			//			}
			//			System.out.println();

			//			System.out.print("cursor: ");
			//			for (int i = 0; i < img.numDimensions(); ++i){
			//				System.out.print(cursor.getLongPosition(i) + " " );
			//			}
			//			System.out.println();

			T maxValue = cursor.next().copy();
			maxValue.set(minValue);
			cursor.reset(); // new 
			// System.out.println("maxValue: " + maxValue);

			final long[] maxCoord = new long[img.numDimensions()];

			long j = 0;

			while(cursor.hasNext()){
				cursor.fwd();

				long[] cCursorPosition = new long[img.numDimensions()];


				cursor.localize(cCursorPosition);
				rImg.localize(cRandomAccessPosition);

				if (!isSame(cCursorPosition, cRandomAccessPosition, img.numDimensions())){ // skip same point 

					// System.out.println(cursor.get().getRealFloat() + " "+ maxValue);

					if(cursor.get().compareTo(maxValue) >= 0){ // not same point
						// TODO: looks redundant
						if (!isSame(cCursorPosition, pRandomAccessPosition, img.numDimensions())){ // not previous point
							rOut.setPosition(cursor);
							if (rOut.get().getRealFloat() < 1){ // not any other point
								maxValue.set(cursor.get().copy());
								cursor.localize(maxCoord);
								//System.out.println("Happens");
							}
						}

					}
				} // do not look like you understand this stuff! :o)
				else{
					// System.out.println("Same point!");
				}
				j++;

				// out of bounds
				for (int i = 0; i < img.numDimensions(); ++i){
					if (maxCoord[i] + 10 >= img.dimension(i)){
						return;
					} 
				}
			}

			//			for (int i = 0; i < img.numDimensions(); ++i){
			//				System.out.println("j: " + j );
			//			}

			//			System.out.print("maxCoord: ");
			//			for (int i = 0; i < img.numDimensions(); ++i){
			//				System.out.print(maxCoord[i] + " " );
			//			}
			//			System.out.println();

			rImg.localize(pRandomAccessPosition);
			rImg.setPosition(maxCoord);
			rOut.setPosition(maxCoord);

			// Something went wrong!
			if (maxCoord[0] == 0){
				System.out.println("wrong maxCoord");
				return;
			}

			if (rOut.get().getRealFloat() == 1.0){
				System.out.println("we are coming back :(");
				return;
			}
			// =======================

			rOut.get().setReal(idx);
			// collect every rem point
			long rem = 10;
			if (idx % rem == 0){
				points.add(new Point(rImg), rImg.get());
			}
			// 
			idx++; 

		}

	}

	public static <T extends RealType<T> & Comparable<T>> boolean isSame(long[] cValue, long[] rValue, int n){
		boolean res = true;

		for (int d = 0; d < n; ++d)
			if (cValue[d] != rValue[d]){
				res = false; 
				break;
			}
		return res;
	}

	public static <T extends RealType <T> & Comparable<T>> void processPoints(PointSampleList<T> points){

	}

	// rotate the vector given by two points
	public static <T extends RealType <T> & Comparable<T>> void rotate(PointSampleList<T> vector){
		// TODO: 
	}

	// translate the vector given by two points
	public static <T extends RealType <T> & Comparable<T>> void translate(PointSampleList<T> vector){
		// TODO: 
	}
	// shear the vector given by two points
	public static <T extends RealType <T> & Comparable<T>> void shear(PointSampleList<T> vector){
		// TODO: 
	}

	// for debugging
	public static void printArray(final long[] a){
		for (int i = 0; i < a.length; ++i)
			System.out.print(a[i] + " ");
		System.out.println();
	} 


	public static <T extends RealType<T>>void main(String[] args){
		if (showImage)
			new ImageJ();
		File file = new File("src/main/resources/distanceImg2.tif");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		if (showImage)
			ImageJFunctions.show(img);

		Img <FloatType> dst = img.factory().create(img, img.firstElement());

		PointSampleList<FloatType> endPoints = new PointSampleList<FloatType>(img.numDimensions());

		PointSampleList<FloatType> dummy = new PointSampleList<FloatType>(img.numDimensions());

		detectEnds(endPoints, dummy, dst, new FloatType(255));
		PointSampleList<FloatType> points = new PointSampleList<FloatType>(img.numDimensions());
		detectCentralCurve(img, endPoints, dst, points, new FloatType(0));
		if (showImage)
			ImageJFunctions.show(dst);

		Img <FloatType> debugPic = img.factory().create(img, img.firstElement());
		Cursor<FloatType> cursor = points.cursor();
		RandomAccess<FloatType> randomAccess = debugPic.randomAccess();
		while(cursor.hasNext()){
			cursor.fwd();
			randomAccess.setPosition(cursor);
			randomAccess.get().set(cursor.get().get());
		}

		if (showImage)
			ImageJFunctions.show(debugPic).setTitle("Every rem Point!");


		// debug 	
		int dim = 3;
		final long[] pPos = new long[]{2, 2, 2};
		final long [] cPos = new long[]{3, 1, 1}; 
		final long [] minmax = new long[2*dim];

		detectDirection(pPos, cPos, minmax, dim);
		printArray(minmax);



		System.out.println("Done!");

	}

}