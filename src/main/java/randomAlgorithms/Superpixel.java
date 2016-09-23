package randomAlgorithms;

import java.io.File;
import java.util.ArrayList;

import com.sun.j3d.utils.universe.ViewerAvatar;

import ij.ImageJ;
import imglib.ops.operator.binary.Max;
import klim.BoundingBox;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.display.projector.RandomAccessibleProjector2D;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelingMapping;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Superpixel {

	public static <T extends RealType<T>> void initializeClusters(
			RandomAccessibleInterval<T> src, 
			RandomAccessibleInterval<T> interval, 
			RandomAccessibleInterval<T> output,
			RandomAccessibleInterval<T> distance,
			RandomAccessibleInterval<T> labels,
			RandomAccessibleInterval<T> tmp, 
			RandomAccessibleInterval<T> kernel,

			ArrayList<Point> clusters
			){

		long S = 20;
		double m = 0.1; // TODO: parameter!

		int numDimensions = src.numDimensions();

		RandomAccess<T> raSrc = src.randomAccess();
		Cursor<T> lc = Views.iterable(src).localizingCursor();

		RandomAccess<T> raDistance = distance.randomAccess();
		RandomAccess<T> raLabels = labels.randomAccess();

		// start at (S,S,S)
		for (int d = 0; d < numDimensions; ++d)
			raSrc.move(S,d);

		ImageJFunctions.show(src);

		// inefficient
		// TODO: change to random access
		while(lc.hasNext()){
			boolean hit = true;

			long[] position = new long[numDimensions];
			lc.localize(position);
			for (int d = 0; d < src.numDimensions(); ++d){
				if(position[d] == 0 || position[d]%S != 0){
					hit = false;
					break;
				}
			} 

			if(hit){
				Point center = new Point(lc);

				long[] minB = new long[numDimensions];
				long[] maxB = new long[numDimensions];
				long neighborhood = 2;

				for (int d = 0; d < numDimensions; ++d){
					minB[d] = center.getLongPosition(d) - neighborhood;
					maxB[d] = center.getLongPosition(d) + neighborhood;

					// TODO: is it possible to adjust this part using views
					//							if (minB[d] < 0) 
					//								minB[d] = 0;
					//							if (maxB[d] >= input.dimension(d)) 
					//								maxB[d] = input.dimension(d) - 1;
				}

				badCopy(Views.interval(src, minB, maxB), tmp); // copied 5x5 block			
				new FFTConvolution<T>(tmp, kernel, new ArrayImgFactory<ComplexFloatType>()).convolve();

				// TODO: change for 3D
				// System.out.println(center. + " " + tmp.dimension(1));
				findMin(tmp, center, minB, maxB);
				// System.out.println(center.toString());
				// find minimum here
				// add minimum to the C

				clusters.add(center); // coordinates
				// @DEBUG:
				// System.out.println(C.get(C.size() - 1).getIntPosition(0) + " " + C.get(C.size() - 1).getIntPosition(1));
			}
			lc.fwd();

		}

	}


	public static <T extends RealType<T>> void assignmentStep(
			RandomAccessibleInterval<T> src, 
			RandomAccessibleInterval<T> interval, 
			RandomAccessibleInterval<T> output,
			RandomAccessibleInterval<T> distance,
			RandomAccessibleInterval<T> labels,
			RandomAccessibleInterval<T> tmp, 
			RandomAccessibleInterval<T> kernel,


			ArrayList<Point> clusters
			){

		long S = 50;
		double m = 0.1; // TODO: parameter!

		int numDimensions = src.numDimensions();


		RandomAccess<T> ra = src.randomAccess();
		Cursor<T> lc = Views.iterable(src).localizingCursor();

		RandomAccess<T> raDistance = distance.randomAccess();
		RandomAccess<T> raLabels = labels.randomAccess();



		for(Point center: clusters){
			long[] minBoundary = new long[numDimensions];
			long[] maxBoundary = new long[numDimensions];

			// TODO: out of bounds check!
			for (int d = 0; d < numDimensions; ++d){
				minBoundary[d] = center.getLongPosition(d) - S;
				maxBoundary[d] = center.getLongPosition(d) + S;

				// TODO: is it possible to adjust this part using views
				if (minBoundary[d] < 0) 
					minBoundary[d] = 0;
				if (maxBoundary[d] >= interval.dimension(d)) 
					maxBoundary[d] = interval.dimension(d) - 1;
			}

			// System.out.println("[" + minBoundary[0] + " " + maxBoundary[0] + "]," + " [" +minBoundary[1] + " " + maxBoundary[1] + "]");


			Interval roi = Views.interval(src, minBoundary, maxBoundary);
			Cursor<T> localCursor = Views.interval(src, roi).cursor();

			// calculate distance and label 

			while(localCursor.hasNext()){
				localCursor.fwd();
				//System.out.println(localCursor.getIntPosition(0));

				raDistance.setPosition(localCursor);
				raLabels.setPosition(localCursor);

				ra.setPosition(center);

				long[] clusterCenter = new long[numDimensions];
				long[] currentPixel = new long[numDimensions];

				ra.localize(currentPixel);
				center.localize(clusterCenter);

				double dc = dc(ra.get().getRealDouble(), localCursor.get().getRealDouble());
				double ds = ds(clusterCenter, currentPixel);
				double D = D(ds, dc, S, m);

				if (D < raDistance.get().getRealDouble()){
					raDistance.get().setReal(D);
					raLabels.get().setReal(clusters.indexOf(center));
				}

			}

		}

	}

	public static <T extends RealType<T>> void runSuperpixel(RandomAccessibleInterval<T> input, RandomAccessibleInterval<T> interval, 
			RandomAccessibleInterval<T> output,
			RandomAccessibleInterval<T> distance,
			RandomAccessibleInterval<T> labels,
			RandomAccessibleInterval<T> tmp, 
			RandomAccessibleInterval<T> kernel
			){
		final RandomAccessibleInterval<T> src = Views.interval(Views.extendMirrorSingle(input), interval);
		// grid step

		long S = 20;
		double m = 0.1; // TODO: parameter!

		int numDimensions = input.numDimensions();

		ArrayList<Point> clusters = new ArrayList<>(); // contains cluster centers
		initializeClusters(src, interval, output, distance, labels, tmp, kernel, clusters);

		// TODO: move cluster center to the lowest gradient position in 3x3 neighborhood
		// set the boundaries for the neighborhood

		// tmp = Views.interval(input, new long[]{}, new long[]{});
		// new FFTConvolution<T>(tmp, kernel, new ArrayImgFactory<ComplexFloatType>()).convolve();
		// // rotate the kernel
		// new FFTConvolution<T>(tmp, kernel, new ArrayImgFactory<ComplexFloatType>()).convolve();

		ImageJFunctions.show(src);

		RandomAccess<T> ra = src.randomAccess();
		Cursor<T> lc = Views.iterable(src).localizingCursor();

		RandomAccess<T> raDistance = distance.randomAccess();
		RandomAccess<T> raLabels = labels.randomAccess();
		

		for (int it = 0; it < 100; ++ it){
			assignmentStep(src, interval, output, distance, labels, tmp, kernel, clusters);

			// TODO: here should come cluster recalculation
			// you should take the mean of all point in the cluster
			// unfortunatly to make it in efficient way you have to use 
			// final ImgLabeling<Integer, ?> labeling = new ImgLabeling<?>(labels);
			// final LabelingMapping<T> labeling = new LabelingMapping<>(new IntType());	
			// final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(labels, new IntType())); 

			
			ArrayList<PointSampleList<T>> objects = new ArrayList<>();
			setPointSampleList(labels, objects);
			
			System.out.print(clusters.size());
			System.out.println(" ?= " + objects.size());

			clusters.clear();

			
			// if this one is constant that is a good sign 
			// System.out.println("# of clusters: " + objects.size());

			// @Debug: print out the size of the clusters
			//			for (int j = 0; j < objects.size(); ++j){
			//				if (j % 50 == 0){
			//					 System.out.println("Size of the cluster: " + objects.get(j).size());
			//				}
			//					
			//			}

			for(PointSampleList<T> one : objects){
				long[] pos = new long[]{0,0};
				// iterate over points in cluster
				Cursor<T> cursor = one.cursor();
				while (cursor.hasNext()){
					cursor.fwd();
					for (int d = 0; d < numDimensions; ++d)
						pos[d] += cursor.getLongPosition(d);			
				}

				// System.out.println(objects.indexOf(one));
				for (int d = 0; d < numDimensions; ++d)
					pos[d] /= one.size();

				// System.out.println(pos[0] + " " + pos[1]);

				clusters.add(new Point(pos));
			}

			int index = 2;
			System.out.println(clusters.get(index).getLongPosition(0) + " " + clusters.get(index).getLongPosition(1));
			// System.out.println(C.size());
			// BoundingBox.setLabeling(thresholdImg, labeling);
		}

	}

	public static <T extends RealType<T>> void setPointSampleList(final RandomAccessibleInterval<T> labels, ArrayList<PointSampleList<T>> objects){
		// this one bellow is the final list you'll need 
		// ArrayList<PointSampleList<T>> objectsList = new ArrayList<>();

		Cursor<T> cursor = Views.iterable(labels).cursor();
		//RandomAccess <T> randomAccess = img.randomAccess();

		// number of objects in the picture
		// calculated dynamically
		int curMax = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			int curElement = (int)cursor.get().getRealFloat();
			// increase the size of the objects list
			if (curElement >= curMax){
				while(curElement >= curMax){
					objects.add(new PointSampleList<T>(labels.numDimensions()));
					++curMax; 
				}
			}
			try{
				objects.get(curElement).add(new Point(cursor), cursor.get().copy()); // .copy here ?
			}
			catch(Exception e){
				System.out.println(objects.get(curElement).size());
			}
		}

	}


	public static <T extends RealType<T> & Comparable<T>> void findMin(RandomAccessibleInterval<T> src, Point min, long[] minB, long[] maxB){
		Cursor<T> c = Views.iterable(Views.interval(src, new long[]{1,1}, new long[]{3,3})).cursor();
		T minVal = c.get();

		while(c.hasNext()){
			c.fwd();
			if (c.get().compareTo(minVal) > 0){
				minVal.set(c.get());
				min.setPosition(c);

				for (int d = 0; d < src.numDimensions(); ++d){
					min.move(minB[d] + 1, d);
				}
			}
		}
	}

	public static  <T extends RealType<T>> void badCopy(RandomAccessibleInterval<T> from, RandomAccessibleInterval<T> to){
		Cursor<T> cursor = Views.iterable(to).cursor();
		RandomAccess<T> randomAccess = from.randomAccess();

		while(cursor.hasNext()){
			cursor.fwd();
			randomAccess.setPosition(cursor);
			cursor.get().set(randomAccess.get());
		}
	}


	// black and white images
	// compare intensities  
	public static double dc(double x, double y){
		return Math.abs(x - y);
	}

	public static double ds(long[] x, long[] y){
		double res = 0;

		for (int i = 0; i < x.length; ++i)
			res += (x[i] - y[i])*(x[i] - y[i]);

		res = Math.sqrt(res);
		return res;
	}

	// m parameter
	public static double D(double ds, double dc, double S, double m){
		return Math.sqrt(dc*dc + m*m*(ds/S)*(ds/S));
	}

	// will be used to initialize labels and distances
	public static <T extends RealType<T>> void initialize(RandomAccessibleInterval<T> img, T val){
		Cursor<T> dc = Views.iterable(img).cursor();	
		while(dc.hasNext()){
			dc.fwd();
			dc.get().set(val);
		}

	} 

	public static void main(String[] args){
		String path = "src/main/resources/LENNA.JPG";
		File file = new File(path);
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());

		final Img<FloatType> distance = img.factory().create(img, img.firstElement());
		final Img<FloatType> labels = img.factory().create(img, img.firstElement());

		Img<FloatType> tmp = img.factory().create(new long[]{5, 5}, img.firstElement()); 
		Img<FloatType> kernel = ArrayImgs.floats( new float[]{-1,0,1}, new long[]{1, 3});


		initialize(distance, new FloatType(Float.MAX_VALUE));
		initialize(labels, new FloatType(-1.0f));

		runSuperpixel(img, img, dst, distance, labels, tmp, kernel);

		new ImageJ();
		ImageJFunctions.show(distance);
		ImageJFunctions.show(labels);


		System.out.println("Doge!");
	}
}
