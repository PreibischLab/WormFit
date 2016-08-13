package randomAlgorithms;

import java.io.File;
import java.util.ArrayList;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Superpixel {

	public static <T extends RealType<T>> void runSuperpixel(RandomAccessibleInterval<T> input, RandomAccessibleInterval<T> interval, 
			RandomAccessibleInterval<T> output,
			RandomAccessibleInterval<T> distance,
			RandomAccessibleInterval<T> labels
			){
		final RandomAccessibleInterval<T> src = Views.interval(Views.extendMirrorSingle(input), interval);
		// grid step
		long S = 20;
		
		ArrayList<Point> C = new ArrayList<>(); // contains cluster centers

		RandomAccess<T> ra = src.randomAccess();
		Cursor<T> lc = Views.iterable(src).localizingCursor();

		// initially start at {S,S,S}
		for (int d = 0; d < src.numDimensions(); ++d)
			ra.move(S,d);

		// inefficient
		while(lc.hasNext()){
			boolean hit = true;
			long[] position = new long[src.numDimensions()];
			lc.localize(position);
			for (int d = 0; d < src.numDimensions(); ++d){
				if(position[d] == 0 || position[d]%S != 0){
					hit = false;
					break;
				}
			}

			if(hit){
				C.add(new Point(lc)); // coordinates
				// @DEBUG:
				System.out.println(C.get(C.size() - 1).getIntPosition(0) + " " + C.get(C.size() - 1).getIntPosition(1));
			}

			lc.fwd();
		}
		
		// TODO: move cluster center to the lowest gradient position in 3x3 neighborhood

		for(Point center: C){
			long[] minBoundary = new long[input.numDimensions()];
			long[] maxBoundary = new long[input.numDimensions()];
			
			for (int d = 0; d < input.numDimensions(); ++d){
				minBoundary[d] = center.getLongPosition(d) - S;
				maxBoundary[d] = center.getLongPosition(d) + S;
			}
				
			Interval roi = Views.interval(src, minBoundary, maxBoundary);
			Cursor<T> localCursor = Views.interval(src, roi).cursor();
			
			// calculate distance and label 
			
		}
		
	}

	
	
	
	// black and white images
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
		File file = new File("src/main/resources/LENNA.JPG");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		
		final Img<FloatType> distance = img.factory().create(img, img.firstElement());
		final Img<FloatType> labels = img.factory().create(img, img.firstElement());
		
		initialize(distance, new FloatType(Float.MAX_VALUE));
		initialize(labels, new FloatType(-1.0f));
		
		runSuperpixel(img, img, dst, distance, labels);

		System.out.println("Doge!");
	}
}
