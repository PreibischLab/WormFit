package klim.deconvolution;

import mpicbg.util.RealSum;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
// used for debug outputs 
import static klim.deconvolution.DeconvolutionTest.debug;


public class Utils {
	// this function adds data to psf
	public static <T extends RealType<T>> void accumulateValues(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf) {
		Cursor<T> c = Views.iterable(img).cursor();
		RandomAccess<T> r = psf.randomAccess();
		while (c.hasNext()) {
			c.fwd();
			r.setPosition(c);
			r.get().add(c.get());
		}
	}

	// this function copies data to psf
	public static <T extends RealType<T>> void copyValues(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf) {
		Cursor<T> c = Views.iterable(img).cursor();
		RandomAccess<T> r = psf.randomAccess();
		while (c.hasNext()) {
			c.fwd();
			r.setPosition(c);
			r.get().set(c.get());
		}
	}

	// averages psf
	public static <T extends RealType<T>> void getAverageValue(RandomAccessibleInterval<T> psf, long total) {
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while (viewCursor.hasNext()) {
			viewCursor.fwd();
			viewCursor.get().setReal(viewCursor.get().getRealDouble() / total);
		}
	}

	// set all pixels to 0
	public static <T extends RealType<T>> void setZero(RandomAccessibleInterval<T> psf) {
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while (viewCursor.hasNext()) {
			viewCursor.fwd();
			viewCursor.get().setZero();
		}
	}
	
	// subtract noise (T val)
	public static <T extends RealType<T> & Comparable<T>> void thresholdNoise(RandomAccessibleInterval<T> psf, T tVal) {
		System.out.println("Noise threshold-value: " + tVal);

		Cursor<T> cursor = Views.iterable(psf).cursor();
		while (cursor.hasNext()) {
			cursor.fwd();
			cursor.get().setReal(Math.max(0, cursor.get().getRealDouble() - tVal.getRealDouble()));
		}
	}
	
	// calculates the sum over all
	// TODO: FIX OVERFLOW
	// TODO: DELETE WITH NEXT COMMIT
	public static <T extends RealType<T>> T sumIntensities(RandomAccessibleInterval<T> img) {
		T sum = img.randomAccess().get().copy();
		sum.setZero();
		Cursor<T> cursor = Views.iterable(img).cursor();

		while (cursor.hasNext()) {
			cursor.fwd();
			sum.add(cursor.get());
		}

		return sum;
	}

	// calculates the sum over all
	public static <T extends RealType<T>> double sumIntensitiesInDouble(RandomAccessibleInterval<T> img) {
		RealSum sumR = new RealSum();

		Cursor<T> cursor = Views.iterable(img).cursor();

		while (cursor.hasNext()) {
			cursor.fwd();
			sumR.add( cursor.get().getRealDouble() );
		}

		return sumR.getSum();
	}
	
	// drop boundaries
	public static <T extends RealType<T>> IntervalView<T> cropImage(RandomAccessibleInterval<T> img,
			RandomAccessibleInterval<T> psf) {

		long[] min = new long[img.numDimensions()];
		long[] max = new long[img.numDimensions()];

		// why is it working ?
		for (int d = 0; d < img.numDimensions(); ++d) {
			min[d] = psf.dimension(d) / 2;
			max[d] = img.dimension(d) - psf.dimension(d) / 2 - 1;
		}

		return Views.interval(img, min, max);
	}
	
	// set the minimal and maximal 
	// TODO: this function is full of hacks and tricks to adjust the size of the interval
	// TODO: DELETE WOTH THIE NEXT COMMIT  
	public static void setMinMax(long[] min, long[] max, double[] parameters){
		int numDimensions = min.length;

		// adjust the last parameter so that it is the radius now
		// b = 1/(2 sigma^2) 
		long b = Math.round(Math.sqrt(1/(2*parameters[parameters.length - 1]))); // parameters.get(parameters.indexOf(element))[element.length - 1];
		// System.out.println("b = " + b);
		for (int d = 0; d < numDimensions; ++d){
			// TODO: the peak intensity should be in the center of the picture !
			min[d] = (long)Math.floor(parameters[d] - 1 - b);
			max[d] = (long)Math.ceil(parameters[d] + b);
			System.out.println(min[d] + " : " + (parameters[d]) + " : " + max[d]);
			// System.out.println(min[d] + " : " + max[d]);
			if ((max[d] - min[d]) % 2 == 1)
				System.out.println("for coordinate d = " + d + ", the psf size is even!");
		}
	}
	
	public static void setRealMinMax( double [] min, double [] max, double[] parameters){
		int numDimensions = min.length;
		// b = 1/(2 sigma^2) 
		double b = Math.sqrt(1/(2*parameters[parameters.length - 1]));
		if (debug)
			System.out.println("Radius = " + b);
		for (int d = 0; d < numDimensions; ++d){
			min[d] = parameters[d] - b;
			max[d] = parameters[d] + b;
			if (debug)
				System.out.println(min[d] + " : " + (parameters[d]) + " : " + max[d]);
			if ((max[d] - min[d]) % 2 == 1)
				System.out.println("For coordinate d = " + d + ", the psf size is even!");
		}
	}
	
	public static <T extends RealType<T>> void setSliceValue(RandomAccessibleInterval<T> img, T val){
		Cursor <T> cursor  = Views.iterable(img).localizingCursor();
		long positionZ = 0; // indicates if the slice index is odd or even

		int lastDim = img.numDimensions() - 1;

		while(cursor.hasNext()){
			cursor.fwd();
			positionZ = cursor.getLongPosition(lastDim);

			if (positionZ % 2 == 1){
				cursor.get().set(val);
			}
		}	
		ImageJFunctions.show(img);
	}
	
}
