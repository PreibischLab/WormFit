package klim.deconvolution;

import java.io.File;
import java.util.ArrayList;

import deconvolution.AdjustInput;
import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.Util;
import net.imglib2.view.Views;
import util.ImgLib2Util;
//used for debug outputs 
import static klim.deconvolution.DeconvolutionTest.debug;

public class Tests {

	// adds gaussian peaks at the specified location with the specified sigma 
	public static void testGaussianFitting(){

		final double[] location = new double[]{50.2, 50, 53.3};
		final double[] sigma = new double[]{7, 5, 9};
		final long [] dimensions = new long[] {300, 300, 300};
		int numDimensions = dimensions.length;

		Img<FloatType> img = new ArrayImgFactory<FloatType>().create(dimensions, new FloatType());
		addGaussian(img, location, sigma);

		addGaussian(img, new double[]{83, 22.4, 42.5}, sigma);
		addGaussian(img, new double[]{200.2, 74.5, 99.48}, sigma);

		ImageJFunctions.show(img);
		Img<BitType> bitImg = new ArrayImgFactory<BitType>().create(dimensions, new BitType());

		ArrayList<double[]> parameters = new ArrayList<>(1);
		DeconvolutionTest.gaussianFitting(img, bitImg, new FloatType(0.98f), parameters);

		for (double[] element : parameters){
			long [] min = new long[numDimensions];
			long [] max = new long[numDimensions];

			// here comes the interpolation check
			RealRandomAccessible<FloatType> realImg = Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<FloatType>());
			// TODO: check if this one is necessary
			Utils.setMinMax(min, max, element);
			FinalInterval interval = new FinalInterval( min, max );

			RandomAccessibleInterval <FloatType> realImgToInt = Views.interval(Views.raster(realImg), interval);
			ImageJFunctions.show(realImgToInt);
		}

	}
	
	// this function used to fix the normalization
	public static <T extends FloatType> void testNormalization() {
		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "psi (deconvolved image)-100.tif"));

		System.out.println("before: " + Utils.sumIntensitiesInDouble(img));

		AdjustInput.normImage(ImgLib2.wrapFloatToImgLib1(img));
		// Normalize.normalize(psf, new FloatType(0.0f), new FloatType(1.0f));

		System.out.println("after: " + Utils.sumIntensitiesInDouble(img));
		ImageJFunctions.show(img);

	}
	
	
	public static void runTestTotalIntensity() {
		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> input = ImgLib2Util.openAs32Bit(new File(path + "DOTs-79-82.tif"));
		Img<FloatType> output = ImgLib2Util.openAs32Bit(new File(path + "it=10-1.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-100.tif"));

		testTotalIntensity(input, output, psf);
	}
	
	// test the total intensity
	public static <T extends RealType<T>> void testTotalIntensity(RandomAccessibleInterval<T> input,
			RandomAccessibleInterval<T> output, RandomAccessibleInterval<T> psf) {
		T iTotal = input.randomAccess().get().copy();
		T oTotal = output.randomAccess().get().copy();
	
		iTotal.setReal(Utils.sumIntensitiesInDouble(input));
		oTotal.setReal(Utils.sumIntensitiesInDouble(Utils.cropImage(output, psf)));

		System.out.println("input = " + iTotal.getRealFloat());
		System.out.println("output = " + oTotal.getRealFloat());

	}
	
	// run the test case in 3D
	public static void test3D() {

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "DOTs-79-82.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-50.tif"));

		Img<FloatType> convImg = img.copy();

		FFTConvolution<FloatType> conv = new FFTConvolution<FloatType>(convImg, psf);
		conv.convolve();

		ImageJFunctions.show(img).setTitle("initial image");
		ImageJFunctions.show(convImg).setTitle("convolved image");

		Run.runDeconvolution(img, psf);
	}
	
	
	final public static void addGaussian( final Img< FloatType > image, final double[] location, final double[] sigma )
	{
		final int numDimensions = image.numDimensions();
		final int[] size = new int[ numDimensions ];
		
		final long[] min = new long[ numDimensions ];
		final long[] max = new long[ numDimensions ];
		
		final double[] two_sq_sigma = new double[ numDimensions ];
		
		for ( int d = 0; d < numDimensions; ++d )
		{
			size[ d ] = Util.getSuggestedKernelDiameter( sigma[ d ] ) * 2;
			min[ d ] = (int)Math.round( location[ d ] ) - size[ d ]/2;
			max[ d ] = min[ d ] + size[ d ] - 1;
			two_sq_sigma[ d ] = 2 * sigma[ d ] * sigma[ d ];
		}

		final RandomAccessible< FloatType > infinite = Views.extendZero( image );
		final RandomAccessibleInterval< FloatType > interval = Views.interval( infinite, min, max );
		final IterableInterval< FloatType > iterable = Views.iterable( interval );
		final Cursor< FloatType > cursor = iterable.localizingCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			
			double value = 1;
			
			for ( int d = 0; d < numDimensions; ++d )
			{
				final double x = location[ d ] - cursor.getIntPosition( d );
				value *= Math.exp( -(x * x) / two_sq_sigma[ d ] );
			}
			
			cursor.get().set( cursor.get().get() + (float)value );
		}
	}
	
}
