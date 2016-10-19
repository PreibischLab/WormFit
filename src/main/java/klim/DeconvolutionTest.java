package klim;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

// import com.sun.glass.ui.View;

import deconvolution.AdjustInput;
import deconvolution.DeconvolveTest;
import deconvolution.LucyRichardson;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.FinalRealInterval;
import net.imglib2.Interval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.RealRandomAccessibleRealInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.localization.DummySolver;
import net.imglib2.algorithm.localization.FitFunction;
import net.imglib2.algorithm.localization.Gaussian;
import net.imglib2.algorithm.localization.MLGaussianEstimator;
import net.imglib2.algorithm.localization.PeakFitter;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.InterpolatorFactory;
import net.imglib2.interpolation.neighborsearch.InverseDistanceWeightingInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.realtransform.AffineGet;
import net.imglib2.realtransform.AffineRandomAccessible;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import test.TestGauss3d;
import klim.Thresholding;
import mpicbg.imglib.algorithm.mirror.MirrorImage;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.wrapper.ImgLib2;
import mpicbg.util.RealSum;
import util.ImgLib2Util;

import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.img.Img;

public class DeconvolutionTest {

	private static final boolean debug = true;

	// calculate the size of the offset
	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long[] offset) {
		for (int d = 0; d < psf.numDimensions(); ++d)
			offset[d] = psf.dimension(d) / 2;
	}

	public static <T extends RealType<T> & NativeType<T>> long getPsf(RandomAccessibleInterval<T> img,
			RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf, T tVal) {
		int numDimensions = img.numDimensions();

		// detect high intensity pixels: candidates for bead centers
		Thresholding.threshold(img, out, tVal);

		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(
				new ArrayImgFactory<IntType>().create(out, new IntType()));
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(numDimensions);

		// TODO: Add check: if two different labels are too close that should be
		// one bead

		ObjectSegmentation.findBeads(img, labeling, beads);
		// offset for beads in all dimensions
		long[] offset = new long[numDimensions];
		getOffset(psf, offset); 

		if (debug){
			long [] min = new long [numDimensions];
			long [] max = new long [numDimensions];

			Cursor<T> beadsCursor = beads.cursor();
			while (beadsCursor.hasNext()){
				beadsCursor.fwd();
				for(int d = 0; d < numDimensions; ++d){
					min[d] = Math.max(0, beadsCursor.getLongPosition(d) - offset[d]);
					max[d] = Math.min(img.dimension(d), beadsCursor.getLongPosition(d) + offset[d]);
				}

				// ImageJFunctions.show(Views.interval(img, min, max));
			}
		}

		long numBeads = 0; // total number of beads
		numBeads = sumPsf(img, psf, beads, offset);
		return numBeads;
	}

	// check if bead is far enough from the boundaries
	public static boolean isFar(long[] min, long[] max, long[] offset, int numDimensions) {
		boolean isBroken = false;
		for (int d = 0; d < numDimensions; d++) {
			if ((max[d] - min[d]) != 2 * offset[d]) {
				isBroken = true;
				break;
			}
		}
		return isBroken;
	}

	// check if the bead is alone in the cropped image
	public static <T extends RealType<T>> boolean isAlone(long[] min, long[] max, long[] offset, long[] position,
			PointSampleList<T> beads, int numDimensions) {
		boolean isAlone = true;
		Cursor<T> jCursor = beads.cursor();

		while (jCursor.hasNext()) {
			jCursor.fwd();
			boolean isSame = true;

			for (int d = 0; d < numDimensions; ++d)
				if (jCursor.getLongPosition(d) != position[d])
					isSame = false;

			if (!isSame) { // not the same bead
				boolean isInside = true;
				for (int d = 0; d < numDimensions; d++) {
					if ((min[d] > jCursor.getLongPosition(d)) || (jCursor.getLongPosition(d) > max[d])) {
						isInside = false;
						break;
					}
				}
				if (isInside) {
					isAlone = false;
					break;
				}
			}
		}

		return isAlone;
	}

	// sum values of the psf (beads) for one image
	public static <T extends RealType<T> & NativeType<T>> long sumPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long[] offset) {
		int numDimensions = img.numDimensions();

		boolean isBroken = false; // bead initially is fine
		long numBrokenBeads = 0; // # of bad bead images
		Cursor<T> beadsCursor = beads.cursor();

		// this two store min/max for current bead window
		long[] min = new long[numDimensions];
		long[] max = new long[numDimensions];

		long[] position = new long[numDimensions];

		while (beadsCursor.hasNext()) {
			beadsCursor.fwd();
			isBroken = false;

			// setting to 0 and img.max(d) ensures the normal program flow
			for (int d = 0; d < numDimensions; d++) {
				min[d] = Math.max(beadsCursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(beadsCursor.getLongPosition(d) + offset[d], img.max(d));
			}

			// if the bead (+ boundary) fits in the image
			isBroken = isFar(min, max, offset, numDimensions);
			if (isBroken) {
				numBrokenBeads++;
			}

			beadsCursor.localize(position);
			// if the bead is alone in the cropped image
			if (!isBroken) {
				isBroken = !isAlone(min, max, offset, position, beads, numDimensions);
				if (isBroken) {
					numBrokenBeads++;
				}
			}
			if (!isBroken) {

				// TODO: Here should come the subpixel fitting routine
				// Guess the min and max can be float if we switch to subpixel accuracy
				// DESCRIPTION: 
				// you have 1 bead here with integer coordinates
				// pipe it to gaussian fitting
				// get back the position + the parameters of the psf
				// accumulateValues for the float values after that
				// continue for the next bead

				// TODO: check this line! maybe erroneous
				// RealRandomAccessible<T> realImg = Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>());


				ArrayList<double[]> parameters = new ArrayList<>(1);
				long [] beadPosition = new long[numDimensions];	
				beadsCursor.localize(beadPosition);

				double typicalSigma = 7.0;
				RandomAccessibleInterval<T> adjustedPsf = new ArrayImgFactory<T>().create(psf, Views.iterable(img).firstElement());
				fitGaussian(img, beadPosition, typicalSigma, adjustedPsf);

				// ImageJFunctions.show(Views.interval(img, min, max));
				// accumulateValues(Views.offset(Views.interval(img, min, max), min), psf);

				System.out.println(Views.iterable(adjustedPsf).size() + " " + Views.iterable(psf).size());
				ImageJFunctions.show(adjustedPsf).setTitle("bead");
				accumulateValues(adjustedPsf, psf);
			}
		}

		return (beads.size() - numBrokenBeads); // return number of clear
		// non-broken beads
	}

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

	// this function adds data to psf
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

	// MAIN FUNCTION! run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf) {
		for ( FloatType t : img )
			t.add( new FloatType( 1 ));
		AdjustInput.adjustImage(ImgLib2.wrapFloatToImgLib1(img), LucyRichardson.minValue, 1);

		System.out.println("before: " + sumIntensities(psf));
		// TODO: Maybe this function is wrong: no it is fine 
		AdjustInput.normImage(ImgLib2.wrapFloatToImgLib1(psf));
		System.out.println("after: " + sumIntensities(psf));

		ImageJFunctions.show(psf).setTitle("normalized psf");
		ImageJFunctions.show(img).setTitle("adjusted image");

		DeconvolveTest.deconvolve(DeconvolveTest.createInput((img), psf));
	}

	// this function used to fix the normalization
	public static <T extends FloatType> void testNormalization() {
		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "psi (deconvolved image)-100.tif"));

		System.out.println("before: " + sumIntensitiesInDouble(img));

		AdjustInput.normImage(ImgLib2.wrapFloatToImgLib1(img));
		// Normalize.normalize(psf, new FloatType(0.0f), new FloatType(1.0f));

		System.out.println("after: " + sumIntensitiesInDouble(img));
		ImageJFunctions.show(img);

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

		runDeconvolution(img, psf);
	}

	// opens multiple files and extracts beads from them
	public static void extractBeadsFromMultipleFiles(Img<FloatType> psf, int numImgs) {
		long totalBeads = 0;

		for (int i = 1; i <= numImgs; ++i) {
			File file = new File("../Desktop/latest_desktop/beads_to_go/beads" + i + ".tif");
			Img<FloatType> beads = ImgLib2Util.openAs32Bit(file);
			Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());
			Normalize.normalize(beads, new FloatType(0), new FloatType(255));
			float tVal = 120;
			totalBeads += getPsf(beads, out, psf, new FloatType(tVal));
		}

		System.out.println("Total number of beads found: " + totalBeads);

		getAverageValue(psf, totalBeads);
		FloatType tVal = new FloatType(22.0f);
		thresholdNoise(psf, tVal);
	}

	public static void runExtractBeads(){
		long totalBeads = 0;		

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/";

		String path = pathMac;
		// size of the beads should be the parameter
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(new long[]{25,25, 51}, new FloatType());
		Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "beads-generated.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		Normalize.normalize(beads, new FloatType(0), new FloatType(255));
		float tVal = 180;
		totalBeads += getPsf(beads, out, psf, new FloatType(tVal));


		System.out.println("Total number of beads found: " + totalBeads);

		getAverageValue(psf, totalBeads);
		thresholdNoise(psf, new FloatType(10.0f));

		ImageJFunctions.show(psf);
	}

	public static void runExtractGeneratedBeads(){
		long totalBeads = 0;		

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/";

		String path = pathMac;
		// size of the beads should be the parameter
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(new long[]{15, 15, 15}, new FloatType());
		Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "beads-generated.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));
		float tVal = 0.98f;
		totalBeads += getPsf(beads, out, psf, new FloatType(tVal));

		System.out.println("Total number of beads found: " + totalBeads);

		getAverageValue(psf, totalBeads);
		// 		thresholdNoise(psf, new FloatType(10.0f));

		ImageJFunctions.show(psf);
	}


	public static void mainDeconvolution() {

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/27_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/27_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "worm-piece.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-100.tif"));

		runDeconvolution(img, psf);
	}

	public static void mainDeconvolution2() {

		String pathUbuntu = "/home/milkyklim/Desktop/input/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "sample-t=1.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "psf-averaged-scaled.tif"));

		runDeconvolution(img, psf);
	}

	// test the total intensity
	public static <T extends RealType<T>> void testTotalIntensity(RandomAccessibleInterval<T> input,
			RandomAccessibleInterval<T> output, RandomAccessibleInterval<T> psf) {
		T iTotal = input.randomAccess().get().copy();
		T oTotal = output.randomAccess().get().copy();

		iTotal.set(sumIntensitiesInDouble(input));
		oTotal.set(sumIntensitiesInDouble(cropImage(output, psf)));

		System.out.println("input = " + iTotal.getRealFloat());
		System.out.println("output = " + oTotal.getRealFloat());

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

	// calculates the sum over all
	// TODO: FIX OVERFLOW
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

	public static void runTestTotalIntensity() {
		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathUbuntu;

		Img<FloatType> input = ImgLib2Util.openAs32Bit(new File(path + "DOTs-79-82.tif"));
		Img<FloatType> output = ImgLib2Util.openAs32Bit(new File(path + "it=10-1.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-100.tif"));

		testTotalIntensity(input, output, psf);
	}

	// this function to try out gaussian fitting 
	public static <T extends RealType<T>> void gaussianFitting(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out,  T tVal, Collection<double[]> results){

		int numDimensions = img.numDimensions();		
		// detect high intensity pixels: candidates for bead centers
		Thresholding.threshold(img, out, tVal);
		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType()));
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(numDimensions);
		// ImageJFunctions.show(img);
		// ImageJFunctions.show(out);


		ObjectSegmentation.findBeads(img, labeling, beads); 

		// System.out.println(beads.size());

		ArrayList<Localizable> peaks = new ArrayList<Localizable>(1);

		long [] position = new long[numDimensions];

		Cursor<T> cursorBeads = beads.cursor();

		while(cursorBeads.hasNext()){
			cursorBeads.fwd();
			cursorBeads.localize(position);
			peaks.add(new Point(position));
		}

		// beads.cursor();
		// part above copied for test 
		// why the fuck this error is here?!
		// TODO: This number below is super important and should be set as the parameter
		// Otherwise fitting won't work
		PeakFitter<T> pf = new PeakFitter(img, peaks, new DummySolver(), new Gaussian(), new MLGaussianEstimator(7.0, numDimensions));
		pf.process();

		// @DEBUG: 
		// looks like parameters are returned in the following order
		// x1, x2, x3, A, b
		// @DEBUG: 
		// maybe you need all elements that belong to the beads, not only central one
		// that means you have to use water-shadding instead of thresholding
		for (double[] element : pf.getResult().values()){
			// System.out.println("At least one element");
			for (int i = 0; i < element.length; ++i){
				System.out.println(i + " : " + element[i]);
			}
		}

		results.addAll(pf.getResult().values());

		// TODO: Remember that the fifth element is b = 1/(2 sigma^2) pf.getResult().values()


		System.out.println(pf.toString());
	}


	public static <T extends RealType<T>> void fitGaussian(RandomAccessibleInterval<T> img, long[] beadPosition, /*ArrayList<Localizable> peaks,*/ /*Collection<double[]> results,*/ double typicalSigma, RandomAccessibleInterval<T> psf){
		// there is always only one element in 
		int numDimensions = img.numDimensions();
		ArrayList<Localizable> peaks = new ArrayList<Localizable>(1);
		peaks.add(new Point(beadPosition));

		PeakFitter<T> pf = new PeakFitter<T>(img, peaks, new DummySolver(), new Gaussian(), new MLGaussianEstimator(typicalSigma, numDimensions));
		pf.process();

		// print out parameters
		//		for (double[] element : pf.getResult().values()){
		//			// System.out.println("At least one element");
		//			for (int i = 0; i < element.length; ++i){
		//				System.out.println(i + " : " + element[i]);
		//			}
		//		}

		// check the results 
		// System.out.println(pf.toString());

		// new part // 		
		double [] element = pf.getResult().values().iterator().next(); // there is only one element in this collection

		long [] min = new long[numDimensions];
		long [] max = new long[numDimensions];

		// here comes the interpolation
		// RealRandomAccessible<T> realImg = Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>());
		// setMinMax(min, max, element);
		// FinalInterval interval = new FinalInterval( min, max );

		// RandomAccessibleInterval <T> adjustedImg = Views.interval(Views.raster(realImg), interval);
		// psf = Views.interval(Views.raster(Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>())), interval);
//		copyValues(Views.offset(
//				Views.interval(
//						Views.raster(
//								Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>())
//								),
//						interval), 
//				min), 
//				psf);
		// ImageJFunctions.show(realImgToInt);
		// ImageJFunctions.show(psf).setTitle("psf");
		
		
		// TODO: incorporate function from below instead of one above
		double [] realMin = new double[numDimensions];
		double [] realMax = new double[numDimensions];
		setRealMinMax(realMin, realMax, element);
		FinalRealInterval realInterval = new FinalRealInterval( realMin, realMax );
		recalculateCoordinates(Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>()), psf, realInterval);

	}
	
	public static void setRealMinMax( double [] min, double [] max, double[] parameters){
		int numDimensions = min.length;
		// adjust the last parameter so that it is the radius now
		// b = 1/(2 sigma^2) 
		double b = Math.sqrt(1/(2*parameters[parameters.length - 1])); // parameters.get(parameters.indexOf(element))[element.length - 1];
		for (int d = 0; d < numDimensions; ++d){
			min[d] = parameters[d] - b;
			max[d] = parameters[d] + b;
			// System.out.println(min[d] + " : " + (parameters[d]) + " : " + max[d]);
			// System.out.println(min[d] + " : " + max[d]);
			if ((max[d] - min[d]) % 2 == 1)
				System.out.println("for coordinate d = " + d + ", the psf size is even!");
		}
	}

	// copies data from float grid to integer one
	public static <T extends RealType<T>> void recalculateCoordinates(RealRandomAccessible<T> realImg, RandomAccessibleInterval<T> img, RealInterval interval){
		int numDimensions = realImg.numDimensions();
		long[] pixelSize = new long[ numDimensions ];
		double[] intervalSize = new double[numDimensions];

		double magnification = 1.0;

		for (int d = 0; d < numDimensions; ++d){
			intervalSize[d] = interval.realMax(d) - interval.realMin(d);
			pixelSize[d] = Math.round(intervalSize[d]*magnification) + 1;
		}

		Cursor<T> cursor = Views.iterable(img).localizingCursor();
		RealRandomAccess<T> realRandomAccess = realImg.realRandomAccess();

		double[] position = new double[ numDimensions ];

		while(cursor.hasNext()){
			cursor.fwd();
			
			for (int d = 0; d < numDimensions; ++d){
				position[d] = cursor.getDoublePosition(d)/img.realMax(d)*intervalSize[d] + interval.realMin(d); 
			}
			
			realRandomAccess.setPosition(position);			
			cursor.get().set(realRandomAccess.get());
		}	
		
	}

	public static void runGaussianFitting(){
		String pathMac = "/Users/kkolyva/Desktop/28_09_16_severin/";

		String path = pathMac;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "bead.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(img, new BitType());
		Normalize.normalize(img, new FloatType(0), new FloatType(255));
		float tVal = 180;

		thresholdNoise(img, new FloatType(10.0f));

		gaussianFitting(img, out, new FloatType(tVal));

		// ImageJFunctions.show(psf);
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

	public static void mainDeconvolutionSliced() {

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathMac;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "worm-piece.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-75.tif"));
		setSliceValue(img, new FloatType(-1));

		runDeconvolution(img, psf);
	}

	// adds gaussian peaks at the specified location with the specified sigma 
	public static void testGaussianFitting(){

		final double[] location = new double[]{50.2, 50, 53.3};
		final double[] sigma = new double[]{7, 5, 9};
		final long [] dimensions = new long[] {300, 300, 300};
		int numDimensions = dimensions.length;

		Img<FloatType> img = new ArrayImgFactory<FloatType>().create(dimensions, new FloatType());
		TestGauss3d.addGaussian(img, location, sigma);

		TestGauss3d.addGaussian(img, new double[]{83, 22.4, 42.5}, sigma);
		TestGauss3d.addGaussian(img, new double[]{200.2, 74.5, 99.48}, sigma);

		ImageJFunctions.show(img);
		Img<BitType> bitImg = new ArrayImgFactory<BitType>().create(dimensions, new BitType());

		ArrayList<double[]> parameters = new ArrayList<>(1);
		gaussianFitting(img, bitImg, new FloatType(0.98f), parameters);

		for (double[] element : parameters){
			long [] min = new long[numDimensions];
			long [] max = new long[numDimensions];

			// here comes the interpolation check
			RealRandomAccessible<FloatType> realImg = Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<FloatType>());
			setMinMax(min, max, element);
			FinalInterval interval = new FinalInterval( min, max );

			RandomAccessibleInterval <FloatType> realImgToInt = Views.interval(Views.raster(realImg), interval);
			ImageJFunctions.show(realImgToInt);
		}

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


	public static void main(String[] args) {
		new ImageJ();

		// testNormalization();
		// test3D();
		// runTestTotalIntensity();
		// mainDeconvolution2();
		runExtractGeneratedBeads();
		// runGaussianFitting();
		// mainDeconvolutionSliced();
		// testGaussianFitting();
		System.out.println("Doge!");

	}
}
