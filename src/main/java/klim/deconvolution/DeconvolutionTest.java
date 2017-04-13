package klim.deconvolution;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;

import ij.ImageJ;
import klim.ObjectSegmentation;
import klim.Thresholding;
import net.imglib2.Cursor;
import net.imglib2.FinalRealInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealInterval;
import net.imglib2.RealRandomAccess;
import net.imglib2.RealRandomAccessible;
import net.imglib2.algorithm.localization.DummySolver;
import net.imglib2.algorithm.localization.Gaussian;
import net.imglib2.algorithm.localization.MLGaussianEstimator;
import net.imglib2.algorithm.localization.PeakFitter;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class DeconvolutionTest {

	protected static final boolean debug = true;

	// calculate the size of the offset
	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long[] offset) {
		for (int d = 0; d < psf.numDimensions(); ++d)
			offset[d] = psf.dimension(d) / 2;
	}

	public static <T extends RealType<T> & NativeType<T>> long getPsf(RandomAccessibleInterval<T> img,
			RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf, T tVal, double typicalSigma) {
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

		if (debug){ // show each bead
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
		numBeads = sumPsf(img, psf, beads, offset, typicalSigma);
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

	// TODO: keep it here ?
	// sum values of the psf (beads) for one image
	public static <T extends RealType<T> & NativeType<T>> long sumPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long[] offset, double typicalSigma) {
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
				ArrayList<double[]> parameters = new ArrayList<>(1);
				long [] beadPosition = new long[numDimensions];	
				beadsCursor.localize(beadPosition);

				RandomAccessibleInterval<T> adjustedPsf = new ArrayImgFactory<T>().create(psf, Views.iterable(img).firstElement());
				fitGaussian(img, beadPosition, typicalSigma, adjustedPsf);

				//				System.out.println(Views.iterable(adjustedPsf).size() + " " + Views.iterable(psf).size());
				ImageJFunctions.show(adjustedPsf).setTitle("bead");
				Utils.accumulateValues(adjustedPsf, psf);
			}
		}

		return (beads.size() - numBrokenBeads); // return number of clear non-broken beads
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
			totalBeads += getPsf(beads, out, psf, new FloatType(tVal),  (long)(psf.dimension(0)/2));
		}

		System.out.println("Total number of beads found: " + totalBeads);

		Utils.getAverageValue(psf, totalBeads);
		FloatType tVal = new FloatType(22.0f);
		Utils.thresholdNoise(psf, tVal);
	}

	public static void mainDeconvolution() {

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/27_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/";

		String path = pathUbuntu;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "worm-piece.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-100.tif"));

		Run.runDeconvolution(img, psf);
	}

	// TODO: test this one; But I think it was working
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
		PeakFitter<T> pf = new PeakFitter<>(img, peaks, new DummySolver(), new Gaussian(), new MLGaussianEstimator(7.0, numDimensions));
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
		int numDimensions = img.numDimensions();
		// there is always only one element in the list
		ArrayList<Localizable> peaks = new ArrayList<Localizable>(1);
		peaks.add(new Point(beadPosition));

		PeakFitter<T> pf = new PeakFitter<T>(img, peaks, new DummySolver(), new Gaussian(), new MLGaussianEstimator(typicalSigma, numDimensions));
		pf.process();

		if (debug){
			// print out parameters
			for (double[] element : pf.getResult().values()){
				for (int i = 0; i < element.length; ++i){
					System.out.println("parameter[" + i + "] : " + element[i]);
				}
			}
			// check the results 
			System.out.println(pf.toString());
		}

		// new part // 		
		double [] element = pf.getResult().values().iterator().next(); // there is only one element in this collection		

		double [] realMin = new double[numDimensions];
		double [] realMax = new double[numDimensions];
		Utils.setRealMinMax(realMin, realMax, element);
		FinalRealInterval realInterval = new FinalRealInterval( realMin, realMax );
		recalculateCoordinates(Views.interpolate( Views.extendMirrorSingle(img), new NLinearInterpolatorFactory<T>()), psf, realInterval);

	}

	// copies data from float grid to integer one
	public static <T extends RealType<T>> void recalculateCoordinates(RealRandomAccessible<T> realImg, RandomAccessibleInterval<T> img, RealInterval interval){
		int numDimensions = realImg.numDimensions();
		long[] pixelSize = new long[ numDimensions ];
		double[] intervalSize = new double[numDimensions];

		double magnification = 1.0; // this one can be deleted later

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

	// run for every second slice set to specific valuesdd
	public static void mainDeconvolutionSliced() {

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/20_09_16_psf_results/";

		String path = pathMac;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "worm-piece.tif"));
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "PSF-done-75.tif"));
		Utils.setSliceValue(img, new FloatType(-1));

		Run.runDeconvolution(img, psf);
	}



	public static void main(String[] args) {
		new ImageJ();

		// testNormalization();
		// test3D();
		// runTestTotalIntensity();
		// mainDeconvolution();
		// runExtractGeneratedBeads();
		Run.runExtractBeads();
		// runGaussianFitting();
		// mainDeconvolutionSliced();
		// testGaussianFitting();
		System.out.println("Doge!");

	}
}
