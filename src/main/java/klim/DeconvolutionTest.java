package klim;

import java.io.File;
import java.util.ArrayList;

// import com.sun.glass.ui.View;

import deconvolution.AdjustInput;
import deconvolution.DeconvolveTest;
import deconvolution.LucyRichardson;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import klim.Thresholding;
import mpicbg.imglib.algorithm.mirror.MirrorImage;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.wrapper.ImgLib2;
import util.ImgLib2Util;

import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.img.Img;

public class DeconvolutionTest {

	// calculate the size of the offset
	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long [] offset){
		for(int d = 0; d < psf.numDimensions(); ++d)
			offset[d] = psf.dimension(d)/2; 
	}


	public static <T extends RealType<T>> long getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf, T tVal){
		int numDimensions = img.numDimensions();		
		
		// detect high intensity pixels: candidates for bead centers
		Thresholding.threshold(img, out, tVal);

		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType())); 		
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(numDimensions);

		// TODO: Add check: if two different labels are too close that should be one bead
		
		ObjectSegmentation.findBeads(img, labeling, beads);
		// offset for beads in all dimensions
		long[] offset = new long[numDimensions];
		getOffset(psf, offset);
		long numBeads = 0; // total number of beads
		numBeads = sumPsf(img, psf, beads, offset);
		return numBeads;
	}

	// check if bead is far enough from the boundaries 
	public static boolean isFar(long[] min, long[] max, long[] offset, int numDimensions){	
		boolean isBroken = false;	
		for (int d = 0; d < numDimensions; d++){
			if ((max[d] - min[d]) != 2*offset[d]){
				isBroken = true;
				break;
			}
		}
		return isBroken;
	}

	// check if the bead is alone in the cropped image
	public static <T extends RealType<T>> boolean isAlone(long[] min, long[] max, long[] offset, long [] position, PointSampleList<T> beads, int numDimensions){	
		boolean isAlone = true;
		Cursor<T> jCursor = beads.cursor();	
		
		while(jCursor.hasNext()){
			jCursor.fwd(); 	
			boolean isSame = true;
			
			for(int d = 0; d < numDimensions; ++d)
				if (jCursor.getLongPosition(d) != position[d])
					isSame = false;
			
			if(!isSame){ // not the same bead
				boolean isInside = true;	
				for (int d = 0; d < numDimensions; d++){
					if ((min[d] > jCursor.getLongPosition(d)) || (jCursor.getLongPosition(d) > max[d])){
						isInside = false;
						break;
					}
				}
				if (isInside){
					isAlone = false;
					break;
				}
			}
		}
				
		return isAlone;
	}
	
	// sum values of the psf (beads) for one image 
	public static <T extends RealType<T>> long sumPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long [] offset){
		int numDimensions = img.numDimensions();
		
		boolean isBroken = false; // bead initially is fine 
		long numBrokenBeads = 0;  // # of bad bead images
		Cursor<T> beadsCursor = beads.cursor();

		// this two store min/max for current bead window
		long[] min = new long[numDimensions];
		long[] max = new long[numDimensions]; 
		
		long[] position = new long[numDimensions];

		while(beadsCursor.hasNext()){
			beadsCursor.fwd();
			isBroken = false;
			
			// setting to 0 and img.max(d) ensures the normal program flow
			for (int d = 0; d < numDimensions; d++){
				min[d] = Math.max(beadsCursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(beadsCursor.getLongPosition(d) + offset[d], img.max(d));	
			}		
			
			// if the bead (+ boundary) fits in the image
			isBroken = isFar(min, max, offset, numDimensions);
			if (isBroken){
				numBrokenBeads++;
			}

			beadsCursor.localize(position);
			// if the bead is alone in the cropped image
			if (!isBroken){
				isBroken = !isAlone(min, max, offset, position, beads, numDimensions);
				if (isBroken){
					numBrokenBeads++;
				}
			}
			if(!isBroken){
				
				// TODO: Here should come the subpixel fitting routine
				
				
				// ImageJFunctions.show(Views.interval(img, min, max));
				accumulateValues(Views.offset(Views.interval(img, min, max), min), psf);
			}
		}
		
		return (beads.size()  - numBrokenBeads); // return number of clear non-broken beads 
	}

	// this function adds data to psf
	public static <T extends RealType<T>> void accumulateValues(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf){
		Cursor<T> c  = Views.iterable(img).cursor();
		RandomAccess<T> r = psf.randomAccess();	
		while(c.hasNext() ){
			c.fwd();
			r.setPosition(c);
			r.get().add(c.get());
		}
	}

	// averages psf
	public static <T extends RealType<T>> void getAverageValue(RandomAccessibleInterval<T> psf, long total){
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while(viewCursor.hasNext()){
			viewCursor.fwd();
			viewCursor.get().setReal(viewCursor.get().getRealDouble()/total);
		}
	}

	// set all pixels to 0
	public static <T extends RealType<T>> void setZero(RandomAccessibleInterval<T> psf){
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while(viewCursor.hasNext()){
			viewCursor.fwd();
			viewCursor.get().setZero();
		}
	}
	
	// subtract noise (T val)
	public static <T extends RealType<T> & Comparable<T>> void thresholdNoise(RandomAccessibleInterval<T> psf, T tVal){		
		System.out.println("Noise threshold-value: " + tVal);

		Cursor<T> cursor = Views.iterable(psf).cursor();
		while(cursor.hasNext()){
			cursor.fwd();			
			cursor.get().setReal(Math.max(0, cursor.get().getRealDouble() - tVal.getRealDouble()));
		}
	}
	
	// MAIN FUNCTION! run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf){
		AdjustInput.adjustImage( ImgLib2.wrapFloatToImgLib1( img ), LucyRichardson.minValue, 1 );
		AdjustInput.normImage( ImgLib2.wrapFloatToImgLib1( psf ) );
		
		ImageJFunctions.show(psf);
		ImageJFunctions.show(img);

		DeconvolveTest.deconvolve( DeconvolveTest.createInput( ( img ), psf ) );
	}


	// this one will be implementation of the subpixel accuracy
	public static <T extends FloatType> void subpixelAccuracy(Img<FloatType> img, Img<FloatType> psf){
		// mpicbg.imglib.algorithm.peak
	}
	
	// run the test case in 3D
	public static void test3D(){
		
		String path = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( path + "DOTs-79-82.tif" ) ); 
		Img <FloatType> psf = ImgLib2Util.openAs32Bit( new File( path + "PSF-done-50.tif" ) );
		
		Img< FloatType > convImg = img.copy();
				
		FFTConvolution< FloatType > conv = new FFTConvolution< FloatType >( convImg, psf );
		conv.convolve();
			
		ImageJFunctions.show( img ).setTitle("initial image");
		ImageJFunctions.show( convImg ).setTitle("convolved image");
			
		 runDeconvolution(img, psf);
	}
	
	//  opens multiple files and extracts beads from them
	public static void extractBeadsFromMultipleFiles(Img<FloatType> psf, int numImgs){
		long totalBeads = 0;
		
		for (int i = 1; i <= numImgs; ++i ){
			File file = new File("../Desktop/latest_desktop/beads_to_go/beads" + i + ".tif");
			Img< FloatType > beads = ImgLib2Util.openAs32Bit(file);
			Img< BitType > out = new ArrayImgFactory<BitType>().create(beads, new BitType());
			Normalize.normalize(beads,  new FloatType(0), new FloatType(255));
			float tVal = 120;
			totalBeads += getPsf(beads, out, psf, new FloatType(tVal));				
		}
		
		System.out.println("Total number of beads found: " + totalBeads);
		
		getAverageValue(psf, totalBeads);
		FloatType tVal = new FloatType(22.0f);
		thresholdNoise(psf, tVal);
	}
	
	public static void mainDeconvolution(){
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "/home/milkyklim/Desktop/latest_desktop/worm-piece.tif" ) ); 
		Img <FloatType> psf = ImgLib2Util.openAs32Bit( new File( "/home/milkyklim/Desktop/latest_desktop/PSF-done-75.tif" ) );
		runDeconvolution(img, psf);
	}
		
	// test the total intensity 
	public static <T extends FloatType> void testTotalIntensity(RandomAccessibleInterval<T> input, RandomAccessibleInterval<T> output){
		T iTotal = input.randomAccess().get();
		T oTotal = output.randomAccess().get();
		iTotal.set(0);
		oTotal.set(0);
		
		Cursor<T> iCursor = Views.iterable(input).cursor();
		RandomAccess<T> oRandomAccess = output.randomAccess();
		
		while(iCursor.hasNext()){
			iCursor.fwd();
			oRandomAccess.setPosition(iCursor);
			
			iTotal.add(iCursor.get());
			oTotal.add(oRandomAccess.get());		
		}
		
		System.out.println("input = " + iTotal.get());
		System.out.println("output = " + iTotal.get());
		
	}
	
	public static void runTestTotalIntensity(){
		
		String path = "/Users/kkolyva/Desktop/latest_desktop/20_09_16_psf_results/";
		
		Img <FloatType> input  = ImgLib2Util.openAs32Bit( new File(path + "DOTs-79-82.tif") ); 
		Img <FloatType> output = ImgLib2Util.openAs32Bit( new File(path + "deconvoled-100.tif" ) );
		
		testTotalIntensity(input, output);			
	}
	
	
	public static void main(String[] args){	
		new ImageJ();
		
		// test3D();
		runTestTotalIntensity();
		System.out.println("Doge!");

	}
}
