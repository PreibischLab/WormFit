package klim;

import java.io.File;
import java.util.ArrayList;

import com.sun.glass.ui.View;

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
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.Type;
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
		long numBeads = 0;
		numBeads = averagePsf(img, psf, beads, offset);
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
	
	public static <T extends RealType<T>> long averagePsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long [] offset){
		int numDimensions = img.numDimensions();
		
		boolean isBroken = false; // bead initially is fine 
		long numBrokenBeads = 0; // # of bad bead images
		Cursor<T> iCursor = beads.cursor();

		// this two store min/max for current bead window
		long[] min = new long[numDimensions];
		long[] max = new long[numDimensions]; 
		
		long[] position = new long[numDimensions];

		while(iCursor.hasNext()){
			iCursor.fwd();
			isBroken = false;

			for (int d = 0; d < numDimensions; d++){
				min[d] = Math.max(iCursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(iCursor.getLongPosition(d) + offset[d], img.max(d));	
			}		
			
			// if the bead (+ boundary) fits in the image
			isBroken = isFar(min, max, offset, numDimensions);
			if (isBroken){
				numBrokenBeads++;
			}

			iCursor.localize(position);
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
				accumulateData(Views.offset(Views.interval(img, min, max), min), psf);
			}
		}
		
		// if not all beads are broken
		if ((beads.size()  - numBrokenBeads) > 0 && false)
			getAverageValue(psf, (beads.size()  - numBrokenBeads));
		
		return (beads.size()  - numBrokenBeads);
	}

	// this function adds data to psf
	public static <T extends RealType<T>> void accumulateData(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf){
		Cursor<T> c  = Views.iterable(img).cursor();
		RandomAccess<T> r = psf.randomAccess();	
		while(c.hasNext() ){
			c.fwd();
			r.setPosition(c);
			r.get().add(c.get());
		}
	}

	public static <T extends RealType<T>> void getAverageValue(RandomAccessibleInterval<T> psf, long total){
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while(viewCursor.hasNext()){
			viewCursor.fwd();
			viewCursor.get().setReal(viewCursor.get().getRealDouble()/total);
		}
	}

	// to be 100% sure 
	public static <T extends RealType<T>> void setZero(RandomAccessibleInterval<T> psf){
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while(viewCursor.hasNext()){
			viewCursor.fwd();
			viewCursor.get().setZero();
		}
	}
	
	// kind of stuff I was doing
	public static <T extends RealType<T> & Comparable<T>> void thresholdNoise(RandomAccessibleInterval<T> psf, T tVal){		
		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		
		while(viewCursor.hasNext()){
			viewCursor.fwd();			
			if (viewCursor.get().compareTo(tVal) < 0)
				tVal.set(viewCursor.get());
		}
		
		System.out.println(tVal);
		tVal.setReal(22);
		
		viewCursor.reset();
	
		while(viewCursor.hasNext()){
			viewCursor.fwd();			
			viewCursor.get().setReal(Math.max(0, viewCursor.get().getRealDouble() - tVal.getRealDouble()));
		}
	}
	

	// run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf){
		AdjustInput.adjustImage( ImgLib2.wrapFloatToImgLib1( img ), LucyRichardson.minValue, 1 );
		AdjustInput.normImage( ImgLib2.wrapFloatToImgLib1( psf ) );
		
		ImageJFunctions.show(psf);
		ImageJFunctions.show(img);

		DeconvolveTest.deconvolve( DeconvolveTest.createInput( ( img ), psf ) );
	}

	public static void openMultipleFiles(Img<FloatType> psf){
		
		// Img< FloatType > psf = new ArrayImgFactory<FloatType>().create(new long[]{191, 191, 121}, new FloatType());
		
		long totalBeads = 0;
		
		for (int i = 1; i <= 4; ++i ){
			File file = new File("../Desktop/latest_desktop/beads_to_go/beads" + i + ".tif");
			Img< FloatType > beads = ImgLib2Util.openAs32Bit(file);
			Img< BitType > out = new ArrayImgFactory<BitType>().create(beads, new BitType());
			Normalize.normalize(beads,  new FloatType(0), new FloatType(255));
			float tVal = 120;
			totalBeads += getPsf(beads, out, psf, new FloatType(tVal));				
		}
		getAverageValue(psf, totalBeads);
		thresholdNoise(psf, psf.firstElement().copy());
	}
	
	public static void testDeconvolution(){
		new ImageJ();

		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/test.tif" ) ); //psi_synthetic.tif" ) );
		AdjustInput.adjustImage( ImgLib2.wrapFloatToImgLib1( img ), LucyRichardson.minValue, 1 );
		Img< FloatType > convImg = img.copy();
		
		Img< FloatType > psf = ImgLib2Util.openAs32Bit( new File( "../Desktop/beads_to_go/psf_next_try_again-small.tif" ) );
		AdjustInput.normImage( ImgLib2.wrapFloatToImgLib1( psf ) );

		FFTConvolution< FloatType > conv = new FFTConvolution< FloatType >( convImg, psf );
		conv.convolve();

		ImageJFunctions.show( img );
		ImageJFunctions.show( convImg );
		//ImageJFunctions.show( cropConvolvedDros( img ) );
		//ImageJFunctions.show( cropConvolvedDros( convImg ) ).setTitle( "input" );
		//ImageJFunctions.show( psf ).setTitle( "psf" );

		DeconvolveTest.deconvolve( DeconvolveTest.createInput( ( convImg ), psf ) );
	}
	
	public static void test2D(){
		// new ImageJ();
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/dot.tif" ) ); 
		Img <FloatType> psf = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/psf-cross-big100-new.tif" ) );
		Img< FloatType > convImg = img.copy();
		
		FFTConvolution< FloatType > conv = new FFTConvolution< FloatType >( convImg, psf );
		conv.convolve();
		
		ImageJFunctions.show( img ).setTitle("initial image");
		ImageJFunctions.show( convImg ).setTitle("convolved image");
		
		runDeconvolution(img, psf);
	}
	
	public static void test3D(){
		// new ImageJ();
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/worm-piece.tif" ) ); 
		Img <FloatType> psf = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/psf-done.tif" ) );
		Img< FloatType > convImg = img.copy();
		
		FFTConvolution< FloatType > conv = new FFTConvolution< FloatType >( convImg, psf );
		conv.convolve();
		
		ImageJFunctions.show( img ).setTitle("initial image");
		ImageJFunctions.show( convImg ).setTitle("convolved image");
		
		runDeconvolution(img, psf);
	}
	
	public static void testFunctionsHere(){

		Image<  mpicbg.imglib.type.numeric.real.FloatType > img = ImgLib2.wrapFloatToImgLib1(ImgLib2Util.openAs32Bit( new File( "../Desktop/worm-piece-one.tif" ) )); 
		
		mpicbg.imglib.image.display.imagej.ImageJFunctions.show(img).setTitle("initital");
		
		for ( int d = 0; d < img.getNumDimensions(); ++d ){
			new MirrorImage< mpicbg.imglib.type.numeric.real.FloatType >( img, d ).process();
			mpicbg.imglib.image.display.imagej.ImageJFunctions.show(img);
		}
	}
	
	public static void main(String[] args){
		// TODO: add the correct pic here
		// Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/worm-piece.tif" ) ); 
		// Img< FloatType > beadsImg = ImgLib2Util.openAs32Bit( new File( "../Desktop/beads-bw-large.tif" ) ); 
		// Img< BitType > out = new ArrayImgFactory<BitType>().create(beadsImg, new BitType());
//		Img< FloatType > psf = new ArrayImgFactory<FloatType>().create(new long[]{191, 191, 121}, new FloatType());

		// Img <FloatType> fakePsf = ImgLib2Util.openAs32Bit( new File( "../Desktop/fakePSF.tif" ) );
		// Img <FloatType> readyPsf = ImgLib2Util.openAs32Bit( new File( "../Desktop/beads_to_go/psf_next_try_scaled.tif" ) );
		
		new ImageJ();

		// TODO: check
		// float minValue =  0;
		// float maxValue = 255;
		// Normalize.normalize(beadsImg,  new FloatType(minValue), new FloatType(maxValue));
		
		// used to adjust parameters
		// float tVal = 14000;
		// Thresholding.threshold(beadsImg, out, new FloatType(tVal));
		// ImageJFunctions.show(out).setTitle("Thresholded");
		// setZero(psf);
		// getPsf(beadsImg, out, psf, new FloatType(tVal));	
		
		// @DEBUG: Noise should not be an issue
		// thresholdNoise(psf, new FloatType(2150));
		// openMultipleFiles(psf);
		
		// ImageJFunctions.show(beadsImg).setTitle("Initial Beads Image");
		// ImageJFunctions.show(readyPsf).setTitle("PSF");
 		// ImageJFunctions.show(psf).setTitle("PSF");
		
		// test2D();
		
		// testFunctionsHere();
		
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/worm-piece.tif" ) ); 
		Img <FloatType> psf = ImgLib2Util.openAs32Bit( new File( "../Desktop/latest_desktop/psf-done.tif" ) );
		runDeconvolution(img, psf);
		//testDeconvolution();
//		
		
		// ImageJFunctions.show(img).setTitle("Deconvolved Image");
		
		
		System.out.println("Doge!");

	}
}
