package klim;

import java.io.File;
import java.util.ArrayList;

import com.sun.glass.ui.View;

import deconvolution.AdjustInput;
import deconvolution.DeconvolveTest;
import deconvolution.LucyRichardson;
import ij.ImageJ;
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
import mpicbg.imglib.wrapper.ImgLib2;
import util.ImgLib2Util;

public class DeconvolutionTest {

	// calculate the size of the offset
	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long [] offset){
		for(int d = 0; d < psf.numDimensions(); ++d)
			offset[d] = psf.dimension(d)/2; 
	}


	public static <T extends RealType<T>> void getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf){
		setZero(psf); // TODO: remove
		
		T tVal = img.randomAccess().get().createVariable();
		tVal.setReal(70);
		Thresholding.threshold(img, out, tVal);

		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType())); 		
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(img.numDimensions());

		ObjectSegmentation.findBeads(img, labeling, beads);
		// offset for beads in all dimensions
		long[] offset = new long[img.numDimensions()];
		getOffset(psf, offset);
		averagePsf(img, psf, beads, offset);

	}

	// check if bead is far enough from the boundaries 
	public static boolean isFar(long[] min, long[] max, long[] offset){	
		boolean isBroken = false;	
		for (int d = 0; d < min.length; d++){
			if ((max[d] - min[d]) != 2*offset[d]){
				isBroken = true;
				break;
			}
		}
		return isBroken;
	}
	
	// check if the bead is alone in the cropped image
	public static <T extends RealType<T>> boolean isAlone(long[] min, long[] max, long[] offset, long [] position, PointSampleList<T> beads){	
		boolean isAlone = true;
		Cursor<T> jCursor = beads.cursor();	
		int numDimensions = min.length;
		
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
	
	public static <T extends RealType<T>> void averagePsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long [] offset){
		boolean isBroken = false; // bead initially is fine 
		long numBrokenBeads = 0; // # of bad bead images
		Cursor<T> iCursor = beads.cursor();

		// this two store min/max for current bead window
		long[] min = new long[img.numDimensions()];
		long[] max = new long[img.numDimensions()]; 
		
		long[] position = new long[img.numDimensions()];

		while(iCursor.hasNext()){
			iCursor.fwd();
			isBroken = false;

			for (int d = 0; d < img.numDimensions(); d++){
				min[d] = Math.max(iCursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(iCursor.getLongPosition(d) + offset[d], img.max(d));	
			}		
			
			// if the bead (+ boundary) fits in the image
			isBroken = isFar(min, max, offset);
			if (isBroken){
				numBrokenBeads++;
			}

			iCursor.localize(position);
			// if the bead is alone in the cropped image
			if (!isBroken){
				isBroken = !isAlone(min, max, offset, position, beads);
				if (isBroken){
					numBrokenBeads++;
				}
			}
			if(!isBroken){
				// ImageJFunctions.show(Views.interval(img, min, max));
				accumulateData(Views.offset(Views.interval(img, min, max), min), psf);
			}
		}
		
		// if not all beads are broken
		if ((beads.size()  - numBrokenBeads) > 0)
			getAverageValue(psf, (beads.size()  - numBrokenBeads));

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

	// run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf){
		AdjustInput.adjustImage( ImgLib2.wrapFloatToImgLib1( img ), LucyRichardson.minValue, 1 );
		AdjustInput.normImage( ImgLib2.wrapFloatToImgLib1( psf ) );

		DeconvolveTest.deconvolve( DeconvolveTest.createInput( ( img ), psf ) );
	}

	public static void main(String[] args){
		// TODO: add the correct pic here
		// Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/inputTest3Dsmall.tif" ) ); 
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/beeds.tif" ) ); 
		Img< BitType > out = new ArrayImgFactory<BitType>().create(img, new BitType());
		Img< FloatType > psf = new ArrayImgFactory<FloatType>().create(new long[]{41, 41, 27}, new FloatType());

		new ImageJ();

		// TODO: check
		float minValue =  0;
		float maxValue = 255;
		Normalize.normalize(img,  new FloatType(minValue), new FloatType(maxValue));

		getPsf(img, out, psf);		
		ImageJFunctions.show(img).setTitle("Initial Image");
		ImageJFunctions.show(Views.rotate(psf, 0, 2)).setTitle("PSF Image");
		// runDeconvolution(img, psf);

		System.out.println("Doge!");

	}
}
