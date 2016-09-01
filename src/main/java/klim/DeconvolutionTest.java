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

	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long [] offset){
		for(int d = 0; d < psf.numDimensions(); ++d){
			offset[d] = psf.dimension(d)/2; 
			// offset[d] = (d <= 1 ? 20 : 13);
		}
	}


	public static <T extends RealType<T>> void getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf){
		setZero(psf);
		T tVal = img.randomAccess().get().createVariable();
		tVal.setReal(70);
		Thresholding.threshold(img, out, tVal);

		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType())); 		
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(img.numDimensions());

		ObjectSegmentation.findBeads(img, labeling, beads);
		// TODO: remove beads with noise in the offset region

		// offset for beads in all dimensions
		long[] offset = new long[img.numDimensions()];
		getOffset(psf, offset);

		// @DEBUG: 
		// for(int i =0; i < img.numDimensions(); ++i)
		// 	System.out.println(offset[i]);

		// not working yet! 
		averagePsf(img, psf, beads, offset);

	}

	// TODO: ? becomes to complicated
	public static <T extends RealType<T>> void setMinMax(long[] position, long[] min, long[] max, long[] offset, long[] imgMax, int n){
		for (int d = 0; d < position.length; d++){
			// TODO: this part should be fixed with mirroring strategy
			min[d] = Math.max(position[d] - offset[d], 0);
			max[d] = Math.min(position[d] + offset[d], imgMax[d]);
			// System.out.print("[" + min[d] + " " + max[d] + "] ");
		}
	}

	public static <T extends RealType<T>> void averagePsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long [] offset){
		Cursor<T> cursor = beads.cursor();

		long numBrokenBeads = 0;
		boolean isFirst = true;
		while(cursor.hasNext()){
			cursor.fwd();
			// this two store min/max for current bead
			long[] min = new long[img.numDimensions()];
			long[] max = new long[img.numDimensions()]; // should this one be excluded ? 

			long extra = 0;

			for (int d = 0; d < img.numDimensions(); d++){
				// TODO: FIXED: this part should be fixed with mirroring strategy
				// To make calculations more precise we do not take into account broken psf's
				min[d] = Math.max(cursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(cursor.getLongPosition(d) + offset[d] + extra, img.max(d));
				System.out.print("[" + min[d] + " " + max[d] + "] ");
			}
			// System.out.println("");

			ImageJFunctions.show(Views.interval(img, min, max));

			// check if the bead is far enough from the boundary 
			boolean isFull = true; 
			for (int d = 0; d < img.numDimensions(); d++){
				if ((max[d] - min[d] - extra) != 2*offset[d]){
					isFull = false;
					numBrokenBeads++;
					break;
				}
			}

			boolean isAlone = true; 
			// TODO: here should come the check that the bead is alone in the pic

			if(isFull && isFirst && isAlone){
				// isFirst = false;
				
				accumulateData(Views.offset(Views.interval(img, min, max), min), psf);
				ImageJFunctions.show(psf);
			}
		}

	}

	// this function copies data from 
	public static <T extends RealType<T>> void accumulateData(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf){

		Cursor<T> c  = Views.iterable(img).cursor();
		RandomAccess<T> r = psf.randomAccess();

		T tmp = Views.iterable(img).firstElement();

		// @Debug: check total number of elements

		long [] total = new long[]{1, 1}; 
		for(int d = 0; d < img.numDimensions(); ++d){
			total[0] *= img.dimension(d);
			total[1] *= psf.dimension(d);

			System.out.println(img.dimension(d) + " ?= " + psf.dimension(d));
		}

		//42820

		System.out.println(total[0] + " ?= " + total[1]);

		long maxIdx = 45386;
		long idx = 0;

		long [] positions = new long[img.numDimensions()];

		while(c.hasNext() ){
			c.fwd();
			c.localize(positions);

			if(idx < 100){
				for (int d = 0; d < img.numDimensions(); ++d)
					System.out.print(positions[d] + " ");
				System.out.println();
			}
			
			
			idx++;
			if (idx >= maxIdx){
				break;
			}
			r.setPosition(c);
			r.localize(positions);
			if(idx < 100){
				for (int d = 0; d < img.numDimensions(); ++d)
					System.out.print(positions[d] + " ");
				System.out.println();
				System.out.println("-------");
			}
			
			
			
			//System.out.println(idx);
			r.setPosition(c);
			// c.get();
			r.get().set(c.get());
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
		// TODO: ?!?!?!?!?!
		ImageJFunctions.show(psf).setTitle("PSF Image");

		// ImageJFunctions.show(out).setTitle("Thresholded Image");

		// runDeconvolution(img, psf);

		System.out.println("Doge!");

	}
}
