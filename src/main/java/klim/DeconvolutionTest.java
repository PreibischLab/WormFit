package klim;

import java.io.File;
import java.util.ArrayList;

import com.sun.glass.ui.View;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
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
import util.ImgLib2Util;

public class DeconvolutionTest {

	public static <T extends RealType<T>> void getOffset(RandomAccessibleInterval<T> psf, long [] offset){
		for(int d = 0; d < psf.numDimensions(); ++d){
			offset[d] = psf.dimension(d)/2; 
			// offset[d] = (d <= 1 ? 20 : 13);
		}
	}


	public static <T extends RealType<T>> void getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf){
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

		for(int i =0; i < img.numDimensions(); ++i)
			System.out.println(offset[i]);


		// TODO: you need boundary check fo each dimension
		//		Cursor<T> cursor = beads.cursor();
		//		while(cursor.hasNext()){
		//			long[] min = new long[3];
		//			long[] max = new long[3];
		//			
		//			cursor.fwd();
		//			for (int d = 0; d < img.numDimensions(); d++){
		//				// System.out.print("[" + (cursor.getLongPosition(d) - offset[d]) + " " + (cursor.getLongPosition(d) + offset[d]) + "]");
		//				min[d] = Math.max(cursor.getLongPosition(d) - offset[d], 0);
		//				max[d] = cursor.getLongPosition(d) + offset[d];
		//			}
		//			// System.out.println();
		//			
		//			// ImageJFunctions.show(Views.rotate(Views.interval(img, min, max), 0, 2));
		//		}

		// not working yet! 
		averagePsf(img, psf, beads, offset);

	}

	public static <T extends RealType<T>> void averagePsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> psf, PointSampleList<T> beads, long [] offset){
		Cursor<T> cursor = beads.cursor();

		long numBrokenBeads = 0;
		while(cursor.hasNext()){
			cursor.fwd();
			// this two store min/max for current bead
			long[] min = new long[img.numDimensions()];
			long[] max = new long[img.numDimensions()]; // should this one be excluded ? 

			for (int d = 0; d < img.numDimensions(); d++){
				// TODO: this part should be fixed with mirroring strategy
				min[d] = Math.max(cursor.getLongPosition(d) - offset[d], 0);
				max[d] = Math.min(cursor.getLongPosition(d) + offset[d], img.max(d));
				// System.out.print("[" + min[d] + " " + max[d] + "] ");
			}
			// System.out.println("");

			boolean isFull = true; 

			for (int d = 0; d < img.numDimensions(); d++){
				if ((max[d] - min[d]) != 2*offset[d]){
					isFull = false;
					numBrokenBeads++;
					break;
				}
			}

			//			if (isFull){
			//				Cursor<T> viewCursor = Views.interval(img, min, max).cursor();
			//				// ImageJFunctions.show(Views.interval(img, min, max));
			//				// ensure that the sizes are correct OKAY!
			////				for (int d = 0; d < img.numDimensions(); d++){
			////					System.out.print(Views.interval(img, min, max).dimension(d) + " == " + psf.dimension(d) + " ? "); 
			////				}
			////				System.out.println();					
			//				RandomAccess<T> ra = psf.randomAccess();
			//				// System.out.println("Hello!2");
			//				long idx =0;
			//				while(viewCursor.hasNext()){
			//					viewCursor.fwd();
			//					ra.setPosition(viewCursor);
			//
			////					for (int d = 0; d < img.numDimensions(); d++){
			////						System.out.print("[" + min[d] + ": " + viewCursor.getLongPosition(d) + " : " + max[d] + "] ");
			////					}
			////					System.out.println();
			////					
			//					try{					
			//						ra.get().add(viewCursor.get());
			//					}
			//					catch(Exception e){
			//
			//					}
			//					//					ra.get().add(viewCursor.get());
			//				}
			//			}
			//
			//			// ImageJFunctions.show(Views.rotate(Views.interval(img, min, max), 0, 2));


			if(isFull){
				Cursor<T> psfCursor = Views.iterable(psf).cursor();
				RandomAccess<T> ra = Views.interval(img, min, max).randomAccess();
				ImageJFunctions.show(Views.interval(img, min, max));
				while(psfCursor.hasNext()){
					psfCursor.fwd();
					ra.setPosition(psfCursor);
					
					for (int d = 0; d < img.numDimensions(); d++){
						System.out.print("[" + "0" + ": " + psfCursor.getLongPosition(d) + " ?= " +  ra.getLongPosition(d) + " : " + psf.dimension(d) + "] ");
					}
					System.out.println();
					
					psfCursor.get().add(ra.get());
				}
				// ImageJFunctions.show(psf);
			}
		}

		Cursor<T> viewCursor = Views.iterable(psf).cursor();
		while(viewCursor.hasNext()){
			viewCursor.fwd();
			viewCursor.get().setReal(viewCursor.get().getRealDouble()/(beads.size() - numBrokenBeads));
		}

		ImageJFunctions.show(psf).setTitle("psf");
	}


	public static void main(String[] args){
		// TODO: add the correct pic here
		// Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/inputTest3Dsmall.tif" ) ); 
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/beeds.tif" ) ); 

		Img< BitType > out = new ArrayImgFactory<BitType>().create(img, new BitType());


		Img< FloatType > psf = new ArrayImgFactory<FloatType>().create(new long[]{41, 41, 27}, new FloatType());


		// TODO: check
		float minValue =  0;
		float maxValue = 255;
		Normalize.normalize(img,  new FloatType(minValue), new FloatType(maxValue));

		getPsf(img, out, psf);

		new ImageJ();
		ImageJFunctions.show(img).setTitle("Initial Image");
		// ImageJFunctions.show(out).setTitle("Thresholded Image");

		System.out.println("Doge!");

	}
}
