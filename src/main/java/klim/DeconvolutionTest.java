package klim;

import java.io.File;
import java.util.ArrayList;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
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
	
	public static <T extends RealType<T>> void getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf){
		T tVal = img.randomAccess().get();
		tVal.setReal(70);
		Thresholding.threshold(img, out, tVal);
		
		// ArrayList<ObjectSegmentation> beadsList = new ArrayList<ObjectSegmentation>();
		
		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType())); 		
		ObjectSegmentation.setLabeling(out, labeling);
		PointSampleList<T> beads = new PointSampleList<T>(img.numDimensions());
		// long numBeads = 0; // -1 ? 
		// System.out.println(labeling.size());
		
		ObjectSegmentation.findBeads(img, labeling, beads);
		// beads 0 is the background!
		
		
		
		
		// ArrayList<ObjectSegmentation> beadsList = ObjectSegmentation.setParameters(out, labeling);
		// @DEBUG: 
		// ImageJFunctions.show(out);
		// 
		// ObjectSegmentation.removeNoise(beads, img);
		
		// PointSampleList<T> beads = new PointSampleList<T>(img.numDimensions());
		// tVal.setReal(100); // if above this value, probably bead center
		// detectBeads(img, beads, tVal);
		
		Cursor<T> cursor = beads.cursor();
		while(cursor.hasNext()){
			long[] min = new long[3];
			long[] max = new long[3];
			
			cursor.fwd();
			for (int d = 0; d < img.numDimensions(); d++){
				System.out.print("[" + (cursor.getLongPosition(d) - 10) + " " + (cursor.getLongPosition(d) + 10) + "]");
				min[d] = cursor.getLongPosition(d) - 10;
				max[d] = cursor.getLongPosition(d) + 10;
			}
			System.out.println();
			
			
			
			ImageJFunctions.show(Views.rotate(Views.interval(img, min, max), 0, 2));
		}
		
	}
	
	public static <T extends RealType<T> & Comparable<T>> void detectBeads(RandomAccessibleInterval<T> img, PointSampleList<T> beads, T tVal){
		Cursor<T> cursor  = Views.iterable(img).cursor();
		
		while(cursor.hasNext()){
			cursor.fwd();
			if (cursor.get().compareTo(tVal) > 0 ){
				beads.add(new Point(cursor), cursor.get());
			}
		}
	}
	
	
	public static void main(String[] args){
		// TODO: add the correct pic here
		// Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/inputTest3Dsmall.tif" ) ); 
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "../Desktop/beeds.tif" ) ); 
		
		Img< BitType > out = new ArrayImgFactory<BitType>().create(img, new BitType());
		// TODO: check
		float minValue =  0;
		float maxValue = 255;
		Normalize.normalize(img,  new FloatType(minValue), new FloatType(maxValue));
		
		getPsf(img, out, img);
		
		new ImageJ();
		ImageJFunctions.show(img).setTitle("Initial Image");
		// ImageJFunctions.show(out).setTitle("Thresholded Image");
		
		System.out.println("Doge!");
		
	}
}
