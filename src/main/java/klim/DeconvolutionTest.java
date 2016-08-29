package klim;

import java.io.File;

import ij.ImageJ;
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
import klim.Thresholding;
import util.ImgLib2Util;

public class DeconvolutionTest {
	
	public static <T extends RealType<T>> void getPsf(RandomAccessibleInterval<T> img, RandomAccessibleInterval<BitType> out, RandomAccessibleInterval<T> psf){
		T tVal = img.randomAccess().get();
		tVal.setReal(128f);
		Thresholding.threshold(img, out, tVal);
		
		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(out, new IntType())); 
		ObjectSegmentation.setLabeling(out, labeling);
		ObjectSegmentation.setParameters(out, labeling);
		// @DEBUG: 
		ImageJFunctions.show(labeling.getIndexImg());
		// 
		
	}
	
	
	
	public static void main(String[] args){
		// TODO: add the correct pic here
		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/inputTest3Dsmall.tif" ) ); 
		Img< BitType > out = new ArrayImgFactory<BitType>().create(img, new BitType());
		// TODO: check
		float minValue =  0;
		float maxValue = 255;
		Normalize.normalize(img,  new FloatType(minValue), new FloatType(maxValue));
		
		getPsf(img, out, img);
		
		new ImageJ();
		// ImageJFunctions.show(img).setTitle("Initial Image");
		// ImageJFunctions.show(out).setTitle("Thresholded Image");
		
		System.out.println("Doge!");
		
	}
}
