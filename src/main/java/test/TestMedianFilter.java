package test;

import java.io.File;
import java.util.concurrent.TimeUnit;

import ij.ImageJ;
import klim.MedianFilter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class TestMedianFilter {
	
	public static void testMedianFilterSliced(){
		new ImageJ();
		File file = new File("/home/milkyklim/Desktop/beads_tifs/cropped/tif-14-c.tif");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		
		ImageJFunctions.show(img);
		int [] kernelDim = new int[]{61, 61};
		long inT = System.nanoTime();
		MedianFilter.medianFilterSliced(img, dst, kernelDim);
		System.out.println("kernel = " + kernelDim[0] + "x" + kernelDim[1] + " : "+ TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - inT)/1000.0);
		ImageJFunctions.show(dst);
	}
	
	public static void main(String [] args){
		testMedianFilterSliced();
	}
}
