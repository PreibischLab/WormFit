package cluster;

import java.io.File;

import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

import ij.io.FileSaver;
import util.ImgLib2Util;

public class ClusterDeconvolution {

	/**
	 * Test for loading and saving ImageJ images on cluster
	 * */
	public static void loadSaveImglib2(String[] args){
		// first arg is the directory + name of the file without extension 
		
		if (args.length != 1) 
			return;
		
		String fullPath = args[0];
		
		Img<FloatType> img =  ImgLib2Util.openAs32Bit(new File(fullPath + ".tif"));		
		new FileSaver(ImageJFunctions.wrap(img, "").duplicate()).saveAsTiffStack(fullPath + "NEW" + ".tif");
	}
	
	
	public static void main(String[] args) {
		loadSaveImglib2(args);
		System.out.println("DOGE!");
	}

}



