package klim.deconvolution;

import java.io.File;

import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class ClusterRun {
	
	public static void main(String [] args){
		
		// parameters via the terminal
		String folderPath;
		String imgFile;
		String psfFile;
		
		// args[] 
		if (args.length != 3){
			System.out.println("You have to provide exactly 3 arguments!");
			// return;
			folderPath = "/home/milkyklim/Desktop/input/";
			imgFile = "N2-DPY23-wdr-5.2-001_C2.tif";
			psfFile = "psf.tif";
		
		}
		else{
			// parameters via the terminal
			folderPath = args[0];
			imgFile = args[1];
			psfFile = args[2];
		}
						
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(folderPath + psfFile));
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(folderPath + imgFile));
		
		String outputFile = folderPath.substring(0, folderPath.length() - 1) + "-result/" + imgFile.substring(0, imgFile.length() - 4);
		
		Run.runDeconvolution(img, psf, 0, outputFile);
						
		System.out.println("Doge!");
	}
}
