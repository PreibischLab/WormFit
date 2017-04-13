package klim.deconvolution;

import java.io.File;

import deconvolution.AdjustInput;
import deconvolution.DeconvolveTest;
import deconvolution.LucyRichardson;
import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class Run {
	
	public static void runGaussianFitting(){
		String pathMac = "/Users/kkolyva/Desktop/28_09_16_severin/";

		String path = pathMac;

		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "bead.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(img, new BitType());
		Normalize.normalize(img, new FloatType(0), new FloatType(255));
		float tVal = 180;

		Utils.thresholdNoise(img, new FloatType(10.0f));

		DeconvolutionTest.gaussianFitting(img, out, new FloatType(tVal));

		// ImageJFunctions.show(psf);
	}
	
	// TODO: This one was used for testing 
	// looks like you can delete this onr
	public static void runExtractGeneratedBeads(){
		long totalBeads = 0;		

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/";
		String pathUbuntu = "/home/milkyklim/Desktop/latest_desktop/";

		String path = pathMac;
		// size of the beads should be the parameter
		long [] psfSize = new long[]{15,15,15}; 
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType());
		Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "beads-generated.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));
		float tVal = 0.98f;
		totalBeads += DeconvolutionTest.getPsf(beads, out, psf, new FloatType(tVal), (long)(psfSize[0]/2));

		System.out.println("Total number of beads found: " + totalBeads);

		Utils.getAverageValue(psf, totalBeads);
		// 		thresholdNoise(psf, new FloatType(10.0f));

		ImageJFunctions.show(psf);
	}
	
	public static void runExtractBeads(){
		long totalBeads = 0;		

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/";
		String pathUbuntu = "/home/milkyklim/Desktop/";

		String path = pathUbuntu;
		// size of the beads should be the parameter
		long [] psfSize = new long[]{15,15,15};
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType());
		Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "beads-bw-large.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		ImageJFunctions.show(beads);

		Normalize.normalize(beads, new FloatType(0), new FloatType(255));

		//thresholdNoise(psf, new FloatType(70.0f));
		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));

		float tVal = 180;
		totalBeads += DeconvolutionTest.getPsf(beads, out, psf, new FloatType(tVal),  (long)(psfSize[0]/2));


		System.out.println("Total number of beads found: " + totalBeads);

		Utils.getAverageValue(psf, totalBeads);
		Utils.thresholdNoise(psf, new FloatType(10.0f));

		ImageJFunctions.show(psf);
	}
	
	// MAIN FUNCTION! run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf) {
		for ( FloatType t : img )
			t.add( new FloatType( 1 ));
		AdjustInput.adjustImage(ImgLib2.wrapFloatToImgLib1(img), LucyRichardson.minValue, 1);

		System.out.println("before: " + Utils.sumIntensities(psf));
		 
		AdjustInput.normImage(ImgLib2.wrapFloatToImgLib1(psf));
		System.out.println("after: " + Utils.sumIntensities(psf));

		ImageJFunctions.show(psf).setTitle("normalized psf");
		ImageJFunctions.show(img).setTitle("adjusted image");

		DeconvolveTest.deconvolve(DeconvolveTest.createInput((img), psf));
	}
	
}
