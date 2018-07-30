package klim.deconvolution;

import java.io.File;
import java.util.ArrayList;

import deconvolution.AdjustInput;
import deconvolution.DeconvolveTest;
import deconvolution.LucyRichardson;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import mpicbg.imglib.wrapper.ImgLib2;
import klim.MedianFilter;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.real.FloatType;
import parameters.GUIParams;
import plane.fit.BackgroundSubtraction;
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

		final ArrayList<double[]> results = new ArrayList<>();
		DeconvolutionTest.gaussianFitting(img, out, new FloatType(tVal), results );

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
		long psfSizeXYZ = 15;
		long [] psfSize = new long[]{psfSizeXYZ,psfSizeXYZ,psfSizeXYZ};
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType());
		Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "beads-bw-large.tif"));
		Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		ImageJFunctions.show(beads);

		Normalize.normalize(beads, new FloatType(0), new FloatType(255));
		
		// TODO: stuff above should be ok
		
		// TODO: Uncomment when done with debugging
		// thresholdNoise(psf, new FloatType(70.0f));
		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));

		float tVal = 180;
		totalBeads += DeconvolutionTest.getPsf(beads, out, psf, new FloatType(tVal), (long)(psfSize[0]/2));


		System.out.println("Total number of beads found: " + totalBeads);
		
		// if (true) return;

		Utils.getAverageValue(psf, totalBeads);
		Utils.thresholdNoise(psf, new FloatType(10.0f));

		ImageJFunctions.show(psf);
	}
	


	public static void runRSExtractBeads(){
		
		// use plane fitting 
		boolean planeFit = false; 
		
		// long totalBeads = 0;		

		String pathMac = "/Users/kkolyva/Desktop/latest_desktop/";
		String pathUbuntu = "/home/milkyklim/Desktop/2017-04-27-beads/";
		// String fileName = "beads-bw-large.tif"; 

		String path = pathUbuntu;
		
		// String[] fileSeries = new String[]{ "3" , "11", "13", "14" }; // 
		int nFiles = 17; // fileSeries.length;
		
		long [] psfSize = new long[]{211, 211, 97};
		long [] psfHalfSize = new long[psfSize.length];
		for(int d = 0; d < psfHalfSize.length; d++)
			psfHalfSize[d] = psfSize[d]/2;
		
		Img<FloatType> rPsf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType()); // global psf over all images and beads
		Utils.setZero(rPsf);
		
		
		// Radial Symmetry parameters must be the same over all images
		// <= because the parameters of the microscope are the same
		final GUIParams params = new GUIParams();
		params.setDefaultValues(); 

		// checked image beforehand (used interactive radial symmetry here)
		params.setSigmaDog(4);
		params.setThresholdDog(0.003f);  // (3.5f); 
		params.setSupportRadius(2);
		params.setMaxError(3);
		params.setInlierRatio(0.9f);
				
		long totalNumBeads = 0;
		Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType()); // local psf over one image and all beads there
		for (int iFile = 1; iFile <= nFiles; ++iFile){
			
			// hack for broken beads
			if (iFile == 14 || iFile == 17 || iFile == 15)
				continue;
			
			Utils.setZero(psf);
			Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + "processed/tif-" + iFile + "-c-m-p.tif"));
			String fullPath = path + "beads/tif-" + iFile + "-c-m-p-b"; // used to save the resulting beads
			if (planeFit){
				
			}
			
			// ArrayList<long[]> center = new ArrayList<>(1);
			// center.add(new long []{0, 0})
			// BackgroundSubtraction.applyBackgroundSubtraction(beads, peaks, offset);
			
			// Img<FloatType> adjustedBeads = ImgLib2Util.openAs32Bit(new File(path + fileSeries[iFile]));	
			Normalize.normalize(beads, new FloatType(0), new FloatType(1));
				
			long curNumBeads = DeconvolutionTest.useRadialSymmetry(beads, psfHalfSize, psf, rPsf, params, fullPath);
			System.out.println("image: " + iFile + " with " + curNumBeads + " beads");
			
			totalNumBeads += curNumBeads;
		}
		
		Utils.getAverageValue(rPsf, totalNumBeads);
		ImageJFunctions.show(rPsf).setTitle("Averaged PSF");
		
		ImagePlus towrite = ImageJFunctions.wrap(psf, "bead").duplicate();
		new FileSaver(towrite).saveAsTiffStack(path + "beads/averaged_bead.tif");
		
		// size of the beads should be the parameter
		// long psfSizeXYZ = 15;
		// long [] psfSize = new long[]{psfSizeXYZ, psfSizeXYZ, psfSizeXYZ};
		// long [] psfSize = new long[]{151, 151, 117};
		// Img<FloatType> psf = new ArrayImgFactory<FloatType>().create(psfSize, new FloatType());
		// Img<FloatType> beads = ImgLib2Util.openAs32Bit(new File(path + fileName));
		// Img<FloatType> adjustedBeads = ImgLib2Util.openAs32Bit(new File(path + fileName));
		// Img<BitType> out = new ArrayImgFactory<BitType>().create(beads, new BitType());

		// ImageJFunctions.show(beads);

		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));
		
		// TODO: Try on the smaller images
		// MedianFilter.medianFilter(beads, adjustedBeads, new int[]{5, 5, 5});
		
		// ImageJFunctions.show(adjustedBeads);
		
		//final GUIParams params = new GUIParams();
		//params.setDefaultValues(); 

		// checked image beforehand (used interactive radial symmetry here)
		// params.setSigmaDog(4);
		// params.setThresholdDoG(3.5f); 
		// params.setSupportRadius(2);
		// params.setMaxError(3);
		// params.setInlierRatio(0.9f);
				
		// DeconvolutionTest.useRadialSymmetry(beads, (long)(psfSize[0]/2), psf, params);
	
		// float thresholdValue = 70;
		// Utils.thresholdNoise(psf, new FloatType(thresholdValue));
		// Normalize.normalize(beads, new FloatType(0), new FloatType(255));

		// float tVal = 180;
		// totalBeads += DeconvolutionTest.getPsf(beads, out, psf, new FloatType(tVal), (long)(psfSize[0]/2));


		// System.out.println("Total number of beads found: " + totalBeads);
		
		// if (true) return;

		// Utils.getAverageValue(psf, totalBeads);
		// Utils.thresholdNoise(psf, new FloatType(10.0f));

		// ImageJFunctions.show(psf);
	}
	
	// MAIN FUNCTION! run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf) {
		runDeconvolution(img, psf, 0, "");		
	}
	
	// MAIN FUNCTION! run deconvolution
	public static <T extends FloatType> void runDeconvolution(Img<FloatType> img, Img<FloatType> psf, float reg, String filename) {
		for ( FloatType t : img )
			t.add( new FloatType( 1 ));
		AdjustInput.adjustImage(ImgLib2.wrapFloatToImgLib1(img), LucyRichardson.minValue, 1);

		System.out.println("before: " + Utils.sumIntensities(psf));
		 
		AdjustInput.normImage(ImgLib2.wrapFloatToImgLib1(psf));
		System.out.println("after: " + Utils.sumIntensities(psf));

		// offline mode
		// ImageJFunctions.show(psf).setTitle("normalized psf");
		// ImageJFunctions.show(img).setTitle("adjusted image");

		DeconvolveTest.deconvolve(DeconvolveTest.createInput((img), psf), reg, filename);
	}
	
	
	
	public static void runDeconvolutionBatch(){
		// 2 different PSF's 
		// with regularization and without one 
		
		String[] psfFiles = new String[]{"psf-0024"};
		float [] regularization = new float[]{0, 0.0006f};
		
		String path = "/home/milkyklim/Desktop/";
				
		for (int j = 1; j <= 2; j++){
			for (int r = 0; r < regularization.length; r++){
				// read img and psf
				Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "2017-04-27-beads/img-" + j + ".tif"));
				Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(path + "2017-04-27-beads/psf/" + psfFiles[0] + ".tif"));
								
				
				String resultPath = path + "2017-04-27-beads/result/" + "img-" + j + "/" + "out-" + regularization[r] + "-";
				
				
				runDeconvolution(img, psf, regularization[r], resultPath);
			}
		}
		
	} 
	
	public static void runDeconvolutionMultipleImages(){
		String[] psfFiles = new String[]{"psf"};
		String folderPath = "/home/milkyklim/Desktop/2017-04-29_smFISH-DAPI-deconvolution"; // /2017-04-29_smFISH-DAPI
		
		Img<FloatType> psf = ImgLib2Util.openAs32Bit(new File(folderPath + "/psf.tif"));
		
		int totalFiles = 40;
		for (int idx = 1; idx <= totalFiles; idx++){
			String inputFile  = folderPath + "/2017-04-29_smFISH-DAPI/" + "N2-DPY23-wdr-5.2-" + String.format("%03d", idx) + "_C2";
			String outputFile = folderPath + "/2017-04-29_smFISH-DAPI-result/" + "N2-DPY23-wdr-5.2-" + String.format("%03d", idx) + "_C2";
			
			Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(inputFile + ".tif"));
			runDeconvolution(img, psf, 0, outputFile);
			
		} 
		
		
	}
	
	public static void main(String [] args){
		runDeconvolutionMultipleImages();
		System.out.println("Doge!");
	}
	
	
}
