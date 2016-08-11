package weka.test;

import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.RandomForest;
import trainableSegmentation.*;
//import hr.irb.fastRandomForest.*;

import java.io.File;

import hr.irb.fastRandomForest.FastRandomForest;
import ij.ImageJ;
import ij.ImagePlus;


public class WekaTest {
	// U is BitType in the beginning 
	// Check how it works with multiple errors
	public static <T extends RealType<T>, U extends RealType<U>> void runMethod(RandomAccessibleInterval<T> trainImage, RandomAccessibleInterval<U> labels, RandomAccessibleInterval<T> testImage){
		// create Weka segmentator 
		WekaSegmentation segmentator = new WekaSegmentation(ImageJFunctions.wrap(trainImage, ""));
		// # of samples to consider 
		int nSamplesToUse = 2000;
		FastRandomForest rf = new FastRandomForest();
		rf.setNumTrees(100);
		rf.setNumFeatures(0);
		rf.setSeed((new java.util.Random()).nextInt());
		
		// set classifier
		segmentator.setClassifier(rf);
		segmentator.setMembranePatchSize(11);
		segmentator.setMaximumSigma(16.0f);
		
		// Selected attributes
		boolean[] enableFeatures = new boolean[]{
		            true,   /* Gaussian_blur */
		            true,   /* Sobel_filter */
		            true,   /* Hessian */
		            true,   /* Difference_of_gaussians */
		            true,   /* Membrane_projections */
		            false,  /* Variance */
		            false,  /* Mean */
		            false,  /* Minimum */
		            false,  /* Maximum */
		            false,  /* Median */
		            false,  /* Anisotropic_diffusion */
		            false,  /* Bilateral */
		            false,  /* Lipschitz */
		            false,  /* Kuwahara */
		            false,  /* Gabor */
		            false,  /* Derivatives */
		            false,  /* Laplacian */
		            false,  /* Structure */
		            false,  /* Entropy */
		            false   /* Neighbors */
		};
		
		segmentator.setEnabledFeatures(enableFeatures);
		segmentator.addRandomBalancedBinaryData(ImageJFunctions.wrap(trainImage, ""), ImageJFunctions.wrap(labels, ""), "class 2", "class 1", nSamplesToUse);
		
		segmentator.trainClassifier();
		
		// Apply classifier to current image
		segmentator.applyClassifier( true );
		// Display classified image
		ImagePlus out1 = segmentator.getClassifiedImage();
		out1.setTitle("train probabilities");
		out1.show();
		
		// Apply classifier to test image
		ImagePlus out2 = segmentator.applyClassifier(ImageJFunctions.wrap(testImage, ""), 0, true);
		// display classified test image
		out2.setTitle("test probabilities");
		out2.show();
		
	}
	

	public static void main(String args[]) throws Exception{
		String file = "src/main/resources/input";
		Img<FloatType> image = util.ImgLib2Util.openAs32Bit(new File(file + ".tif"));
		Img<FloatType> labels = util.ImgLib2Util.openAs32Bit(new File(file + "Labeled.tif"));
		Img<FloatType> testImage = util.ImgLib2Util.openAs32Bit(new File(file + "Test.tif"));
		
		runMethod(image, labels, testImage);
		// Done!
		System.out.println("Doge!");
	}
}
