package weka.test;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.RandomForest;
import trainableSegmentation.*;
//import hr.irb.fastRandomForest.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import hr.irb.fastRandomForest.FastRandomForest;
import ij.ImageJ;
import ij.ImagePlus;



public class WekaTest {

	public static class processThread<T extends RealType<T>> implements Callable<Void>{
		final ImagePortion task;
		final RandomAccessibleInterval<T> src;
		final RandomAccessibleInterval< T > dst;
		final int n;
		final WekaSegmentation segmentator;

		public processThread( final ImagePortion task, final RandomAccessibleInterval<T> src, final RandomAccessibleInterval< T > dst, WekaSegmentation segmentator)
		{
			this.task = task;
			this.src = src;
			this.dst = dst;
			this.n = src.numDimensions();
			this.segmentator = segmentator;
		}

		@Override
		public Void call() throws Exception{
			
			ImagePlus imp = ImageJFunctions.wrap(src, "");
			
			
			// final Cursor<T> cSrc = Views.iterable(src).localizingCursor();
			// final RandomAccess<T> rDst = dst.randomAccess();

			// cSrc.jumpFwd(task.getStartPosition());

			// do something with image 			
			for (int j = 0; j < task.loopSize; ++j){
				cSrc.fwd();
				
				// TODO: the segmentation should come here 
				
				
				
				
				rDst.setPosition(cSrc);
				rDst.get().set(cSrc.get());
			}

			return null;
		} 		
	}

	public static class ImagePortion
	{
		public ImagePortion( final long startPosition, long loopSize )
		{
			this.startPosition = startPosition;
			this.loopSize = loopSize;
		}

		public long getStartPosition() { return startPosition; }
		public long getLoopSize() { return loopSize; }

		protected long startPosition;
		protected long loopSize;

		@Override
		public String toString() { return "Portion [" + getStartPosition() + " ... " + ( getStartPosition() + getLoopSize() - 1 ) + " ]"; }
	}

	public static final Vector<ImagePortion> divideIntoPortions( final long imageSize, final int numPortions )
	{
		final long threadChunkSize = imageSize / numPortions;
		final long threadChunkMod = imageSize % numPortions;

		final Vector<ImagePortion> portions = new Vector<ImagePortion>();

		for ( int portionID = 0; portionID < numPortions; ++portionID )
		{
			// move to the starting position of the current thread
			final long startPosition = portionID * threadChunkSize;

			// the last thread may has to run longer if the number of pixels cannot be divided by the number of threads
			final long loopSize;
			if ( portionID == numPortions - 1 )
				loopSize = threadChunkSize + threadChunkMod;
			else
				loopSize = threadChunkSize;

			portions.add( new ImagePortion( startPosition, loopSize ) );
		}

		return portions;
	}

	public static <T extends RealType<T>> void runParallelMethod(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> dst){
		int numThreads = 4;
		int numTasks = 16; 
		// TODO: change numTasks to an array of chunks 
		
		// Weka code // Sequential
		ImagePlus imp = ImageJFunctions.wrap(img, "");
		WekaSegmentation segmentator = new WekaSegmentation(); 
		
		boolean isLoaded = segmentator.loadClassifier("src/main/resources/classifier.model"); // OK! 
		if (!isLoaded){
			System.out.println("Problem loading classifier");
			return;
		}
		
		
		// @Parallel : 
		final ExecutorService taskExecutor = Executors.newFixedThreadPool( numThreads );
		final ArrayList< Callable<Void> > taskList = new ArrayList< Callable< Void > >(); 

		// @Parallel :
		for (final ImagePortion task : divideIntoPortions(Views.iterable(img).size(), numTasks) )
			taskList.add(new processThread< T >( task, img, dst, segmentator));

		try
		{
			// invokeAll() returns when all tasks are complete
			// synchronization point
			taskExecutor.invokeAll(taskList);
		}
		catch(final Exception e){
			System.out.println( "Failed to invoke all tasks" + e );
			e.printStackTrace();
		}
		taskExecutor.shutdown();

	}

	public static <T extends RealType<T> & Comparable<T>> void compareRA(RandomAccessibleInterval<T> img1, RandomAccessibleInterval<T> img2){
		Cursor<T> c = Views.iterable(img1).cursor();
		RandomAccess<T> r = img2.randomAccess();
		while(c.hasNext()){
			c.fwd();
			r.setPosition(c);
			if (c.get().compareTo(r.get()) != 0){
				System.out.println("Wrong?");
				return;
			}
		}
		System.out.println("Correct!");
	}

	public static <T extends RealType<T>> void runMethod(RandomAccessibleInterval<T> img){
		// create Weka segmentator 

		ImagePlus ip = ImageJFunctions.wrap(img, "");
		ip.show();
		WekaSegmentation segmentator = new WekaSegmentation(ip); 

		boolean isLoaded = segmentator.loadClassifier("src/main/resources/classifier.model"); // OK! 
		if (!isLoaded){
			System.out.println("Problem loading classifier");
			return;
		}


		new ImageJ();

		// Apply classifier to current image
		segmentator.applyClassifier(true);
		ImagePlus out1 = segmentator.getClassifiedImage();
		// Display classified image
		out1.setTitle("test probabilities");
		out1.show();

	}


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
				true   /* Neighbors */
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
		Img<FloatType> dst = image.factory().create(image, image.firstElement());

		// Normalize.normalize(image, new FloatType(0),  new FloatType(255));
		// Normalize.normalize(testImage, new FloatType(0),  new FloatType(255));

		// runMethod(image, labels, testImage);
		// runMethod(testImage);

		runParallelMethod(image, dst);

		ImageJFunctions.show(image);
		ImageJFunctions.show(dst);

		compareRA(image, dst);

		// Done!
		System.out.println("Doge!");
	}
}
