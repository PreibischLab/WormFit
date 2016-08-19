package weka.test;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;
import spim.process.cuda.Block;
import spim.process.cuda.BlockGeneratorFixedSizePrecise;
import spim.process.fusion.FusionHelper;
import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.classifiers.trees.RandomForest;
import trainableSegmentation.*;
//import hr.irb.fastRandomForest.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import hr.irb.fastRandomForest.FastRandomForest;
import ij.ImageJ;
import ij.ImagePlus;
import mpicbg.imglib.image.ImageFactory;

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
			try
			{
				ImagePlus out1 = segmentator.applyClassifier(ImageJFunctions.wrap(Views.interval(src, task.min, task.max), ""), 0, false);
				// Display classified image
				out1.setTitle("");
				out1.show();

				// testFunction();
			}
			catch (Exception e ) {
				System.out.println(task.toString());
				e.printStackTrace();
			}

			return null;
		} 		

		// simply the test function 
		// сhecks whether all threads are running
		public void testFunction(){
			final Cursor<T> cSrc = Views.interval(src, task.min, task.max).localizingCursor();
			final RandomAccess<T> rDst = dst.randomAccess();

			System.out.println("Called!");
			// cSrc.jumpFwd(task.getStartPosition());
			while(cSrc.hasNext()){
				cSrc.fwd();
				rDst.setPosition(cSrc);
				rDst.get().set(cSrc.get());
			}
		}
	}

	public static class processThreadBlock implements Callable<Void>{
		final Block block;
		final RandomAccessibleInterval<FloatType> src; 	
		final RandomAccessibleInterval<FloatType> dst;
		final WekaSegmentation segmentator;

		public processThreadBlock(Block block, RandomAccessibleInterval<FloatType> src, final RandomAccessibleInterval< FloatType > dst, WekaSegmentation segmentator)
		{
			this.block = block;
			this.src = src;
			this.dst = dst;
			this.segmentator = segmentator;
		}

		@Override
		public Void call() throws Exception{
			try
			{
				final ImagePlus out1 = segmentator.applyClassifier(ImageJFunctions.wrap(src, ""), 0, false);
				block.pasteBlock(dst, ImageJFunctions.wrap(out1));
			}
			catch (Exception e ) {
				e.printStackTrace();
			}
			return null;
		} 		

	}


	public static class ImagePortion
	{

		public ImagePortion(final long[] min, final long[] max){

			// allocate memory
			this.min = new long[min.length];
			this.max = new long[max.length];

			for (int i = 0; i < min.length; ++i){
				this.min[i] = min[i];
				this.max[i] = max[i];
			}
		}

		// write the min values to result
		public void getMin(long [] result){
			for (int i = 0; i < min.length; ++i){
				result[i] = this.min[i];
			}
		}

		// write the min values to result
		public void getMax(long [] result){
			for (int i = 0; i < max.length; ++i){
				result[i] = this.max[i];
			}
		}

		// set min value for the provided dimension
		public void setMin(int d, long val){
			min[d] = val;
		}

		// set min value for the provided dimension
		public void setMax(int d, long val){
			max[d] = val;
		}		


		protected long[] min;
		protected long[] max;

		@Override
		public String toString() { return "Portion [" + min[0] + " ... " + max[0] + " ]" +  "[" + min[1] + " ... " + max[1] + " ]"; }
	}

	public static final Vector<ImagePortion> divideIntoPortions( final long[] dimensions, final int [] numPortions){

		final long[] threadChunkSize = new long[dimensions.length];
		final long[] threadChunkMod = new long[dimensions.length];

		for (int d = 0; d < dimensions.length; ++d){
			// TODO: quick fix to make the code scalable / rewrite
			threadChunkSize[d] = (d == 0 ? dimensions[d] / numPortions[d] : dimensions[d]);
			threadChunkMod[d] =  (d == 0 ? dimensions[d] % numPortions[d] : 0            );			
		}

		final Vector<ImagePortion> portions = new Vector<ImagePortion>();

		long [] min = new long[dimensions.length];
		long [] max = new long[dimensions.length];	

		for (int d = 1; d < dimensions.length; ++d){
			min[d] = 0;
			max[d] = dimensions[d] - 1;
		}

		// create all portions without coordinates
		// TODO: useful
		//		for (int d = 0; d < dimensions.length; ++d){
		//			for(int portionID = 0; portionID < numPortions[d]; ++portionID){
		//				portions.add( new ImagePortion( new long[dimensions.length], new long[dimensions.length] ) );
		//			}
		//		}

		// iterate only over x-axis (we currently take the full slice of the image)
		for(int portionID = 0; portionID < numPortions[0]; ++portionID){
			// move to the starting position of the current thread
			min[0] = portionID * threadChunkSize[0];
			// the last thread may has to run longer if the number of pixels cannot be divided by the number of threads
			if ( portionID == numPortions[0] - 1 )
				max[0] = threadChunkSize[0]*(portionID + 1) + threadChunkMod[0] - 1;
			else
				max[0] = threadChunkSize[0]*(portionID + 1) - 1;

			portions.addElement(new ImagePortion(min, max));
		}


		return portions;
	}

	public static <T extends RealType<T>> void runParallelMethod(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> dst){
		int numThreads = 4;
		int [] numTasks = new int[img.numDimensions()];

		// Weka code // Sequential
		// TODO:  remove? 
		ImagePlus imp = ImageJFunctions.wrap(img, "");
		WekaSegmentation segmentator = new WekaSegmentation(imp); 

		segmentator.getMembranePatchSize(); // this one should be included in the size of the image

		boolean isLoaded = segmentator.loadClassifier("src/main/resources/classifier.model"); // OK! 
		if (!isLoaded){
			System.out.println("Problem loading classifier");
			// return;
		}
		else{
			System.out.println("Loaded!");
		}

		// @Parallel : 
		final ExecutorService taskExecutor = Executors.newFixedThreadPool( numThreads );
		final ArrayList< Callable<Void> > taskList = new ArrayList< Callable< Void > >(); 

		long [] dimensions = new long[img.numDimensions()];
		for (int d = 0; d < img.numDimensions(); ++d){
			dimensions[d] = img.dimension(d);
			numTasks[d] = 1; 
		}

		// split only x axis
		numTasks[0] = 3; // this one should be user defined

		// @Parallel :
		for (final ImagePortion task : divideIntoPortions(dimensions, numTasks) ){
			taskList.add(new processThread< T >( task, img, dst, segmentator));
			// System.out.println("tasks created");
		}

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

	// this version uses Block implementation from spim reconstruction
	// bacause Block is optimized for floats we use them here 
	public static void runBlockMethod(RandomAccessibleInterval<FloatType> img, RandomAccessibleInterval<FloatType> dst, RandomAccessibleInterval<FloatType> testImg, long[] blockSize){
		int numDimensions = testImg.numDimensions();
		final long[] testImgSize = new long[numDimensions];
		testImg.dimensions( testImgSize );

		// Weka code
		final ImagePlus imp = ImageJFunctions.wrap(img, "");
		final WekaSegmentation segmentator = new WekaSegmentation(imp); 
		boolean isLoaded = segmentator.loadClassifier("src/main/resources/classifier.model"); 

		if (!isLoaded){
			System.out.println("Problem loading classifier");
			// return;
		}
		else{
			System.out.println("Loaded!");
		}

		final long[] kernelSize = new long[numDimensions]; // used by classifier
		long kernel = (long) Math.max(segmentator.getMembranePatchSize(), segmentator.getMaximumSigma()); 
		// TODO: maybe for 3D z-axis offset is not necessary		
		for (int d = 0; d < numDimensions; ++d)
			kernelSize[d] = kernel;	
		
		final BlockGeneratorFixedSizePrecise blockGenerator = new BlockGeneratorFixedSizePrecise(blockSize);
		
		final Block[] blocks = blockGenerator.divideIntoBlocks( testImgSize, kernelSize ); 
		final long [] extendedBlockSize = new long[numDimensions + 1]; // used as temporary storage
		for (int d = 0; d < numDimensions; ++d)
			extendedBlockSize[d] = blockSize[d];
		extendedBlockSize[numDimensions] = blocks.length;
		
		final RandomAccessibleInterval< FloatType > block = ArrayImgs.floats(extendedBlockSize);
		
		int numThreads = Runtime.getRuntime().availableProcessors(); // either predefine or give as a parameter

		// @Parallel : 
		final ExecutorService taskExecutor = Executors.newFixedThreadPool( numThreads );
		final ArrayList< Callable<Void> > taskList = new ArrayList< Callable< Void > >(); 
		
		for (int idx = 0; idx < blocks.length; ++idx){
			final IntervalView<FloatType> slice  = Views.hyperSlice(block, numDimensions, idx);	
			blocks[idx].copyBlock(Views.extendMirrorDouble(testImg), slice);			
			taskList.add(new processThreadBlock(blocks[idx], slice, dst, segmentator));			
		}
		ImageJFunctions.show(block);

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
		Img<FloatType> dst = image.factory().create(testImage, testImage.firstElement());

		// Normalize.normalize(image, new FloatType(0),  new FloatType(255));
		// Normalize.normalize(testImage, new FloatType(0),  new FloatType(255));

		// runMethod(image, labels, testImage);
		// runMethod(testImage);

		// runParallelMethod(image, dst);
		// TODO: move this one to the function

		new ImageJ();
		
		long[] blockSize = new long[]{384, 384};				
		runBlockMethod(image, dst, testImage, blockSize);

		ImageJFunctions.show(testImage);
		ImageJFunctions.show(dst).setTitle("This one should be the correct image");

		compareRA(image, dst);

		// Done!
		System.out.println("Doge!");
	}
}
