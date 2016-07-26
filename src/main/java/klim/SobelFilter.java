package klim;

import java.io.File;
import java.util.concurrent.Executors;

import stephan.MeanFilter;
import fiji.threshold.Auto_Threshold;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.plugin.frame.MemoryMonitor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.wrapper.ImgLib2;
import mpicbg.stitching.PairWiseStitchingImgLib;
import mpicbg.stitching.PairWiseStitchingResult;
import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.Interval;
import net.imglib2.IterableInterval;
import net.imglib2.KDTree;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.exception.ImgLibException;
import net.imglib2.exception.IncompatibleTypeException;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.FloatImagePlus;
import net.imglib2.img.imageplus.ImagePlusImgFactory;
import net.imglib2.multithreading.SimpleMultiThreading;
import net.imglib2.neighborsearch.NearestNeighborSearchOnKDTree;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class SobelFilter {

	private static final boolean debug = false;
	private static final boolean showImage = true;

	// this class applies Sobel Filter
	public static <T extends RealType<T>> void sobelFilter(
			final RandomAccessibleInterval< T > src, final RandomAccessibleInterval< T > dst, final int[] kernelDim){
		sobelFilter(Views.extendMirrorSingle(src), src, dst, kernelDim);
	}

	public static < T extends RealType<T>> void sobelFilter(
			final RandomAccessible< T > infSrc, final Interval srcInterval, final RandomAccessibleInterval< T > dst, final int[] kernelDim){

		final RandomAccessibleInterval<T> src = Views.interval(infSrc, srcInterval);
		sobelFilter(src, dst, kernelDim); // TODO: check if this one is corrects

	}

	public static <T extends RealType<  T > > void getKernel(final int dim, final RandomAccessibleInterval<T> kernel){
		if (dim == 2){

		}
	}

	// set the values for the Sobel kernel
	// the values are set only for one axis
	// to get other kernels one should rotate 
	// the initial kernel 
	// the normalization is provided but 
	// it is not crucial
	public static float[] getKernelValues(int dim, int kType){
		float[] k = new float[9];
		if (dim == 2){
			if (kType == 0){
				k[0] = 0;
				k[1] = 1; 
				k[2] = 2;
				k[3] = -1;
				k[4] = 0;
				k[5] = 1;
				k[6] = -2;
				k[7] = 1;
				k[8] = 0;
			}
			else{
				k[0] = 1;
				k[1] = 2; 
				k[2] = 1;
				k[3] = 0;
				k[4] = 0;
				k[5] = 0;
				k[6] = -1;
				k[7] = -2;
				k[8] = -1;
			}

			for (int i = 0; i < k.length; ++i) {
				k[i] /= 4;
			}
			return k;			
		} 
		if (dim == 3){
			float[] kernel = new float[]{
					1, 2, 1, 	2, 4, 2, 	1, 2, 1,
					0, 0, 0, 	0, 0, 0,	0, 0, 0,
					-1,-2,-1,  -2,-4,-2,   -1,-2,-1};
			for (int i = 0; i < kernel.length; ++i) {
				kernel[i] /= 16;
			}

			return kernel;


		}
		// otherwise we don't know what to do
		// exit program
		System.out.println("dimensionality is wrong");
		System.exit(1);
		return new float[]{};
	}

	// set the values for the Sobel kernel
	// the values are set only for one axis
	// to get other kernels one should rotate 
	// the initial kernel 
	// the normalization is provided but 
	// it is not crucial
	public static void getKernelValues(int dim, int kType, float[] kValues){
		if (dim == 2){
			if (kType == 0){ // pi/4
				kValues[0] = 0;
				kValues[1] = 1; 
				kValues[2] = 2;
				kValues[3] = -1;
				kValues[4] = 0;
				kValues[5] = 1;
				kValues[6] = -2;
				kValues[7] = 1;
				kValues[8] = 0;
			}
			else{ // pi/2
				kValues[0] = 1;
				kValues[1] = 2; 
				kValues[2] = 1;
				kValues[3] = 0;
				kValues[4] = 0;
				kValues[5] = 0;
				kValues[6] = -1;
				kValues[7] = -2;
				kValues[8] = -1;
			}
			// normalization
			for (int i = 0; i < kValues.length; ++i) {
				kValues[i] /= 4;
			}		
		} 
		if (dim == 3){
			float[] kernel = new float[]{
					1,2,1, 		2,4,2, 		1,2,1,
					0,0,0, 		0,0,0,		0,0,0,
					-1,-2,-1, -2,-4,-2, -1,-2,-1};
			for (int i = 0; i < kValues.length; ++i) {
				kValues[i] = kernel[i];
				kValues[i] /= 16;
			}
		}
		// otherwise we don't know what to do
		// exit program
		if (dim != 2 && dim != 3){
			System.out.println("dimensionality is wrong");
			System.exit(1);
		}
		// return new float[]{};
	}

	public static  < T extends RealType< T >> void applyMedianFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2){
		final int n = img.numDimensions();
		final int[] mD = new int[n];
		// here we consider symmetric kernel
		for (int d = 0; d < n; ++d)
			mD[d] = 3;
		MedianFilter.medianFilter(img, img2, mD);
	}

	public static < T extends RealType< T >> void  applyGaussianFilter(RandomAccessibleInterval <T> img, RandomAccessibleInterval <T> img2) throws IncompatibleTypeException{
		double[] sigma = new double[ img.numDimensions() ];
		for (int d = 0; d < img.numDimensions(); ++d)
			sigma[d] = 2; // size of the radius
		Gauss3.gauss(sigma, Views.extendMirrorSingle(img), img2);
	}

	public static Img <FloatType> setKernel(int n, int kType){
		// int kType = 0; // defines the direction of the sobel filter
		// fill in the kernel with proper values
		float[] kernelValues = new float [(int) Math.pow(3, n)]; 

		//System.out.println(kernelValues.length);

		getKernelValues(n, kType, kernelValues); 

		long [] kernelDimensions = new long[n];

		for (int d = 0; d < n; d++) 
			kernelDimensions[d] = 3; // the value is always set to 3 = size of the stencil

		// TODO: think which directions you need for a 3D case 

		// convert kernel to image 
		return ArrayImgs.floats( kernelValues, kernelDimensions);
	}

	// naive copy function
	public static <T extends RealType<T>> void copy(IterableInterval<T> in, RandomAccessibleInterval<T> out){

		Cursor<T> cursor = in.cursor();
		RandomAccess<T> randomAccess = out.randomAccess();

		while (cursor.hasNext()){
			cursor.fwd();
			randomAccess.setPosition(cursor);
			randomAccess.get().set(cursor.get().copy());		
		}
	}


	public static <T extends RealType<T>, U extends RealType<U>> void copyBitToFloat(IterableInterval<U> in, RandomAccessibleInterval<T> out){

		Cursor<U> cursor = in.cursor();
		RandomAccess<T> randomAccess = out.randomAccess();

		while (cursor.hasNext()){
			cursor.fwd();
			randomAccess.setPosition(cursor);
			randomAccess.get().setReal(cursor.get().getRealFloat());		
		}
	}

	public static <T extends RealType<T>> void applySobelFilter(Img<T> img, RandomAccessibleInterval<T> img2, RandomAccessibleInterval<T> kernel, T minValue, T maxValue){
		applySobelFilter(img, img2, kernel, minValue, maxValue, img.factory());
	}

	public static <T extends RealType<T>> void applySobelFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2, RandomAccessibleInterval<T> kernel, T minValue, T maxValue, ImgFactory<T> factory){

		applySobelFilter(img, img2, kernel, minValue, maxValue,  factory.create(img, minValue ));
	}

	public static <T extends RealType<T>> void applySobelFilter(RandomAccessibleInterval<T> img, RandomAccessibleInterval<T> img2, RandomAccessibleInterval<T> kernel, T minValue, T maxValue, RandomAccessibleInterval<T> tmp){


		// int[] mapping = new int []{0, 0, 2, 1, 1, 2};

		// apply Sobel filter # of dimensions times
		for (int d = 0; d < img.numDimensions(); d++) {
			// if (d == 2) continue;

			copy(Views.iterable(img), tmp);
			// convolve with each kernel
			// DEBUG: This one only works for 2D
			if (img.numDimensions() == 2)
				new FFTConvolution<T>(tmp, Views.rotate(kernel, 0, d), new ArrayImgFactory<ComplexFloatType>()).convolve();
			else{
				// for 3D you want 0-1 0-2 1-2 
				int[] mapping = new int []{0, 0, 1, 1, 2, 2};
				new FFTConvolution<T>(tmp, Views.rotate(kernel, mapping[d], mapping[d + img.numDimensions()]), new ArrayImgFactory<ComplexFloatType>()).convolve();
				// ImageJFunctions.show(Views.rotate(kernel, mapping[d], mapping[d + img.numDimensions()]));
			}

			// here we copy data to the destination image
			Cursor<T> tmpCursor = Views.iterable(tmp).cursor();
			RandomAccess<T> dstRandomAccess = img2.randomAccess();

			while (tmpCursor.hasNext()){
				tmpCursor.fwd();
				dstRandomAccess.setPosition(tmpCursor);
				T val = tmpCursor.get().copy();
				//System.out.println(val);
				val.mul(val);
				//System.out.println(val);
				dstRandomAccess.get().add(val);
			}
		}

		Cursor<T> dstCursor = Views.iterable(img2).cursor();
		while (dstCursor.hasNext()){
			dstCursor.fwd();	
			dstCursor.get().setReal(Math.sqrt(dstCursor.get().getRealFloat()));
		}

		Normalize.normalize(Views.iterable(img2), minValue, maxValue);

		// ImageJFunctions.show(img2);

	}

	public static <T extends RealType<T>> void distanceTransformKDTree(PointSampleList<T> worm, final RandomAccessibleInterval< BitType > in, final RandomAccessibleInterval< FloatType > out, PointSampleList<BitType> wormOutline){

		final PointSampleList< BitType > oneList = new PointSampleList< BitType >( in.numDimensions() );
		Cursor<T> wormCursor = worm.cursor(); 
		final RandomAccess<BitType> rIn = in.randomAccess();

		while(wormCursor.hasNext()){
			wormCursor.fwd();
			rIn.setPosition(wormCursor);
			if (rIn.get().get()) {
				oneList.add(new Point(wormCursor), new BitType(true));
			}
		}

		final KDTree<BitType> tree = new KDTree<BitType>(oneList);
		final NearestNeighborSearchOnKDTree< BitType > search = new NearestNeighborSearchOnKDTree< BitType >(tree);


		final RandomAccess<FloatType> rOut = out.randomAccess();
		wormCursor.reset();

		while(wormCursor.hasNext()){
			wormCursor.fwd();
			rIn.setPosition(wormCursor);
			rOut.setPosition(wormCursor);
			// System.out.println(wormCursor.get().getRealFloat());

			if(rIn.get().getInteger() == 0){  
				search.search(rIn);

				rOut.get().setReal(search.getDistance());
			}
			else{
				rOut.get().setZero();
			}

			// TODO: check this one: adjust the brightness of the boundary, if necessary
			if(rIn.get().getInteger() == 1){
				rOut.get().set(80);
				// TODO: this is the additional parameter! necessary for tail detection
				wormOutline.add(new Point(wormCursor), new BitType(true));
			}

		}

	}

	// process the worm and return the distance transform of it 
	// fix the declaration
	public static <T extends RealType<T> & NativeType<T>, U extends RealType<U>> void processWorm(RandomAccessibleInterval<T> initialImg, RandomAccessibleInterval<T> filterImg, RandomAccessibleInterval<T> edgeImg, RandomAccessibleInterval<BitType> thresholdImg, Img <T> distanceImg,
			T minValue, T maxValue,
			T tVal,
			U minVal, U maxVal,
			PointSampleList<BitType> wormOutline
			) throws IncompatibleTypeException{		

		Normalize.normalize(Views.iterable(initialImg), minValue, maxValue);
		System.out.println("Normalization: done!");
		final int n = initialImg.numDimensions();
		// ImageJFunctions.show(initialImg);	
		applyGaussianFilter(initialImg, filterImg);	
		System.out.println("Gaussian filter: done!");
		// applyMedianFilter(initialImg, filterImg);
		// System.out.println("Median filter: done!");


		Img<T> kernel = (Img<T>) setKernel(n, 0);

		final Img<T> edgeTmpImg = new ArrayImgFactory<T>().create(edgeImg, minValue);

		applySobelFilter(filterImg, edgeTmpImg, kernel, minValue, maxValue, new ArrayImgFactory<T>());
		System.out.println("Sobel filter: done!");
		if (debug){
			ImageJFunctions.show(filterImg).setTitle("DEBUG: filterImg: in function");	
			ImageJFunctions.show(edgeTmpImg).setTitle("DEBUG: edgeTmpImg: in function");

		}


		Thresholding.threshold(edgeTmpImg, thresholdImg, tVal);
		System.out.println("Thresholding: done!");

		if (debug)
			ImageJFunctions.show(thresholdImg).setTitle("DEBUG: thresholdImg: in function");

		// if (true) return;

		copyBitToFloat(Views.iterable(thresholdImg), edgeTmpImg);
		Normalize.normalize(Views.iterable(edgeTmpImg), minValue, maxValue);
		System.out.println("Normalization: done!");
		//if (debug)
		// 	ImageJFunctions.show(edgeTmpImg).setTitle("DEBUG: edgeImg: in function");



		kernel = (Img<T>) setKernel(n, 1);
		applySobelFilter(edgeTmpImg, edgeImg, kernel, minValue, maxValue, new ArrayImgFactory<T>());
		System.out.println("Sobel filter: done!");
		ImageJFunctions.show(edgeImg).setTitle("edgeImg: in function");

		if (false) return;


		final ImgLabeling<Integer, IntType> labeling = new ImgLabeling<Integer, IntType>(new ArrayImgFactory<IntType>().create(thresholdImg, new IntType())); 
		BoundingBox.setLabeling(thresholdImg, labeling);
		System.out.println("Labeling: done!");
		if (debug)
			ImageJFunctions.show(labeling.getIndexImg()).setTitle("DEBUG: labeling: in function");

		PointSampleList<T> worm = new PointSampleList<T>(initialImg.numDimensions());
		BoundingBox.setPointSampleList(labeling, (RandomAccessible<T>)initialImg, worm);
		Thresholding.threshold(edgeImg, thresholdImg, tVal);
		System.out.println("Thresholding: done!");

		if (debug)
			ImageJFunctions.show(thresholdImg).setTitle("DEBUG: thresholdImg: in function");

		// PointSampleList<BitType> wormOutline = new PointSampleList<BitType>(initialImg.numDimensions());
		distanceTransformKDTree(worm, thresholdImg, (RandomAccessibleInterval<FloatType>)distanceImg, wormOutline);
		System.out.println("Distance transform: done!");
		Normalize.normalize(distanceImg, minValue, maxValue);
		System.out.println("Normalization: done!");
		if (debug)
			ImageJFunctions.show(distanceImg).setTitle("DEBUG: distanceImg: in function");	
	}


	public static <T extends RealType<T>> void subtractImg(RandomAccessibleInterval<T> first, RandomAccessibleInterval<T> second){
		Cursor<T> f = Views.iterable(first).cursor();
		RandomAccess<T> s = second.randomAccess();

		while(f.hasNext()){
			f.fwd();
			s.setPosition(f);
			f.get().sub(s.get());
		}

	}

	public static <T extends RealType<T>>void main(String[] args) throws IncompatibleTypeException{
		new ImageJ();

		// to have memory monitor
		Executors.newFixedThreadPool( 1 ).execute(new Runnable()
		{
			@Override
			public void run() {new MemoryMonitor();}
		});


		System.out.println("Max memory: " + java.lang.Runtime.getRuntime().maxMemory());

		for (int num = 0; num < 1; ++num){

			// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
			// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one-1.tif");
			File file = new File("../Documents/Useful/initial_worms_pics/stack/1001-yellow-one0" + num +".tif");

			//Auto_Threshold.
			//IJ.run(new ImagePlus(file.getAbsolutePath()), "Auto Threshold", "method=Default ignore_white white stack");
			//ImageJFunctions.convertFloat(imp)
			// File file = new File("../Documents/Useful/initial_worms_pics/1003_first_3D.tif");
			// File file = new File("../Documents/Useful/initial_worms_pics/1001-red-one-1.tif");
			// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");

			//ImagePlusImgFactory<FloatType> >> directly compatible with ImagePlus
			Img<FloatType> initialImg = ImgLib2Util.openAs32Bit(file);		
			//Img<FloatType> im2 = initialImg.copy();
			
			// wrap them to ImgLib1 for old phase correlation code
//			Image<mpicbg.imglib.type.numeric.real.FloatType> i1 = ImgLib2.wrapFloatToImgLib1( initialImg );
//			Image<mpicbg.imglib.type.numeric.real.FloatType> i2 = ImgLib2.wrapFloatToImgLib1( im2 );
//			PairWiseStitchingResult r = PairWiseStitchingImgLib.computePhaseCorrelation(i1, i2, 5, false);
//			System.out.println( r.getCrossCorrelation() );
//			for ( int d = 0; d < initialImg.numDimensions(); ++ d )
//				System.out.println( r.getOffset( d )  );
//			
//			SimpleMultiThreading.threadHaltUnClean();
			
			ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >();
			ImgFactory< BitType > bitFactory = new ArrayImgFactory< BitType >();

//			ImagePlus i = ImageJFunctions.wrapFloat(initialImg, "bdfs");
//			FileSaver saver = new FileSaver(i);
//			saver.saveAsTiffStack("out.tiff" ); //3d
//			saver.saveAsTiff("out.tiff" ); //2d
						
			// lets have a look if this preprocessing works! 
			Img < FloatType > preprocessedImg = imgFactory.create( initialImg, new FloatType());

			// MedianFilter.medianFilter(initialImg, preprocessedImg, new int []{51, 51});

			final int n = initialImg.numDimensions();
			final long[] min = new long[ n ];
			final long[] max = new long[ n ];

			for ( int d = 0; d < n; ++d )
			{
				min[ d ] = -25;
				max[ d ] = 25;
			}

			System.out.println("before mean!");
			stephan.MeanFilter.filterRect(initialImg, preprocessedImg, new FinalInterval(min, max));
			System.out.println("mean done!");
			// subtractImg(initialImg, preprocessedImg);

			subtractImg(preprocessedImg, initialImg);
			
			Img< FloatType > filterImg = imgFactory.create( initialImg, new FloatType() );
			Img< BitType > thresholdImg = bitFactory.create( initialImg, new BitType() );
			Img< FloatType > edgeImg = imgFactory.create( initialImg, new FloatType() );
			Img< FloatType > distanceImg = imgFactory.create( initialImg, new FloatType() );


			// applyGaussianFilter(initialImg, filterImg);	
			// System.out.println("Gaussian filter: done!");


			PointSampleList<BitType> wormOutline = new PointSampleList<BitType>(initialImg.numDimensions());
			processWorm(preprocessedImg, filterImg, edgeImg, thresholdImg, distanceImg,
					new FloatType((float) 0), new FloatType((float) 255), new FloatType((float) 110.0), // 52
					new BitType(false), new BitType(true), wormOutline);

			if (showImage){
				ImageJFunctions.show(initialImg).setTitle("initialImg");
				ImageJFunctions.show(filterImg).setTitle("filterImg");
				ImageJFunctions.show(edgeImg).setTitle("edgeImg");
				ImageJFunctions.show(thresholdImg).setTitle("thresholdImg");
				ImageJFunctions.show(distanceImg).setTitle("distanceImg");
			}

			// ImageJFunctions.show(distanceImg).setTitle("distanceImg " + num);

			// debugging + sandbox 

			//		Img< FloatType > tmpImg = imgFactory.create( initialImg, new FloatType() );
			//		Cursor<BitType> wormOulineCursor = wormOutline.cursor(); 
			//		RandomAccess<FloatType> ra = tmpImg.randomAccess();
			//		long idx = 0;
			//		while(wormOulineCursor.hasNext()){
			//			wormOulineCursor.fwd();
			//			ra.setPosition(wormOulineCursor);
			//			ra.get().set(idx++);
			//
			//			for (int d = 0; d < wormOutline.numDimensions(); ++d){
			//				System.out.print(wormOulineCursor.getIntPosition(d) + " ");
			//			}
			//			System.out.println();
			//
			//		}
			//		ImageJFunctions.show(tmpImg).setTitle("tmpImg");
			//		// wormOutline.
			//		PointSampleList<FloatType> endPoints = new PointSampleList<FloatType>(initialImg.numDimensions());
			//
			//		CentralCurveDetection.detectEnds(endPoints, wormOutline, initialImg, new FloatType(255));

			
			
			System.out.println("DONE!");
		}
	}
}
