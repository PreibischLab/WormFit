package klim;

import java.io.File;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane.SystemMenuBar;

import fiji.threshold.Auto_Threshold;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.Undo;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackConverter;
import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

public class MyAutoThreshold {

	public static <T extends RealType<T>> void myAutoThreshold(RandomAccessibleInterval<T> iImg, RandomAccessibleInterval<BitType> bImg){

		final int zDim = 2; // stands for z-axis
		final boolean doIstack = true; 		// is picture a stack 
		final boolean doIstackHistogram = false; 		// use histogram of the whole image
		// String myMethod = "Default"; 
		boolean noWhite = false;
		boolean noBlack = true; 
		// boolean doIwhite = false; // SetThreshold instead if Threshold (single images)
		// boolean doIset = false;  
		// boolean doIlog = false; 

		RandomAccess <T> rImg = iImg.randomAccess();
		rImg.setPosition(new long[]{0, 0, -1}); // initialize random access
		final long stackSize = iImg.max(zDim); 
		long pos = rImg.getLongPosition(zDim);
		
		T tVal = rImg.get(); 

		if (stackSize>1 &&( doIstack || doIstackHistogram) ) { //whole stack
			// slice by slice
			while (pos < stackSize){
				rImg.fwd(zDim);
				++pos;
				System.out.println(pos);
				
				int nBins = 256;
				int [] histogram = new int [nBins];
				int [] temp = new int [nBins];
				
				getHistogram(iImg, histogram);

				if (noBlack) histogram[0] = 0;
				if (noWhite) histogram[nBins - 1] = 0;

				// bracket the histogram to the range that holds data to make it quicker
				int minbin = -1, maxbin = -1;
				for (int i = 0; i < nBins; i++){
					if (histogram[i] > 0) maxbin = i;
				}
				for (int i = nBins - 1; i >= 0; i--){
					if (histogram[i] > 0) minbin = i;
				}

				int [] data2 = new int [(maxbin - minbin)+1];
				for (int i = minbin; i <= maxbin; i++){
					data2[i - minbin]= histogram[i];
				}

				tVal.setReal(minbin);
				
				IJDefault(data2, tVal);
				Thresholding.threshold(Views.hyperSlice(iImg, zDim, pos), Views.hyperSlice(bImg, zDim, pos), tVal);

			} 
		}
	}

	// from fiji plugin
	// to use it add tVal as a parameter minimal req to make this thing work
	public static <T extends RealType<T>> void IJDefault(int [] data, T tVal) {	
		// Original IJ implementation for compatibility.
		int level;
		int maxValue = data.length - 1;
		double result, sum1, sum2, sum3, sum4;
		
		boolean isLevelSet = false;
		
		int min = 0;
		while ((data[min]==0) && (min<maxValue))
			min++;
		int max = maxValue;
		while ((data[max]==0) && (max>0))
			max--;
		if (min>=max) {
			level = data.length/2;
			tVal.setReal(level + tVal.getRealFloat());
			isLevelSet = true;
		}
		
		if (!isLevelSet){
			int movingIndex = min;
			int inc = Math.max(max/40, 1);
			do {
				sum1=sum2=sum3=sum4=0.0;
				for (int i=min; i<=movingIndex; i++) {
					sum1 += i*data[i];
					sum2 += data[i];
				}
				for (int i=(movingIndex+1); i<=max; i++) {
					sum3 += i*data[i];
					sum4 += data[i];
				}			
				result = (sum1/sum2 + sum3/sum4)/2.0;
				movingIndex++;
			} while ((movingIndex+1)<=result && movingIndex<max-1);

			//.showProgress(1.0);
			level = (int)Math.round(result);

			// TODO:
			tVal.setReal(level + tVal.getRealFloat());
		}

	}

	// TODO: switch to longs?  
	public static <T extends RealType<T>> void getHistogram(RandomAccessibleInterval<T> img, int[] histogram){
		int nBins = histogram.length;
		
		Cursor<T> cursor = Views.iterable(img).cursor();
		
		while(cursor.hasNext()){
			cursor.fwd();
			int idx = (int)cursor.get().getRealFloat();
			if (idx >= nBins) 
				System.out.println("sasaj");
			++histogram[idx];		
		}
		
	}


	public static void main(String args[]){
		new ImageJ();
		// File file = new File("../Documents/Useful/initial_worms_pics/1001-yellow-one.tif");
		
		for (int i = 1; i <= 4; ++i){
			File file = new File("../Documents/Useful/initial_worms_pics/2016-04-20_N2wt_lea-1/100" + i +"_fool_one.tif");
			Img<FloatType> iImg = util.ImgLib2Util.openAs32Bit(file);
			Img<BitType> bImg = new ArrayImgFactory<BitType>().create(iImg, new BitType());
			Normalize.normalize(iImg, new FloatType(0), new FloatType(255));
			myAutoThreshold(iImg, bImg);
			ImageJFunctions.show(bImg);
		}
		


		// copy (.duplicate()) any RandomAccessibleInterval into an ImageJ ImagePlus
		// ImagePlus imp = ImageJFunctions.wrapFloat(iImg, "").duplicate();

		// autothreshold...
		// imp.show();

		// autoThresholdCall(iImg);

		// wrap the ImageJ ImagePlus into an Img< FloatType > (no copy)
		// but the ImgFactory is now an ImagePlusImgFactory< FloatType >
		// Img< FloatType > i2 = ImageJFunctions.wrap(imp);

	}
}
