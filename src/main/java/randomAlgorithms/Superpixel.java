package randomAlgorithms;

import java.io.File;
import java.util.ArrayList;

import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class Superpixel {

	
	
	public static <T extends RealType<T>> void runSuperpixel(RandomAccessibleInterval<T> input, RandomAccessibleInterval<T> interval, RandomAccessibleInterval<T> output){
		final RandomAccessibleInterval<T> src = Views.interval(Views.extendMirrorSingle(input), interval);
		// grid step
		long S = 20;
		
		ArrayList<T> C = new ArrayList<>();
		
		
	}
	
	// black and white images
	public static double dc(double x, double y){
		return Math.abs(x - y);
	}
	
	public static double ds(long[] x, long[] y){
		double res = 0;
		
		for (int i = 0; i < x.length; ++i)
			res += (x[i] - y[i])*(x[i] - y[i]);
		
		res = Math.sqrt(res);
		return res;
	}
	
	// m parameter
	public static double D(double ds, double dc, double S, double m){
		return Math.sqrt(dc*dc + m*m*(ds/S)*(ds/S));
	}
	
	
	public static void main(String[] args){
		File file = new File("src/main/resources/LENNA.JPG");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		runSuperpixel(img, img, dst);
		
		System.out.println("Done!");
	}
}
