package randomAlgorithms;

import java.io.File;

import net.imglib2.RandomAccessible;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class Superpixel {

	public static <T extends NumericType<T>> void runSuperpixel(RandomAccessible<T> input, RandomAccessible<T> putput){
		
	}
	
	public static void main(String[] args){
		File file = new File("../resources/LENNA.JPG");
		final Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final Img<FloatType> dst = img.factory().create(img, img.firstElement());
		runSuperpixel(img, dst);
	}
}
