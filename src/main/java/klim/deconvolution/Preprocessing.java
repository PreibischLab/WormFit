package klim.deconvolution;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;

import deconvolution.AdjustInput;
import mpicbg.imglib.wrapper.ImgLib2;
import util.ImgLib2Util;

public class Preprocessing {
	public static void preprocessImage(){
		// grab initial image 
		// grab bg image 
		
		// make a z projection of the median filter and get an average value 
		// median - average (get zero mean)
		
		// initial image - (median - average)
		String path = "";
		String file = "";
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + file + ".tif"));
		Img<FloatType> bg  = ImgLib2Util.openAs32Bit(new File(path + "bg/" + file + "-m.tif"));
		
		subtractAverage(bg);
		subtractImgs(img, bg);
		
	}
	
	
	public static void subtractAverage(Img<FloatType> bg){		
		// compute average over all pixels
		double sum = AdjustInput.sumImage(ImgLib2.wrapFloatToImgLib1(bg));
		for (int d = 0; d < bg.numDimensions(); ++d)
			sum /= bg.dimension(d);
		
		Cursor<FloatType> cursor = bg.cursor();
		
		while(cursor.hasNext()){
			cursor.fwd();		
			cursor.get().sub(new FloatType((float)sum));
		}
	}
	
	public static void subtractImgs(Img<FloatType> img, Img<FloatType> bg){
		Cursor<FloatType> cursor = img.cursor();
		RandomAccess<FloatType> ra = bg.randomAccess();
		
		while(cursor.hasNext()){
			cursor.fwd();		
			ra.setPosition(cursor);		
			cursor.get().sub(ra.get());
		}
	}
	
	public static void main(String [] args){
		preprocessImage();
	}
	
}