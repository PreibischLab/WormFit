package klim.deconvolution;

import java.io.File;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import mpicbg.util.RealSum;
import util.ImgLib2Util;

public class Preprocessing {
	public static void preprocessImage(String path, int i){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File(path + "/unprocessed/tif-" + i + "-c.tif"));//path + file + ".tif"));
		Img<FloatType> bg  = ImgLib2Util.openAs32Bit(new File(path + "/bg/tif-" + i + "-c-m.tif"));//path + "bg/" + file + "-m.tif"));
		
		double average = getImageAverage(bg);	
		System.out.println(average);
		
		subtractValue(bg, average);
		subtractImg(img, bg);	
		
		// ImageJFunctions.show(img);
		
		ImagePlus toWrite = ImageJFunctions.wrap(img, "").duplicate();
		new FileSaver(toWrite).saveAsTiffStack(path + "/processed/tif-" + i + "-c-m-p.tif");
	}
	
	
	public static void batchProcess(){
		new ImageJ();
		String path = "/home/milkyklim/Desktop/2017-04-27-beads";		
		int numFiles = 17;
		
		for (int i = 1; i <= numFiles; ++i){
			preprocessImage(path, i);
			System.out.println("image " + i + " processed!");
		}				
	}
	
		
	public static void subtractValue(Img<FloatType> bg, double value){			
		Cursor<FloatType> cursor = bg.cursor();
		
		while(cursor.hasNext()){
			cursor.fwd();
			float val = (float) (cursor.get().get() - value);
			cursor.get().set(val);
		}
	}
	
	// subtract value and cut negative values
	public static void subtractValueNonNegative(Img<FloatType> img, double value){			
		Cursor<FloatType> cursor = img.cursor();		
		while(cursor.hasNext()){
			cursor.fwd();
			float val = Math.max((float) (cursor.get().get() - value), 0);
			cursor.get().set(val);
		}
	}
	
	public static double getImageAverage(Img<FloatType> bg){
		// compute average over all pixels
		double sum = sumImage(bg);
		for (int d = 0; d < bg.numDimensions(); ++d)
			sum /= bg.dimension(d);
		return sum;
	}
	
	public static double sumImage(Img<FloatType> img)
	{
		final RealSum sum = new RealSum();		
		for ( final FloatType t : img )
			sum.add( t.get() );
		return sum.getSum();
	}
	
	public static void subtractImg(Img<FloatType> img, Img<FloatType> bg){
		Cursor<FloatType> cursor = img.cursor();
		RandomAccess<FloatType> ra = bg.randomAccess();
		
		while(cursor.hasNext()){
			cursor.fwd();		
			ra.setPosition(cursor);			
			float val = cursor.get().get() - ra.get().get();			
			cursor.get().set(val);
		}
	}
	
	public static void main(String [] args){
		Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("/home/milkyklim/Desktop/2017-04-27-beads/beads/averaged_bead.tif"));//path + file + ".tif"));
		subtractValueNonNegative(img,  0.024);
		
		ImagePlus toWrite = ImageJFunctions.wrap(img, "").duplicate();
		new FileSaver(toWrite).saveAsTiffStack("/home/milkyklim/Desktop/2017-04-27-beads/beads/averaged_bead_done.tif");
		
		// batchProcess();		
		System.out.println("DOGE!");
	}
	
}