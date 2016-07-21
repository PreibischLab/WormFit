package klim;

import java.io.File;

import ij.ImageJ;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.wrapper.ImgLib2;
import mpicbg.stitching.PairWiseStitchingImgLib;
import mpicbg.stitching.PairWiseStitchingResult;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;
import net.imglib2.type.numeric.real.FloatType;

// example of stitching via imglib
public class MyStitching {

	public static void myStitching(){
		File file = new File("../Documents/Useful/initial_worms_pics/2016-04-20_N2wt_lea-1/1001_fool_one.tif");
		Img<net.imglib2.type.numeric.real.FloatType> img = ImgLib2Util.openAs32Bit(file);
		Img<net.imglib2.type.numeric.real.FloatType> img2 = img.copy();
		
		Image<mpicbg.imglib.type.numeric.real.FloatType> i1 = ImgLib2.wrapFloatToImgLib1( img );
		Image<mpicbg.imglib.type.numeric.real.FloatType> i2 = ImgLib2.wrapFloatToImgLib1( img2 );
		PairWiseStitchingResult r = PairWiseStitchingImgLib.computePhaseCorrelation(i1, i2, 5, false);
		System.out.println( r.getCrossCorrelation() );
		for ( int d = 0; d < img.numDimensions(); ++ d )
			System.out.println( r.getOffset( d )  );	
	}
	
	
	public static void main(String args[]){
		new ImageJ();
		myStitching();
	}
}
