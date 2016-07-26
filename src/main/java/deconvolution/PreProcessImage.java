package deconvolution;

import ij.ImageJ;

import java.io.File;

import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;

public class PreProcessImage
{
	public static Img< FloatType > cropConvolvedDros( final Img< FloatType > conv )
	{
		// crop a smaller version of the convolved image
		// 66 >> 274+66

		final Img< FloatType > croppedConv = new ArrayImgFactory< FloatType >().create( new int[]{ 615, 274 }, new FloatType() );

		final Cursor< FloatType  > c = croppedConv.cursor();
		final RandomAccess< FloatType > r = conv.randomAccess();

		while ( c.hasNext() )
		{
			c.fwd();
			r.setPosition( c );
			r.move( 66, 1 );
			c.get().set( r.get() );
		}

		return croppedConv;
	}
	
	public static void main( String[] args )
	{
		File f = new File( "." );
		for ( final String s : f.list() )
			System.out.println( s );

		new ImageJ();

		Img< FloatType > img = ImgLib2Util.openAs32Bit( new File( "src/main/resources/decon/dros-1-pb.tif" ) ); //psi_synthetic.tif" ) );
		AdjustInput.adjustImage( ImgLib2.wrapFloatToImgLib1( img ), LucyRichardson.minValue, 1 );
		Img< FloatType > convImg = img.copy();
		Img< FloatType > psf = ImgLib2Util.openAs32Bit( new File( "src/main/resources/decon/psf-1b.tif" ) );
		AdjustInput.normImage( ImgLib2.wrapFloatToImgLib1( psf ) );

		FFTConvolution< FloatType > conv = new FFTConvolution< FloatType >( convImg, psf );
		conv.convolve();

		//ImageJFunctions.show( img );
		//ImageJFunctions.show( convImg );
		//ImageJFunctions.show( cropConvolvedDros( img ) );
		ImageJFunctions.show( cropConvolvedDros( convImg ) ).setTitle( "input" );
		ImageJFunctions.show( psf ).setTitle( "psf" );

		DeconvolveTest.deconvolve( DeconvolveTest.createInput( cropConvolvedDros( convImg ), psf ) );
	}
}
