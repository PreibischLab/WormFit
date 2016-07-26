package deconvolution;

import ij.ImagePlus;

import java.util.Date;

import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.img.Img;

public class DeconvolveTest
{
	public static LRFFT_Test createInput( final Img< net.imglib2.type.numeric.real.FloatType > im, final Img< net.imglib2.type.numeric.real.FloatType > psf )
	{
		return new LRFFT_Test( ImgLib2.wrapFloatToImgLib1( im.copy() ), ImgLib2.wrapFloatToImgLib1( psf.copy() ) );
	}

	public static void deconvolve( final LRFFT_Test deconvolutionData )
	{
		final Deconvolver d = new LucyRichardson( deconvolutionData.init(), 0.0, "independent (sequential)" );

		for ( int i = 0; i < 1000; ++ i )
		{
			System.out.println( new Date( System.currentTimeMillis() ) + " " +  i );
			d.runIteration();
			
			if ( i > 0 && ( i == 10 ||  i == 50 || i % 100 == 0 ) )
			{
				ImagePlus imp = ImageJFunctions.copyToImagePlus( d.getPsi() );
				imp.setTitle( "it=" + i );
				imp.resetDisplayRange();
				imp.show();
			}
		}
		
		ImageJFunctions.show( d.getPsi() );
	}
}
