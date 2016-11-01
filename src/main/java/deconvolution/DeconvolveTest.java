package deconvolution;

import ij.ImagePlus;
import ij.ImageStack;
import klim.DeconvolutionTest;

import java.util.Date;

import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.wrapper.ImgLib1;
import mpicbg.imglib.wrapper.ImgLib2;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.FloatType;

public class DeconvolveTest
{
	public static LRFFT_Test createInput( final Img< net.imglib2.type.numeric.real.FloatType > im, final Img< net.imglib2.type.numeric.real.FloatType > psf )
	{
		return new LRFFT_Test( ImgLib2.wrapFloatToImgLib1( im.copy() ), ImgLib2.wrapFloatToImgLib1( psf.copy() ) );
	}

	public static void deconvolve( final LRFFT_Test deconvolutionData )
	{
		final LucyRichardson d = new LucyRichardson( deconvolutionData, 1.0e-4 );

		// will be used to show the result
		ImageStack stack = null;
		
		// value 
		double totalEnergy = AdjustInput.sumImage(d.getPsi())/(double)d.getPsi().getNumPixels();
		double currentEnergy = AdjustInput.sumImage(d.getPsi())/(double)d.getPsi().getNumPixels();
		// this coefficent will be used to adjust the total image intensity
		double ratio = totalEnergy/currentEnergy;

		// number of deconvolve iterations
		for ( int i = 0; i < 200; ++ i )
		{
			System.out.println( new Date( System.currentTimeMillis() ) + " " +  i );
			d.runIteration();

			// every, every 10th and every 100th 
			
			if ( i < 10 || ( i < 100 && i % 10 == 0 ) || ( i >= 100 && i % 100 == 0 ) )
			{
				ImagePlus imp = ImageJFunctions.copyToImagePlus( d.getPsi() );
	
				if ( stack == null )
					stack = new ImageStack( imp.getWidth(), imp.getHeight() );
		
				imp.setSlice(0);
				
				stack.addSlice( "iteration_"  + i, imp.getProcessor());
			}
			
			if ( i == 500 )
			{
				ImagePlus impStack = new ImagePlus( "decon", stack );
				impStack.show();
				// impStack.setRoi(15,39,615,274);
			}
			
			
			if ( i > 0 && ( i == 10 ||  i == 50 || i % 99 == 0 ) )
			{
				ImagePlus imp = ImageJFunctions.copyToImagePlus( d.getPsi() );
				imp.setTitle( "it=" + i );
				imp.resetDisplayRange();
				imp.show();
			}					
		}
		
		// ImagePlus impStack = new ImagePlus( "decon", stack );
		// impStack.show();
		// impStack.setRoi(15,39,615,274);
		
		System.out.println("image before: " + DeconvolutionTest.sumIntensitiesInDouble(ImgLib1.wrapArrayFloatToImgLib2(d.getPsi()))/(double)d.getPsi().getNumPixels());
		currentEnergy = AdjustInput.sumImage(d.getPsi())/(double)d.getPsi().getNumPixels();
		ratio = totalEnergy/currentEnergy;
		
		// adjust the Image intensities so that they sum up to 1
		AdjustInput.adjustImage(d.getPsi(), ratio);
		System.out.println("image after : " + DeconvolutionTest.sumIntensitiesInDouble(ImgLib1.wrapArrayFloatToImgLib2(d.getPsi()))/(double)d.getPsi().getNumPixels());
		
		ImageJFunctions.show( d.getPsi() );
	}
}
