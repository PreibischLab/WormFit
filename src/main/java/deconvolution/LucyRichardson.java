package deconvolution;

import ij.IJ;

import java.util.ArrayList;

import mpicbg.imglib.algorithm.fft.FourierConvolution;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LucyRichardson implements Deconvolver
{
	final int numViews, numDimensions;
	final public static float minValue = 0.0001f;
    final float avg;
    final double lambda;
    
    // current iteration
    int i = 0;
    
	// the multi-view deconvolved image
	public Image<FloatType> psi;
	
	// the input data
	final LRInput views;
	final ArrayList<LRFFT_Test> data;
	final String name;

	
	public LucyRichardson( final LRInput views, final double lambda, String name )
	{
		this.name = name;
		this.data = views.getViews();
		this.views = views;
		this.numViews = data.size();
		this.numDimensions = data.get( 0 ).getImage().getNumDimensions();
		this.lambda = lambda;
		
		this.psi = data.get( 0 ).getImage().createNewImage( "psi (deconvolved image)" );
		
		this.avg = (float)AdjustInput.normAllImages( data );
		
		//System.out.println( "Average intensity in overlapping area: " + avg );        
		
		//
		// the real data image psi is initialized with the average 
		//	
		for ( final FloatType f : psi )
			f.set( avg );		
	}

	public LRInput getData() { return views; }
	public String getName() { return name; }
	public double getAvg() { return avg; }
	
	public Image<FloatType> getPsi() { return psi; }	
	public int getCurrentIteration() { return i; }
	
	public DeconvolveRuntimeStatistics runIteration() 
	{
		//if ( kernelsUpdate && i > 1 && i < 300  )
		//	updateKernelIteration( psi, data, lambda, minValue, i );
		//else
		DeconvolveRuntimeStatistics d = runIteration( psi, data, lambda > 0, lambda, minValue, i );
		i++;
		return d;
	}

	final private static DeconvolveRuntimeStatistics runIteration( final Image< FloatType> psi, final ArrayList<LRFFT_Test> data, 
			final boolean tikhonov, final double lambda, final float minValue, final int iteration )
	{
		final int numViews = data.size();

		// a different view at each iteration
		//final int view = iteration % numViews;
		
		double sumChange = 0;
		double maxChange = -1;

		for ( int view = 0; view < numViews; ++view )
		{
			final LRFFT_Test processingData = data.get( view );
			
			// convolve psi (current guess of the image) with the PSF of the current view
			final FourierConvolution<FloatType, FloatType> fftConvolution = processingData.getFFTConvolution1();
			fftConvolution.replaceImage( psi );
			//fftConvolution.replaceKernel( processingData.getKernel() );			
			fftConvolution.process();
			
			final Image<FloatType> psiBlurred = fftConvolution.getResult();
			
			// compute quotient img/psiBlurred
			final Cursor<FloatType> cursorImg = processingData.getImage().createCursor();
			final Cursor<FloatType> cursorPsiBlurred = psiBlurred.createCursor();
			
			while ( cursorImg.hasNext() )
			{
				cursorImg.fwd();
				cursorPsiBlurred.fwd();
				
				final float imgValue = cursorImg.getType().get();
				final float psiBlurredValue = cursorPsiBlurred.getType().get();
				
				if ( psiBlurredValue == 0 )
					cursorPsiBlurred.getType().set( 0 );
				else
					cursorPsiBlurred.getType().set( imgValue / psiBlurredValue );
			}
			
			cursorImg.close();
			cursorPsiBlurred.close();
			
			// blur the residuals image with the kernel
			final FourierConvolution<FloatType, FloatType> invFFConvolution = processingData.getFFTConvolution2();
			invFFConvolution.replaceImage( psiBlurred );
			invFFConvolution.process();
	
			// the result from the previous iteration
			final Cursor< FloatType > cursorPsi = psi.createCursor();
			final Cursor< FloatType > cursorIntegral = invFFConvolution.getResult().createCursor();
			final Cursor< FloatType > cursorWeight = processingData.getWeight().createCursor();
			
			while ( cursorPsi.hasNext() )
			{
				//cursorPsi.fwd();
				//cursorIntegral.fwd();
				final FloatType lastPsi = cursorPsi.next();
				final float lastPsiValue = lastPsi.get();				
				double value = lastPsiValue * cursorIntegral.next().get();
				
				///
				/// new
				///
				// compute the difference between old and new
				double changed = value - lastPsiValue;				
				
				// apply the apropriate amount
				changed *= cursorWeight.next().get();
				value = lastPsiValue + changed;
				
				if ( value > 0 )
				{
					//
					// perform Tikhonov regularization if desired
					//		
					if ( tikhonov )
						value = ( (float)( (Math.sqrt( 1.0 + 2.0*lambda*value ) - 1.0) / lambda ) );
				}
				else
				{
					value = minValue;
				}
				
				//
				// get the final value and some statistics
				//
				final float nextPsiValue;
				
				if ( Double.isNaN( value ) )
					nextPsiValue = (float)minValue;
				else
					nextPsiValue = (float)Math.max( minValue, value );
				
				final float change = Math.abs( lastPsiValue - nextPsiValue );				
				sumChange += change;
				maxChange = Math.max( maxChange, change );

				// store the new value
				lastPsi.set( (float)nextPsiValue );
			}
		}
		
		IJ.log("iteration: " + iteration + " --- sum change: " + sumChange + " --- max change per pixel: " + maxChange );
		
		return new DeconvolveRuntimeStatistics( maxChange, sumChange );
	}

	@Override
	public void setDebug(boolean debug) {
		// TODO Auto-generated method stub
		
	}
}
