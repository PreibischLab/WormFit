package deconvolution;

import ij.IJ;

import java.util.ArrayList;

import mpicbg.imglib.algorithm.fft.FourierConvolution;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LucyRichardsonSeq implements Deconvolver
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
	final boolean kernelsUpdate;
	
	public LucyRichardsonSeq( final LRInput views, final boolean kernelsUpdate, final double lambda, String name )
	{
		this.name = name;
		this.kernelsUpdate = kernelsUpdate;
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

	public boolean kernelsUpdated() { return kernelsUpdate; }
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

	final private static void updateKernelIteration( final Image< FloatType> psi, final ArrayList<LRFFT_Test> data, 
			final double lambda, final float minValue, final int iteration )
	{		
		final int numViews = data.size();

		double sumChange = 0;
		double maxChange = -1;

		final Image< FloatType > psiNew = psi;// LRFFT.computeExponentialKernel( psi, numViews );
		//LRInput.normImage( psiNew );
		
		//int view = k;
		for ( int view = 0; view < numViews; ++view )
		{
			final LRFFT_Test processingData = data.get( view );
			
			// convolve kernel with psi (current guess of the image)
			// it does not matter which way around, this way we have no boundary artifacts (which we get when we convolve the kernel and use mirroring)
			final FourierConvolution<FloatType, FloatType> fftConvolution = new FourierConvolution<FloatType, FloatType>( psi, processingData.getKernel1() );
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
				
				cursorPsiBlurred.getType().set( imgValue / psiBlurredValue );
			}
	
			//ImageJFunctions.show( psiBlurred );
			
			cursorImg.close();
			cursorPsiBlurred.close();
			
			// blur the residuals image with the kernel
			final FourierConvolution<FloatType, FloatType> fftConvolution2 = new FourierConvolution<FloatType, FloatType>( psiBlurred, LRFFT_Test.computeInvertedKernel( psiNew ) );
			fftConvolution2.process();
	
			//ImageJFunctions.show( fftConvolution2.getResult() );
						
			// the result from the previous iteration
			final Cursor< FloatType > cursorKernel1 = processingData.getKernel1().createCursor();
			final Cursor< FloatType > cursorIntegral = fftConvolution2.getResult().createCursor();
			
			//LRInput.normImage( fftConvolution2.getResult() );
			
			while ( cursorKernel1.hasNext() )
			{
				cursorKernel1.fwd();
				cursorIntegral.fwd();
				final float lastKernelValue = cursorKernel1.getType().get();
				
				double value = lastKernelValue * cursorIntegral.getType().get();
	
				if ( value > 0 )
				{
					//
					// perform Tikhonov regularization if desired
					//		
					if ( lambda > 0 )
						value = ( (float)( (Math.sqrt( 1.0 + 2.0*lambda*value ) - 1.0) / lambda ) );
				}
				else
				{
					value = minValue;
				}
				//
				// get the final value and some statistics
				//
				final float nextKernelValue;
				
				if ( Double.isNaN( value ) )
					nextKernelValue = (float)minValue;
				else
					nextKernelValue = (float)Math.max( minValue, value );
				
				final float change = Math.abs( lastKernelValue - nextKernelValue );				
				sumChange += change;
				maxChange = Math.max( maxChange, change );
				
				// store the new value
				cursorKernel1.getType().set( (float)nextKernelValue );
			}
			
			//
			// update the other fourier convolutions with the new kernel
			//
			
			// replace the fft for kernel1, the kernel itself is already updated
			AdjustInput.normImage( processingData.getKernel1() );

			//ImageJFunctions.show( processingData.getKernel1() );
			//SimpleMultiThreading.threadHaltUnClean();

			processingData.getFFTConvolution1().replaceKernel( processingData.getKernel1() );
			
			// compute the exponential kernel
			final Image< FloatType > exponentialKernel = LRFFT_Test.computeExponentialKernel( processingData.getKernel1(), numViews );
			
			// norm the exponential kernel
			AdjustInput.normImage( exponentialKernel );
			
			// replace it in the second fourier convolution
			final Image< FloatType > kernel2 = LRFFT_Test.computeInvertedKernel( processingData.getKernel2() );
			processingData.getFFTConvolution2().replaceKernel( kernel2 );
			
			// update the "original kernel2"
			final Cursor< FloatType > c = processingData.getKernel2().createCursor();
			
			for ( final FloatType f : kernel2 )
			{
				c.fwd();
				c.getType().set( f );
			}
		}
		
		IJ.log("------------------------------------------------");
		IJ.log(" Kernel-Iteration: " + iteration );
		IJ.log(" Kernel-Change: " + sumChange );
		IJ.log(" Kernel-Max Change per Pixel: " + maxChange );
		IJ.log("------------------------------------------------");			
	}

	@Override
	public void setDebug(boolean debug) {
		// TODO Auto-generated method stub
		
	}
}
