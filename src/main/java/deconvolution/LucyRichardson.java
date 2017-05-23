package deconvolution;

import net.imglib2.multithreading.SimpleMultiThreading;
import ij.IJ;
import mpicbg.imglib.algorithm.fft.FourierConvolution;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LucyRichardson
{
	final public static float minValue = 0.0001f;
    //final float avg;
    final double lambda;
    
    // current iteration
    int i = 0;
    
	// the multi-view deconvolved image
	public Image<FloatType> psi;
	
	// the input data
	final LRFFT_Test data;

	final int n;
	final int[] psfHalf;

	private final FourierConvolution<FloatType, FloatType> fftConvolution1, fftConvolution2;

	public LucyRichardson( final LRFFT_Test data, final double lambda )
	{
		this.data = data;
		this.lambda = lambda;

		this.n = data.getImage().getNumDimensions();
		this.psfHalf = new int[ n ];

		for ( int d = 0; d < n; ++d )
			psfHalf[ d ] = data.getKernel1().getDimension( d ) / 2;

		// psi is 1/2 psf bigger in all dimensions
		this.psi = createPSI( data.getImage(), data.getKernel1() );

		//this.avg = (float)AdjustInput.normImage( data );
		//System.out.println( "Average intensity in overlapping area: " + avg );

		//
		// the real data image psi is initialized with the average
		//
		final LocalizableCursor< FloatType > cursor = psi.createLocalizableCursor();
		final LocalizableByDimCursor< FloatType > randomAccess = data.getImage().createLocalizableByDimCursor();
		final int[] location = new int[ n ];
		final int[] size = data.getImage().getDimensions();

		while ( cursor.hasNext() )
		{
			cursor.next().get();
			cursor.getPosition( location );

			boolean isInside = true;

			for ( int d = 0; d < n; ++d )
			{
				location[ d ] -= psfHalf[ d ];
				if ( location[ d ] < 0 || location[ d ] >= size[ d ] )
				{
					isInside = false;
					break;
				}
			}

			if ( isInside )
			{
				randomAccess.setPosition( location );
				cursor.getType().set( randomAccess.getType().get() );
			}
			else
				cursor.getType().set( LucyRichardson.minValue );
			
			cursor.getType().set( 1.0f ); // ???
		}

		//for ( final FloatType f : psi )
		//	f.set( avg );

		// set FFT convolutions right
		this.fftConvolution1 = new FourierConvolution<FloatType, FloatType>( psi, data.getKernel1() );
		this.fftConvolution1.setNumThreads();
		this.fftConvolution1.setKeepImgFFT( false );
		this.fftConvolution1.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>( new FloatType( LucyRichardson.minValue ) ) ); // values outside are minimal values, but should be irrelevant

		this.fftConvolution2 = new FourierConvolution<FloatType, FloatType>( psi, data.getKernel2() );
		this.fftConvolution2.setNumThreads();
		this.fftConvolution2.setKeepImgFFT( false );
		this.fftConvolution2.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>( new FloatType( 1.0f ) ) ); // ratio is 1

		// blur psi twice as the starting point
		//fftConvolution1.process();
		//fftConvolution1.replaceImage( fftConvolution1.getResult() );
		//this.psi = fftConvolution1.getResult();
	}

	public LRFFT_Test getData() { return data; }

	public Image<FloatType> getPsi() { return psi; }
	public int getCurrentIteration() { return i; }

	public DeconvolveRuntimeStatistics runIteration() 
	{
		DeconvolveRuntimeStatistics d = runIteration( psi, data, fftConvolution1, fftConvolution2, n, psfHalf, lambda > 0, lambda, minValue, i );
		i++;
		return d;
	}

	public static boolean debug = false;

	final private static DeconvolveRuntimeStatistics runIteration(
			final Image< FloatType> psi, final LRFFT_Test processingData,
			final FourierConvolution<FloatType, FloatType> fftConvolution1,
			final FourierConvolution<FloatType, FloatType> fftConvolution2,
			final int n, final int[] psfHalf,
			final boolean tikhonov, final double lambda, final float minValue, final int iteration )
	{
		double sumChange = 0;
		double maxChange = -1;

		if ( debug )
			ImageJFunctions.show( psi.clone() ).setTitle( "psi_" + iteration );

		//
		// convolve psi (current guess of the image) with the PSF of the current view (outofbounds = MIN_VALUE!)
		//
		fftConvolution1.replaceImage( psi );
		fftConvolution1.process();
		final Image<FloatType> psiBlurred = fftConvolution1.getResult();

		if ( debug )
			ImageJFunctions.show( psiBlurred.clone() ).setTitle( "psiBlurred_" + iteration );

		//
		// compute ratio img/psiBlurred (only for the part where we observed the image, the rest we set to 1.0 - the ratio saying that everything is alright, this will be only influenced later during convolution)
		//
		final LocalizableCursor<FloatType> cursorPsiBlurred = psiBlurred.createLocalizableCursor();
		final LocalizableByDimCursor< FloatType > randomAccess = processingData.getImage().createLocalizableByDimCursor();
		final int[] location = new int[ n ];
		final int[] size = processingData.getImage().getDimensions();

		while ( cursorPsiBlurred.hasNext() )
		{
			final float psiBlurredValue = cursorPsiBlurred.next().get();
			cursorPsiBlurred.getPosition( location );

			boolean isInside = true;

			for ( int d = 0; d < n; ++d )
			{
				location[ d ] -= psfHalf[ d ];
				if ( location[ d ] < 0 || location[ d ] >= size[ d ] )
				{
					isInside = false;
					break;
				}
			}

			if ( isInside )
			{
				randomAccess.setPosition( location );
				final float imgValue = randomAccess.getType().get();
				
				if ( psiBlurredValue <= 0 )
					cursorPsiBlurred.getType().set( 1.0f );
				else
					cursorPsiBlurred.getType().set( imgValue / psiBlurredValue );
			}
			else
			{
				cursorPsiBlurred.getType().set( 1.0f );
			}
		}

		randomAccess.close();
		cursorPsiBlurred.close();

		if ( debug ) // || iteration == 0 )
			ImageJFunctions.show( psiBlurred.clone() ).setTitle( "ratio_" + iteration );

		//
		// blur the residuals image with the kernel (outofbounds = 1, ratio!)
		//
		fftConvolution2.replaceImage( psiBlurred );
		fftConvolution2.process();

		if ( debug )
			ImageJFunctions.show( fftConvolution2.getResult().clone() ).setTitle( "ratioBlurred_" + iteration );

		// the result from the previous iteration
		final Cursor< FloatType > cursorPsi = psi.createCursor();
		final Cursor< FloatType > cursorIntegral = fftConvolution2.getResult().createCursor();
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

		if ( debug )
			ImageJFunctions.show( psi.clone() ).setTitle( "psiNew_" + iteration );

		IJ.log("iteration: " + iteration + " --- sum change: " + sumChange + " --- max change per pixel: " + maxChange );

		return new DeconvolveRuntimeStatistics( maxChange, sumChange );
	}

	/**
	 * psi is 1/2 psf [int] * 2 bigger in all dimensions
	 * 
	 * @param observed - the observed image
	 * @param psf - the psf
	 * @return the initial image
	 */
	protected static Image< FloatType > createPSI( final Image< FloatType > observed, final Image< FloatType > psf )
	{
		final int n = observed.getNumDimensions();
		final int[] psfHalf = new int[ n ];
		final int[] psiSize = new int[ n ];

		for ( int d = 0; d < n; ++d )
		{
			psfHalf[ d ] = psf.getDimension( d ) / 2;
			psiSize[ d ] = observed.getDimension( d ) + 2 * psfHalf[ d ];
		}

		return observed.createNewImage( psiSize, "psi (deconvolved image)" );
	}

}
