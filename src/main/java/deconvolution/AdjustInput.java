package deconvolution;

import java.util.ArrayList;
import java.util.Random;

import org.uncommons.maths.random.PoissonGenerator;

import net.imglib2.util.Util;

import mpicbg.imglib.algorithm.gauss.GaussianConvolution;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.interpolation.Interpolator;
import mpicbg.imglib.interpolation.linear.LinearInterpolatorFactory;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.util.RealSum;

public class AdjustInput 
{	
	public static Random rnd = new Random( 14235235 );

	/**
	 * Norms an image so that the sum over all pixels is 1.
	 * 
	 * @param img - the {@link Image} to normalize
	 */
	final public static void normImage( final Image<FloatType> img )
	{
		final double sum = sumImage( img );	

		for ( final FloatType t : img )
			t.set( (float) ((double)t.get() / sum) );
	}
	
	/**
	 * @param img - the input {@link Image}
	 * @return - the sum of all pixels using {@link RealSum}
	 */
	final public static double sumImage( final Image<FloatType> img )
	{
		final RealSum sum = new RealSum();		

		for ( final FloatType t : img )
			sum.add( t.get() );

		return sum.getSum();
	}	

	public static double normImage( final LRFFT_Test data )
	{
		final ArrayList< LRFFT_Test > l = new ArrayList< LRFFT_Test >();
		l.add(  data );
		return normAllImages( l );
	}
	public static double normAllImages( final ArrayList<LRFFT_Test> data )
	{
		// the individual sums of the overlapping area
		//final double[] sums = new double[ data.size() ];
		
		final RealSum sum = new RealSum();
		// the number of overlapping pixels
		long count = 0;
		
		final ArrayList<Cursor<FloatType>> cursorsImage = new ArrayList<Cursor<FloatType>>();
		final ArrayList<Cursor<FloatType>> cursorsWeight = new ArrayList<Cursor<FloatType>>();
		
		for ( final LRFFT_Test fft : data )
		{
			cursorsImage.add( fft.getImage().createCursor() );
			if ( fft.getWeight() != null )
				cursorsWeight.add( fft.getWeight().createCursor() );
		}
		
		final Cursor<FloatType> cursor = cursorsImage.get( 0 );

		// sum overlapping area individually
		while ( cursor.hasNext() )
		{
			for ( final Cursor<FloatType> c : cursorsImage )
				c.fwd();
			
			for ( final Cursor<FloatType> c : cursorsWeight )
				c.fwd();
			
			// sum up individual intensities
			double sumLocal = 0;
			int countLocal = 0;
			
			for ( int i = 0; i < cursorsImage.size(); ++i )
			{
				if ( cursorsWeight.get( i ).getType().get() != 0 )
				{
					sumLocal += cursorsImage.get( i ).getType().get();
					countLocal++;
				}
			}

			// at least two overlap
			if ( countLocal > 1 )
			{
				sum.add( sumLocal );
				count += countLocal;
			}
		}

		if ( count == 0 )
			return 1;
		
		// compute the average sum
		final double avg = sum.getSum() / (double)count;
		
		// return the average intensity in the overlapping area
		return avg;
	}

	/**
	 * Adds additive gaussian noise: i = i + gauss(x, sigma)
	 * 
	 * @param amount - how many times sigma
	 * @return the signal-to-noise ratio (measured)
	 */
	public static double addGaussianNoise( final Image< FloatType > img, final Random rnd, final float sigma, boolean onlyPositive )
	{
		double[] values = new double[ (int)img.getNumPixels() ];
		int i = 0;
		
		for ( final FloatType f : img )
		{
			final double gauss = rnd.nextGaussian() * sigma;
			values[ i++ ] = gauss;
			
			float newValue = f.get() + (float)( gauss );
			
			if ( onlyPositive )
				newValue = Math.max( 0, newValue );
			
			f.set( newValue );
		}
		
		//System.out.print( "\t" + computeStDev( values ) );
		
		return 1;
	}
	
	public static double computeStDev( final double[] values )
	{
		double avg = 0;

		for ( final double d : values )
			avg += d;

		avg /= (double)values.length;

		double stDev = 0;
		
		for ( final double value : values )
			stDev += (avg - value) * (avg - value );
		
		stDev /= (double) values.length;
		
		return Math.sqrt( stDev );
	}

	/**
	 * Adds additive and multiplicative gaussian noise: i = i*gauss(x,sigma) + gauss(x, sigma)
	 * 
	 * @param amount - how many times sigma
	 * @return the signal-to-noise ratio (measured)
	 */
	public static double addGaussianNoiseAddMul( final Image< FloatType > img, final Random rnd, final float sigma, boolean onlyPositive )
	{
		for ( final FloatType f : img )
		{
			final float value = f.get();
			float newValue = value*(1+(float)( rnd.nextGaussian() * sigma/3 )) + (float)( Math.abs( rnd.nextGaussian() ) * sigma );
			
			if ( onlyPositive )
				newValue = Math.max( 0, newValue );
			
			f.set( newValue );
		}
		
		return 1;
	}

	public static void translate( final Image< FloatType > img, final float[] vector )
	{
		final Image< FloatType > tmp = img.clone();
		
		final LocalizableCursor< FloatType > cursor1 = img.createLocalizableCursor();		
		final Interpolator< FloatType > interpolator = tmp.createInterpolator( new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
		
		final int numDimensions = img.getNumDimensions();
		final float[] pos = new float[ numDimensions ];
		
		while ( cursor1.hasNext() )
		{
			cursor1.fwd();
			
			for( int d = 0; d < numDimensions; ++d )
				pos[ d ] = cursor1.getPosition( d ) - vector[ d ];
			
			interpolator.setPosition( pos );
			cursor1.getType().set( interpolator.getType() );
		}
		
		cursor1.close();
		interpolator.close();
	}
	
	/**
	 * Adjusts an image so that the minimal intensity is minValue and the average is average
	 * 
	 * @param image - the image to norm
	 * @param minValue - the minimal value
	 * @param targetAverage - the average that we want to have
	 * 
	 * @return by which number all intensities were multiplied
	 */
	public static double adjustImage( final Image<FloatType> image, final float minValue, final float targetAverage )
	{
		// first norm the image to an average of (targetAverage - minValue)
		final double avg = sumImage( image )/(double)image.getNumPixels();
		final double correction = ( targetAverage - minValue ) / avg;

		// correct 
		for ( final FloatType t : image )
			t.set( (float)( t.get() * correction ) );
			
		// now add minValue to all pixels
		for ( final FloatType t : image )
			t.set( t.get() + minValue );
		
		//System.out.println( sumImage( image )/(double)image.getNumPixels() );
		return correction;
	}

	public static void addPoissonNoise( final Image<FloatType> img, final double SNR )
	{
		// based on an average intensity of 5, a multiplicator of 1 corresponds to a SNR of 2.23 = sqrt( 5 );	
		final double mul = Math.pow( SNR / Math.sqrt( 5 ), 2 );
		
		final NumberGeneratorImage< FloatType> ng = new NumberGeneratorImage< FloatType>( img, mul );
		final PoissonGenerator pg = new PoissonGenerator( ng, rnd );
		
		for ( final FloatType v : img )
		{
			ng.fwd();
			//v.set( Math.max( LucyRichardsonSeq.minValue, pg.nextValue().floatValue() ) );
			v.set( pg.nextValue().floatValue() );
		}
	}

	public static int poissonValue( final double pixVal )
	{
		double L = Math.exp(-(pixVal));
		int k = 0;
		double p = 1;
		do
		{
			++k;
			p *= rnd.nextDouble();
		}
		while ( p >= L );
		
	    return ( k - 1 );
	}

	/**
	 * 
	 * @param img
	 * @param percentdeadpixels - number of dead pixels in percent
	 */
	public static void addDeadPixels( final Image<FloatType> img, final float percentdeadpixels, final float minValue )
	{
		for ( final FloatType f : img )
			if ( rnd.nextDouble() <= percentdeadpixels )
				f.set( minValue );		
	}
}
