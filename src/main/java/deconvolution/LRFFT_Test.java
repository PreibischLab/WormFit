package deconvolution;

import java.util.ArrayList;

import mpicbg.imglib.algorithm.fft.FourierConvolution;
import mpicbg.imglib.algorithm.mirror.MirrorImage;
import mpicbg.imglib.container.constant.ConstantContainer;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LRFFT_Test 
{
	public static enum PSFTYPE { SIMPLE, EXPONENT, CONDITIONAL, CONDITIONAL_NEW, MAPG };
	
	private Image<FloatType> image, weight, kernel1, kernel2;
	Image<FloatType> viewContribution = null;
	FourierConvolution<FloatType, FloatType> fftConvolution1, fftConvolution2;
	protected int numViews = 0;
	
	/**
	 * The rotation angle of the PSF, just for internal use
	 */
	public int angle;
	
	PSFTYPE type = PSFTYPE.SIMPLE;
	ArrayList< LRFFT_Test > views;
	
	/**
	 * Used to determine if the Convolutions already have been computed for the current iteration
	 */
	int i = -1;

	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> weight, final Image<FloatType> kernel )
	{
		this( image, weight, kernel, 0 );
	}
	
	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> weight, final Image<FloatType> kernel, final int angle )
	{
		this.image = image;
		this.kernel1 = kernel;
		this.weight = weight;
		this.angle = angle;
	}

	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> kernel )
	{
		this( image, kernel, 0 );
	}
	
	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> kernel, final int angle )
	{
		this( image, new Image< FloatType > ( new ConstantContainer< FloatType >( image.getDimensions(), new FloatType( 1 ) ), new FloatType() ), kernel, angle );
	}

	/**
	 * @param numViews - the number of views in the acquisition, determines the exponential of the kernel
	 */
	protected void setNumViews( final int numViews ) { this.numViews = numViews; }
	
	/**
	 * This method is called once all views are added to the {@link LRInput}
	 */
	protected void init( final PSFTYPE type, final ArrayList< LRFFT_Test > views ) 
	{
		this.type = type;
		this.views = views;
		
		// normalize kernel so that sum of all pixels == 1
		AdjustInput.normImage( kernel1 );

		if ( numViews == 0 )
		{
			System.out.println( "Warning, numViews was not set." );
			numViews = 1;
		}
		
		if ( type == PSFTYPE.MAPG )
		{	
			this.kernel2 = null;
		}
		else if ( numViews == 1 || type == PSFTYPE.SIMPLE )
		{
			// compute the inverted kernel (switch dimensions)
			this.kernel2 = computeInvertedKernel( this.kernel1 );			
		}
		else if ( type == PSFTYPE.CONDITIONAL )
		{
			// compute the kernel using conditional probabilities
			// if this is x1, then compute
			// P(x1|psi) * P(x2|x1) * P(x3|x1) * ... * P(xn|x1)
			// where
			// P(xi|x1) = P(x1|psi) convolved with P(xi|psi)
			
			// we first get P(x1|psi)
			final Image< FloatType > tmp = computeInvertedKernel( this.kernel1.clone() );

			// now convolve P(x1|psi) with all other kernels 
			for ( final LRFFT_Test view : views )
			{
				if ( view != this )
				{
					final FourierConvolution<FloatType, FloatType> conv = new FourierConvolution<FloatType, FloatType>( computeInvertedKernel( this.kernel1.clone() ), ( view.kernel1 ) );
					conv.setNumThreads();
					conv.setKeepImgFFT( false );
					conv.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
					conv.process();
					
					// multiply with the kernel
					final Cursor<FloatType> cursor = tmp.createCursor();
					for ( final FloatType t : ( conv.getResult() ) )
					{
						cursor.fwd();
						cursor.getType().set( t.get() * cursor.getType().get() );
					}
				}
			}
			
			// norm the compound kernel
			AdjustInput.normImage( tmp );
						
			// compute the inverted kernel
			this.kernel2 = ( tmp );
			
			// close the temp image
			//tmp.close();
		}
		else if ( type == PSFTYPE.CONDITIONAL_NEW )
		{
			// compute the compound kernel P_v^compound of the efficient bayesian multi-view deconvolution
			// for the current view \phi_v(x_v)
			//
			// P_v^compound = P_v^{*} prod{w \in W_v} P_v^{*} \ast P_w \ast P_w^{*}
			//              = P_v^{*} prod{w \in W_v} P_{W_v}^compound
			
			// we first get P_v^{*} -> * refers to the inverted coordinates
			final Image< FloatType > tmp = computeInvertedKernel( this.kernel1.clone() );

			// now for each view: w \in W_v
			for ( final LRFFT_Test view : views )
			{
				if ( view != this )
				{
					// convolve first P_v^{*} with P_w
					final FourierConvolution<FloatType, FloatType> conv1 = new FourierConvolution<FloatType, FloatType>( computeInvertedKernel( this.kernel1 ), view.kernel1 );
					conv1.setNumThreads();
					conv1.setKeepImgFFT( false );
					conv1.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
					conv1.process();
					
					//ImageJFunctions.show( conv1.getResult() );
		
					// and now convolve the result with P_w^{*}
					final FourierConvolution<FloatType, FloatType> conv2 = new FourierConvolution<FloatType, FloatType>( conv1.getResult(), computeInvertedKernel( view.kernel1 ) );
					conv2.setNumThreads();
					conv2.setKeepImgFFT( false );
					conv2.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
					conv2.process();
					
					// multiply the result with P_v^{*} yielding the compound kernel
					final Cursor<FloatType> cursor = tmp.createCursor();
					for ( final FloatType t : ( conv2.getResult() ) )
					{
						cursor.fwd();
						cursor.getType().set( t.get() * cursor.getType().get() );
					}					
				}
			}
			
			// norm the compound kernel
			AdjustInput.normImage( tmp );
						
			// set it as kernel2 of the deconvolution
			this.kernel2 = ( tmp );									
		}
		else if ( type == PSFTYPE.EXPONENT )
		{
			// compute the squared kernel and its inverse
			final Image< FloatType > exponentialKernel = computeExponentialKernel( this.kernel1, numViews );
			
			// norm the squared kernel
			AdjustInput.normImage( exponentialKernel );
			
			// compute the inverted squared kernel
			this.kernel2 = computeInvertedKernel( exponentialKernel );
		}
		else
		{
			throw new RuntimeException( "Kernel type not supported: " + type );
		}
		
		//ImageJFunctions.show( this.kernel2 ).setTitle( "kernel2 of " + image.getName() );
		
		// the first kernel is to compute the ratio between original image and current guess,
		// so we need no exponent here
		this.fftConvolution1 = new FourierConvolution<FloatType, FloatType>( this.image, this.kernel1 );	
		this.fftConvolution1.setNumThreads();
		this.fftConvolution1.setKeepImgFFT( false );
		//this.fftConvolution1.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
		
		this.fftConvolution2 = new FourierConvolution<FloatType, FloatType>( this.image, this.kernel2 );	
		this.fftConvolution2.setNumThreads();
		this.fftConvolution2.setKeepImgFFT( false );
		//this.fftConvolution2.setImageOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
	}
	
	public static Image<FloatType> computeExponentialKernel( final Image<FloatType> kernel, final int numViews )
	{
		final Image<FloatType> exponentialKernel = kernel.clone();
		
		for ( final FloatType f : exponentialKernel )
			f.set( pow( f.get(), numViews ) );
		
		return exponentialKernel;
	}

	public static Image< FloatType > computeInvertedKernel( final Image< FloatType > kernel )
	{
		final Image< FloatType > invKernel = kernel.clone();
		
		for ( int d = 0; d < invKernel.getNumDimensions(); ++d )
			new MirrorImage< FloatType >( invKernel, d ).process();
		
		return invKernel;
	}

	final private static float pow( final float value, final int power )
	{
		float result = value;
		
		for ( int i = 1; i < power; ++i )
			result *= value;
		
		return result;
	}


	public void setImage( final Image<FloatType> image ) 
	{
		this.image = image;
		setCurrentIteration( -1 );
	}
	public void setWeight( final Image<FloatType> weight ) { this.weight = weight; }
	public void setKernel( final Image<FloatType> kernel ) 
	{
		this.kernel1 = kernel;
		
		init( type, views );

		setCurrentIteration( -1 );
	}

	public Image<FloatType> getImage() { return image; }
	public Image<FloatType> getWeight() { return weight; }
	public Image<FloatType> getKernel1() { return kernel1; }
	public Image<FloatType> getKernel2() { return kernel2; }
	public Image<FloatType> getViewContribution() { return viewContribution; }
	
	public void setCurrentIteration( final int i ) { this.i = i; }
	public int getCurrentIteration() { return i; }
	
	public FourierConvolution<FloatType, FloatType> getFFTConvolution1() { return fftConvolution1; }
	public FourierConvolution<FloatType, FloatType> getFFTConvolution2() { return fftConvolution2; }
	
	public void setViewContribution( final Image<FloatType> viewContribution )
	{
		if ( this.viewContribution != null )
			this.viewContribution.close();
		
		this.viewContribution = viewContribution;
	}
	
	@Override
	public LRFFT_Test clone()
	{
		final LRFFT_Test viewClone = new LRFFT_Test( this.image.clone(), this.weight.clone(), this.kernel1.clone(), angle );
	
		viewClone.numViews = numViews;
		viewClone.type = type;
		viewClone.views = views;
		viewClone.i = i;
		
		viewClone.image.setName( this.image.getName() + " cloned" );
		viewClone.weight.setName( this.weight.getName() + " cloned" );
		viewClone.kernel1.setName( this.kernel1.getName() + " cloned" );
		
		if ( this.kernel2 != null )
			viewClone.kernel2 = kernel2.clone();
		
		if ( this.viewContribution != null )
			viewClone.viewContribution = this.viewContribution.clone();
		
		if ( this.fftConvolution1 != null )
		{
			viewClone.fftConvolution1 = new FourierConvolution<FloatType, FloatType>( fftConvolution1.getImage(), fftConvolution1.getKernel() );
			viewClone.fftConvolution1.process();
		}

		if ( this.fftConvolution2 != null )
		{
			viewClone.fftConvolution2 = new FourierConvolution<FloatType, FloatType>( fftConvolution2.getImage(), fftConvolution2.getKernel() );
			viewClone.fftConvolution2.process();
		}

		return viewClone;
	}	
}
