package deconvolution;

import mpicbg.imglib.algorithm.mirror.MirrorImage;
import mpicbg.imglib.container.constant.ConstantContainer;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LRFFT_Test 
{
	final private Image< FloatType > image, weight, kernel1, kernel2;

	/**
	 * Used to determine if the Convolutions already have been computed for the current iteration
	 */
	int i = -1;
	
	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> weight, final Image<FloatType> kernel )
	{
		this.image = image;
		this.kernel1 = kernel;
		this.weight = weight;

		// normalize kernel so that sum of all pixels == 1
		AdjustInput.normImage( kernel1 );

		// compute the inverted kernel (switch dimensions)
		this.kernel2 = computeInvertedKernel( this.kernel1 );
	}

	public LRFFT_Test( final Image<FloatType> image, final Image<FloatType> kernel )
	{
		this( image, new Image< FloatType > ( new ConstantContainer< FloatType >( image.getDimensions(), new FloatType( 1 ) ), new FloatType() ), kernel );
	}

	
	public static Image< FloatType > computeInvertedKernel( final Image< FloatType > kernel )
	{
		final Image< FloatType > invKernel = kernel.clone();
		
		for ( int d = 0; d < invKernel.getNumDimensions(); ++d )
			new MirrorImage< FloatType >( invKernel, d ).process();
		
		return invKernel;
	}

	public Image<FloatType> getImage() { return image; }
	public Image<FloatType> getWeight() { return weight; }
	public Image<FloatType> getKernel1() { return kernel1; }
	public Image<FloatType> getKernel2() { return kernel2; }
	
	public void setCurrentIteration( final int i ) { this.i = i; }
	public int getCurrentIteration() { return i; }
}
