package deconvolution;

import ij.ImagePlus;
import ij.ImageStack;
import ij.io.Opener;
import ij.process.ImageProcessor;

import java.io.IOException;
import java.util.ArrayList;

import deconvolution.LRFFT_Test.PSFTYPE;
import mpicbg.imglib.container.ContainerFactory;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.cursor.LocalizablePlaneCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.ImageFactory;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.type.numeric.real.FloatType;

public class LRInput 
{
	final ArrayList< LRFFT_Test > views = new ArrayList<LRFFT_Test>();
	
	public void add( final LRFFT_Test view )
	{
		views.add( view );
		
		for ( final LRFFT_Test v : views )
			v.setNumViews( getNumViews() );
	}
		
	/**
	 * init all views
	 * 
	 * @param exponentialKernel - use exponential kernel?
	 * 
	 * @return the same instance again for convinience
	 */
	public LRInput init( final PSFTYPE type )
	{
		for ( final LRFFT_Test view : views )
			view.init( type, views );

		return this;
	}
	
	/**
	 * @return - the image data
	 */
	public ArrayList< LRFFT_Test > getViews() { return views; }
	
	/**
	 * The number of views for this deconvolution
	 * @return
	 */
	public int getNumViews() { return views.size(); }
	
	@Override
	public LRInput clone()
	{
		LRInput clone = new LRInput();
		
		for ( final LRFFT_Test view : views )
			clone.add( view.clone() );
		
		return clone;
	}

	public static Image<FloatType> open( final String name, final ContainerFactory containerFactory )
	{
		final Opener io = new Opener();
		ImagePlus imp = io.openImage( name );
		
		if ( imp.getStack().getSize() > 1 )
		{
			final int depth = imp.getStack().getSize();
			
			ImageFactory<FloatType> factory = new ImageFactory<FloatType>( new FloatType(), containerFactory );
			Image<FloatType> img = factory.createImage( new int[]{ imp.getWidth(), imp.getHeight(), depth } );
			
			for ( int i = 0; i < depth; ++i )
			{
				final ImageProcessor ip = imp.getStack().getProcessor( i + 1 );
				
				final LocalizablePlaneCursor<FloatType> cursorOut = img.createLocalizablePlaneCursor();
				cursorOut.reset( 0, 1, new int[]{ 0, 0, i } );
				
				while ( cursorOut.hasNext() )
				{
					cursorOut.fwd();
					cursorOut.getType().set( ip.getPixelValue( cursorOut.getPosition( 0 ), cursorOut.getPosition( 1 ) ) );
				}
			}			

			return img;
		}
		else
		{
			ImageFactory<FloatType> factory = new ImageFactory<FloatType>( new FloatType(), containerFactory );
			Image<FloatType> img = factory.createImage( new int[]{ imp.getWidth(), imp.getHeight() } );
			
			final ImageProcessor ip = imp.getProcessor();
			
			final LocalizableCursor<FloatType> cursorOut = img.createLocalizableCursor();
			
			while ( cursorOut.hasNext() )
			{
				cursorOut.fwd();
				cursorOut.getType().set( ip.getPixelValue( cursorOut.getPosition( 0 ), cursorOut.getPosition( 1 ) ) );
			}
			
			return img;
		}
	}

}
