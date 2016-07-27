package deconvolution;

import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.RealType;

import org.uncommons.maths.number.NumberGenerator;

public class NumberGeneratorImage< T extends RealType< T > > implements NumberGenerator< Double >
{
	final Cursor< T > cursor;
	final double multiplicativeFactor;
	
	public NumberGeneratorImage( final Image<T> image, final double multiplicativeFactor )
	{
		this.cursor = image.createCursor();
		this.multiplicativeFactor = multiplicativeFactor;
	}
	
	/**
	 * Otherwise it gets out of sync for some reason
	 */
	public void fwd()
	{
		cursor.fwd();
	}
	
	public void reset()
	{
		cursor.reset();
	}
	
	@Override
	public Double nextValue()
	{
		return cursor.getType().getRealDouble() * multiplicativeFactor;
	}
}
