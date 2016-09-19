package stephan;

import java.io.File;

import ij.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;

public class AvgPSFs
{
	public static void main( String[] args )
	{
		Img< FloatType > img = null;

		for ( int i = 0; i < 10; ++i )
		{
			File f = new File( "/Users/spreibi/Downloads/PSFs", "Image " + i + "-1.tif" );
			if ( f.exists() )
			{
				System.out.println( "Opening: " + f );
				Img< FloatType > psf = ImgLib2Util.openAs32Bit( f );

				if ( img == null )
				{
					// first one is just copied
					img = psf.copy();
				}
				else
				{
					// otherwise averaged
					Cursor<FloatType> psfCursor = Views.iterable(psf).cursor();
					RandomAccess<FloatType> ra = img.randomAccess();

					while(psfCursor.hasNext())
					{
						psfCursor.fwd();
						ra.setPosition(psfCursor);
						psfCursor.get().add(ra.get());
					}
				}
			}
		}
		
		new ImageJ();
		ImageJFunctions.show( img );
	}
}
