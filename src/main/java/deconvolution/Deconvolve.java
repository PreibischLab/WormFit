package deconvolution;

import ij.ImageJ;
import ij.ImagePlus;

import java.util.ArrayList;
import java.util.Date;

import deconvolution.LRFFT_Test.PSFTYPE;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.type.numeric.real.FloatType;

public class Deconvolve
{
	public final static float avgIntensity = 1;
	public final float minValue = LucyRichardsonSeq.minValue;

	public Deconvolve( final String[] images, final String[] psfs )
	{
		final ArrayList< Deconvolver > decons = new ArrayList<Deconvolver>() ;
		
		final LRInput deconvolutionData = new LRInput();

		for ( int i = 0; i < images.length; ++i )
		{
			final Image< FloatType > im = LRInput.open( images[ i ], new ArrayContainerFactory() );
			final Image< FloatType > psf = LRInput.open( psfs[ i ], new ArrayContainerFactory() );
		
			AdjustInput.normImage( psf );
			AdjustInput.adjustImage( im, LucyRichardsonSeq.minValue, avgIntensity );
			
			ImageJFunctions.show( im );
			ImageJFunctions.show( psf );
			
			deconvolutionData.add( new LRFFT_Test( im, psf, i ) );
		}

		final Deconvolver d = new LucyRichardsonSeq( deconvolutionData.clone().init( PSFTYPE.SIMPLE ), false, 6.0E-4, "independent (sequential)" );
		//new MAPG( deconvolutionData.clone().init( PSFTYPE.SIMPLE ) );
		/*
		final Image< FloatType > it = TestLR.open( "/Volumes/ExternalBay/SPIM/SPIMfixed/patrizia/2p_decon/it=100.tif", new ArrayContainerFactory() );
		final Cursor< FloatType > c1 = d.getPsi().createCursor();
		final Cursor< FloatType > c2 = it.createCursor();
		
		while ( c1.hasNext() )
		{
			c1.next().set( c2.next().get() );
		}
		*/
		for ( int i = 0; i < 110; ++ i )
		{
			System.out.println( new Date( System.currentTimeMillis() ) + " " +  i );
			d.runIteration();
			
			if ( i > 0 && ( i == 10 ||  i == 50 || i % 100 == 0 ) )
			{
				ImagePlus imp = ImageJFunctions.copyToImagePlus( d.getPsi() );
				imp.setTitle( "it=" + i );
				imp.resetDisplayRange();
				imp.show();
			}
		}
		
		ImageJFunctions.show( d.getPsi() );
	}
	
	public static void main( String[] args )
	{
		final String dir = "/Users/spreibi/Downloads/";///Volumes/ExternalBay/Multiview Imaging/l1/";
		//"/Volumes/ExternalBay/SPIM/SPIMfixed/patrizia/";//"
		final String[] images = new String[]{ dir + "Channel 2.tif" };//"spim_TL1_Angle999-crop.tif" };
		final String[] psfs = new String[]{ dir + "claire_psf.tif" }; //"psf_2p.tif" };

		new ImageJ();
		new Deconvolve( images, psfs );
	}
}
