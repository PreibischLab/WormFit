package klim.matching;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.models.RigidModel2D;
import mpicbg.models.TranslationModel2D;
import net.imglib2.util.Util;

public class RANSAC
{

	public static void fitPointCloudsTranslationOutlierRemoval()
	{
		RigidModel2D m2 = new RigidModel2D();
		m2.set( Math.toRadians( 45), 400, 120 );
		
		// l[] - local coordinates
		Point p1a = new Point( new double[]{ 0.0, 0.0 } );
		Point p2a = new Point( new double[]{ 100.0, 0.0 } );
		Point p3a = new Point( new double[]{ 0.0, 200.0 } );
		Point p4a = new Point( new double[]{ 100.0, 200.0 } );
		Point p5a = new Point( new double[]{ 0.0, 0.0 } );

		// l[] - local coordinates
		Point p1b = new Point( new double[]{ 10.0, 20.0 } );
		Point p2b = new Point( new double[]{ 110.0, 20.0 } );
		Point p3b = new Point( new double[]{ 10.0, 220.0 } );
		Point p4b = new Point( new double[]{ 110.0, 220.0 } );
		Point p5b = new Point( new double[]{ 1000.0, 2000.0 } );

		// potentially corresponding
		final List< PointMatch > candidates = new ArrayList<>();
		
		candidates.add( new PointMatch( p1a, p1b ) );
		candidates.add( new PointMatch( p2a, p2b ) );
		candidates.add( new PointMatch( p3a, p3b ) );
		candidates.add( new PointMatch( p4a, p4b ) );
		candidates.add( new PointMatch( p5a, p5b ) );

		for ( PointMatch pm : candidates )
			System.out.println( "p1(l[]): " + Util.printCoordinates( pm.getP1().getL() ) + " >> p2(l[]): " + Util.printCoordinates( pm.getP2().getL() ) );
		System.out.println();
		
		// maps the w[] - world coordinates of p1 to l[] of p2
		TranslationModel2D m = new TranslationModel2D();
		//RigidModel2D m;
		
		System.out.println( m.getClass().getSimpleName() + " needs " + m.getMinNumMatches() + " matches, we have: " + candidates.size() );

		// the true ones go in here
		final ArrayList< PointMatch > inliers = new ArrayList<>();
		
		try
		{
			// now it finds the transformation mapping p1(w[]) >> p2(l[])
			m.ransac(candidates, inliers, 10000, 0.01, 0.7 );

			if ( inliers.size() == 0 )
			{
				System.out.println( "No model found." );
				return;
			}
			// from the inlies compute the actual model
			m.fit( inliers );

			for ( PointMatch pm : inliers )
				pm.apply( m );

			m.setCost( PointMatch.meanDistance( inliers ) );
		}
		catch (NotEnoughDataPointsException e)
		{
			e.printStackTrace();
			return;
		}

		System.out.println( "Remaining inliers: " + inliers.size() + " (" + (100*(double)inliers.size()/candidates.size()) + "%)" );
		System.out.println( m );
		System.out.println();
				
		for ( PointMatch pm : inliers )
		{
			// apply the transformation to p1 of the pointmatch
			//pm.apply( m );
			
			// alternatively (identical)
			pm.getP1().apply( m );
			
			// alternatively (identical)
			//double[] w = m.apply( pm.getP1().getL() );
			//for ( int d = 0; d < pm.getP1().getW().length; ++ d )
			//	pm.getP1().getW()[ d ] = w[ d ];
			
			System.out.println( "p1(w[]): " + Util.printCoordinates( pm.getP1().getW() ) + " >> p2(l[]): " + Util.printCoordinates( pm.getP2().getW() ) + " error: " + pm.getDistance() );
		}
	}

	public static void main( String[] args )
	{
		fitPointCloudsTranslationOutlierRemoval();
	}

}
