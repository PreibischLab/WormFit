package klim.matching;

import java.util.ArrayList;
import java.util.List;

import mpicbg.models.CoordinateTransform;
import mpicbg.models.CoordinateTransformMap2D;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.models.RigidModel2D;
import mpicbg.models.TranslationModel2D;
import net.imglib2.util.Util; 

public class RansacRigidExample {

	public static long count = 0;

	public static void setPoints(List< Point > iPoints, List< Point > fPoints){
		// l[] - local coordinates
		// set initial points
		iPoints.add(new Point( new double[]{ 0.0, 0.0 } ));
		iPoints.add(new Point( new double[]{ 2.0, 0.0 } ));
		iPoints.add(new Point( new double[]{ 4.0, 0.0 } ));
		iPoints.add(new Point( new double[]{ 6.0, 0.0 } ));
		iPoints.add(new Point( new double[]{ 0.0, 2.0 } ));
		iPoints.add(new Point( new double[]{ 0.0, 4.0 } ));
		
		// TODO: add some noise 

		// the transformation that is _not_ known in general 
		// used here as an example

		// rotate by 90 degrees and translate by (1, 1)
		for(Point p : iPoints){
			double x = p.getL()[0]* Math.cos(Math.toRadians(90)) - p.getL()[1]* Math.sin(Math.toRadians(90));
			double y = p.getL()[0]* Math.sin(Math.toRadians(90)) + p.getL()[1]* Math.cos(Math.toRadians(90));
			x += 1.0;
			y += 1.0;
			fPoints.add(new Point(new double[]{x, y}));
		}
	}


	public static void getPermutation(){

	}

	public static void getPermutation(int in[], final List< Point > iPoints, final List< Point > fPoints) {
		int n = in.length;
		int[] a = new int[n];
		for (int i = 0; i < n; i++)
			a[i] = in[i];
		getPermutation(a, n, iPoints, fPoints);
	}
	
	public static void getPermutation(int [] a, int n, final List< Point > iPoints, final List< Point > fPoints){
		if (n == 1) {
			
			for (int i : a)
				System.out.print(i);
			System.out.println();
			//  System.out.println(a);

			// here is the place where you set matching of the elements
			// and run the algorithm

			runAlgo(a, iPoints, fPoints);

			return;
		}
		for (int i = 0; i < n; i++) {
			swap(a, i, n-1);
			getPermutation(a, n - 1, iPoints, fPoints);
			swap(a, i, n-1);
		}
	}
	
	
	public static void swap(int [] a, int i, int j){
		int c = a[i];
		a[i] = a[j];
		a[j] = c;
	}

	public static void runAlgo(int[] a, final List< Point > iPoints, final List< Point > fPoints){

		RigidModel2D m = new RigidModel2D();
		// System.out.println( m.getClass().getSimpleName() + " needs " + m.getMinNumMatches() + " matches.");

		final ArrayList< PointMatch > candidates = new ArrayList<>();

		for(int i = 0; i < a.length; ++i){
			// candidates.add( new PointMatch(iPoints.get(i), fPoints.get(Character.getNumericValue(a[i]))));
			candidates.add( new PointMatch(iPoints.get(i), fPoints.get(a[i])));
			
		}

		final ArrayList< PointMatch > inliers = new ArrayList<>();

		try{
			m.ransac(candidates, inliers, 10000, 0.01, 0.4 );
			if ( inliers.size() == 0 )
			{
				// System.out.println( "No model found." );
				return;
			}
			// from the inliers compute the actual model
			m.fit( inliers );

			for ( PointMatch pm : inliers )
				pm.apply( m );

			m.setCost( PointMatch.meanDistance( inliers ) );
		}
		catch(NotEnoughDataPointsException e){
			e.printStackTrace();
			return;
		}

		// System.out.println( "Remaining inliers: " + inliers.size() + " (" + (100*(double)inliers.size()/candidates.size()) + "%)" );
		// System.out.println( m );
		// System.out.println();

		for ( PointMatch pm : inliers )
		{			
			// alternatively (identical)
			pm.getP1().apply( m );	
			int idx1 = candidates.indexOf(pm);
			int idx2 = candidates.indexOf(pm);
			// System.out.println( "p1(w[" + idx1 + "]):"  + Util.printCoordinates( pm.getP1().getW() ) + " >> " + "p2(l[" + idx2 + "]):" + Util.printCoordinates( pm.getP2().getW() ) + " error: " + pm.getDistance() );
		}
	}


	public static void fitPointCloudsRigidOutlierRemoval()
	{

		final List< Point > iPoints = new ArrayList<>();
		final List< Point > fPoints = new ArrayList<>();

		setPoints(iPoints, fPoints);
		
		int n = 6; // number of points 
		int[] initial = new int[n];
		for(int i = 0; i < n; ++i)
			initial[i] = n - i - 1;
		getPermutation(initial, iPoints, fPoints);
	}

	public static void main( String[] args )
	{
		fitPointCloudsRigidOutlierRemoval();
	}

}
