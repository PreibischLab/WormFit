package varun;


import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.Localizable;
import net.imglib2.Point;
import net.imglib2.PointSampleList;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessible;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.RealCursor;
import net.imglib2.img.Img;
import net.imglib2.img.cell.CellImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import util.ImgLib2Util;


public class TwoDtree {
	
		
		
		public static <T extends RealType<T>> PointSampleList<T> getList(RandomAccessibleInterval<T> img) {

			final RandomAccessible<T> infinite = Views.extendZero(img);

			final int n = img.numDimensions();
			long min[] = new long[n];
			long max[] = new long[n];

			for (int d = 0; d < n; ++d) {

				min[d] = img.min(d);
				max[d] = img.max(d);

			}

			FinalInterval interval = new FinalInterval(min, max);

			final IterableInterval<T> imgav = Views.interval(infinite, interval);

			final Cursor<T> first = imgav.cursor();

			// A point sample list with coordinates declared and initialized.
			PointSampleList<T> parent = new PointSampleList<T>(n);

			while (first.hasNext()) {
				first.fwd();
				Point cord = new Point(n);

				cord.setPosition(first);

				parent.add(cord, first.get().copy());

			}

			return parent;

		}
		
		
		
		public static <T extends RealType<T>> void splitbyCoordinate(PointSampleList<T> list, int direction ) {

			int n = list.numDimensions();
			/****
			 * To ward against running over the dimensionality, creating some
			 * local restrictions on the global variable direction
			 ****/
			if (direction == list.numDimensions())
				direction = 0;
			if (list.dimension(direction) <= 1)
				return;

			else {

				/****
				 * To ward against running over the dimensionality, creating some
				 * local restrictions on the global variable direction
				 ****/
				if (direction == list.numDimensions())
					direction = 0;

				
				
				
				
				// the first element belonging to the right list childB
				final int splitIndex = (int)list.dimension(direction) / 2;

				System.out.println( splitIndex );
				
				final PointSampleList<T> childA = new PointSampleList<T>(n);
				final PointSampleList<T> childB = new PointSampleList<T>(n);

				final Cursor<T> listCursor = list.localizingCursor();
				
				
				final ArrayList<  Long  > values = new ArrayList< Long  >( ( int ) list.dimension(direction) );

				while ( listCursor.hasNext() )
				{
					listCursor.next();
					values.add( listCursor.getLongPosition(direction) );
				}
				
				
		//	sort(values, direction);
				
				int i = 0;
				while (listCursor.hasNext()) {

					listCursor.fwd();

					Point cord = new Point(listCursor);


					if ( i < splitIndex )
					{

						childA.add(cord, listCursor.get().copy());

					//	System.out.println("childA: "+listCursor.get());

					} else

					{

						childB.add(cord, listCursor.get().copy());
					//	System.out.println("childB: "+listCursor.get());
					}
					i++;
				}

				
				splitbyCoordinate(childA, direction+1);

				splitbyCoordinate(childB, direction+1);

				mergeList(list, childA, childB, direction);
			}

		}
		
		
		
		
		public static <T extends RealType<T>> Localizable firstLocation(final IterableInterval<T> interval) {
			Cursor<T> c = interval.localizingCursor();
			c.fwd();
			return c;
		}
		
		public static <T extends RealType<T>> void split(ArrayList< Long > coordinateList, int direction ) {

			

			if (coordinateList.size() <= 1)
				return;

			else {

				
				// the first element belonging to the right list childB
				final int splitIndex = (int)coordinateList.size()/ 2;

				
				
				final ArrayList<Long > childA = new ArrayList<Long >(( int ) coordinateList.size()/2 );
			
				final ArrayList<Long> childB = new ArrayList<Long >(( int ) coordinateList.size()/2+coordinateList.size()%2 );
				
				
				
				

				
for(int xindex=0;xindex<splitIndex; ++xindex){
						childA.add( coordinateList.get(xindex));
}
for(int xindex=0;xindex<splitIndex+coordinateList.size()%2; ++xindex){

						childB.add(coordinateList.get(xindex));
					//	System.out.println("childB: "+listCursor.get());
					}
					
				

				
				split(childA, direction);

				split(childB, direction);

				mergeListValue(coordinateList, childA, childB);
			}

		}
		
	
		
		///*****       Returns a sorted list *********////
		public static <T extends RealType<T>> void mergeListValue(ArrayList<Long> sortedlist, ArrayList<Long> listA,
				ArrayList<Long> listB) {
/*
			final Cursor<T> cursorA = listA.localizingCursor();
			final Cursor<T> cursorB = listB.localizingCursor();
			final Cursor<T> cursor = list.localizingCursor();

			
			
			
			cursorA.fwd();
			cursorB.fwd();

		//	System.out.println("listA : " + cursorA.get());
		//	System.out.println("listB : " + cursorB.get());

			boolean cannotMoveOn = false;
			
			do
			{
				// here is where you decide what you sort after
				if (cursorA.getLongPosition(direction)<cursorB.getLongPosition(direction)) {

					cursor.fwd();
					cursor.get().set( cursorA.get() );
					if ( cursorA.hasNext() )
						cursorA.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorB until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorB.get() );
							if ( cursorB.hasNext() )
								cursorB.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
					
			//		System.out.println("In here");
				}

				else

				{

					cursor.fwd();
					cursor.get().set( cursorB.get() );
					if ( cursorB.hasNext() )
						cursorB.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorA until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorA.get() );
							if ( cursorA.hasNext() )
								cursorA.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
			//		System.out.println("Out here");
				}

			}
			while ( cannotMoveOn == false );
			
			*/
		}
		
		
		
		
		
		
		
		
		
		
		///*****       Returns a sorted list *********////
		public static <T extends RealType<T>> void mergeList(PointSampleList<T> list, PointSampleList<T> listA,
				PointSampleList<T> listB, int direction) {

			final Cursor<T> cursorA = listA.localizingCursor();
			final Cursor<T> cursorB = listB.localizingCursor();
			final Cursor<T> cursor = list.localizingCursor();

			
			
			
			cursorA.fwd();
			cursorB.fwd();

		//	System.out.println("listA : " + cursorA.get());
		//	System.out.println("listB : " + cursorB.get());

			boolean cannotMoveOn = false;
			
			do
			{
				// here is where you decide what you sort after
				if (cursorA.getLongPosition(direction)<cursorB.getLongPosition(direction)) {

					cursor.fwd();
					cursor.get().set( cursorA.get() );
					if ( cursorA.hasNext() )
						cursorA.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorB until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorB.get() );
							if ( cursorB.hasNext() )
								cursorB.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
					
			//		System.out.println("In here");
				}

				else

				{

					cursor.fwd();
					cursor.get().set( cursorB.get() );
					if ( cursorB.hasNext() )
						cursorB.fwd();
					else
					{
						cannotMoveOn = true;
						
						// move cursorA until the end
						boolean stopped = false;
						do
						{
							cursor.fwd();
							cursor.get().set( cursorA.get() );
							if ( cursorA.hasNext() )
								cursorA.fwd();
							else
								stopped = true;					
						}
						while ( stopped == false );
					}
			//		System.out.println("Out here");
				}

			}
			while ( cannotMoveOn == false );
		}
		
		
		
		public static void main(String[] args) {

			final Img<FloatType> img = ImgLib2Util.openAs32Bit(new File("src/main/resources/bridge.png"));

			

			PointSampleList<FloatType> list = new PointSampleList<FloatType>(img.numDimensions());
			
			
			//Make a 1D list along the X direction by setting an appropriate interval on the image. 

			IterableInterval<FloatType> view = Views.interval(img, new long[] { 0, 0 }, new long[] { 2, 2 });

			final Cursor<FloatType> first = view.cursor();

			while (first.hasNext()) {
				first.fwd();
				Point cord = new Point(img.numDimensions());

				cord.setPosition(first);

				list.add(cord, first.get().copy());
				 System.out.println("Set of x co-ordinates Initial List : " +
				 cord.getDoublePosition(0));
				 System.out.println("Set of y co-ordinates Initial List : " +
				 cord.getDoublePosition(1));
				System.out.println("Values Initial list : " + first.get());

			}

			

		//	split(list, 0); // Split list along X direction
			

			Cursor<FloatType> testtwo = list.cursor();

			while (testtwo.hasNext()) {
				testtwo.fwd();
				Point newpoint = new Point(img.numDimensions());

				newpoint.setPosition(testtwo);

				 System.out.println("Set of x co-ordinates sorted List : " +
				 newpoint.getDoublePosition(0));
				 System.out.println("Set of y co-ordinates sorted List : " +
				 newpoint.getDoublePosition(1));
				System.out.println("Values sorted list : " + testtwo.get());

			}

		}

	}

	
	


