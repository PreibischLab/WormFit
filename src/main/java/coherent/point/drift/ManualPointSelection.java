package coherent.point.drift;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.Roi;
import mpicbg.imglib.util.Util;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;
import util.opencsv.CSVWriter;

/**
 * ManualPointSelection -- fast way to get the points manually 
 * and write them to cvs after that
 */

public class ManualPointSelection {

	// selected points
	final ArrayList<double []> points =  new ArrayList<double[]>(1);
	final Img <FloatType> img;
	final ImagePlus imp; // wrapper
	final int numDimensions;
	
	// used for detection
	final MyMouseListener mouseListener;
	final MyKeyListener keyListener;
	
	// file names
	String name =  "worm-folded";
	String dir = "/Users/kkolyva/Desktop/cpd_ex/";
	File output = new File(dir + name + ".csv");
	File input = new File(dir + name + ".tif");

	public ManualPointSelection(){
		img = ImgLib2Util.openAs32Bit(input);
		numDimensions = img.numDimensions();		
		imp = ImageJFunctions.wrapFloat(img, "");
		
		imp.show();
		mouseListener = new MyMouseListener();
		keyListener = new MyKeyListener();
		imp.getCanvas().addMouseListener(mouseListener);
		imp.getCanvas().addKeyListener(keyListener);	
	}

	public void addOverlay(ImagePlus imp){
		Overlay overlay = imp.getOverlay();
		if (overlay == null){
			overlay = new Overlay();
			imp.setOverlay( overlay );
		}		
		overlay.clear();
		
		// here should be the T coordinate of the fitted points
		double [] location = new double[numDimensions];
		for (double[] point : points){
			location = point.clone();
			// special 
			final PointRoi or = new PointRoi(location[0], location[1]);
			or.setStrokeColor(Color.RED);
			overlay.add(or);
		}
	}
	

	protected class MyMouseListener implements MouseListener
	{
		@Override
		public void mouseClicked(MouseEvent e) {}

		@Override
		public void mouseEntered(MouseEvent e) {}

		@Override
		public void mouseExited(MouseEvent e) {}

		@Override
		public void mousePressed(MouseEvent e) {}

		@Override
		public void mouseReleased( final MouseEvent e )
		{
			// convert the coordinates to "proper" coordinates
			final int xm = imp.getCanvas().offScreenX( e.getX() );
			final int ym = imp.getCanvas().offScreenY( e.getY() );
		
			points.add(new double[]{xm, ym});
			System.out.println("Added:\n" + "x: " + xm + " y: " + ym);
			addOverlay(imp);
		}

	}

	protected class MyKeyListener implements KeyListener
	{
		public void keyTyped(KeyEvent e) {}

		public void keyPressed(KeyEvent e) {}

		public void keyReleased(KeyEvent e) {
			int id = e.getID();	
			if (id == KeyEvent.KEY_RELEASED) {
				char c = e.getKeyChar();
				switch (c){
				case 's':
				case 'f':
				case 'q':
					savePointsToFile(output);
					System.out.println("Saved:\n" + points.size() + " Points");
					break;
				case 'r':
				case 'd':	
					if (points.size() - 1 > 0){
						System.out.println("Removed:\n" + "x: " + points.get(points.size() - 1)[0] + " y: " + points.get(points.size() - 1)[1]);
						points.remove(points.size() - 1);
						addOverlay(imp);
					}
					else{
						// do nothing 
					}
					break;
				default:
					System.out.println("s: save file; \nr: delete last point");
				}
			} 
		}
	}

	public void savePointsToFile(File file){
		CSVWriter writer = null;
		String [] nextLine;

		try {
			writer = new CSVWriter(new FileWriter(file), '\t');
		} catch (Exception e) {
			e.printStackTrace();
		}

		try {	
			for (double[] point : points){
				nextLine = Arrays.toString(point).replaceAll("\\[|\\]|\\s", "").split(",");
				writer.writeNext(nextLine);
			}
			writer.close();			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args){
		new ImageJ();
		new ManualPointSelection();
	}

}
