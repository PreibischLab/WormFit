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
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import util.ImgLib2Util;
import util.opencsv.CSVWriter;

/**
 * ManualPointSelection -- fast way to get the points manually 
 * and write them to txt after that
 * 
 */

public class ManualPointSelection {

	final ArrayList<double []> points =  new ArrayList<double[]>(1);

	Img <FloatType> img;
	MyMouseListener mouseListener;
	MyKeyListener keyListener;
	String dir = "/Users/kkolyva/Desktop/latest_desktop/";
	File txtOutput = new File(dir + "init-1.csv");

	public ManualPointSelection(){
		File file = new File(dir + "init-1.tif");
		Img<FloatType> img = ImgLib2Util.openAs32Bit(file);
		final ImagePlus imp = ImageJFunctions.wrapFloat(img, "");	

		imp.show();
		mouseListener = new MyMouseListener(imp);
		keyListener = new MyKeyListener();
		imp.getCanvas().addMouseListener(mouseListener);
		imp.getCanvas().addKeyListener(keyListener);	
	}

	public void addOverlay(ImagePlus imp){
		Overlay overlay = imp.getOverlay();
		if (overlay == null){
			// System.out.println("addOverlay: overlay is null!");	
			overlay = new Overlay();
			imp.setOverlay( overlay );
		}		
		overlay.clear();
		
		// here should be the T coordinate of the fitted points
		int numDimensions = 2;
		double [] sigma = new double[]{5, 5};
		double [] location = new double[numDimensions];
		for (double[] point : points){
			for (int d = 0; d < numDimensions; ++d){
				location[d] = point[d];
			} 
			final PointRoi or = new PointRoi(location[0], location[1]);
			or.setSize(100);
			overlay.add(or);
		}

//		imp.setOverlay(overlay);
//		imp.updateAndDraw();	
	}
	

	protected class MyMouseListener implements MouseListener
	{
		ImagePlus imp;
		
		public MyMouseListener(ImagePlus imp){
			this.imp = imp;
		}
		
		@Override
		public void mouseClicked(MouseEvent e) {
		}

		@Override
		public void mouseEntered(MouseEvent e) {}

		@Override
		public void mouseExited(MouseEvent e) {}

		@Override
		public void mousePressed(MouseEvent e) {}

		@Override
		public void mouseReleased( final MouseEvent e )
		{
			final int xm = imp.getCanvas().offScreenX( e.getX() );
			final int ym = imp.getCanvas().offScreenY( e.getY() );
			// imp.getCanvas().getMagnification()
			
			double [] location = new double[]{xm, ym};
			points.add(location);
			System.out.println("Added: ");
			System.out.println("x: " + location[0] + " y: " + location[1]);
			addOverlay(imp);
		}

	}

	protected class MyKeyListener implements KeyListener
	{
		public void keyTyped(KeyEvent e) {
			// nothing
		}

		public void keyPressed(KeyEvent e) {
			// nothing
		}

		public void keyReleased(KeyEvent e) {
			int id = e.getID();	
			if (id == KeyEvent.KEY_RELEASED) {
				char c = e.getKeyChar();
				switch (c){
				case 's':
					savePointsToFile(txtOutput);
					System.out.println("Saved: ");
					System.out.println(points.size() + " Points");
					break;
				case 'f':
					break;
				case 'r':
				case 'd':	
					if (points.size() - 1 > 0){
						System.out.println("Removed: ");
						System.out.println("x: " + points.get(points.size() - 1)[0] + " y: " + points.get(points.size() - 1)[1]);
						points.remove(points.size() - 1);
					}
					else{
						// do nothing 
					}
					break;
				default:
					System.out.print("s: save file; \nr: delete last point");
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
