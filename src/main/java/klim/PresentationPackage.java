package klim;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import net.imglib2.Cursor;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.region.hypersphere.HyperSphere;
import net.imglib2.algorithm.stats.Normalize;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

import ij.ImageJ;
import util.ImgLib2Util;
import util.opencsv.CSVReader;

public class PresentationPackage {

	// merges two images with the split
	public static void mergeImages(Img<FloatType> iImg, Img<FloatType> dImg,  Img<FloatType> rImg){ // , long I, long J, long K){

		int numDimensions = iImg.numDimensions(); 

		Cursor<FloatType> cursor = iImg.cursor();
		RandomAccess<FloatType> dRa = dImg.randomAccess();
		RandomAccess<FloatType> rRa = rImg.randomAccess();



		// long diagonal = 4*(iImg.dimension(0) + iImg.dimension(1))/6;

		while (cursor.hasNext()){
			cursor.fwd();

			dRa.setPosition(cursor);
			rRa.setPosition(cursor);

			long [] position = new long [numDimensions];
			cursor.localize(position);

			if (position[0] <= iImg.dimension(0)/2){
				rRa.get().set(cursor.get());
			}		
			else{
				rRa.get().set(dRa.get());

			}
		}

		ImageJFunctions.show(rImg).setTitle("Fine?");

	}

	public static void createModel(Img<FloatType> rImg, ArrayList<double[]> nucleiPos, ArrayList<double[]> mRnaPos, long scale, boolean showNuclei, boolean showMRna){
		// debug
		for (int j = 0; j < nucleiPos.size(); j++){
			for (int d = 0; d < 3; d++){
				System.out.print(nucleiPos.get(j)[d] + " ");
			}
			System.out.println();
		}	


		// TODO: pass the color and the size of the spheres

		if (showMRna){
			long radius = 2*scale;
			int intensity = 128;
			drawSpheres(rImg, mRnaPos, intensity, radius, scale);
		}		

		if (showNuclei){
			long radius = 8*scale;
			int intensity = 10;
			drawSpheres(rImg, nucleiPos, intensity, radius, scale);
		}	
		ImageJFunctions.show(rImg);
	}

	// reads the data from the csv file 
	public static ArrayList<double[]> readCSV(String from, long [] dimensions) {
		CSVReader reader = null;
		String[] nextLine;
		ArrayList<double[]> positions = new ArrayList<>();

		int numDimensions = 3;

		try {
			reader = new CSVReader(new FileReader(from), ',');
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		try {
			int spotId = 0;
			// read the header
			if ((nextLine = reader.readNext()) != null){
				for (int d = 0; d < numDimensions; d++){
					System.out.println(nextLine[d]);
					dimensions[d] = Long.parseLong(nextLine[d]);
				}
			}	
			while ((nextLine = reader.readNext()) != null) {
				positions.add(new double[numDimensions]);
				for (int d = 0; d < numDimensions; d++){
					positions.get(spotId)[d] = Double.parseDouble(nextLine[d]);
				}
				spotId++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return positions;

	}

	/**
	 * Draws a sphere that contains lots of small spheres into the center of the interval
	 *
	 * @param randomAccessible - the image data to write to
	 * @param minValue - the minimal intensity of one of the small spheres
	 * @param maxValue - the maximal intensity of one of the small spheres
	 */
	public static void drawSpheres(
		final RandomAccessibleInterval< FloatType > rImg, ArrayList<double[]> positions, float intensity, long radius, long scale)
	{
		// the number of dimensions
		int numDimensions = rImg.numDimensions();

		RandomAccess< FloatType> ra = Views.extendMirrorDouble(rImg).randomAccess();

		long total = 0;

		for (int idx = 0; idx < positions.size(); ++idx){
			long [] pos = new long [numDimensions]; 

			double ratioXZ = 0.2/0.13; // ratio between x and z scaling

			for (int d = 0; d < numDimensions; d++){
				pos[d] = scale*Math.round(positions.get(idx)[d]);

				if (d == numDimensions - 1) 
					pos[d] *= ratioXZ;

			}

			// a model drop wrong detections
			if (pos[0] + pos[1] <= 370) continue;
			if (pos[0] > 817 && pos[1] < 170) continue;


			// int radius = 10;
			ra.setPosition(pos);      	
			HyperSphere< FloatType> sphere =
					new HyperSphere< FloatType >( Views.extendMirrorDouble(rImg), ra, radius );

			for ( final FloatType value : sphere ){     	
				value.setReal( intensity);
			}

			total++;
		}
		System.out.println(total + "!");
	}


	public static void startCreateModel(){

		new ImageJ();
		String fileNuclei = "/Users/kkolyva/Dropbox/PhD/2017-07-13-lab-meeting/data/dog/points/d.csv";
		String fileFISH = "/Users/kkolyva/Dropbox/PhD/2017-07-13-lab-meeting/data/smFISH/points/a.csv";
		long [] dimensions  = new long [3];


		boolean showNuclei = false; 
		boolean showMRna = true;

		ArrayList<double[]> nucleiPos = new ArrayList<>();
		ArrayList<double[]> mRnaPos = new ArrayList<>();

		if (showNuclei)	{
			nucleiPos = readCSV(fileNuclei, dimensions);
		}
		if (showMRna){
			mRnaPos = readCSV(fileFISH, dimensions);
		}
		// Img<FloatType> rImg = ImgLib2Util.openAs32Bit(new File("/Users/kkolyva/Dropbox/PhD/2017-07-13-lab-meeting/data/b-model-full.tif")); 

		long scale = 2;

		scaleImage(dimensions, scale);
		dimensions[2] *= 0.2/0.13; // ratio between x and z scaling

		Img<FloatType> rImg = new ArrayImgFactory<FloatType>().create(dimensions, new FloatType());
		createModel(rImg, nucleiPos, mRnaPos, scale, showNuclei, showMRna);

	}

	public static void scaleImage(long[] dimensions, long scale){
		for (int j = 0; j < dimensions.length; ++j)
			dimensions[j] *= scale; 
	}	


	public static void startMergeImages(){
		new ImageJ();

		File [] paths = new File[2];
		paths[0] = new File( "/Users/kkolyva/Dropbox/PhD/2017-07-13-lab-meeting/data/SEA12_dpy23_wdr52_mdh1_005.nd2 - SEA12_dpy23_wdr52_mdh1_005.nd2 (series 10) - C=4.tif" );
		paths[1] = new File( "/Users/kkolyva/Dropbox/PhD/2017-07-13-lab-meeting/data/SEA12_dpy23_wdr52_mdh1_005.nd2 - SEA12_dpy23_wdr52_mdh1_005.nd2 (series 10) - C=4 2.tif" );

		Img<FloatType> iImg = ImgLib2Util.openAs32Bit(paths[0]);
		Img<FloatType> dImg = ImgLib2Util.openAs32Bit(paths[1]);

		Img<FloatType> rImg = new ArrayImgFactory<FloatType>().create(iImg, iImg.cursor().get());

		FloatType minT = new FloatType(0);
		FloatType maxT = new FloatType(1);

		Normalize.normalize(iImg, minT, maxT);
		Normalize.normalize(dImg, minT, maxT);

		mergeImages(iImg, dImg, rImg);
	}


	public static void main(String [] args){
		// startMergeImages();
		startCreateModel();


		System.out.println("Doge!");
	}
}
