package coherent.point.drift.refactored;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import util.opencsv.CSVReader;
import util.opencsv.CSVWriter;

public class IOUtils {
	// TESTED!
	public static int getNumLines(String filepath, char delimiter) {
		int numLines = 0;
		try {
			CSVReader reader = new CSVReader(new FileReader(filepath), delimiter); // we don't care about the delimiter
			for (;reader.readNext() != null; numLines++);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return numLines;
	}

	// TESTED!
	public static int getNumDimensions(String filepath, char delimiter) {
		int numDimensions = 0;
		try {
			CSVReader reader = new CSVReader(new FileReader(filepath), delimiter); // we don't care about the delimiter
			String[] nextLine = reader.readNext();
			numDimensions = nextLine.length; 
		} catch (IOException e) {
			e.printStackTrace();
		}
		return numDimensions;
	}

	// TESTED!
	public static RealMatrix readPositionsFromCSV(String filepath, char delimiter) {
		String[] nextLine;

		int numLines = getNumLines(filepath, delimiter); 
		int numDimensions = getNumDimensions(filepath, delimiter);
		RealMatrix mA = MatrixUtils.createRealMatrix(numLines, numDimensions);
		try {
			CSVReader reader = new CSVReader(new FileReader(filepath), delimiter);
			int i = 0;
			while ((nextLine = reader.readNext()) != null) {
				for (int d = 0; d < numDimensions; d++)
					mA.setEntry(i, d, Double.parseDouble(nextLine[d]));
				i++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		return mA;
	}

	// TESTED!
	public static void writePositionsToCSV(RealMatrix mA, String filepath, char delimiter) {
		try {
			CSVWriter writer = new CSVWriter(new FileWriter(filepath), delimiter);
			String[] nextLine = new String[mA.getColumnDimension()];
			
			for (int i = 0; i < mA.getRowDimension(); ++i) {
				for (int d = 0; d < mA.getColumnDimension(); ++d)
					nextLine[d] = Double.toString(mA.getEntry(i, d));
				writer.writeNext(nextLine);
			}
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// testWritePositionsToCSV();
	}
	
}