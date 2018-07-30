package cpd;

import static org.junit.Assert.*;

import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

public class CPDTest {

	RealMatrix resMatlab;
	RealMatrix mX;
	RealMatrix mY;
	
	@Before
	public void readData(RealMatrix yMatlab, RealMatrix X, RealMatrix Y){
//		String from;
//		String to;
//		new ApacheCPD().readCSV(from, to);
		
	}
	
	@Before
	public void useAllData(){
//		String []
	}
	
/*	
	@Test
	public void nonRigidRegistrationTest(){
		ApacheCPD cpd = new ApacheCPD();// change to parameters and load the data 
		
		// compare result for the= matlab and you implementation
		assertEquals( 1.0, 1.0 );
	}
*/
}
