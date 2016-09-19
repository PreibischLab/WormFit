package klim;

public class Psf {
	private final long[] min; 
	private final long[] max; 
	private final long[] size;
	private long index;
	private int numDimensions;
	
	public Psf(long[] min, long[] max, long[] size, long index){
		this.numDimensions = min.length;
		this.min = new long[min.length];
		this.max = new long[max.length];
		this.size = new long[size.length];
		this.index = index;
		
		for (int d = 0; d < min.length; d++){
			this.min[d] = min[d];
			this.max[d] = max[d];
			this.size[d] = size[d];
		}
	}
}
