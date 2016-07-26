package deconvolution;

public class DeconvolveRuntimeStatistics
{
	public double maxChange = -1;
	public double sumChange = -1;
	
	public DeconvolveRuntimeStatistics( final double maxChange, final double sumChange )
	{
		this.maxChange = maxChange;
		this.sumChange = sumChange;
	}
}
