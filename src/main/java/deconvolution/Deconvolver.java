package deconvolution;

import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public interface Deconvolver 
{	
	public boolean kernelsUpdated();
	public String getName();
	public double getAvg();
	public LRInput getData();
	
	public Image<FloatType> getPsi();
	public DeconvolveRuntimeStatistics runIteration();
	public void setDebug( boolean debug );
}
