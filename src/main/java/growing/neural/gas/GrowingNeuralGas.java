package growing.neural.gas;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import net.imglib2.Point;

public class GrowingNeuralGas {

	long N = 50;

	long MaxIt = 20;

	long L = 50;

	double epsilon_b = 0.2;
	double epsilon_n = 0.005;

	double alpha = 0.5;
	double delta = 0.995;

	long T = 50;
	
	// constructor
	public GrowingNeuralGas(){
		// initialise parameters here
	}
	
	public double dist(double[] x, double [] y){
		double res = 0;
		for (int d = 0; d < x.length; ++d){
			res += (x[d] - y[d])*(x[d] - y[d]);
		}
		res = Math.sqrt(res);
		return res;
	}
	
	public void runGNG(List<Point> x){
		long nData = x.size(); // number of points you have
		int nDim = x.get(0).numDimensions(); // TODO: check for zero length
		
		Collections.shuffle(x);
		
		Point xMin = Collections.min(x);
		Point xMax = Collections.max(x);
		
	    int Ni = 2;    
	    double[] w = new double [Ni*nDim];
	    
	    for (int i = 0; i < Ni; ++i){
	    	for (int d = 0; d < nDim; ++d){
	    		double val = ThreadLocalRandom.current().nextInt(xMin.getIntPosition(d), xMax.getIntPosition(d) + 1);
		    	w[i*nDim + d] = val; // TODO: ?
	    	}
	    }
	    
	    // TODO: figure out what is what 
	    double [] E = new double[Ni]; // error
	    double [] C = new double[Ni*Ni]; // links
	    double [] t = new double[Ni*Ni]; // TODO: ?
	    		
	    long nx = 0;

	    for (int it = 0; it < MaxIt; ++it){
	    	// Fucking dickheads use l as a variable!!!!
	    	for (int l = 0; l < nData; ++l ){
	    		++nx; 
	    		// here we work with l'th row; 
	    		
	    		double[] d = new double[];
	    		
	    		for (int i = 0; i < x.get(l).numDimensions(); ++i){
	    			for (int j = 0; j < nDim; ++j){	
	    				// each with each multiplication
	    				
	    			}
	    		}
	    		
	    	}
	    }

	    
//
//	    for it = 1:MaxIt
//	        for l = 1:nData
//	            % Select Input
//	            nx = nx + 1;
//	            x = X(l,:);
//
//	            % Competion and Ranking
//	            d = pdist2(x, w);
//	            [~, SortOrder] = sort(d);
//	            s1 = SortOrder(1);
//	            s2 = SortOrder(2);
//
//	            % Aging
//	            t(s1, :) = t(s1, :) + 1;
//	            t(:, s1) = t(:, s1) + 1;
//
//	            % Add Error
//	            E(s1) = E(s1) + d(s1)^2;
//
//	            % Adaptation
//	            w(s1,:) = w(s1,:) + epsilon_b*(x-w(s1,:));
//	            Ns1 = find(C(s1,:)==1);
//	            for j=Ns1
//	                w(j,:) = w(j,:) + epsilon_n*(x-w(j,:));
//	            end
//
//	            % Create Link
//	            C(s1,s2) = 1;
//	            C(s2,s1) = 1;
//	            t(s1,s2) = 0;
//	            t(s2,s1) = 0;
//
//	            % Remove Old Links
//	            C(t>T) = 0;
//	            nNeighbor = sum(C);
//	            AloneNodes = (nNeighbor==0);
//	            C(AloneNodes, :) = [];
//	            C(:, AloneNodes) = [];
//	            t(AloneNodes, :) = [];
//	            t(:, AloneNodes) = [];
//	            w(AloneNodes, :) = [];
//	            E(AloneNodes) = [];
//
//	            % Add New Nodes
//	            if mod(nx, L) == 0 && size(w,1) < N
//	                [~, q] = max(E);
//	                [~, f] = max(C(:,q).*E);
//	                r = size(w,1) + 1;
//	                w(r,:) = (w(q,:) + w(f,:))/2;
//	                C(q,f) = 0;
//	                C(f,q) = 0;
//	                C(q,r) = 1;
//	                C(r,q) = 1;
//	                C(r,f) = 1;
//	                C(f,r) = 1;
//	                t(r,:) = 0;
//	                t(:, r) = 0;
//	                E(q) = alpha*E(q);
//	                E(f) = alpha*E(f);
//	                E(r) = E(q);
//	            end
//
//	            % Decrease Errors
//	            E = delta*E;
//	        end
//
//	        % Plot Results
//	        if PlotFlag
//	            figure(1);
//	            PlotResults(X, w, C);
//	            pause(0.01);
//	        end
//	    end
//
//	    %% Export Results
//	    net.w = w;
//	    net.E = E;
//	    net.C = C;
//	    net.t = t;
	}
	
	public static void main(String[] args) {
		

	}

}
