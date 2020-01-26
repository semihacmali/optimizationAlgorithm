package BigData.OptimizationAlgoritms;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args )
    {
    	int iter = 10000;
    	int UpLevel = 100;
    	int LowLevel = -100; 
    	int searchAgent = 50;
    	int dimension = 10;
    	
    	//Grey Wolf Optimizer
       GreyWolfOptimizer GO = new GreyWolfOptimizer(iter, UpLevel, LowLevel, searchAgent, dimension);
		double[][] best = GO.solution();
		System.out.println("Optimized value = " + best[0][0]);
		for(int i = 0; i < dimension; i++) {
			System.out.println("x["+i+"] = " + best[1][i]);
		}
	
    	//Artificial Algae Algorithm
    	//cutting = 2.0;
		//energyLoss = 0.3;
		//adaptationConstant = 0.2;
		double cutting = 2.0;
		double energyLoss = 0.3;
		double AC = 0.5;
    	
    	ArtificialAlgaeAlgorithm AAA = new ArtificialAlgaeAlgorithm(iter, UpLevel, LowLevel, searchAgent, dimension, cutting, energyLoss, AC);
    	double[][] best2 = AAA.solution();
    	System.out.println("Optimized value = " + best2[0][0]);
		for(int i = 0; i < dimension; i++) {
			System.out.println("x["+i+"] = " + best2[1][i]);
		}
		
		
    	//Harmony Search Algoritm
		//HMCR = 0.9;
		//PAR = 0.4;
		//BW = 0.2;
		double HMCRValue = 0.99;
		double PARValue = 0.5;
		double BWValue = 0.5;
	    HarmonySearchAlgoritms HSA = new HarmonySearchAlgoritms(iter, UpLevel, LowLevel, searchAgent, dimension, HMCRValue, PARValue, BWValue);
		double[][] best1 = HSA.solution();
		System.out.println("Optimized value = " + best1[0][0]);
		for(int i = 0; i < dimension; i++) {
			System.out.println("x["+i+"] = " + best1[1][i]);
		}
    }
}
