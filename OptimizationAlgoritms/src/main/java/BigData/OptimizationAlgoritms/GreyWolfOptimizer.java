package BigData.OptimizationAlgoritms;



public class GreyWolfOptimizer {
	
	static int N; //Population
	static int D; //Dimension
	static int maxiter; // İteration
	static double a; // between 0 - 2 by iteration number
	static double r1; //random number
	static double r2; //random number
	static double alfa[]; // Best Position
	static double beta[]; // Second Position
	static double delta[]; // Third Position
	static double X1, X2, X3;
	static double A1, C1; //Alfa update value
	static double A2, C2; //Beta update value
	static double A3, C3; // delta uptade value
	static double Lower; // Min Value
	static double Upper; // Max Value
	static double positions[][]; // Population Position
	static double BestVal[];
	static double fitness, alfaScore, betaScore, deltaScore  = Double.MAX_VALUE;
	
	public GreyWolfOptimizer(int iter, int UpLevel, int LowLevel, int searchAgent, int dimension) {
		
		maxiter = iter;
		Lower = LowLevel;
		Upper = UpLevel;
		N = searchAgent;
		D = dimension;
		fitness = Double.MAX_VALUE;
		alfaScore = Double.MAX_VALUE;
		betaScore = Double.MAX_VALUE;
		deltaScore = Double.MAX_VALUE;
		positions = new double[N][D];
		alfa =  new double[D];
		beta =  new double[D];
		delta = new double[D];
		BestVal = new double[maxiter];
		
	}
	
	//Bechmark Function
	static double bechmark(double[] position) {
		double result = 0;
		for(int i = 0; i < D; i++) {
			result += Math.pow(position[i], 2);
		}
		return result;
	}
	// update best solution
	static void sort_and_index(double[][] position) {
		double score;
		for(int i = 0; i < N; i++) {
			score = bechmark(position[i]);
			if(score < alfaScore) {
				alfaScore = score;
				for(int j = 0; j < D; j++) {
					alfa[j] = position[i][j];
				}
			}
			if(score > alfaScore && score < betaScore) {
				betaScore = score;
				for(int j = 0; j < D; j++) {
					beta[j] = position[i][j];
				}
			}
			if(score > alfaScore && score > betaScore && score < deltaScore) {
				deltaScore = score;
				for(int j = 0; j < D; j++) {
					delta[j] = position[i][j];
				}
			}
		}
	}		
	//first iterition and initialization		
	static void init() {
		for(int i =0; i < N; i++) {
			for(int j = 0; j < D; j++) {
				positions[i][j] = Lower + (Upper - Lower) * Math.random();
			}
		}
		// Position Sort
		sort_and_index(positions);
		
		BestVal[0] = bechmark(alfa);
	}
	//update population positions
	static double[][] solution(){
		init();
		//int iter = 1;
		for(int iter = 1; iter < maxiter ; iter++) {
			a = 2.0 - ((double)iter * (2.0 / (double) maxiter));
			for( int i = 0; i < N; i++) {
				for(int j = 0; j < D; j++) {
					//Update Values for Alfa
					r1 = Math.random();
					r2 = Math.random();
					A1 = (2.0 * a * r1) - a;
					C1 = 2.0 * r2;
					
					//Update position by Alfa
					X1 = alfa[j] - A1 * (Math.abs(C1 * alfa[j] - positions[i][j]));
					
					//Update VAlues for beta
					r1 = Math.random();
					r2 = Math.random();
					A2 = (2.0 * a * r1) - a;
					C2 = 2.0 * r2;
					
					//Update position by Beta
					X2 = beta[j] - A2 * (Math.abs(C2 * beta[j] - positions[i][j]));
					
					//Update VAlues for Delta
					r1 = Math.random();
					r2 = Math.random();
					A3 = (2.0 * a * r1) - a;
					C3 = 2.0 * r2;
					
					//Update position by Delta
					X3 = delta[j] - A3 * (Math.abs(C3 * delta[j] - positions[i][j]));
					
					//update Population Positions
					positions[i][j] = (X1 + X2 + X3) / 3.0;
					positions[i][j] = simplebounds(positions[i][j]);
				}
			}
			sort_and_index(positions);

			BestVal[iter] = bechmark(alfa);
		}
		double[][] out = new double[2][D];
		for(int i = 0; i < D; i++) {
			out[1][i] = alfa[i];
		}
		out[0][0] = bechmark(alfa);
		return out;
	}
	
	static double simplebounds(double s){
		if(s < Lower) {
			s = Lower;
		}
		if(s > Upper) {
			s = Upper;
		}
		return s;
	}
}
