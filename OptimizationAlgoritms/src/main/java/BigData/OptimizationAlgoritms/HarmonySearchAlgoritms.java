package BigData.OptimizationAlgoritms;

public class HarmonySearchAlgoritms {
	static int maxiter;
	static int iter;
	static int Lower;
	static int Upper;
	static int D; // harmony dimension
	static int N; // harmony size
	static double bestScore;
	static double[][] harmonyMemory;
	static double[] harmonyScore;
	static double[] bestVal;
	static double[] worstVal;
	static double[] bestHarmony;
	static double HMCR; // harmony Memory Considering Rate
	static double PAR; // Pitch Adjusting Rate
	static double BW; // Band Width
	
	public HarmonySearchAlgoritms(int itera, int UpLevel, int LowLevel, int searchAgent, int dimension, double HMCRValue, double PARValue, double BWValue) {
		maxiter = itera;
		Lower = LowLevel;
		Upper = UpLevel;
		N = searchAgent;
		D = dimension;
		bestScore = Double.MAX_VALUE;
		harmonyMemory = new double[N][D];
		harmonyScore = new double[N];
		bestVal = new double[maxiter];
		worstVal = new double[maxiter];
		bestHarmony = new double[D];
		HMCR = HMCRValue;
		PAR = PARValue;
		BW = BWValue;
		
	}
	// uygunluk fonksiyonu
		static double bechmark(double[] position) {
			double result = 0;
			for(int i = 0; i < D; i++) {
				result += Math.pow(position[i], 2);
			}
			return result;
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
		
		//Find Min Max İndex
			static int[] minMaxIndex(double[] iArray) {
				int maxIndex = 0;
				int minIndex = 0;
				for(int i = 1; i < iArray.length; i++) {
					if(iArray[maxIndex] < iArray[i]) {
						maxIndex = i;
					}
					if(iArray[minIndex] > iArray[i]) {
						minIndex = i;
					}
				}
				int[] indexes = new int[] {minIndex, maxIndex};
				return indexes;
			}
			
		static void init() {
			for(int i =0; i < N; i++) {
				for(int j = 0; j < D; j++) {
					harmonyMemory[i][j] = Lower + (Upper - Lower) * Math.random();
				}
				harmonyScore[i] = bechmark(harmonyMemory[i]);
			}
			
			
			bestVal[0] = bestScore = harmonyScore[minMaxIndex(harmonyScore)[0]];
			worstVal[0] = harmonyScore[minMaxIndex(harmonyScore)[1]];
		}
		static double[][] solution() {
			init();
			for(iter = 1; iter < maxiter; iter++) {
				double[] newHarmony = new double[D]; 
				
				for(int i = 0; i < D; i++) {
					if(Math.random() < HMCR) {
						newHarmony[i] =  harmonyMemory[(int)(N * Math.random())][i];
						if(Math.random() < PAR) {
							newHarmony[i] = pitchAdjust(newHarmony[i]);
						}
					}
					else {
						newHarmony[i] = Lower + (Upper - Lower) * Math.random();
					}
				}
				double newScore = bechmark(newHarmony);
				//updateHarmonyMemory(newScore, newHarmony);
				//find worst harmony
				int[] minMaxIn = minMaxIndex(harmonyScore);
				int worstIndex = minMaxIn[1];
				double worst = harmonyScore[worstIndex];
				worstVal[iter] = worst;
				
				if(newScore < worst) {
					harmonyMemory[worstIndex] = newHarmony.clone();
					harmonyScore[worstIndex] = newScore;
				}
				//find best harmony
				int bestIndex = minMaxIn[0];
				double best = harmonyScore[bestIndex];
				if(best != bestVal[iter -1 ]) {
					bestHarmony = harmonyMemory[bestIndex].clone();
					bestVal[iter] = best;
				}
				
			}
			double[][] out = new double[2][D];
			for(int i = 0; i < D; i++) {
				out[1][i] = bestHarmony[i];
			}
			out[0][0] = bechmark(bestHarmony);
			return out;
			
		}
		
		static double pitchAdjust(double index) {
			double newIndex;
			if(Math.random() < 0.5) {
				newIndex = index + Math.random() * BW;
			}
			else {
				newIndex = index - Math.random() * BW;
			}
			return(simplebounds(newIndex));
		}
		static void updateHarmonyMemory(double newFitness, double[] newHar) {
			//find worst harmony
			int[] minMaxIn = minMaxIndex(harmonyScore);
			int worstIndex = minMaxIn[1];
			double worst = harmonyScore[worstIndex];
			worstVal[iter] = worst;
			
			if(newFitness < worst) {
				harmonyMemory[worstIndex] = newHar.clone();
				harmonyScore[worstIndex] = newFitness;
			}
			//find best harmony
			int bestIndex = minMaxIn[0];
			double best = harmonyScore[bestIndex];
			if(iter > 0 && best != bestVal[iter -1 ]) {
				bestHarmony = harmonyMemory[bestIndex].clone();
				bestVal[iter] = best;
			}
		}
}
