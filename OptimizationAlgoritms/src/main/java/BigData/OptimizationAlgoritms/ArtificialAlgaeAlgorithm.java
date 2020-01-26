package BigData.OptimizationAlgoritms;

import java.util.ArrayList;
import java.util.Collections;

public class ArtificialAlgaeAlgorithm {

	static int N; //Population
	static int D; //Dimension
	static int maxiter; // İteration
	static int iter;
	static double cutting;
	static double energyLoss;
	static double adaptationConstant;
	static double score;
	
	static double Lower; // Min Value
	static double Upper; // Max Value
	static double positions[][]; // Population Position
	static double positionsScore[];
	static double colonySize[];
	static double colonyStarving[];
	static double colonyEnergy[];
	static double bestVal;
	static double bestPosition[];
	static double convergenceCurve[];
	static double colonyFriction[];
	
	public ArtificialAlgaeAlgorithm(int iter, int UpLevel, int LowLevel, int searchAgent, int dimension, double cuttingValue, double energyLossValue, double AC) {
		
		maxiter = iter;
		N = searchAgent;
		D = dimension;
		Lower = LowLevel;
		Upper = UpLevel;
		cutting = cuttingValue;
		energyLoss = energyLossValue;
		adaptationConstant = AC;
		score = Double.MAX_VALUE;
		positions = new double[N][D];
		positionsScore = new double[N];
		colonySize = new double[N];
		colonyStarving = new double[N];
		colonyEnergy = new double[N];
		colonyFriction = new double[N]; 
		bestPosition = new double[D];
		convergenceCurve = new double[maxiter];
	}
	// uygunluk fonksiyonu
	static double bechmark(double[] position) {
		double result = 0;
		for(int i = 0; i < D; i++) {
			result += Math.pow(position[i], 2);
		}
		return result;
	}
	static void init() {
		for(int i =0; i < N; i++) {
			for(int j = 0; j < D; j++) {
				positions[i][j] = Lower + (Upper - Lower) * Math.random();
			}
			positionsScore[i] = bechmark(positions[i]);
			colonySize[i] = 1;
			colonyStarving[i] = 0;
		}
	}
	static double[][] solution(){
		init();		
		bestVal = positionsScore[minMaxIndex(positionsScore)[0]];
		if(bestVal < score) {
			score = bestVal;
			bestPosition = positions[minMaxIndex(positionsScore)[0]].clone();
		}
		convergenceCurve[0] = score;
		colonySize = calculateGreatness(colonySize, positionsScore);
		for(int iter = 1; iter < maxiter; iter++) {
			colonyEnergy = calculate_energy(colonySize);
			colonyFriction = frictionForce(colonySize);
			for(int i = 0; i < N; i++) {
				int iStarve = 0;
				while(colonyEnergy[i] >= 0) {
					int neighbor = tournamentMethod(positionsScore);
					while(neighbor == i) {
						neighbor = tournamentMethod(positionsScore);
					}
					int dim1 = (int)(D * Math.random());
					int dim2 = (int)(D * Math.random());
					int dim3 = (int)(D * Math.random());
					while(dim1 == dim2 || dim1 == dim3 || dim2 == dim3) {
						dim2 = (int)(D * Math.random());
						dim3 = (int)(D * Math.random());
					}
					double[] newColony = new double[D];
					newColony = positions[i].clone();
					double p = -1 + (2 * Math.random());
					int degree1 = (int)(360 * Math.random()); 
					int degree2 = (int)(360 * Math.random()); 
					newColony[dim1] = newColony[dim1] + (positions[neighbor][dim1] - newColony[dim1]) * (cutting - colonyFriction[i]) * p;
					newColony[dim1] = simplebounds(newColony[dim1]);
					newColony[dim2] = newColony[dim2] + (positions[neighbor][dim2] - newColony[dim2]) * (cutting - colonyFriction[i]) * Math.cos(Math.toRadians(degree1));
					newColony[dim2] = simplebounds(newColony[dim2]);
					newColony[dim3] = newColony[dim3] + (positions[neighbor][dim3] - newColony[dim3]) * (cutting - colonyFriction[i]) * Math.sin(Math.toRadians(degree2));
					newColony[dim3] = simplebounds(newColony[dim3]);
					
					double newScore = bechmark(newColony);
					colonyEnergy[i] = colonyEnergy[i] - (energyLoss /2);
					if(newScore < positionsScore[i]) {
						positions[i] = newColony.clone();
						positionsScore[i] = newScore;
						iStarve = 1;
					}else {
						colonyEnergy[i] = colonyEnergy[i] - (energyLoss /2);
					}
				}
				bestVal = positionsScore[minMaxIndex(positionsScore)[0]];
				if(bestVal < score) {
					score = bestVal;
					bestPosition = positions[minMaxIndex(positionsScore)[0]].clone();
				}
				if(iStarve == 0) {
					colonyStarving[i] = colonyStarving[i] + 1;
				}
			}
			//Evrimsel Süreç
			colonySize = calculateGreatness(colonySize, positionsScore);
			int dim = (int)(D * Math.random());
			int[] index = new int[2];
			index = minMaxIndex(colonySize);
			positions[index[0]][dim] = positions[index[1]][dim];
			
			//Adaptanson İşlemi
			int maxStarving = minMaxIndex(colonyStarving)[1];
			if(Math.random() < adaptationConstant) {
				for(int i = 0; i < D; i++) {
					positions[maxStarving][i] = positions[maxStarving][i] + (bestPosition[i] - positions[maxStarving][i]) * Math.random(); 
				}
			}
			bestVal = positionsScore[minMaxIndex(positionsScore)[0]];
			if(bestVal < score) {
				score = bestVal;
				bestPosition = positions[minMaxIndex(positionsScore)[0]].clone();
			}
			convergenceCurve[iter] = score;
		}
		double[][] out = new double[2][D];
		for(int i = 0; i < D; i++) {
			out[1][i] = bestPosition[i];
		}
		out[0][0] = bechmark(bestPosition);
		return out;
	}
	//Herbir koloninin enerjisinin hesaplanması
	static double[] normalization(double[] size) {
		ArrayList<Double> sortedScores = new ArrayList<Double>();
		for(int i = 0; i < size.length; i++) {
			sortedScores.add(size[i]);
		}
		Collections.sort(sortedScores);
		double maxValue = sortedScores.get(sortedScores.size() - 1);
		double minValue = sortedScores.get(0);
		double energy[] = new double[N];
		for(int i = 0; i < size.length; i++) {
			energy[i] = (size[i] - minValue) / (maxValue - minValue);
		}
		return energy;
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
	static double[] calculateGreatness(double[] greatness, double[] scores) {
		double[] newGreatness = new double[N];
		ArrayList<Double> sortedScores = new ArrayList<Double>();
		for(int i = 0; i < scores.length; i++) {
			sortedScores.add(scores[i]);
		}
		Collections.sort(sortedScores);
		double maxValue = sortedScores.get(sortedScores.size() - 1);
		double minValue = sortedScores.get(0);
		double normalScore[] = new double[N];
		for(int i = 0; i < scores.length; i++) {
			normalScore[i] = (scores[i] - minValue) / (maxValue - minValue);
			normalScore[i] = 1 - normalScore[i]; //min bulmaya çalışıldığı için
		}
		for(int i = 0; i < scores.length; i++) {
			double fKs = Math.abs(greatness[i] / 2); //half saturation constant
			double M = normalScore[i] / (fKs + normalScore[i]);
			double dX = M * greatness[i]; //greatness rate
			newGreatness[i] = greatness[i] + dX;
		}
		return newGreatness;
	}
	static double[] calculate_energy(double[] greatness) {
		int[] sort = new int[greatness.length];
		double[] fGreatSurface = new double[greatness.length];
		for(int i = 0; i < greatness.length; i++) {
			sort[i] = i;
		}
		for(int i = 0; i < (greatness.length - 1); i++) {
			for(int j = (i + 1); j < greatness.length; j++) {
				if(greatness[sort[i]] > greatness[sort[j]]) {
					int value;
					value = sort[i];
					sort[i] = sort[j];
					sort[j] = value;
				}
			}
			fGreatSurface[sort[i]] = Math.pow(i, 2);
		}
		fGreatSurface[sort[(greatness.length - 1)]] = Math.pow((greatness.length - 1), 2);
		int[] indis = new int[2];
		indis = minMaxIndex(fGreatSurface);
		double minValue = fGreatSurface[indis[0]];
		double maxValue = fGreatSurface[indis[1]]; 
		for(int i = 0; i < fGreatSurface.length; i++) {
			fGreatSurface[i] = (fGreatSurface[i] - minValue) / (maxValue - minValue);
		}
		return fGreatSurface;
	}
	static double[] frictionForce(double[] size) {
		double[] friction = new double[N];
		for(int i = 0; i < N; i++) {
			double r = (double)(size[i] * 3) / (double)(4 * Math.PI);
			r = Math.pow(r,  (1.0/3));
			friction[i] = 2 * Math.PI * Math.pow(r, 2);
		}
		friction = normalization(friction); //normalizasyon amacli
		return friction;
	}
	static int tournamentMethod(double[] score) {
		int colonyOne = (int)((N - 1) * Math.random());
		int colonyTwo = (int)((N - 1) * Math.random());
		while(colonyOne == colonyTwo) {
			colonyTwo = (int)((N - 1) * Math.random());
		}
		if(score[colonyOne] < score[colonyTwo]) {
			return colonyOne;
		}else {
			return colonyTwo;
		}
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
