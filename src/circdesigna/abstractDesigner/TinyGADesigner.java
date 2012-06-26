/*
  Part of the CircDesigNA Project - http://cssb.utexas.edu/circdesigna
  
  Copyright (c) 2010-11 Ben Braun
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation, version 2.1.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General
  Public License along with this library; if not, write to the
  Free Software Foundation, Inc., 59 Temple Place, Suite 330,
  Boston, MA  02111-1307  USA
*/
package circdesigna.abstractDesigner;

import circdesigna.CircDesigNA;


/**
 * An improved version of GADesigner, inspired by Tiny-GA in
 * "A field guide to genetic programming" Riccardo Poli
 * 
 * Implements MOGA.
 */
public class TinyGADesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	//Parameters used in Tiny-GA
	private int POPSIZE;
	private int TSIZE;
	private int NUM_CHILDREN_PER_GENERATION; 
	
	public TinyGADesigner(SingleMemberDesigner<T> SingleDesigner, ParetoSort<T> psort, boolean isMOGA) {
		super(SingleDesigner);
		this.psort = psort;
		this.isMOGA = isMOGA;
		
		//Settings used in the sin(x) example from the same book.
		TSIZE = 2; //Arbitrary!
		CROSSOVER_PROB = .9;
		SingleDesigner.setMutationProbabilities(1, 1);
	}
	private double CROSSOVER_PROB;
	private boolean isMOGA;
	private ParetoSort psort;
	private T tempOffspring;
	private T persistentBest;
	
	@Override
	public void initialize(T init, int numCopies) {
		POPSIZE = numCopies;
		NUM_CHILDREN_PER_GENERATION = POPSIZE; //Arbitrary
		tempOffspring = init.designerCopyConstructor(numCopies+2);
		super.initialize(init, numCopies+1);
		persistentBest = population_mutable[numCopies];
	}
	private int rnegtourn(int tsize) {
		return rtourn_(tsize,true);
	}
	private int rtourn(int tsize) {
		return rtourn_(tsize,false);
	}
	private int rtourn_(int tsize, boolean negativeTournament){
		int toRet = (int) (Math.random()*POPSIZE);
		double bestScore = SingleDesigner.getOverallScore(population_mutable[toRet]);
		for(int i = 0; i < tsize; i++){
			int competitor = (int) (Math.random()*POPSIZE);
			double competitorScore = SingleDesigner.getOverallScore(population_mutable[competitor]);
			if ((competitorScore < bestScore) ^ negativeTournament){
				toRet = competitor;
				bestScore = competitorScore;
			}
		}
		return toRet;
	}
	/**
	 * Returns how many strings were evaluated in this iteration.
	 */
	public void runBlockIteration_ (CircDesigNA runner, double endThreshold) {
		updateBestChild(); //Make sure bestChild is current.
		
		for(int i = 0; i < NUM_CHILDREN_PER_GENERATION; i++){
			//A new child will be born, by either SEXUAL or ASEXUAL reproduction.
			//This child will immediately knock out a randomly selected member of the population.
			if (Math.random() < CROSSOVER_PROB){ //Strange use of the word "Probability"
				//SEXUAL
				int parent1 = rtourn(TSIZE);
				int parent2 = rtourn(TSIZE);
				//Crossover:
				SingleDesigner.fourPtCrossoverAndEval(population_mutable[parent1],population_mutable[parent2],tempOffspring);
				runner.evaluated_strings++;
			} else {
				//ASEXUAL
				int parent1 = rtourn(TSIZE);
				SingleDesigner.mutateAndEval(population_mutable[parent1], tempOffspring);
				runner.evaluated_strings++;
			}
			
			//Kick out the result of a negative tournament with tempoffspring
			AddToPopulation(rnegtourn(TSIZE));
			setProgress((i+1), NUM_CHILDREN_PER_GENERATION);
		}
		/*
		//Adjust for pareto fitness
		if (isMOGA){
			psort.adjustFitness(redist,0,populationSize*2,maxFitness+1,1);
		}
		//Redistribute members to mutable and backups
		Arrays.sort(redist);
		setBestChild(redist[0].myKey);
		for(int i = 0; i < populationSize; i++){
			population_mutable[i] = redist[i].myKey;
			population_backups[i] = redist[i+populationSize].myKey;
		}
		*/
	}
	private void AddToPopulation(int victim) {

		//Dedicate all penalties:
		//Knock out a random member. Note that the mutated member ALWAYS enters the population!
		int knockout = victim;

		double nScore = SingleDesigner.getOverallScore(tempOffspring);
		double bestScore = SingleDesigner.getOverallScore(getBestPerformingChild());
		
		//Don't kick out the best population member unless the new guy is better
		if (population_mutable[knockout] == getBestPerformingChild() 
				&& nScore >= bestScore){
			return;
		}
		
		//Swap it out.
		T temp = tempOffspring;
		tempOffspring = population_mutable[knockout];
		population_mutable[knockout] = temp;
		
		//Does the introduced member beat the old best?
		if (nScore < bestScore){
			setBestChild(population_mutable[knockout]);
		}

	}
	private void updateBestChild() {
		T best = null; double bestScore = Double.MAX_VALUE;
		for(int i = 0; i < population_mutable.length; i++){
			double score = SingleDesigner.getOverallScore(population_mutable[i]);
			//System.out.println(i+" "+score+" "+population_mutable[i]);
			if (score < bestScore){
				best = population_mutable[i];
				bestScore = score;
			}
		}
		
		setBestChild(best);
	}
}