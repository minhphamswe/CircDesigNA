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
public class StochasticOptDesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	
	public StochasticOptDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
		T_MAX = 10;
		CHILDREN_PER_ITERATION = 100;
	}
	private T tempOffspring;
	private int CHILDREN_PER_ITERATION;
	private double T_MAX;
	
	public void initialize(T init, int numCopies) {
		//clone once, to ID numCopies + 1.
		tempOffspring = init.designerCopyConstructor(numCopies+1);
		super.initialize(init, numCopies);
	}
	public void runBlockIteration_ (CircDesigNA runner, double endThreshold) {
		for(int child = 0; child < populationSize; child++){
			//Linearly cool down temperature, from T_MAX to zero.
			for(int i = 0; i < CHILDREN_PER_ITERATION; i++)
			{
				SingleDesigner.mutateAndEval(population_mutable[child], tempOffspring);
				runner.evaluated_strings++;

				double nScore = SingleDesigner.getOverallScore(tempOffspring);
				double oScore = SingleDesigner.getOverallScore(population_mutable[child]);
				double T = (CHILDREN_PER_ITERATION-1-i) / (double)CHILDREN_PER_ITERATION * T_MAX;

				if (Math.random() < P(oScore, nScore, T)){
					//Swap it out.
					T temp = tempOffspring;
					tempOffspring = population_mutable[child];
					population_mutable[child] = temp;
				}

			}
			setProgress(child+1, populationSize);
		}
		updateBestChild();
	}

	private void updateBestChild() {
		int best = 0; double bestScore = Double.MAX_VALUE;
		for(int i = 0; i < population_mutable.length; i++){
			double score = SingleDesigner.getOverallScore(population_mutable[i]);
			//System.out.println(i+" "+score+" "+population_mutable[i]);
			if (score < bestScore){
				best = i;
				bestScore = score;
			}
		}
		setBestChild(population_mutable[best]);
	}
	private double P(double e, double eprime, double T) {
		if (eprime < e){
			return 1;
		}
		if (T <= 0){
			return 0;
		}
		return Math.exp((e - eprime)/T);
	}
}