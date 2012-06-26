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
 * Implementation of hill climbing
 */
public class HillClimbingDesigner <T extends PopulationDesignMember<T>> extends BlockDesigner<T>{
	public HillClimbingDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
		SingleDesigner.setMutationProbabilities(.1, .01);
		TIME_PER_ITERATION = (long)1e9;
	}
	public long TIME_PER_ITERATION;
	
	public void runBlockIteration_(CircDesigNA runner, double endThreshold) {
		long now = System.nanoTime();
		while(true){
			for(int i = 0; i < populationSize; i++){
				boolean mutationSuccessful = SingleDesigner.mutateAndTestAndBackup(population_mutable[i]);
				runner.evaluated_strings++;
				//System.out.println(mutationSuccessful);
				if(runner!=null && runner.abort){
					return; //Abort
				}
			}

			long dt = System.nanoTime() - now;
			int p = (int)(100 * dt / (double)TIME_PER_ITERATION);
			setProgress((p+1), 100);
			if (dt > TIME_PER_ITERATION){
				break;
			}
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
}
