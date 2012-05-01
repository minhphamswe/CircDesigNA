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
 * Implementation of an standard tournament designer, where each member is allowed
 * to consume resources (mutation tests) for a finite amount of cycles. 
 * Then the fittest members "reproduce" at some probability, P. To keep the population size constant,
 * for every reproduction, the least fit member is dropped.
 * 
 * @author Benjamin
 */
public class StandardTournament <T extends PopulationDesignMember<T>>  extends TournamentDesigner <T> {
	public StandardTournament(SingleMemberDesigner<T> SingleDesigner, double d) {
		super(SingleDesigner);
		designTime = (long) (d*1e9);
	}
	private int numElites = 1;
	private long designTime;
	public void runBlockIteration_(CircDesigNA runner, double endThreshold) {
		long now = System.nanoTime();
		while(true){
			for(int i = 0; i < populationSize; i++){
				boolean mutationSuccessful = SingleDesigner.mutateAndTestAndBackup(population_mutable[i]);
				//System.out.println(mutationSuccessful);
				if(runner!=null && runner.abort){
					return; //Abort
				}	

				setProgress((i+1), populationSize);
			}
			tournamentSelect(numElites);
			if (System.nanoTime()-now > designTime){
				break; //Timeup
			}
		}
	}
}
