package DnaDesign.AbstractDesigner;

import DnaDesign.DomainDesigner;

/**
 * Implementation of an standard tournament designer, where each member is allowed
 * to consume resources (mutation tests) for a finite amount of cycles. 
 * Then the fittest members "reproduce" at some probability, P. To keep the population size constant,
 * for every reproduction, the least fit member is dropped.
 * 
 * @author Benjamin
 */
public class StandardTournament <T extends PopulationDesignMember<T>>  extends TournamentDesigner <T> {
	public StandardTournament(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}
	private int numElites = 1;
	public void runBlockIteration_(DomainDesigner runner, double endThreshold) {
		long designTime = (long).1e9;
		for(int i = 0; i < populationSize; i++){
			long now = System.nanoTime();
			while(true){
				boolean mutationSuccessful = SingleDesigner.mutateAndTestAndBackup(population_mutable[i]);
				if(runner.abort){
					return; //OUT OUT OUT
				}	
				if (System.nanoTime()-now > designTime){
					break; //Timeup
				}
			}
		}

		tournamentSelect(numElites);
	}
}
