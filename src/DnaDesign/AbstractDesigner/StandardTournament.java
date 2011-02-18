package DnaDesign.AbstractDesigner;

import DnaDesign.DomainDesigner;
import DnaDesign.DomainDesigner.ScorePenalty;
import DnaDesign.impl.DomainDesignPMemberImpl;

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
	public void runBlockIteration_(DomainDesigner runner, double endThreshold) {
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
