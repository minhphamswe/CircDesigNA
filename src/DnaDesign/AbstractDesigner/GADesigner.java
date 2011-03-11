package DnaDesign.AbstractDesigner;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Comparator;

import DnaDesign.DomainDesigner;


/**
 * This is a Best N designer
 * (N mutated products, select down to N.)
 * 
 * Implements MOGA.
 * 
 * No sharing. No crossover.
 */
public class GADesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {

	public GADesigner(SingleMemberDesigner<T> SingleDesigner, ParetoSort<T> psort, boolean isMOGA) {
		super(SingleDesigner);
		this.psort = psort;
		this.isMOGA = isMOGA;
	}
	private boolean isMOGA;
	private ParetoSort psort;
	private T[] population_backups;
	private FitnessPopulationDesignMember<T>[] redist;
	
	//override initialize to also produce the N added members each iteration.
	@Override
	public void initialize(T init, int numCopies) {
		super.initialize(init, numCopies);

		population_backups = (T[])Array.newInstance(init.getClass(),populationSize);
		redist = (FitnessPopulationDesignMember<T>[])Array.newInstance(new FitnessPopulationDesignMember<T>().getClass(),populationSize+populationSize);
		for(int i = 0; i < populationSize; i++){
			population_backups[i] = init.designerCopyConstructor(numCopies+1+i);
		}
		for(int i = 0; i < redist.length; i++){
			redist[i] = new FitnessPopulationDesignMember<T>();
		}
	}
	
	public void runBlockIteration_ (DomainDesigner runner, double endThreshold) {
		double maxFitness = 0;
		for(int i = 0; i < populationSize; i++){
			T toMutate = population_mutable[i];
			T backup = population_backups[i];
			SingleDesigner.mutateAndTest(toMutate,backup);
			redist[i*2].myKey=toMutate;
			redist[i*2].myScore=SingleDesigner.getOverallScore(toMutate);
			redist[i*2+1].myKey=backup;
			redist[i*2+1].myScore=SingleDesigner.getOverallScore(backup);
			if (runner!=null && runner.abort){
				return;
			}
		}
		for(int i = 0; i < redist.length; i++){
			maxFitness = Math.max(maxFitness,redist[i].myScore);
		}
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
	}
}
	