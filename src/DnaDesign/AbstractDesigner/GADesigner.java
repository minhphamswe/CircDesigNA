package DnaDesign.AbstractDesigner;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Comparator;

import DnaDesign.DomainDesigner;

public class GADesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	public GADesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}
	private T[] population_backups;
	private T[] redist;
	
	//override initialize to also produce the N added members each iteration.
	@Override
	public void initialize(T init, int numCopies) {
		super.initialize(init, numCopies);

		population_backups = (T[])Array.newInstance((Class<T>)init.getClass(),populationSize);
		redist = (T[])Array.newInstance((Class<T>)init.getClass(),populationSize+populationSize);
		for(int i = 0; i < populationSize; i++){
			population_backups[i] = init.designerCopyConstructor(numCopies+1+i);
		}
	}
	
	public void runBlockIteration_ (DomainDesigner runner, double endThreshold) {
		for(int i = 0; i < populationSize; i++){
			T toMutate = population_mutable[i];
			T backup = population_backups[i];
			SingleDesigner.mutateAndTest(toMutate,backup);
			redist[i*2]=toMutate;
			redist[i*2+1]=backup;
		}
		//Redistribute members to mutable and backups
		Arrays.sort(redist, new Comparator<T>(){
			public int compare(T o1, T o2) {
				return ((Double)SingleDesigner.getOverallScore(o1)).compareTo((Double)SingleDesigner.getOverallScore(o2));
			}
		});
		setBestChild(redist[0]);
		for(int i = 0; i < populationSize; i++){
			population_mutable[i] = redist[i];
			population_backups[i] = redist[i+populationSize];
		}
	}
}
