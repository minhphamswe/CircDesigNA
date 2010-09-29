package DnaDesign.AbstractDesigner;

import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Possible other names: Equal Opportunity Designer, or Patient Designer
 * 
 * @author Benjamin
 */
public abstract class BlockDesigner <T extends PopulationDesignMember<T>> {
	
	private PopulationDesignMember<T>[] population_mutable;
	private PopulationDesignMember<T>[] population_backup;
	private int populationSize = 0;
	/**
	 * Initializes this designer with one member, and numCopies-1 number of newly created members with that same seed.
	 * Necessary before designing.
	 * 
	 * In actuality, numCopies * 2 copies are made, because the designer will keep a backup copy for you.
	 */
	public void initialize(T init, int numCopies){
		populationSize = numCopies;
		population_backup = new PopulationDesignMember[numCopies];
		population_mutable = new PopulationDesignMember[numCopies];
		for(int k = 0; k < numCopies; k++){
			if (k!=0){
				population_mutable[k] = init.designerCopyConstructor(k);
			} else {
				population_mutable[k] = init;
			}
			population_backup[k] = init.designerCopyConstructor(k);
		}
	}
	public void runBlockIteration(){
		TreeSet<PopulationDesignMember<T>> blockIterationLevel = new TreeSet();
		for(int k = 0; k < populationSize; k++){
			blockIterationLevel.add(population_mutable[k]);
		}
		while(true){
			Iterator<PopulationDesignMember<T>> qb = blockIterationLevel.iterator();
			for(int k = 0; k < 5; k++){ //multiple times for random chance.
				while(qb.hasNext()){
					PopulationDesignMember<T> q = qb.next();
					boolean mutationSuccessful = mutateAndTest(q);
					if (mutationSuccessful){
						//commit
						population_backup[q.myID].seedFromOther(q);
						qb.remove();
					} else {
						//need to backup.
						q.seedFromOther(population_backup[q.myID]);
					}
				}
			}
			if (blockIterationLevel.size()<populationSize*.5){
				break;
			}
		}
		//Seed the fittest
		TreeMap<Double, PopulationDesignMember<T>> populationView = new TreeMap();
		for(PopulationDesignMember<T> q : population_mutable){
			double score = getOverallScore(q);
			populationView.put(-score, q); //sort descending
		}
		PopulationDesignMember<T> top = populationView.remove(populationView.lastKey());
		int bottomToSeed = (int) (populationView.size()*.4f);
		for(int k = 0; k < bottomToSeed; k++){
			PopulationDesignMember<T> bottom = populationView.remove(populationView.firstKey());
			bottom.seedFromOther(top);
		}
	}
	public abstract double getOverallScore(PopulationDesignMember<T> q);
	public abstract boolean mutateAndTest(PopulationDesignMember<T> q);
}
