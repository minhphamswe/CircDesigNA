package DnaDesign.AbstractDesigner;

import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import DnaDesign.DomainDesigner;

/**
 * Possible other names: Equal Opportunity Designer, or Patient Designer
 * 
 * Optimizes scores to be LOW. The fittest individual has the smallest score.
 *  
 * 
 * @author Benjamin
 */
public abstract class BlockDesigner <T extends PopulationDesignMember<T>> {
	private T[] population_mutable;
	private T fittest;
	private int populationSize = 0;
	private int iterations = 0;
	private int param_iterationShortcut = 1000;
	/**
	 * Initializes this designer with one member, and numCopies-1 number of newly created members with that same seed.
	 * Necessary before designing.
	 * 
	 * In actuality, numCopies * 2 copies are made, because the designer will keep a backup copy for you.
	 */
	public void initialize(T init, int numCopies){
		populationSize = numCopies;
		population_mutable = (T[]) new PopulationDesignMember[numCopies];
		for(int k = 0; k < numCopies; k++){
			if (k!=0){
				population_mutable[k] = init.designerCopyConstructor(k);
				population_mutable[k].seedFromOther(init);
			} else {
				population_mutable[k] = init;
			}
		}
	}
	/**
	 * Break if any child is at least as optimal as endThreshold.
	 */
	public void runBlockIteration(DomainDesigner runner, double endThreshold){
		TreeSet<T> blockIterationLevel = new TreeSet();
		for(int k = 0; k < populationSize; k++){
			blockIterationLevel.add(population_mutable[k]);
		}
		iterations++;
		blockItr: for(int itrCount = 0;; itrCount++){
			Iterator<T> qb = blockIterationLevel.iterator();
			//for(int k = 0; k < 1; k++){ //multiple times for random chance.
				while(qb.hasNext()){
					T q = qb.next();
					boolean mutationSuccessful = mutateAndTestAndBackup(q);
					if (mutationSuccessful){
						//commit. Did we solve the problem entirely?
						double newScore = getOverallScore(q);
						if (newScore <= endThreshold){
							break blockItr;
						}
						qb.remove();
					} else {
						//need to backup. It should have
					}
				}
			//}
			if (blockIterationLevel.size()<populationSize*.5){
				break;
			}
			if(runner.abort){
				break;
			}
			if (param_iterationShortcut>=0){
				if (itrCount>param_iterationShortcut && populationSize-blockIterationLevel.size()>0){
					//Someone made it through. Break;
					System.err.println("Breaking iteration early.");
					break;
				}
			}
		}
		//Seed the fittest
		TreeMap<Double, T> populationView = new TreeMap();
		for(T q : population_mutable){
			double score = getOverallScore(q);
			populationView.put(-score, q); //sort descending
		}
		fittest = populationView.remove(populationView.lastKey());
		System.out.println("Iteration "+iterations+" Score "+getOverallScore(fittest));
		int bottomToSeed = (int) (populationView.size()*.4f);
		for(int k = 0; k < bottomToSeed; k++){
			T bottom = populationView.remove(populationView.firstKey());
			bottom.seedFromOther(fittest);
		}
	}
	/**
	 * Returns null before blockIteration is called.
	 */
	public final T getBestPerformingChild(){
		return fittest;
	}
	public final PopulationDesignMember<T>[] getPopulation(){
		return population_mutable;
	}
	public abstract double getOverallScore(T q);
	public abstract boolean mutateAndTestAndBackup(T q);
}
