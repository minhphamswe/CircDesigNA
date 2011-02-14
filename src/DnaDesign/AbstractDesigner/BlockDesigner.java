package DnaDesign.AbstractDesigner;

import DnaDesign.DomainDesigner;

/**
 * Possible other names: Equal Opportunity Designer, or Patient Designer
 * 
 * Optimizes scores to be LOW. The fittest individual has the smallest score.
 *  An implementation of a genetic algorithm inspired designer.
 * 
 * @author Benjamin
 */
public abstract class BlockDesigner <T extends PopulationDesignMember<T>> {
	public BlockDesigner(SingleMemberDesigner<T> SingleDesigner){
		this.SingleDesigner = SingleDesigner;
	}
	public SingleMemberDesigner<T> SingleDesigner;
	public T[] population_mutable;
	public int populationSize = 0;
	private T fittest;
	private int iterations = 0;
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
				population_mutable[k] = null;
				population_mutable[k] = init.designerCopyConstructor(k);
				population_mutable[k].seedFromOther(init);
			} else {
				population_mutable[k] = init;
			}
		}
	}
	/**
	 * Wraps 
	 */
	public void runBlockIteration(DomainDesigner runner, double endThreshold){
		iterations++;
		System.out.print("Iteration "+iterations);
		
		runBlockIteration_(runner,endThreshold);
		
		System.out.println(" Score "+SingleDesigner.getOverallScore(getBestPerformingChild()));
	}
	/**
	 * Overridden by design implementations.
	 * Runs whatever an "iteration" means.
	 */
	public abstract void runBlockIteration_(DomainDesigner runner, double endThreshold);
	/**
	 * Returns null before blockIteration is called.
	 */
	public T getBestPerformingChild(){
		return fittest;
	}
	public void setBestChild(T child){
		fittest = child;
	}
	public final PopulationDesignMember<T>[] getPopulation(){
		return population_mutable;
	}
	public final int getIterationCount(){
		return iterations;
	}
}
