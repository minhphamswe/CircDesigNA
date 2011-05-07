package DnaDesign.AbstractDesigner;

import DnaDesign.DomainDesigner;
import DnaDesign.DomainDesigner.ScorePenalty;
import DnaDesign.impl.DomainDesignPMemberImpl;

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
	//set to true to add a dirty test that the scorefunctions are indeed isolated by mutation domains
	//This checks both the "affectedBy" methods and also the "clone" methods of ScorePenalties.
	public static final boolean ASSERT_SCOREFUNC_ISOLATION = false;
	
	public SingleMemberDesigner<T> SingleDesigner;
	public T[] population_mutable;
	public int populationSize = 0;
	private T fittest;
	private int iterations = 0;
	/**
	 * Initializes this designer with one member, and numCopies-1 number of newly created members with that same seed.
	 * Necessary before designing.
	 */
	public void initialize(T init, int numCopies){
		populationSize = numCopies;
		population_mutable = (T[]) new PopulationDesignMember[numCopies];
		for(int k = 0; k < numCopies; k++){
			if (k!=0){
				population_mutable[k] = null;
				population_mutable[k] = init.designerCopyConstructor(k);
				//Which calls population_mutable[k].seedFromOther(init);
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
		System.out.print("Iteration "+iterations+" ");
		
		runBlockIteration_(runner,endThreshold);
	
		if (ASSERT_SCOREFUNC_ISOLATION){
			for(int i = 0; i < population_mutable.length; i++){
				DomainDesignPMemberImpl q = (DomainDesignPMemberImpl)(Object)population_mutable[i];
				//Check: before mutation, all penalties should have change in score of 0.
				for(ScorePenalty s : q.penalties){
					if (s.evalScore(q.domain,q.domain_markings)!=0){
						throw new RuntimeException("FAIL");
					}
					s.dedicate();
				}
			}
		}

		
		final T bestPerformingChild = getBestPerformingChild();
		if (bestPerformingChild==null){
			return;
		}
		System.out.println("Score "+SingleDesigner.getOverallScore(bestPerformingChild));
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
