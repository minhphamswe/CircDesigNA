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

import java.lang.reflect.Array;

import circdesigna.CircDesigNA;
import circdesigna.CircDesigNA.ScorePenalty;
import circdesigna.impl.CircDesigNAPMemberImpl;

/**
 * Framework for sequence designers that operate through iterative improvement of a population
 * of candidate solutions. Population may have size 1. 
 * 
 * @author Benjamin
 */
public abstract class BlockDesigner <T extends PopulationDesignMember<T>>{
	public BlockDesigner(SingleMemberDesigner<T> SingleDesigner){
		this.SingleDesigner = SingleDesigner;
	}
	//set to true to add an expensive test that the scorefunctions are indeed isolated by mutation domains
	//Technical note: Coverage of both the "affectedBy" methods and also the "clone" methods of ScorePenalties.
	public static final boolean ASSERT_SCOREFUNC_ISOLATION = false;
	
	public SingleMemberDesigner<T> SingleDesigner;
	
	protected T[] population_mutable;
	protected int populationSize = 0;
	private T fittest;
	private int iterations = 0;
	private double progressInIteration = 0;
	/**
	 * Initializes this designer with one member, and numCopies-1 number of newly created members with that same seed.
	 * Necessary before designing.
	 */
	public void initialize(T init, int popSize){
		populationSize = popSize;
		population_mutable = (T[])Array.newInstance(init.getClass(),popSize);
		for(int k = 0; k < populationSize; k++){
			if (k!=0){
				population_mutable[k] = init.designerCopyConstructor(k);
			} else {
				population_mutable[k] = init;
			}
		}
		progressInIteration = 0;
	}
	public void setProgress(int a, int max){
		progressInIteration = a / (double)max;
	}
	public double getProgress(){
		return progressInIteration;
	}
	/**
	 * Wraps runBlockIteration_  
	 */
	public void runBlockIteration(CircDesigNA runner, double endThreshold){
		iterations++;
		System.out.print("Iteration "+iterations+" ");
		
		runBlockIteration_(runner,endThreshold);
	
		if (ASSERT_SCOREFUNC_ISOLATION){
			for(int i = 0; i < population_mutable.length; i++){
				CircDesigNAPMemberImpl q = (CircDesigNAPMemberImpl)(Object)population_mutable[i];
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
	 */
	public abstract void runBlockIteration_(CircDesigNA runner, double endThreshold);
	/**
	 * Returns null before blockIteration is called.
	 */
	public T getBestPerformingChild(){
		return fittest;
	}
	public void setBestChild(T child){
		fittest = child;
	}
	public final T[] getPopulation(){
		return population_mutable;
	}
	public final int getIterationCount(){
		return iterations;
	}
}
