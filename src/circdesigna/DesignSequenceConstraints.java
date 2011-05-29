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
package circdesigna;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math.optimization.GoalType;
import org.apache.commons.math.optimization.RealPointValuePair;
import org.apache.commons.math.optimization.linear.LinearConstraint;
import org.apache.commons.math.optimization.linear.LinearObjectiveFunction;
import org.apache.commons.math.optimization.linear.Relationship;
import org.apache.commons.math.optimization.linear.SimplexSolver;

import circdesigna.abstractpolymer.MonomerDefinition;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;


/**
 * Handles mutations on constrained sequence space. 
 */
public class DesignSequenceConstraints extends CircDesigNASystemElement{
	//-1 means no upper bound
	//Format: int[b0...bn-1, constraint]
	//So, int[{1 in A,T position},-1] in minConstituents means no lower bound on A+T
	private class Constraint{
		private boolean[] regulates;
		private int constraintValue;
		public Constraint() {
			regulates = new boolean[Std.monomer.getNumMonomers()];
			Arrays.fill(regulates,false);
			constraintValue = -1;
		}
		/**
		 * Two constraints "overlap" if the regulated set of one is a subset of the other.
		 */
		public boolean equals(Object otherO){
			if (!(otherO instanceof Constraint)){
				return false;
			}
			Constraint other = (Constraint) otherO;
			return isSubset(regulates,other.regulates)&&isSubset(other.regulates,regulates);
		}
		private boolean isSubset(boolean[] regulates2, boolean[] regulates3) {
			for(int k = 0; k < regulates2.length; k++){
				if (regulates2[k]&&!regulates3[k]){
					return false;
				}
			}
			return true;
		}
	}
	private ArrayList<Constraint> maxConstituents;
	private ArrayList<Constraint> minConstituents;
	private double[] simplexSolution;
	//-1 means no lower bound
	public DesignSequenceConstraints(CircDesigNAConfig system){
		super(system);
		maxConstituents = new ArrayList<Constraint>();
		minConstituents = new ArrayList<Constraint>();
	}
	/**
	 * -1 to set unconstrained
	 * Returns true if a conflicting constraint existed, (and was overwritten)
	 */
	public boolean setMaxConstraint(int maxVal, int ... base){
		return addConstraint(maxConstituents,maxVal,base);
	}
	/**
	 * -1 to set unconstrained
	 * Returns true if a conflicting constraint existed, (and was overwritten)
	 */
	public boolean setMinConstraint(int minVal, int ... base){
		return addConstraint(minConstituents,minVal,base);
	}
	/**
	 * Returns true if an old entry was removed to add this constraint.
	 */
	private boolean addConstraint(ArrayList<Constraint> toSet, int maxVal, int[] base) {
		if (base.length==0){
			throw new RuntimeException("Invalid constraint, no bases");
		}
		Constraint made = new Constraint();
		for(int i = 0; i < base.length; i++){
			base[i] = Std.monomer.noFlags(base[i]);
			made.regulates[base[i]] = true;
		}
		made.constraintValue = maxVal;
		boolean removed = toSet.remove(made);
		toSet.add(made);
		
		//Check consistency. If constraints have no solution, an exception is thrown. 
		solveSimplex();
		
		return removed;
	}
	private void solveSimplex(){
		//Closest-To-Origin objective
		double[] ones = new double[Std.monomer.getNumMonomers()];
		for(int i = 0; i < ones.length; i++){
			ones[i] = 1;
		}
		LinearObjectiveFunction f = new LinearObjectiveFunction(ones, 0);
		
		List<LinearConstraint> constraints = new ArrayList();
		for(Constraint d : maxConstituents){
			if (d.constraintValue==-1){
				continue;
			}
			double[] ei = new double[Std.monomer.getNumMonomers()];
			for(int i = 0; i < ei.length; i++){
				if (d.regulates[i]){
					ei[i] = 1;
				}
			}
			constraints.add(new LinearConstraint(ei, Relationship.LEQ, d.constraintValue));
		}
		for(Constraint d : minConstituents){
			if (d.constraintValue==-1){
				continue;
			}
			double[] ei = new double[Std.monomer.getNumMonomers()];
			for(int i = 0; i < ei.length; i++){
				if (d.regulates[i]){
					ei[i] = 1;
				}
			}
			constraints.add(new LinearConstraint(ei, Relationship.GEQ, d.constraintValue));
		}
		try {
			RealPointValuePair optimize = new SimplexSolver().optimize(f, constraints, GoalType.MINIMIZE, true);
			simplexSolution = optimize.getPoint();
			//System.out.println(Arrays.toString(simplexSolution));
		} catch (Throwable e) {
			throw new RuntimeException("Constraints are too strict: "+e.getMessage());
		}
		
	}
	/**
	 * Not multithreaded.
	 */
	public boolean isValid(int[] domain){
		if (getBaseCounts(domain)){
			for(int k = 0; k < Std.monomer.getNumMonomers(); k++){
				if (isOutOfValidRange(k)){
					return false;
				}
			}
		}
		return true;
	}
	/**
	 * Modifies getBaseCounts_shared.
	 */
	private boolean getBaseCounts(int[] domain) {
		return getBaseCounts(domain,-1);
	}
	private int[] getBaseCounts_shared;
	private boolean getBaseCounts(int[] domain, int ignoreIndex) {
		//Ensure correct length of share vector
		if (getBaseCounts_shared==null || getBaseCounts_shared.length != Std.monomer.getNumMonomers()){
			getBaseCounts_shared = new int[Std.monomer.getNumMonomers()];
		}
		Arrays.fill(getBaseCounts_shared,0);
		for(int k = 0; k < domain.length; k++){
			if (k==ignoreIndex) continue;
			getBaseCounts_shared[Std.monomer.noFlags(domain[k])]++;
		}
		return true;
	}
	/**
	 * Tests whether given flags allow a certain base.
	 */
	private boolean isAllowableBaseforFlags(int constraintBase, int testBase) {
		if (Std.isNAmode()){
			//Make sure testBase is pure base
			testBase = Std.monomer.noFlags(testBase);
			//If we are allowed to make this mutation, return true.
			return Std.monomer.allowBase(constraintBase,testBase);
		}
		return true;
	}
	

	private boolean isOverValidMax(int base){
		return isOverValidMax(base, 0);
	}
	private boolean isOverValidMax(int base, int addToComp){
		return checkInvalidating(maxConstituents,base,true,addToComp);
	}
	private boolean isUnderValidMin(int base){
		return isUnderValidMin(base, 0);
	}
	private boolean isUnderValidMin(int base, int addToComp){
		return checkInvalidating(minConstituents,base,false,addToComp);
	}

	private boolean checkInvalidating(ArrayList<Constraint> constraints, int base, boolean isMaxConstraint, int addToComp) {
		base = Std.monomer.noFlags(base);
		for(Constraint q : constraints){
			if (q.constraintValue!=-1){
				if (q.regulates[base]){
					int sum = addToComp;
					for(int k = 0; k < Std.monomer.getNumMonomers(); k++){
						if (q.regulates[k]){
							sum += getBaseCounts_shared[k];
						}
					}
					if (isMaxConstraint){
						if (sum > q.constraintValue){
							return true;
						}
					} else {
						if (sum < q.constraintValue){
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	private boolean isOutOfValidRange(int base) {
		return isUnderValidMin(base) || isOverValidMax(base);
	}
	/**
	 * Enumerates the possible mutations that could be made, without invalidating the sequence constraints.
	 * @author Benjamin
	 */
	private class ConstraintChoiceIterator{
		/*
		private final int[] baseOrders = new int[]{
				//No specific order.
			A,T,C,G,D,P,H,Z	
		};
		*/
		
		private int[] mut_new;
		private int j;
		//Index in the baseOrders array.
		private int i_inBaseOrders;
		//Base baseOrders[i_inBaseOrders], or overriden
		private int b_inBaseOrders;
		private int lastSwapIndex;
		private boolean isDirectMutation = false;
		private final int UnderQuota_Not = -2;
		private final int UnderQuotaButImpossibleToChange = -3;
		private int isUnderQuotaInBase = UnderQuota_Not;
		/**
		 * Returns false if no more valid mutations are possible.
		 * @return
		 */
		public boolean nextChoice(){
			if (isUnderQuotaInBase!=UnderQuota_Not){
				if (isUnderQuotaInBase==UnderQuotaButImpossibleToChange){
					return false; 
				}
				boolean alreadyReturnedBypass = isUnderQuotaInBase==b_inBaseOrders;
				b_inBaseOrders = isUnderQuotaInBase;
				isDirectMutation = true;
				return !alreadyReturnedBypass;
			}
			int oldBase = mut_new[j];
			int oldBase_pure = Std.monomer.noFlags(oldBase);
			int oldBase_flag = oldBase - oldBase_pure;
			if (oldBase_flag == Std.monomer.LOCK_FLAG()){
				//Locked bases are not mutable.
				return false;
			}
			while(true){
				//State machine tick.
				if (i_inBaseOrders+1>=Std.monomer.getMonomers().length){
					break;
				}
				++i_inBaseOrders;
				b_inBaseOrders = Std.monomer.getMonomers()[i_inBaseOrders];
				int testBase = b_inBaseOrders;
				if (testBase==oldBase_pure){
					continue; //This wouldn't be a changing mutation.
				}
				isDirectMutation = false;
				if (isDirectMutation = canMutateBaseDirectly()){
					return isDirectMutation;
				}
				if (canMutateBySwapping()){
					return true;
				}
			}
			return false;
		}
		private boolean canMutateBySwapping() {
			int oldBase = mut_new[j];
			int oldBase_pure = Std.monomer.noFlags(oldBase);
			if (oldBase_pure==MonomerDefinition.NOBASE){
				//Special case: don't swap "uninitialized" bases around (code 0)
				return false; 
			}
			
			int oldBase_flag = oldBase - oldBase_pure;
			//Currently, this is always possible as long as we have some other base which is of 
			//base type "testBase". This is because swapping doesn't change composition.
			//If we allow more complicated constraints, this assumption will change.
			int testBase = b_inBaseOrders;
			int i_off = (int) (Math.random()*mut_new.length);
			for(int i_true = 0; i_true < mut_new.length; i_true++){
				int i = (i_off + i_true) % mut_new.length;
				
				//Can we swap i with j, to achieve testBase at j?
				if (i==j){ //No, same base.
					continue;
				}
				int swapBase = Std.monomer.noFlags(mut_new[i]);
				int swapBase_flag = mut_new[i] - swapBase;
				//Does this swap achieve the goal of changing mut_new[j] to testBase?
				if (!(swapBase==testBase)){
					continue;
				}
				//so, is the swap allowed?
				if (!(isAllowableBaseforFlags(oldBase_flag,swapBase) && isAllowableBaseforFlags(swapBase_flag,oldBase_pure))){
					continue;
				}
				//Then we should be able to swap only the BASE parts of the number. (keep constraint flag intact)
				lastSwapIndex = i;
				return true;
			}
			return false;
		}
		/**
		 * This method assumes that the values in getBaseCounts_shared do NOT include the addition
		 * of base j. (In other words, remove its contribution before calling this method)
		 */
		private boolean canMutateBaseDirectly(){
			int oldBase = mut_new[j];
			int oldBase_pure = Std.monomer.noFlags(oldBase);
			int oldBase_flag = oldBase - oldBase_pure;
			//Will replacing oldBase cause us to go under a minimum quota?
			if (isUnderValidMin(oldBase_pure)){
				//We can't remove it, it's keeping us above a quota.
				return false;
			}
			int testBase = b_inBaseOrders;
			//Count the number of occurrences of testBase, which is not the same as oldBase, so we can assume
			//that the index j is irrelevent, and see if the incremented case is out of range
			int old_test_ct = getBaseCounts_shared[testBase]; 
			getBaseCounts_shared[testBase]++;
			boolean isOverValidMax = isOverValidMax(testBase);
			getBaseCounts_shared[testBase] = old_test_ct;
			if (isOverValidMax){
				return false;
			}
			if (!isAllowableBaseforFlags(oldBase_flag,testBase)){
				return false;
			}
			//Ok! this worked!
			return true;
		}
		/**
		 * Resets (or initializes) the iterator on a certain domain, at position j.
		 */
		public void reset(int[] mut_new, int j) {
			this.mut_new = mut_new;
			this.j = j;
			i_inBaseOrders = -1;
			isDirectMutation = true;
			b_inBaseOrders = -1;
			//Get base counts.
			getBaseCounts(mut_new,j);
			
			//This call provides the ability to move "into" valid sequence space from outside.
			//If isUnderQuotaInBase is equal to "UnderQuotaButImpossibleToChange", then
			//this mutation must fail, because none of the bases which are under quota can be applied
			//here.
			//Otherwise, if isUnderQuota is positive, then set this mutation to be base isUnderQuota.
			//The target point is calculated by evaluating the simplex solution of the system of linear equations.
			//It may not be attainable, because the base-specific sequence constraints may further restrict
			//the solution space.
			isUnderQuotaInBase = checkUnderQuota(j);
		}

		private int checkUnderQuota(int j) {
			int[] underQuota = new int[Std.monomer.getNumMonomers()];
			int underQuota_Used = 0;
			boolean hadUnderValid = false;
			for(int test_base : Std.monomer.getMonomers()){
				if (isUnderValidMin(test_base) && !isOverValidMax(test_base,1)){
					if (isAllowableBaseforFlags(mut_new[j],test_base)){
						//Ok, we MUST add more of this kind of base, to make our quota.
						hadUnderValid = true;
						//The following inequality returns false only if there are multiple
						//ways to make a quota (say, a G+C limit). It chooses the base which
						//is still under the simple solution, which is a valid tradeoff.
						//Leaving out this additional check leads to the possibility of us putting
						//too many of a single base in a multiple-base quota, and then not having
						//any way of fixing thir problem. It's similar to resource allocation.
						if (getBaseCounts_shared[test_base] < simplexSolution[test_base]){
							underQuota[underQuota_Used++] = test_base;	
						}
					}
				}
			}
			if (underQuota_Used>0){
				return underQuota[(int) (Math.random()*underQuota_Used)]; 
			}
			return hadUnderValid?UnderQuotaButImpossibleToChange:UnderQuota_Not;
		}
		public int getMutationBase() {
			if (b_inBaseOrders <0){
				throw new RuntimeException("Assertion error: possible mutation < 0");
			}
			return b_inBaseOrders;
		}
		public boolean isDirectMutationChoice() {
			return isDirectMutation;
		}
		public int getSwapIndex() {
			return lastSwapIndex;
		}
	}
	ConstraintChoiceIterator itr = new ConstraintChoiceIterator();
	
	public int countAvailableMutations(int[] mut_new, int j) {
		itr.reset(mut_new,j);
		int count = 0;
		while(itr.nextChoice()){
			count++;
		}
		return count;
	}
	public int getMutationNumberForNewBase(int[] mut_new, int j, int new_base){
		new_base = Std.monomer.noFlags(new_base);
		itr.reset(mut_new,j);
		int count = 0;
		while(itr.nextChoice()){
			if (itr.getMutationBase()==new_base){
				return count;
			}
			count++;
		}
		return -1;
		
	}
	public void makeAvailableMutation(int choice, int[] mut_new, int j) {
		itr.reset(mut_new,j);
		for(int k = 0; k <= choice; k++){
			itr.nextChoice();
		}
		if (itr.isDirectMutationChoice()){
			int newBase = itr.getMutationBase();
			//new Base < getNumMonomers
			mut_new[j] = mut_new[j] - Std.monomer.noFlags(mut_new[j]) + newBase;
		} else {
			int swapIndex = itr.getSwapIndex();
			int tmp = mut_new[j];
			mut_new[j] = mut_new[j] - Std.monomer.noFlags(mut_new[j]) + Std.monomer.noFlags(mut_new[swapIndex]);
			mut_new[swapIndex] = mut_new[swapIndex] - Std.monomer.noFlags(mut_new[swapIndex]) + Std.monomer.noFlags(tmp);
		}
	}
}