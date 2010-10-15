package DnaDesign;

import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.D;
import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.H;
import static DnaDesign.DnaDefinition.P;
import static DnaDesign.DnaDefinition.T;
import static DnaDesign.DnaDefinition.Z;
import static DnaDesign.DnaDefinition.noFlags;

import java.util.ArrayList;
import java.util.Arrays;
public class DesignSequenceConstraints {

	public static final int GCL_FLAG =  0 + DNAFLAG_ADD;
	public static final int LOCK_FLAG = GCL_FLAG + DNAFLAG_ADD;

	public static DesignSequenceConstraints getDefaultConstraints() {
		DesignSequenceConstraints defaults = new DesignSequenceConstraints();
		defaults.setMaxConstraint(0, H);
		defaults.setMaxConstraint(0, P);
		defaults.setMaxConstraint(0, Z);
		defaults.setMaxConstraint(0, D);
		return defaults;
	}
	//-1 means no upper bound
	//Format: int[b0...bn-1, constraint]
	//So, int[{1 in A,T position},-1] in minConstituents means no lower bound on A+T
	private static class Constraint{
		private boolean[] regulates;
		private int constraintValue;
		public Constraint() {
			regulates = new boolean[DNAFLAG_ADD];
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
	//-1 means no lower bound
	public DesignSequenceConstraints(){
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
			base[i] = DnaDefinition.noFlags(base[i]);
			made.regulates[base[i]] = true;
		}
		made.constraintValue = maxVal;
		boolean removed = toSet.remove(made);
		toSet.add(made);
		return removed;
	}
	/**
	 * Not multithreaded.
	 */
	public boolean isValid(int[] domain){
		if (getBaseCounts(domain)){
			for(int k = 0; k < DNAFLAG_ADD; k++){
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
	private int[] getBaseCounts_shared = new int[DNAFLAG_ADD];
	private boolean getBaseCounts(int[] domain, int ignoreIndex) {
		Arrays.fill(getBaseCounts_shared,0);
		for(int k = 0; k < domain.length; k++){
			if (k==ignoreIndex) continue;
			getBaseCounts_shared[DnaDefinition.noFlags(domain[k])]++;
		}
		return true;
	}
	/**
	 * Tests whether a given flags are acceptable with a certain base.
	 */
	private boolean isAllowableBaseforFlags(int oldBase_flag, int testBase) {
		//Make sure oldBase_flag is pure flag,
		oldBase_flag -= noFlags(oldBase_flag);
		//Make sure testBase is pure base
		testBase = noFlags(testBase);
		//If we are allowed to make this mutation, return true.
		if (oldBase_flag == GCL_FLAG){
			//Then, only G / C mutations valid.
			return testBase==G || testBase==C;
		}
		return true;
	}
	

	private boolean isOverValidMax(int base){
		return checkInvalidating(maxConstituents,base,true);
	}
	private boolean isUnderValidMin(int base){
		return checkInvalidating(minConstituents,base,false);
	}
	/**
	 * Assumes isValid_shared is filled with the base counts.
	 */
	private boolean checkInvalidating(ArrayList<Constraint> constraints, int base, boolean isMaxConstraint) {
		base = noFlags(base);
		for(Constraint q : constraints){
			if (q.constraintValue!=-1){
				if (q.regulates[base]){
					int sum = 0;
					for(int k = 0; k < DNAFLAG_ADD; k++){
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
		private final int[] baseOrders = new int[]{
				//No specific order.
			A,T,C,G,D,P,H,Z	
		};
		
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
			int oldBase_pure = (oldBase % DNAFLAG_ADD);
			int oldBase_flag = oldBase - oldBase_pure;
			if (oldBase_flag == LOCK_FLAG){
				//Locked bases are not mutable.
				return false;
			}
			while(true){
				if (i_inBaseOrders+1>=baseOrders.length){
					break;
				}
				++i_inBaseOrders;
				b_inBaseOrders = baseOrders[i_inBaseOrders];
				int testBase = b_inBaseOrders;
				if (testBase==oldBase_pure){
					continue; //This wouldn't be a changing mutation.
				}
				if (isDirectMutation = canMutateBaseDirectly()){
					return isDirectMutation;
				}
				isDirectMutation = false;
				if (canMutateBySwapping()){
					return true;
				}
			}
			return false;
		}
		private boolean canMutateBySwapping() {
			int oldBase = mut_new[j];
			int oldBase_pure = (oldBase % DNAFLAG_ADD);
			if (oldBase_pure==DnaDefinition.NOBASE){
				//Special case: don't swap "uninitialized" bases around (code 0)
				return false; 
			}
			
			int oldBase_flag = oldBase - oldBase_pure;
			//Currently, this is always possible as long as we have some other base which is of 
			//base type "testBase". This is because swapping doesn't change composition.
			//If we allow more complicated constraints, this assumption will change.
			int testBase = b_inBaseOrders;
			for(int i = 0; i < mut_new.length; i++){
				if (i==j){
					continue;
				}
				int swapBase = noFlags(mut_new[i]);
				int swapBase_flag = mut_new[i] - swapBase;
				//Does this swap achieve the goal of changing mut_new[j] to testBase?
				if (!(swapBase==testBase)){
					continue;
				}
				//ok. so, is the swap allowed?
				if (!(isAllowableBaseforFlags(oldBase_flag,swapBase) && isAllowableBaseforFlags(swapBase_flag,oldBase_pure))){
					continue;
				}
				//Then we should be able to swap only the BASE parts of the number. (keep constraint flag intact)
				lastSwapIndex = i;
				return true;
			}
			return false;
		}
		private boolean canMutateBaseDirectly(){
			int oldBase = mut_new[j];
			int oldBase_pure = (oldBase % DNAFLAG_ADD);
			int oldBase_flag = oldBase - oldBase_pure;
			//Will replacing oldBase cause us to go under a minimum quota?
			if (getBaseCounts(mut_new,j)){ //Counts all BUT j
				if (isUnderValidMin(oldBase_pure)){
					//We can't remove it, it's keeping us above a quota.
					return false;
				}
			}
			int testBase = b_inBaseOrders;
			//Count the number of occurrences of testBase, which is not the same as oldBase, so we can assume
			//that the index j is irrelevent, and see if the incremented case is out of range
			if (getBaseCounts(mut_new,j)){
				getBaseCounts_shared[testBase]++;
				if (isOverValidMax(testBase)){
					return false;
				}
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
			isUnderQuotaInBase = checkUnderQuota(j);
		}

		private int checkUnderQuota(int j) {
			boolean hadUnderValid = false;
			//Is mut_new[j] underquota?
			if (getBaseCounts(mut_new)){
				if (isUnderValidMin(noFlags(mut_new[j]))){
					return UnderQuotaButImpossibleToChange;
				}
				for(int test_base : baseOrders){
					if (noFlags(mut_new[j])==test_base){
						continue; //These don't count.
					}
					if (isUnderValidMin(test_base)){
						hadUnderValid = true;
						if (isAllowableBaseforFlags(mut_new[j],test_base)){
							return test_base;
						}
					}
				}
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
	public void makeAvailableMutationNo(int choice, int[] mut_new, int j) {
		itr.reset(mut_new,j);
		for(int k = 0; k <= choice; k++){
			itr.nextChoice();
		}
		if (itr.isDirectMutationChoice()){
			int newBase = itr.getMutationBase();
			//new Base < DNA_FLAGADD
			mut_new[j] = mut_new[j] - (mut_new[j]%DNAFLAG_ADD) + newBase;
		} else {
			int swapIndex = itr.getSwapIndex();
			int tmp = mut_new[j];
			mut_new[j] = mut_new[j] - noFlags(mut_new[j]) + noFlags(mut_new[swapIndex]);
			mut_new[swapIndex] = mut_new[swapIndex] - noFlags(mut_new[swapIndex]) + noFlags(tmp);
		}
	}	
}