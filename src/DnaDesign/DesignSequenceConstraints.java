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

import java.util.Arrays;
public class DesignSequenceConstraints {

	public static final int GCL_FLAG =  0 + DNAFLAG_ADD;
	public static final int LOCK_FLAG = GCL_FLAG + DNAFLAG_ADD;

	public static DesignSequenceConstraints getDefaultConstraints() {
		DesignSequenceConstraints defaults = new DesignSequenceConstraints();
		defaults.setMaxConstraint(H, 0);
		defaults.setMaxConstraint(P, 0);
		defaults.setMaxConstraint(Z, 0);
		defaults.setMaxConstraint(D, 0);
		return defaults;
	}
	//-1 means no upper bound
	private int[] maxConstituents;
	private int[] minConstituents;
	//-1 means no lower bound
	public DesignSequenceConstraints(){
		maxConstituents = new int[DNAFLAG_ADD];
		minConstituents = new int[DNAFLAG_ADD];
		Arrays.fill(maxConstituents,-1);
		Arrays.fill(minConstituents,-1);
	}
	/**
	 * -1 to set unconstrained
	 */
	public void setMaxConstraint(int base, int maxVal){
		base = DnaDefinition.noFlags(base);
		maxConstituents[base] = maxVal;
	}
	/**
	 * -1 to set unconstrained
	 */
	public void setMinConstraint(int base, int minVal){
		base = DnaDefinition.noFlags(base);
		minConstituents[base] = minVal;
	}
	private int[] isValid_shared = new int[DNAFLAG_ADD];
	/**
	 * Not multithreaded.
	 */
	public boolean isValid(int[] domain){
		Arrays.fill(isValid_shared,0);
		for(int k = 0; k < domain.length; k++){
			isValid_shared[DnaDefinition.noFlags(domain[k])]++;
		}
		for(int k = 0; k < DNAFLAG_ADD; k++){
			if (isOutOfValidRange(k,isValid_shared[k])){
				return false;
			}
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
	
	private boolean isOutOfValidRange(int base, int basecount) {
		if (minConstituents[base]!=-1){
			if (basecount < minConstituents[base]){
				return true;
			}
		}
		if (maxConstituents[base]!=-1){
			if (basecount > maxConstituents[base]){
				return true;
			}
		}
		return false;
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
		private int lastSwapIndex;
		private boolean isDirectMutation = false;
		/**
		 * Returns false if no more valid mutations are possible.
		 * @return
		 */
		public boolean nextChoice(){
			int oldBase = mut_new[j];
			int oldBase_pure = (oldBase % DNAFLAG_ADD);
			int oldBase_flag = oldBase - oldBase_pure;
			if (oldBase_flag == LOCK_FLAG){
				//Locked bases are not mutable.
				return false;
			}
			while(++i_inBaseOrders < baseOrders.length){
				int testBase = baseOrders[i_inBaseOrders];
				if (testBase==oldBase_pure){
					continue; //This wouldn't be a changing mutation.
				}
				if (canMutateBaseDirectly()){
					return isDirectMutation = true;
				}
				isDirectMutation = false;
				if (canMutateBySwapping()){
					return true;
				}
			}
			return false;
		}
		private boolean canMutateBySwapping() {
			//TODO: check if a mutation is possible by swapping with some other base.
			int oldBase = mut_new[j];
			int oldBase_pure = (oldBase % DNAFLAG_ADD);
			if (oldBase_pure==0){
				//Special case: don't swap "uninitialized" bases around (code 0)
				return false; 
			}
			
			int oldBase_flag = oldBase - oldBase_pure;
			//Currently, this is always possible as long as we have some other base which is of 
			//base type "testBase". This is because swapping doesn't change composition.
			//If we allow more complicated constraints, this assumption will change.
			int testBase = baseOrders[i_inBaseOrders];
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
			if (isOutOfValidRange(oldBase_pure,countBaseOccurrences(oldBase_pure)-1)){
				return false;
				//We can't remove it, it's keeping us above a quota.
			}
			int testBase = baseOrders[i_inBaseOrders];
			//Count the number of occurrences of testBase, which is not the same as oldBase, so we can assume
			//that the index j is irrelevent, and see if that count + 1 is out of range
			int count = countBaseOccurrences(testBase);
			if (isOutOfValidRange(testBase,count+1)){
				return false;
			}
			if (!isAllowableBaseforFlags(oldBase_flag,testBase)){
				return false;
			}
			//Ok! this worked!
			return true;
		}
		/**
		 * Counts the occurrences of base in the domain.
		 */
		private int countBaseOccurrences(int base) {
			base = DnaDefinition.noFlags(base);
			int count = 0;
			for(int i = 0; i < mut_new.length; i++){
				if (DnaDefinition.noFlags(mut_new[i])==base){
					count++;
				}
			}
			return count;
		}

		/**
		 * Resets (or initializes) the iterator on a certain domain, at position j.
		 */
		public void reset(int[] mut_new, int j) {
			this.mut_new = mut_new;
			this.j = j;
			i_inBaseOrders = -1;
		}

		public int getMutationBase() {
			return baseOrders[i_inBaseOrders];
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