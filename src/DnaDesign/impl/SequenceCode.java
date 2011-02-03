package DnaDesign.impl;

import DnaDesign.DesignSequenceConstraints;
import DnaDesign.DesignerCode;

public class SequenceCode implements DesignerCode{
	private DesignSequenceConstraints dsc;
	public SequenceCode(){
	}
	public void setConstraints(DesignSequenceConstraints dsc){
		this.dsc = dsc;
	}
	/**
	 * Mutates a single base, maintaining constraint conditions. 
	 */
	public boolean mutateToOther(int[][] domain, int mut_domain, int j) {
		int[] mut_new = domain[mut_domain];
		int oldJ = mut_new[j];
		//if oldJ is 0, change to any.
		
		//Available does NOT count the current value - a mutation must change the base.
		int choices = dsc.countAvailableMutations(mut_new,j);
		if (choices==0){
			return false;
		}
		int choice = (int) (Math.random()*choices);
		dsc.makeAvailableMutationNo(choice,mut_new,j);
		return true;
		
		/*
		if (forceBaseConservingMutations){
			if (len > 1 && isGClessAT(mut_new[j])!=0){
				//Handle the new creation, by swapping another base to conserve the GC / AT balance.
				int otherJ;
				int candidateL3 = -1;
				int candidateL2 = -1;
				int candidateL1 = -1;
				int otherJ_startCheck = int_urn(0,len-1);
				searchOtherJ: for(int otherJi = 0; otherJi < len; otherJi++){
					otherJ = (otherJ_startCheck+otherJi)%len;
					if (otherJ != j){
						if (is[otherJ] < DNAFLAG_ADD){
							//Necessary conditions over.
							if (isGC(is[otherJ])!=isGC(oldJ)){
								if (domain_markings[mut_domain][otherJ]!=DNAMARKER_DONTMUTATE){
									candidateL1 = otherJ; //All conditions (L1) succeeded.
									break searchOtherJ;
								} else {
									candidateL2 = otherJ;
								}
							} else {
								candidateL3 = otherJ;
							}
						}
					}
				}
				otherJ = candidateL1;
				if (otherJ==-1){
					otherJ = candidateL2;
				}
				if (otherJ==-1){
					otherJ = candidateL3;
				}
				if (otherJ==-1){
					//It absolutely was not possible to find a partner. Oh well.
				} else {
					int newJ = is[j]%DNAFLAG_ADD;
					if ((oldJ==G || oldJ==C) && (newJ==A || newJ==T)){
						is[otherJ] = int_urn(0, 1)==0?G:C;
					}
					if ((oldJ==A || oldJ==T) && (newJ==G || newJ==C)){
						is[otherJ] = int_urn(0, 1)==0?A:T;
					}
				}
			}
		}
		*/
	}
	public boolean isValid(int[][] domain, int whichDomain) {
		return dsc.isValid(domain[whichDomain]);
	}
}
