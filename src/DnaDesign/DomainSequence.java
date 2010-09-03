package DnaDesign;

import java.util.ArrayList;

import DnaDesign.DomainStructureData.DomainStructure;


public class DomainSequence {
	public static final int DNA_COMPLEMENT_FLAG = 0x8000;
	public static final int DNA_SEQ_FLAGSINVERSE = ~(DNA_COMPLEMENT_FLAG);

	/**
	 * & 0x8000 - complement.
	 */
	protected int[] domainList = new int[2];
	protected int numDomains = 0;
	private String moleculeName;
	public String getMoleculeName(){
		return moleculeName;
	}
	private void setFromDSD(DomainStructureData dsd){
		if (dsd==null){
			this.moleculeName = "???";
			return;
		}
		this.moleculeName = dsd.moleculeName;
	}
	public void appendMoleculeNames(DomainSequence seq) {
		ArrayList<String> current = new ArrayList();
		for(String q : moleculeName.split(",")){
			current.add(q.trim());
		}
		for(String q : seq.getMoleculeName().split(",")){
			if (!current.contains(q)){
				current.add(q);
			}
		}
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i < current.size(); i++){
			sb.append(current.get(i));
			if (i+1 < current.size()){
				sb.append(",");
			}
		}
		moleculeName = sb.toString();
	}
	public void setDomains(int a, DomainStructureData dsd) {
		setFromDSD(dsd);
		numDomains = 1;
		domainList[0] = a;
	}
	public void setDomains(int a, int b, DomainStructureData dsd) {
		setFromDSD(dsd);
		numDomains = 2;
		domainList[0] = a;
		domainList[1] = b;
	}
	public void setDomains(ArrayList<Integer> freeList, DomainStructureData dsd) {
		setFromDSD(dsd);
		domainList = new int[freeList.size()];
		numDomains = freeList.size();
		for(int k = 0 ; k < numDomains; k++){
			domainList[k] = freeList.get(k);
		}
	}
	public void setDomains(String subStrand, DomainStructureData dsd) {
		setFromDSD(dsd);
		domainList = DomainDesigner_SharedUtils.utilReadSequence(subStrand,dsd);
		numDomains = domainList.length;
	}
	public void setDomains(DomainStructure ds,DomainStructureData dsd) {
		setFromDSD(dsd);
		domainList = new int[ds.sequencePartsInvolved.length];
		for(int i = 0; i < domainList.length; i++){
			int k = ds.sequencePartsInvolved[i];
			domainList[i] = dsd.domains[k];
		}
		numDomains = domainList.length;
	}
	//GETTERS
	public int length(int[][] domain){
		int length = 0, seq;
		for(int i = 0; i < numDomains; i++){
			//Extract which domain
			seq = domainList[i] & DNA_SEQ_FLAGSINVERSE;
			length += domain[seq].length;
		}
		return length;
	}
	public void mark(int q, int i, int[][] domain, int[][] domain_markings) {
		if(domain_markings==null) return;
		if (i < 0){
			for(int j = q+i+1; j <= q; j++){
				mark(j,domain,domain_markings);
			}
		} else {
			for(int j = q; j < q+i; j++){
				mark(j,domain,domain_markings);
			}
		}
	}
	public void mark(int i, int[][] domain, int[][] domain_markings) {
		mark(i,domain,domain_markings,1);
	}
	public void mark(int i, int[][] domain, int[][] domain_markings, int markerValue) {
		if(domain_markings==null) return;
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < numDomains; r++){
			int dNum = domainList[r] & DNA_SEQ_FLAGSINVERSE;
			d = domain[dNum];
			if (q < d.length){
				if ((domainList[r]&DNA_COMPLEMENT_FLAG)!=0){
					domain_markings[dNum][d.length-1-q] = markerValue;
					return;
				} else{
					domain_markings[dNum][q] = markerValue;
					return;
				}
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	public int base(int i, int[][] domain){
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < numDomains; r++){
			d = domain[domainList[r] & DNA_SEQ_FLAGSINVERSE];
			if (q < d.length){
				if ((domainList[r] & DNA_COMPLEMENT_FLAG)!=0){
					return DnaDefinition.complement(d[d.length-1-q]);
				} else{
					return DnaDefinition.noFlags(d[q]);
				}
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	/**
	 * Returns true if this sequence contains domain i, or its complement.
	 */
	public boolean contains(int i) {
		for(int k = 0; k < numDomains; k++){
			if ((domainList[k] & DNA_SEQ_FLAGSINVERSE) == i){
				return true;
			}
		}
		return false;
	}
	public String toString(DomainStructureData dsd){
		StringBuffer sb = new StringBuffer();
		sb.append("(Partial) "+moleculeName+": ");
		for(int i = 0; i < numDomains; i++){
			sb.append(dsd.getDomainName(domainList[i] & DNA_SEQ_FLAGSINVERSE));
			if ((domainList[i]&DNA_COMPLEMENT_FLAG)!=0){
				sb.append("*");
			}
			if (i+1<numDomains){
				sb.append("|");
			}
		}
		return sb.toString();
	}
	public void makeReverseComplement(DomainSequence out) {
		out.numDomains = numDomains;
		for(int i = 0; i < numDomains; i++){
			out.domainList[numDomains-1-i] = domainList[i];
			out.domainList[numDomains-1-i] ^= DNA_COMPLEMENT_FLAG;
		}
	}
	public boolean equals(Object other){
		if (!(other instanceof DomainSequence)){
			return super.equals(other);
		}
		int len1 = numDomains;
		int len2 = ((DomainSequence)other).numDomains;
		if (len1!=len2){
			return false;
		}
		for(int i = 0; i < len1; i++){
			if (domainList[i]!=((DomainSequence)other).domainList[i]){
				return false;
			}
		}
		return true;
	}
}
