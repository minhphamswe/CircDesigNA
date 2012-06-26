package circdesigna;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import circdesigna.DomainStructureBNFTree.DomainStructure;
import circdesigna.abstractpolymer.MonomerDefinition;

/**
 * An abstraction of a interactive region, consisting of continuous unpaired domains with energetically relevant "connector" base pairs
 */
public class GeneralizedInteractiveRegion {
	public static final int NA_COMPLEMENT_FLAG = 0x8000;
	public static final int NA_COMPLEMENT_FLAGINV = ~(NA_COMPLEMENT_FLAG);
	public static final int DNAMARKER_DONTMUTATE = -1;

	/**
	 * | DNA_COMPLEMENT_FLAG - complement.
	 */
	//domainList.length >= 2.
	public int[] domainList;
	/**
	 * connectorsToPrior[0] contains the connectors, in 5' to 3' order, which are on the 5' side of domain 0
	 * connectorsToPrior[domainList.length] contains the connectors, in 5' to 3' order, which are on the 3' side of domain domainList.length-1
	 * 
	 * If the molecule is circular, then connectorsToPrior[0] is equal to connectorsToPrior[domainList.length].
	 */
	private ArrayList<Connector>[] connectorsToPrior; 
	private String moleculeName;
	private boolean circular = false;

	/**
	 * Returns connectors which are connected to the 5' end of base i
	 */
	public List<Connector> getConnectorsTo5Of(int i, int[][] domain){
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			d = domain[domainList[r] & NA_COMPLEMENT_FLAGINV];
			if (q == 0){
				return connectorsToPrior[r];
			}
			if (q < d.length){
				return null;
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	/**
	 * Returns connectors which are connected to the 3' end of base i
	 */
	public List<Connector> getConnectorsTo3Of(int i, int[][] domain){
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			d = domain[domainList[r] & NA_COMPLEMENT_FLAGINV];
			if (q == d.length - 1){
				return connectorsToPrior[r+1];
			}
			if (q < d.length){
				return null;
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	
	public String getMoleculeName(){
		return moleculeName;
	}
	public void setMoleculeName(String newName){
		moleculeName = newName;
	}
	/**
	 * Marks the "circular" flag on this domain sequence.
	 */
	public void setCircular(boolean circular){
		this.circular = circular;
	}
	public boolean isCircular(){
		return circular;
	}
	public void appendMoleculeNames(GeneralizedInteractiveRegion seq) {
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
	public void setDomains(int a, AbstractComplex dsd) {
		domainList = new int[]{a};
		init_state(dsd);
	}
	public void setDomains(int a, int b, AbstractComplex dsd) {
		domainList = new int[]{a, b};
		init_state(dsd);
	}	
	public void setDomains(List<Integer> freeList, AbstractComplex dsd) {
		domainList = new int[freeList.size()];
		for(int k = 0 ; k < domainList.length; k++){
			domainList[k] = freeList.get(k);
		}
		init_state(dsd);
	}
	public void setDomains(String subStrand, DomainDefinitions dd, AbstractComplex dsd) {
		domainList = CircDesigNA_SharedUtils.utilReadSequence(subStrand,dd);
		init_state(dsd);
	}
	public void setDomains(DomainStructure ds, DomainStructureBNFTree dsd) {
		domainList = new int[ds.sequencePartsInvolved.length];
		for(int i = 0; i < domainList.length; i++){
			int k = ds.sequencePartsInvolved[i];
			domainList[i] = dsd.domains[k];
		}
		init_state(dsd);
	}
	public void setDomainsAndConnectors(List<Object> sequence, AbstractComplex dsd){
		//Integers are domain numbers, Connectors are connector objects.
		ArrayList<Integer> domainListNeu = new ArrayList();
		ArrayList<ArrayList<Connector>> connectorsNeu = new ArrayList();
		connectorsNeu.add(new ArrayList()); //Initial list (3' of most 3' domain)
		
		for(Object q : sequence){
			if (q instanceof Integer){
				domainListNeu.add((Integer)q);
				connectorsNeu.add(new ArrayList()); //New connector list
			} else {
				connectorsNeu.get(connectorsNeu.size() - 1).add((Connector)q);
			}
		}
		
		domainList = new int[domainListNeu.size()];
		connectorsToPrior = (ArrayList<Connector>[]) new ArrayList[domainList.length + 1];
		for(int i = 0; i < domainList.length; i++){
			domainList[i] = domainListNeu.get(i);
			connectorsToPrior[i] = connectorsNeu.get(i);
		}
		connectorsToPrior[domainList.length] = connectorsNeu.get(domainList.length);
		init_state(dsd);
	}
	/**
	 * Build lookup tables and register molecule name
	 */
	private void init_state(AbstractComplex dsd){
		if (dsd==null){
			this.moleculeName = "???";
		} else {
			this.moleculeName = dsd.getMoleculeName();
		}
		if (connectorsToPrior == null){
			connectorsToPrior = new ArrayList[domainList.length+1];
		}
	}
	//GETTERS
	public int length(int[][] domain){
		int length = 0, seq;
		for(int i = 0; i < domainList.length; i++){
			//Extract which domain
			seq = domainList[i] & NA_COMPLEMENT_FLAGINV;
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
	private void mark(int i, int[][] domain, int[][] domain_markings, int markerValue) {
		if(domain_markings==null) return;
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			int dNum = domainList[r] & NA_COMPLEMENT_FLAGINV;
			d = domain[dNum];
			if (q < d.length){
				if ((domainList[r]&NA_COMPLEMENT_FLAG)!=0){
					int old = domain_markings[dNum][d.length-1-q];
					domain_markings[dNum][d.length-1-q] = old==DNAMARKER_DONTMUTATE?markerValue:old+markerValue;
					return;
				} else {
					int old = domain_markings[dNum][q];
					domain_markings[dNum][q] = old==DNAMARKER_DONTMUTATE?markerValue:old+markerValue;
					return;
				}
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	public int base(int i, int[][] domain, MonomerDefinition monomer){
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			d = domain[domainList[r] & NA_COMPLEMENT_FLAGINV];
			if (q < d.length){
				if ((domainList[r] & NA_COMPLEMENT_FLAG)!=0){
					return monomer.complement(d[d.length-1-q]);
				} else{
					return monomer.noFlags(d[q]);
				}
			}
			q -= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	
	public int domainAt(int i, int[][] domain) {
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			d = domain[domainList[r] & NA_COMPLEMENT_FLAGINV];
			if (q < d.length){
				return domainList[r];
			}
			q-= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	/**
	 * @param i
	 * @param domain
	 * @param offsetIntoUncomplemented - if true, the offset is into the uncomplemented version of the base.
	 * If false, the offset is calculated into a domain of the form in which it appears in this domainsequence.
	 */
	public int offsetInto(int i, int[][] domain, boolean offsetIntoUncomplemented) {
		int q = i, r = 0;
		int[] d;
		for(r = 0; r < domainList.length; r++){
			d = domain[domainList[r] & NA_COMPLEMENT_FLAGINV];
			if (q < d.length){
				if (offsetIntoUncomplemented){
					if ((domainList[r] & NA_COMPLEMENT_FLAG) != 0){
						return d.length-1-q;
					} else {
						return q;
					}
				} else {
					return q;
				}
			}
			q-= d.length;
		}
		throw new ArrayIndexOutOfBoundsException(i);
	}
	/**
	 * Returns true if this sequence contains domain i, or its complement.
	 */
	public boolean contains(int i) {
		i &= NA_COMPLEMENT_FLAGINV;
		for(int k = 0; k < domainList.length; k++){
			if ((domainList[k] & NA_COMPLEMENT_FLAGINV) == i){
				return true;
			}
		}
		for(int k = 0; k < connectorsToPrior.length; k++){
			if (connectorsToPrior[k] != null){
				for(Connector c : connectorsToPrior[k]){
					if (c.contains(i)){
						return true;
					}
				}
			}
		}
		return false;
	}
	public void makeReverseComplement(DomainSequence out) {
		out.domainList = new int[domainList.length];
		for(int i = 0; i < domainList.length; i++){
			out.domainList[domainList.length-1-i] = domainList[i];
			out.domainList[domainList.length-1-i] ^= NA_COMPLEMENT_FLAG;
		}
		out.setMoleculeName(getMoleculeName());
	}
	public boolean isSubsequenceOf(DomainSequence q) {
		int len1 = domainList.length;
		int len2 = q.domainList.length;
		for(int i = 0; i < len2; i++){
			boolean isSubsequence = true;
			for(int j = 0; j < len1; j++){
				if (i+j>=len2){
					isSubsequence = false;
					break;
				}
				if (q.domainList[i+j]!=domainList[j]){
					isSubsequence = false;
					break;
				}
			}
			if (isSubsequence){
				return true;
			}
		}
		return false;
	}
	public String toString(DomainDefinitions dsd){
		StringBuffer sb = new StringBuffer();
		sb.append("(Unpaired domains) "+getMoleculeName()+": ");
		for(int i = 0; i < domainList.length; i++){
			sb.append(dsd.getDomainName(domainList[i] & NA_COMPLEMENT_FLAGINV));
			if ((domainList[i]&NA_COMPLEMENT_FLAG)!=0){
				sb.append("*");
			}
			if (i+1<domainList.length){
				sb.append("|");
			}
		}
		if (isCircular()){
			sb.append(" (circular)");
		}
		return sb.toString();
	}

	public boolean equals(Object other){
		if (!(other instanceof GeneralizedInteractiveRegion)){
			return super.equals(other);
		}
		return Arrays.equals(connectorsToPrior, ((GeneralizedInteractiveRegion)other).connectorsToPrior) && 
				Arrays.equals(domainList, ((GeneralizedInteractiveRegion)other).domainList) &&
				circular == ((GeneralizedInteractiveRegion)other).circular; 
	}
}
