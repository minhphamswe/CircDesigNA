package circdesigna;

import circdesigna.abstractpolymer.MonomerDefinition;
import static circdesigna.DomainSequence.*;

/**
 * Connectors consist of two bases A and B, arranged 5' to 3', which are paired.
 */
public class Connector {
	/**
	 * Two bases of two domains, the bases are cyclicly numbered in the complemented form of the domain
	 * 
	 * Negative base numbering is allowed, for instance, -1 means the last base of the domain
	 */
	public Connector(int domainA, int baseA, int domainB, int baseB){
		this.domainA = domainA;
		this.baseA = baseA;
		this.domainB = domainB;
		this.baseB = baseB;
	}
	private int domainA;
	private int baseA;
	private int domainB;
	private int baseB;
	public int base(int i, int[][] domain, MonomerDefinition monomer) {
		if (i < 0 || i >= 2){
			throw new ArrayIndexOutOfBoundsException(i);
		}
		
		int domainNum = i==0?domainA : domainB;
		int baseNum = i==0?baseA : baseB;
		int[] d = domain[domainNum & NA_COMPLEMENT_FLAGINV];
		
		while (baseNum < 0){
			baseNum += d.length;
		}
		baseNum %= d.length;
		
		if ((domainNum & NA_COMPLEMENT_FLAG)!=0){
			return monomer.complement(d[d.length-1-baseNum]);
		} else{
			return monomer.noFlags(d[baseNum]);
		}
	}

	public boolean contains(int i) {
		i &= NA_COMPLEMENT_FLAGINV;
		if ((domainA & NA_COMPLEMENT_FLAGINV) == i){
			return true;
		}
		if ((domainB & NA_COMPLEMENT_FLAGINV) == i){
			return true;
		}
		return false;
	}
}