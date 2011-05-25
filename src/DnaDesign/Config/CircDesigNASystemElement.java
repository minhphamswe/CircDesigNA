package DnaDesign.Config;

import edu.utexas.cssb.circdesigna.DomainSequence;

public class CircDesigNASystemElement extends SystemElement<CircDesigNAConfig>{
	/**
	 * Convenience function ,returns base i of domain sequence ds
	 */
	public final int base(DomainSequence ds, int i, int[][] domain){
		return ds.base(i,domain,Std.monomer);
	}
	public CircDesigNASystemElement(CircDesigNAConfig System) {
		super(System);
	}
}
