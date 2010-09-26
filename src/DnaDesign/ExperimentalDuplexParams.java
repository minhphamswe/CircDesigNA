package DnaDesign;

/**
 * <pre>
 * Contains experimentally backed parameters determining the deltaG for 
 * 1) Nearest Neighbor duplexes (WY
 *                               XZ)
 * 2) "Terminal" (Final pair before hairpin loop or internal loop opens) scores.
 * 
 * DNA bases are marked by their identifiers in DnaDefinition.java
 *  </pre>
 */
public interface ExperimentalDuplexParams {
	/**
	 * Returns the delta G for a (W,X) pair neighboring a (Y,Z) pair. 
	 * Nearest Neighbor pairing need not be symmetric; such that 
	 * d((W,X),(Y,Z))!=d((X,W),(Z,Y)).
	 */
	public double getNNdeltaG(int W, int X, int Y, int Z);
	/**
	 * Returns the delta G for a (W,X) pair neighboring a (Y,Z) pair. 
	 * Specifically, Y,Z must be a MISMATCH and X,W must be a PAIR.
	 */
	public double getNNdeltaGterm(int W, int X, int Y, int Z);
	/**
	 * Returns the dangle penalty for base D, on 3primeEnd, where X,Y are the terminal pair of the helix.
	 */
	public double getDanglePenalty(int X, int Y, int D, boolean PrimeEnd);
}
