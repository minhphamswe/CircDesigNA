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
package circdesigna.energy;

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
public interface NAExperimentDatabase {
	/**
	 * Returns the delta G for a (W,X) pair neighboring a (Y,Z) pair. 
	 * Nearest Neighbor pairing need not be symmetric; such that 
	 * d((W,X),(Y,Z))!=d((X,W),(Z,Y)).
	 */
	public double getNNdeltaG(int W, int X, int Y, int Z);
	public int getNNdeltaG_deci(int W, int X, int Y, int Z);
	/**
	 * Returns the delta G for a (W,X) pair neighboring a (Y,Z) pair. 
	 * Specifically, Y,Z must be a MISMATCH and X,W must be a PAIR.
	 */
	public double getNNdeltaGterm(int W, int X, int Y, int Z);
	public int getNNdeltaGterm_deci(int W, int X, int Y, int Z);
	/**
	 * Returns the dangle penalty for base D, on 3primeEnd, where X,Y are the terminal pair of the helix.
	 */
	public double getDanglePenalty(int X, int Y, int D, boolean PrimeEnd);
	public int getDanglePenalty_deci(int X, int Y, int D, boolean PrimeEnd);
	
	/**
	 * Returns the delta G bonus for the association of numStrands strands at T kelvin.
	 */
	public double getDeltaGAssoc(int numStrands, double T);
}
