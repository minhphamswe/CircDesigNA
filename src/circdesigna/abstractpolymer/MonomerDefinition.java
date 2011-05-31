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
package circdesigna.abstractpolymer;

import static circdesigna.abstractpolymer.DnaDefinition.*;
/**
 * MonomerDefinitions handle the notational aspects of representing linear polymers as integer arrays.
 * They can convert from this format back and forth to text Strings. Constraint information is also combined
 * together, because separating the design flags from sequence information was too cumbersome.
 * 
 * REQUIREMENTS:
 * 
 * Monomers must be encoded as integers in the space 0 (reserved for no monomer) ... up to getNumMonomers()-1.
 * There may be no "holes" in this encoding
 * 
 * (Classes that mutate strings of monomers will assume that incrementing a value and wrapping at getNumMonomers
 * is a valid transformation)
 * 
 * Flags may be added by adding in multiples of getNumMonomers().
 */
public abstract class MonomerDefinition {
	public static final int NO_MONOMER = 0, NOBASE=NO_MONOMER, NOACID=NO_MONOMER;

	public final int noFlags(int oldBase){
		return oldBase % getNumMonomers();
	}
	/**
	 * Equivalent to the 'S' degenerate basepair
	 */
	public int GCL_FLAG(){
		return getNumMonomers()*S;
	}
	public int LOCK_FLAG(){
		return getNumMonomers();
	}
	
	/**
	 * IUPAC's definitions:
	 * 
	 * Flag	Char	Allows
	 * 0	N		A,C,G,U
	 * 1	(lock) 
	 * 2	R		A,G
	 *  	Y		C,U
	 *  	M		A,C
	 *  	K		G,U
	 *  	S		C,G
	 *  	W		A,U
	 *  	V		A,C,G
	 *  	H		A,C,U
	 *  	B		C,G,U
	 *  	D		A,G,U
	 * 
	 */
	private static final int N = 0, R = 2, Y = R+1, M = Y+1, K = M+1, S = K+1, W = S+1,
	V = W+1, H = V+1, B = H+1, D = B+1;
	private final int[][] flagsAllowsMap = new int[12][getNumMonomers()];
	{
		registerAllowRule(N, A,T,C,G);
		registerAllowRule(R, A, G);
		registerAllowRule(Y, C, T);
		registerAllowRule(M, A, C);
		registerAllowRule(K, G, T);
		registerAllowRule(S, G, C);
		registerAllowRule(W, A, T);
		registerAllowRule(V, A, C, G);
		registerAllowRule(H, A, C, T);
		registerAllowRule(B, C, G, T);
		registerAllowRule(D, A, G, T);
	}
	private void registerAllowRule(int index, int ... bases){
		for(int q : bases){
			flagsAllowsMap[index][getNormalBaseFromZero(q)] = 1;
		}
	}
	public int decodeConstraintChar(char charAt){
		int flagMult = getNumMonomers();
		int degenerate = -1;
		switch(Character.toUpperCase(charAt)){
		case 'N':
			degenerate = N*flagMult; break;
		case 'R':
			degenerate = R*flagMult; break;
		case 'Y':
			degenerate = Y*flagMult; break;
		case 'M':
			degenerate = M*flagMult; break;
		case 'K':
			degenerate = K*flagMult; break;
		case 'S':
			degenerate = S*flagMult; break;
		case 'W':
			degenerate = W*flagMult; break;
		case 'V':
			degenerate = V*flagMult; break;
		case 'H':
			degenerate = H*flagMult; break;
		case 'B':
			degenerate = B*flagMult; break;
		case 'D':
			degenerate = D*flagMult; break;
		}
		if (degenerate!=-1){
			//Ok, initial degenerate. Not locked.
			return degenerate;
		}
		//Locked flag.
		int ret = decodeBaseChar(charAt);
		ret += LOCK_FLAG();
		return ret;
	}
	public int decodeInitializationChar(char u) {
		switch(Character.toUpperCase(u)){
		case 'N':
		case '-':
			return NOBASE;
		}
		return decodeBaseChar(u);
	}
	/**
	 * Assumes both 
	 */
	public boolean allowBase(int oldBase_wflag, int testBase) {
		if (oldBase_wflag - noFlags(oldBase_wflag) == LOCK_FLAG()){
			return noFlags(oldBase_wflag) == testBase;
		}
		oldBase_wflag /= getNumMonomers();
		return flagsAllowsMap[oldBase_wflag][getNormalBaseFromZero(testBase)]==1;
	}

	public abstract int decodeBaseChar(char charAt);
	
	public abstract String displayBase(int i);

	/**
	 * Include the "No Monomer" in this.
	 * 
	 * Must be == getMonomers().length + 1
	 */
	public abstract int getNumMonomers();

	/**
	 * Do not include the "no monomer" in this.
	 */
	public abstract int[] getMonomers();

	public abstract int complement(int i);
	
	public abstract int getNormalBaseFromZero(int nonnormalBase);
	
	public abstract int bindScore(int base, int base2);
}
