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


/**
 * This class defines the syntax / encoding of DNA, and provides some basic properties (complement)
 */
public class DnaDefinition extends MonomerDefinition{
	//Explicit definition of DNA: Mapping from a subset of Z to all of DNA
	public static final int A = 1, T = 2, G = 3, C = 4, P = 5, Z = 6, MAX_DNA_BASE = Z+1;
	private final int[] baseOrdering = new int[]{A,T,G,C,P,Z}; 
	
	//Used by bindScore
	private final int GCstr = -20;
	private final int ATstr = -10;
	private final int DHstr = GCstr;
	private final int PZstr = GCstr;
	private final int GTstr = -1;
	/**
	 * Really rough scoring function, returns some representative negative value if two bases
	 * are complementary (wobble gets less, etc.)
	 * 
	 * A 0 value should be seen as "don't even consider pairing these".
	 */
	public final int bindScore(int a, int b){
		a = noFlags(a);
		b = noFlags(b);
		if ((a==A && b==T) || (a==T && b==A)){
			return ATstr;
		} else
			if ((a==G && b==C) || (a==C && b==G)){
				return GCstr;
			} else
				//if ((a==D && b==H) || (a==H && b==D)){
				//	return DHstr;
				//} else 
					if ((a==P && b==Z) || (a==Z && b==P)){
						return PZstr;
					} else
						if ((a==G && b==T) || (a==T && b==G)){
							return GTstr;
						} else {
							return 0;
						}
	}
	/**
	 * Returns a string representation of a base int.
	 * @param base
	 * @return
	 */
	public String displayBase(int base) {
		base = noFlags(base);
		switch(base){
		case G:
			return "G";
		case A:
			return "A";
		case C:
			return "C";
		case T:
			return "T";
//		case D:
//			return "D";
//		case H:
//			return "H";
		case P:
			return "P";
		case Z:
			return "Z";
		default:
			throw new IllegalArgumentException("Unrecognized Base: "+base);
		}
	}
	/**
	 * See getNormalBase, but the number returned will be in [0,3].
	 */
	public int getNormalBaseFromZero(int nonnormalBase) {
		return getNormalBase(nonnormalBase)-1;
	}
	/**
	 * Returns a base from {A,C,T,G} which is most like 'x'.
	 */
	public int getNormalBase(int x){
		switch(x){
		case G:
			return G;
		case A:
			return A;
		case C:
			return C;
		case T:
			return T;
//		case D:
//			return C;
//		case H:
//			return G;
		case P:
			return C;
		case Z:
			return G;
		default:
			throw new IllegalArgumentException("Unrecognized Base: "+x);
		}
	} 

	/**
	 * Inverse of displayBase: turns characters into numbers.
	 * @param charAt
	 * @return
	 */
	public int decodeBaseChar(char charAt){
		charAt = Character.toUpperCase(charAt);
		switch(charAt){
		case 'G':
			return G;
		case 'A':
			return A;
		case 'C':
			return C;
		case 'T':
			return T;
//		case 'D':
//			return D;
//		case 'H':
//			return H;
		case 'P':
			return P;
		case 'Z':
			return Z;
		default:
			throw new RuntimeException("Invalid char: "+charAt);
		}
	}

	/**
	 * Returns the unflagged complement of base i.
	 * Needs to know the numeric definitions of the bases.
	 */
	public int complement(int i) {
		i = noFlags(i);
		if ((i&1)!=0){
			return i+1;
		}
		return i-1;
	}
	public int getNumMonomers() {
		return MAX_DNA_BASE;
	}
	public int[] getMonomers() {
		return baseOrdering;
	}
}
