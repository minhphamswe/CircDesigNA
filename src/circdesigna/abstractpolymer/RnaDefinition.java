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
 * This class defines the syntax / encoding of RNA
 */
public class RnaDefinition extends DnaDefinition{
	//Explicit definition of DNA: Mapping from a subset of Z to all of DNA
	public static final int A = DnaDefinition.A, U = DnaDefinition.T, 
		G = DnaDefinition.G, C = DnaDefinition.C, 
		//D = DnaDefinition.D, H = DnaDefinition.H, 
		P = DnaDefinition.P, Z = DnaDefinition.Z;

	public int decodeBaseChar(char charAt) {
		charAt = Character.toUpperCase(charAt);
		if (charAt=='U' || charAt=='T'){
			return U;
		}
		return super.decodeBaseChar(charAt);
	}
	public String displayBase(int base) {
		if (base==U){
			return "U";
		} else {
			return super.displayBase(base);
		}
	}
}
