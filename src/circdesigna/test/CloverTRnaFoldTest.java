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
package circdesigna.test;

import circdesigna.DomainDefinitions;
import circdesigna.CircDesigNA;
import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.ConstraintsNAFoldingImpl;
import circdesigna.energy.NAFolding;

public class CloverTRnaFoldTest {
	//The following sequence should fold into a clover.
	//AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT
	//Correct duplexes:
	//1-4x
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		NAFolding na = new ConstraintsNAFoldingImpl(config);
		String seq = "AAATGGCCAAACAGGCCGGCGCCGAACGCCCGGGAGCAGCCCGATTT";
		int[][] domain = new int[1][seq.length()];
		int[][] domainMark= new int[1][seq.length()];
		for(int k = 0; k < seq.length(); k++){
			domain[0][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
		}
		DomainDefinitions dsd = new DomainDefinitions(config);
		DomainSequence ds = new DomainSequence();
		ds.setDomains(0, null);
		System.out.println(na.mfe(ds, domain, domainMark));
	}
}
