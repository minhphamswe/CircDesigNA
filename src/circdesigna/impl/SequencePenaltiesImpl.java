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
package circdesigna.impl;

import static circdesigna.abstractpolymer.DnaDefinition.A;
import static circdesigna.abstractpolymer.DnaDefinition.C;
import static circdesigna.abstractpolymer.DnaDefinition.G;
import static circdesigna.abstractpolymer.DnaDefinition.T;
import circdesigna.DomainSequence;
import circdesigna.SequencePenalties;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

public class SequencePenaltiesImpl extends CircDesigNASystemElement implements SequencePenalties{
	public SequencePenaltiesImpl(CircDesigNAConfig System) {
		super(System);
	}

	/**
	 * This routine checks for potentially problematic (hard to synthesize) DNA sequences.
	 * 
	 * Amounts to poly-N checking, and uses the same routines as David Zhang's Domain Designer
	 * Penalizes:
	 *    +1 for each GGGG or CCCC
	 *    +1 for each run of As and Ts of length 6
	 *    +1 for each run of Gs and Cs of length 6
	 * 
	 * @param domain_markings 
	 */
	public double getSynthesizabilityScore(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		int n = seq.length(domain_markings);
		double sumResult = 0;
		int[] baseCounts4 = new int[Std.monomer.getNumMonomers()];
		int[] baseCounts6 = new int[Std.monomer.getNumMonomers()];
		for(int i = 0; i < n; i++){
			//Remove the counts 4/6 spaces back
			if (i >= 4){
				int prior = base(seq,i-4,domain);
				baseCounts4[prior]--;
			}
			if (i >= 6){
				int prior = base(seq,i-6,domain);
				baseCounts6[prior]--;
			}
			//Add to the current count
			int now = base(seq,i,domain);
			baseCounts4[now]++;
			baseCounts6[now]++;
			//System.out.println(now+" "+Arrays.toString(baseCounts4)+" "+Arrays.toString(baseCounts6));
			
			if (baseCounts4[G]==4 || baseCounts4[C]==4){
				sumResult++;
				seq.mark(i, -4, domain, domain_markings);
			}
			
			if (baseCounts6[A]+baseCounts6[T]==6 || baseCounts6[G]+baseCounts6[C]==6){
				sumResult++;
				seq.mark(i, -6, domain, domain_markings);
			}
		}
		
		return sumResult;
	}
}
