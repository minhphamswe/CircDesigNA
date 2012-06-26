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
import circdesigna.BannedPatterns;
import circdesigna.DomainSequence;
import circdesigna.SequencePenalties;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

public class SequencePenaltiesImpl extends CircDesigNASystemElement implements SequencePenalties{
	public SequencePenaltiesImpl(CircDesigNAConfig System) {
		super(System);
		banlist = new BannedPatterns(Std.bannedPatternsList, Std);
	}
	private BannedPatterns banlist;

	/**
	 * This score adds a penalty for each match in the "banned patterns" list, according
	 * to the weights specified there. 
	 * 
	 * @param domain_markings 
	 */
	public double getSequenceScore(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		int n = seq.length(domain_markings);
		double sumResult = 0;
		for(int i = 0; i < n; i++){
			for(int pattern = 0; pattern < banlist.patternSize(); pattern++){
				if (banlist.matches(pattern, seq, domain, domain_markings, i)){
					sumResult += banlist.patternWeight(pattern);
				}
			}
		}
		
		return sumResult;
	}
}
