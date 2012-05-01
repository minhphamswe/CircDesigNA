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

import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;

/**
 * Run this after modifying the interaction scoring function; run a linear regression on results
 * to find a baseline for self-similarity scoring
 * @author Benjamin
 */
public class SelfSimilarityBaseline {
	public static void main(String[] args){
		CircDesigNAConfig cfg = new CircDesigNAConfig();
		cfg.setMode(CircDesigNAConfig.RNA_MODE);
		ConstraintsNAFoldingImpl fli = new ConstraintsNAFoldingImpl(cfg);
		for(int u = 0; u < 1000; u++){
			int len = (int) (Math.random()*8000+10);
			int[][] domain = new int[1][len];
			int[][] nullMark = new int[1][len];
			for(int k = 0; k < domain.length; k++){
				for(int y = 0; y < domain[k].length; y++){
					domain[k][y] = (int) (Math.random()*4 + 1);
				}
			}
			DomainSequence ds1 = new DomainSequence();
			DomainSequence ds2 = new DomainSequence();
			ds1.setDomains(0,null);
			ds2.setDomains(0 | DomainSequence.NA_COMPLEMENT_FLAG,null);
			//System.out.println(len+" "+fli.mfeNoDiag_NoBaseline(ds1, ds2, domain, nullMark));
		}
	}
}
