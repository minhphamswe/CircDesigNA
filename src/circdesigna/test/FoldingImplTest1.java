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

import java.util.Arrays;

import circdesigna.DomainDefinitions;
import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.CircDesigNAMCSFolder;

public class FoldingImplTest1 {
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		config.setMode(CircDesigNAConfig.RNA_MODE);
		for(int i = 0; i < 1; i++){
			DomainSequence seqS = new DomainSequence();
			DomainDefinitions dsd = new DomainDefinitions(config);
			seqS.setDomains(0, null);
			int[][] domain = new int[1][];
			String seq = "GTGATAGACAC";
			int seqLength;
			if (seq==null){
				seqLength = (int)(Math.random()*10+10);
			} else {
				seqLength = seq.length();
			}
			domain[0] = new int[seqLength];
			int[][] domainMark= new int[1][seqLength];
			if (seq!=null){
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
				}
			} else {
				for(int k = 0; k < seqLength; k++){
					domain[0][k] = randomChoice(1,config.monomer.getNumMonomers()-1);
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			CircDesigNAMCSFolder fl = new CircDesigNAMCSFolder(config);
			final double viaMatrix = fl.mfe(seqS, domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			System.out.println(viaMatrix);
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			//final double viaUnafold = fl.foldSingleStranded_viaUnafold(seqS, domain, domainMark);
			//System.out.println(Arrays.deepToString(domainMark));
			//System.out.println(seqLength+" "+-viaMatrix+" "+-viaUnafold);
		}
	}

	private static int randomChoice(int i, int j) {
		return ((int) (Math.random()*(j-i)))+i;
	}
}
