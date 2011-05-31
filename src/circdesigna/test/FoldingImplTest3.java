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

import static circdesigna.abstractpolymer.DnaDefinition.A;
import static circdesigna.abstractpolymer.DnaDefinition.C;
import static circdesigna.abstractpolymer.DnaDefinition.G;
import static circdesigna.abstractpolymer.DnaDefinition.T;

import java.util.Arrays;


import circdesigna.AbstractDomainDesignTarget;
import circdesigna.DomainDefinitions;
import circdesigna.DomainSequence;
import circdesigna.AbstractDomainDesignTarget.HairpinClosingTarget;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.CircDesigNAMCSFolder;
import circdesigna.impl.CircDesigNAImpl;
import circdesigna.impl.SequencePenaltiesImpl;
import circdesigna.impl.CircDesigNAImpl.HairpinOpening;


public class FoldingImplTest3 {
	public static void main(String[] args){
		CircDesigNAConfig config = new CircDesigNAConfig();
		for(int i = 0; i < 1; i++){
			CircDesigNAImpl impl = new CircDesigNAImpl(new CircDesigNAMCSFolder(config),new SequencePenaltiesImpl(config), config);
			DomainDefinitions dsd = new DomainDefinitions(config);
			HairpinClosingTarget hairpin = new AbstractDomainDesignTarget(dsd,config).new HairpinClosingTarget(1,0,0|DomainSequence.NA_COMPLEMENT_FLAG,2,true,null);
			HairpinOpening hairpinOpening = impl.new HairpinOpening(hairpin, null);
			DomainSequence seqS = new DomainSequence();
			seqS.setDomains(0, null);
			String[] seq = new String[3];
			int[][] domain = new int[seq.length][];
			int[][] domainMark= new int[seq.length][];
			seq[0] = "GTGATAGACAC";
			seq[1] = "CTAT";
			seq[2] = "ATAGCATAG";
			for(int j = 0; j < seq.length; j++){
				int seqLength;
				if (seq[j]==null){
					seqLength = (int)(Math.random()*10+10);
				} else {
					seqLength = seq[j].length();
				}
				domain[j] = new int[seqLength];
				domainMark[j] = new int[seqLength];
				if (seq!=null){
					for(int k = 0; k < seqLength; k++){
						domain[j][k] = config.monomer.decodeConstraintChar(seq[j].charAt(k));
					}
				} else {
					for(int k = 0; k < seqLength; k++){
						domain[j][k] = randomChoice(A,C,G,T);
					}
				}
			}
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			final double viaMatrix = hairpinOpening.evalScoreSub(domain, domainMark);
			System.out.println(Arrays.deepToString(domainMark));
			System.out.println(viaMatrix);
			for(int k = 0; k < domainMark.length; k++)Arrays.fill(domainMark[k],0);
			
			//final double viaUnafold = fl.foldSingleStranded_viaUnafold(seqS, domain, domainMark);
			//System.out.println(Arrays.deepToString(domainMark));
			//System.out.println(seqLength+" "+-viaMatrix+" "+-viaUnafold);
		}
	}

	private static int randomChoice(int a, int b, int c, int d) {
		switch((int)(Math.random()*4)){
		case 0:return a;
		case 1:return b;
		case 2:return c;
		case 3:return d;
		}
		throw new RuntimeException("Assertion error: randomChoice");
	}
}
