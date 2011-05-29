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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import circdesigna.CircDesigNA.ScorePenalty;
import circdesigna.abstractDesigner.PopulationDesignMember;


/**
 * Concrete population member.
 */
public class CircDesigNAPMemberImpl extends PopulationDesignMember<CircDesigNAPMemberImpl>{
	protected CircDesigNAPMemberImpl designerCopyConstructor() {
		return new CircDesigNAPMemberImpl(null, null, null, null);
	}
	public void seedFromOther(CircDesigNAPMemberImpl pdm) {
		if (penalties==null){
			penalties = new ArrayList(); 
		}
		penalties.clear(); //Penalties are lightweight!
		for(ScorePenalty q : pdm.penalties){
			penalties.add((ScorePenalty)q.clone());
		}
		scoredElements = copy2dInt(pdm.scoredElements,scoredElements);
		domain = copy2dInt(pdm.domain,domain);
		domain_markings = copy2dInt(pdm.domain_markings,domain_markings);
		//totalBases = pdm.totalBases;
	}
	private static int[][] copy2dInt(int[][] se, int[][] toRet) {
		if (toRet==null || toRet.length < se.length){
			toRet = new int[se.length][];
			for(int k = 0; k < se.length; k++){
				int[] toCopy = se[k];
				toRet[k] = new int[toCopy.length];
			}
		}
		for(int k = 0; k < se.length; k++){
			int[] toCopy = se[k];
			System.arraycopy(toCopy, 0, toRet[k],0,toRet[k].length);
		}
		return toRet;
	}
	public List<ScorePenalty> penalties; 
	public int[][] scoredElements; 
	public int[][] domain; 
	public int[][] domain_markings;
	public int totalBases;
	public CircDesigNAPMemberImpl(List<ScorePenalty> penalties, int[][] scoredElements, int[][] domain, int[][] domain_markings){
		this.penalties = penalties;
		this.scoredElements = scoredElements;
		this.domain = domain; 
		this.domain_markings = domain_markings;
	}
}
