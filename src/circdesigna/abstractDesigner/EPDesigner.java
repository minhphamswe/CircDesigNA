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
package circdesigna.abstractDesigner;

import java.util.TreeMap;

import circdesigna.CircDesigNA;


/**
 * A (senseless?) simple design algorithm which merely mutates each population member, if the improvements are deleterious, revert them.
 *   
 * @author Benjamin
 */
public class EPDesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	public EPDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}
	
	public void runBlockIteration_ (CircDesigNA runner, double endThreshold) {
		TreeMap<Double, T> sorted = new TreeMap();
		for(int i = 0; i < populationSize; i++){
			T toMutate = population_mutable[i];
			SingleDesigner.mutateAndTestAndBackup(toMutate);
			sorted.put(SingleDesigner.getOverallScore(toMutate),toMutate);
		};
		setBestChild(sorted.get(sorted.firstKey()));
	}
}
