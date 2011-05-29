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

import circdesigna.CircDesigNA.ScorePenalty;
import circdesigna.impl.CircDesigNAPMemberImpl;



public abstract class TournamentDesigner <T extends PopulationDesignMember<T>>  extends BlockDesigner <T> {
	public TournamentDesigner(SingleMemberDesigner<T> SingleDesigner) {
		super(SingleDesigner);
	}

	public void tournamentSelect(int numReproduce){
		//Seed the fittest
		TreeMap<Double, T> populationView = new TreeMap();
		for(T q : population_mutable){
			double score = SingleDesigner.getOverallScore(q);
			populationView.put(-score, q); //sort descending
		}
		T fittest = populationView.get(populationView.lastKey());
		for(int k = 0; k < numReproduce && populationView.size()>=2; k++){
			T bottom = populationView.remove(populationView.firstKey());
			T replace = populationView.remove(populationView.lastKey());
			if (replace==bottom){
				throw new RuntimeException("Assertion error");
			}
			bottom.seedFromOther(replace);
		}
		setBestChild(fittest);
	}
}
