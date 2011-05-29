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


/**
 * Self documenting class name!
 */
public class FitnessPopulationDesignMember <T extends PopulationDesignMember<T>> implements Comparable <FitnessPopulationDesignMember<T>>{
	public T myKey;
	public double myScore;
	public int compareTo(FitnessPopulationDesignMember<T> o) {
		double diff = myScore - o.myScore;
		if (diff==0){
			return 0;
		} else if (diff < 0){
			return -1;
		} else {
			return 1;
		}
	}

}
