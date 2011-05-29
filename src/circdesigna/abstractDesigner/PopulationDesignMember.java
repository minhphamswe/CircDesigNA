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
 * A population member for a genetic algorithm inspired designer.
 */
public abstract class PopulationDesignMember<T extends PopulationDesignMember> implements Comparable<T> {
	//Population design members are sorted by number in the population.
	private int myID = 0;
	public final int compareTo(T o) {
		return myID - o.myID;
	}
	public int getID(){
		return myID;
	}
	public T designerCopyConstructor(int myID){
		T toRet = designerCopyConstructor();
		toRet.myID = myID;
		toRet.seedFromOther(this);
		return toRet;
	}
	/**
	 * Creates a new instance of this design member. A deep copy is not required, as
	 * this call will always be followed up with a call to "seed".
	 */
	protected abstract T designerCopyConstructor();
	public abstract void seedFromOther(T pdm);
}
