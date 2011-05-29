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
package circdesigna;

/**
 * Determines a set of mutation rules.
 */
public interface DesignerCode {
	/**
	 * By some arbitrary set of constraints, determine whether the sequence of domain #"whichDomain" is
	 * "valid". @see mutateToOther
	 */
	public boolean isValid(int[][] domain, int whichDomain);
	/**
	 * Given that domain no. whichDomain is valid, returns true if it is possible to change base i of that domain
	 * to another base while still maintaining validity. If so, the post condition is that domain no. whichDomain 
	 * has been mutated to a domain which is valid (conforms to isValid) but has a different base at location i 
	 * that it had before. In a sense, false is only returned if base i cannot be mutated without producing an invalid
	 * domain.
	 * 
	 * If domain no. whichDomain is invalid, then this method returns true if it is possible to generate any new domain which
	 * is valid, this becomes domain #whichDomain. A case that this returns false is when the sequence constraints
	 * are too strict to allow for any satisfying domain. 
	 * 
	 * Remark: While this is the specification of the API, a weaker set of conditions are probably actually implemented.
	 * Namely, the method may only return true if it is able to detect a way of mutating base i and producing a valid
	 * domain. Even if a mutation exists, if that method is not detected by the algorithm implemented, than false will
	 * be returned. A good mutation scheme should minimize this problem.
	 * 
	 * Remark: do not assume that base i was the only base that was mutated.
	 */
	public boolean mutateToOther(int[][] domain, int whichDomain, int i);
	/**
	 * If possible, make base i of whichdomain newbase. Returns false if this would invalidate
	 * the sequence constraints.
	 */
	public boolean mutateToOther(int[][] domain, int whichDomain, int i, int newbase);
}
