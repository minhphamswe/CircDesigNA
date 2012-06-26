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
package circdesigna.config;

import circdesigna.Connector;
import circdesigna.GeneralizedInteractiveRegion;

public class CircDesigNASystemElement extends SystemElement<CircDesigNAConfig>{
	/**
	 * Convenience function returns base i of the sequence of a GIR
	 */
	public final int base(GeneralizedInteractiveRegion ds, int i, int[][] domain){
		return ds.base(i,domain,Std.monomer);
	}
	public final int base(Connector con, int i, int[][] domain){
		return con.base(i, domain, Std.monomer);
	}
	public final int Nbase(GeneralizedInteractiveRegion ds, int i, int[][] domain){
		return Std.monomer.getNormalBaseFromZero(ds.base(i,domain,Std.monomer));
	}
	public CircDesigNASystemElement(CircDesigNAConfig System) {
		super(System);
	}
}
