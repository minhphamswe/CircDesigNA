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
package circdesigna.batch;

import java.util.Scanner;

public class StringSplitter {
	/**
	 * Used for using NUPACK as a designer.
	 */
	public static void main(String[] args){
		Scanner in = new Scanner(System.in);
		String line = in.nextLine();
		int numWrap = new Integer(in.nextLine().split("\\s+")[0]);
		for(int k = 0; k < line.length();){
			for(int y = 0; y < numWrap; y++, k++){
				System.out.print(line.charAt(k));
			}
			if (k < line.length()){
				System.out.println();
			}
		}
	}
}
