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

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import circdesigna.DomainDefinitions;
import circdesigna.config.CircDesigNAConfig;


import static circdesigna.batch.DesignMultipleTimes.*;

public class ExtractScoresFromRun {
	public static void main(String[] args) throws IOException{
		Scanner in = new Scanner(System.in);
		String file = in.nextLine();

		String Domains = readToEnd(in);
		String Molecules = readToEnd(in);
		
		DomainDefinitions dsd = new DomainDefinitions(new CircDesigNAConfig());
		dsd.readDomainDefs(Domains, dsd);
		
		DesignMultipleTimes.RunEvaluation(new File(file), Molecules, dsd, -1, false);
	}
}
