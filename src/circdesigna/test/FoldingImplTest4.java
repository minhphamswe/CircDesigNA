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

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.NAFolding;
import circdesigna.energy.UnafoldFolder;
import circdesigna.energy.UnpseudoknottedFolder;

public class FoldingImplTest4 {
	public static void main(String[] args) throws FileNotFoundException{
		File outFile = new File(System.getProperty("inFile")+".out");
		CircDesigNAConfig config = new CircDesigNAConfig();
		System.out.println(outFile.getAbsolutePath());
		//System.setOut(new PrintStream(new FileOutputStream(outFile)));
		//Error in this file. third column represents pairing of seq1 with seq1! Embarrasing.
		//Scanner in = new Scanner(new File("C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\Circ_DesigNAPaper\\AssemblaRepo\\circdesignapaper_w_figures\\scoresComparison_SS.txt"));
		Scanner in = new Scanner(new File(System.getProperty("inFile")));
		in.nextLine();
		UnpseudoknottedFolder fl = new UnpseudoknottedFolder(config);
		NAFolding fl_unafold = new UnafoldFolder(config);
		int[][] domain = new int[2][];
		int[][] domain_markings = new int[2][];

		DomainSequence ds1 = new DomainSequence();
		ds1.setDomains(0, null);
		DomainSequence ds2 = new DomainSequence();
		ds2.setDomains(1, null);
		while(in.hasNextLine()){
			String line2 = in.nextLine();
			String[] line = line2.trim().split(" ");
			for(int j = 0; j < 2; j++){
				int seqLength = line[j+1].length();
				domain[j] = new int[seqLength];
				domain_markings[j] = new int[seqLength];
				for(int k = 0; k < seqLength; k++){
					domain[j][k] = config.monomer.decodeConstraintChar(line[j+1].charAt(k));
				}
			}
			//double resultLocal = fl.mfe(ds1, domain, domain_markings);
			long now = System.nanoTime();
			fl.setOrder(3);
			double resultLocal3 = fl.mfe(ds1, domain, domain_markings);
			double dt3 = (System.nanoTime()-now)/1e9;
			now = System.nanoTime();
			fl.setOrder(2);
			double resultLocal2 = fl.mfe(ds1, domain, domain_markings);
			double dt2 = (System.nanoTime()-now)/1e9;
			//System.out.println(String.format("%d %.3e %.3e %.3e %.3e",ds1.length(domain), resultLocal2, resultLocal3, dt2, dt3));
			if (line.length > 3){
				double error3 = Math.abs(new Double(line[3]) - resultLocal3);
				double error2 = Math.abs(new Double(line[3]) - resultLocal2);
				System.out.printf("Nupack: %.2f errN^3: %.2f errN^2: %.2f\n", new Double(line[3]), error3, error2);
			} else {
				//System.out.println(line2+" "+resultLocal3);
			}
			//System.out.println(Arrays.toString(domain_markings[0]));
			Arrays.fill(domain_markings[0],0);
			//double result = fl.mfeHybridDeltaG_viaUnafold(ds1, ds1, domain, domain_markings);
			double result = 0; //fl_unafold.mfe(ds1, domain, domain_markings);
			//System.out.println(Arrays.toString(domain_markings[0]));
		}
		System.out.flush();
	}
}
