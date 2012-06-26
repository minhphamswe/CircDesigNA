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
import circdesigna.energy.ConstraintsNAFoldingImpl;
import circdesigna.energy.NAFolding;
import circdesigna.energy.UnafoldFolder;

public class FoldingImplTest4 {
	public static void main(String[] args) throws FileNotFoundException{
		boolean doCorrelation = true;
		
		File outFile = new File(System.getProperty("inFile")+".out");
		CircDesigNAConfig config = new CircDesigNAConfig();
		System.out.println(outFile.getAbsolutePath());
		//System.setOut(new PrintStream(new FileOutputStream(outFile)));
		//Error in this file. third column represents pairing of seq1 with seq1! Embarrasing.
		//Scanner in = new Scanner(new File("C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\Circ_DesigNAPaper\\AssemblaRepo\\circdesignapaper_w_figures\\scoresComparison_SS.txt"));
		Scanner in = new Scanner(new File(System.getProperty("inFile")));
		in.nextLine();
		ConstraintsNAFoldingImpl fl = new ConstraintsNAFoldingImpl(config);
		NAFolding fl_unafold = new UnafoldFolder(config);
		int[][] domain = new int[2][];
		int[][] domain_markings = new int[2][];

		DomainSequence ds1 = new DomainSequence();
		ds1.setDomains(0, null);
		DomainSequence ds2 = new DomainSequence();
		ds2.setDomains(1, null);
		if (doCorrelation){
			System.out.printf("Nupack\tFoldNoLoops\n");
		}
		
		
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
			fl.setScoringModel(3);
			double resultLocal3 = 0;
			try {
				resultLocal3 = fl.mfe(ds1, domain, domain_markings, true);
			} catch (Throwable e){
				//e.printStackTrace();
			}
			double dt3 = (System.nanoTime()-now)/1e9;
			now = System.nanoTime();
			fl.setScoringModel(2);
			double resultLocal2 = 0;
			try {
				resultLocal2 = fl.mfe(ds1, domain, domain_markings, true);
			} catch (Throwable e){
				//e.printStackTrace();
			}
			double dt2 = (System.nanoTime()-now)/1e9;
			now = System.nanoTime();
			fl.setScoringModel(1);
			double resultLocal1 = 0;
			try {
				resultLocal1 = fl.mfe(ds1, domain, domain_markings, true);
			} catch (Throwable e){
				//e.printStackTrace();
			}
			double dt1 = (System.nanoTime()-now)/1e9;
			//System.out.println(String.format("%d M3 %.2f %.3e M2 %.2f %.3e M1 %.2f %.3e",ds1.length(domain), resultLocal3, dt3, resultLocal2, dt2, resultLocal1, dt1));
			if (line.length > 3){
				Double nupackResult = new Double(line[3]);
				double error3 = (resultLocal3 - nupackResult);
				double error2 = (resultLocal2 - nupackResult);
				double error1 = (resultLocal1 - nupackResult);
				//System.out.printf("Nupack: %.2f errN^3: %.2f errN^2: %.2f errN^1: %.2f\n", new Double(line[3]), error3, error2, error1);
				if (doCorrelation){
					if (error1 < 0){
						throw new RuntimeException("Assertion error");
					}
					System.out.printf("%.3f\t%.3f\n",nupackResult, resultLocal1);
				}
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
