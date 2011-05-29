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
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import circdesigna.CircDesigNAOptions;
import circdesigna.DesignScoreBreakdown;
import circdesigna.CircDesigNA;
import circdesigna.CircDesigNA_SharedUtils;
import circdesigna.SequenceDesigner;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.impl.CircDesigNAImpl;



/**
 * Takes a file in the format
 * Domain Defs
 * Blank Line
 * Molecules
 * 
 * reproduces it N times
 * 
 * Runs the system for K iterations P times
 * 
 * Outputs 4 files: 
 *  - Net score
 *  - Interactions score
 *  - Hairpins
 *  - Selffolding score
 *
 * where each file has P elements on K lines corresponding to the score of each run of the system at every iteration.
 * 
 * All default options are used (population size is 30, domains begin/end with G/C, etc.)
 */
public class DomainDesignerTest2 {
	public static void main(String[] args) throws FileNotFoundException{
		String inFile = System.getProperty("infile")+"/";
		String outFile = System.getProperty("outfile")+"/";
		//runTest(inFile+"Annihilator.txt",outFile+"Annihilator.txt", true);
		//runTest(inFile+"Transducer.txt",outFile+"Transducer.txt", false);
		System.out.println("WARNING!");
		runTest(inFile+"8000kb.txt",outFile+"8000kb.txt", false);
		//CircDesigNA.ENABLE_MARKINGS = false;
		//runTest(inFile+"8000kb.txt",outFile+"8000kb-NoMarkings.txt", false);
		//runTest(inFile+"Transducer.txt",outFile+"Transducer-NoMarkings.txt", false);
		//runTest(inFile+"TwoWayJoin.txt",outFile+"TwoWayJoin.txt");
	}
	public static void runTest(String infile, String outfile, boolean append) throws FileNotFoundException{
		int N = 1;
		int K = 20;
		int P = 3;

		int startRunNumber = 0;
		OutputHandler[] oh = makeOutputHandlers(new String[]{"net_score","cross_interaction","helixbreathe","selffold"},K,P);
		if (append){
			for(OutputHandler o : oh){
				File outFile = new File(outfile+o.name+".txt");
				Scanner in = new Scanner(outFile);
				int itr = 0;
				while(in.hasNextLine()){
					String[] line = in.nextLine().trim().split("\\s+");
					if (line.length<=1){
						continue;
					}
					startRunNumber = line.length;
					for(int p = 0; p < startRunNumber; p++){
						o.addScore(itr, new Double(line[p]));
					}
					itr++;
				}
				//System.out.println(o.toString());
			}
		}
		
		if (startRunNumber > 0){
			System.out.printf("Found first iterations %d. Continuing.\n", startRunNumber);
		}
		//System.exit(1);
		
		File inFile = new File(infile);
		{
			File outFile = new File(outfile+"desc"+".txt");
			outFile.getParentFile().mkdir();
			PrintStream out = new PrintStream(new FileOutputStream(outFile));
			out.printf("%s %d %d %d\n",inFile.toString(),N,K,P);
			out.close();
		}
		try {
			runFile(inFile,N,K,P,oh,startRunNumber,outfile);
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}
	private static class OutputHandler{
		private ArrayList<Double>[] rows;
		private String name;
		public OutputHandler(String name, int rownum){
			rows = new ArrayList[rownum];
			for(int i = 0; i < rows.length; i++){
				rows[i] = new ArrayList();
			}
			this.name = name;
		}
		public void addScore(int itr, double x) {
			rows[itr].add(x);
		}
		public String toString(){
			StringBuffer toRet = new StringBuffer();
			for(ArrayList<Double> k : rows){
				for(Double d : k){
					toRet.append(String.format("%.3f ",d));
				}
				toRet.append("\n");
			}
			return toRet.toString();
		}
	}
	private static OutputHandler[] makeOutputHandlers(String[] names, int k, int p) {
		OutputHandler[] toRet = new OutputHandler[names.length];
		for(int i = 0; i < toRet.length; i++){
			toRet[i] = new OutputHandler(names[i],k);
		}
		return toRet;
	}
	private static void runFile(File inFile, int n, int k, int p, OutputHandler[] oh, int startRunNumber, String outfile) throws Throwable {
		String domainDefs, moleculeDefs;
		{
			Scanner in = new Scanner(inFile);
			StringBuffer[] out = new StringBuffer[]{new StringBuffer(), new StringBuffer()};
			int ind = 0;
			while(in.hasNextLine() && ind < out.length){
				String line = in.nextLine().trim();
				if (line.equals("")){
					ind++;
				} else {
					out[ind].append(line+"\n");
				}
			}
			domainDefs = out[0].toString();
			moleculeDefs = out[1].toString();
		}
		for(int i = 0; i < n-1; i++){
			String[] result = CircDesigNA_SharedUtils.duplicateSystem(domainDefs, moleculeDefs);
			domainDefs = result[0];
			moleculeDefs = result[1];
		}

		System.out.println("Read file: ");
		System.out.println("Domain definitions:");
		System.out.println(domainDefs);
		System.out.println("Molecules:");
		System.out.println(moleculeDefs);
		
		//Run the system.
		for(int kp = startRunNumber; kp < p; kp++){
			CircDesigNAConfig Std = new CircDesigNAConfig();
			SequenceDesigner<CircDesigNAOptions> dd = CircDesigNA.getDefaultDesigner(moleculeDefs, domainDefs, Std);
			CircDesigNAOptions options = dd.getOptions();
			
			do {
				dd.runIteration();
				if (dd.isEndConditionError()){
					System.err.println("Broke on error:");
					System.err.println(dd.getResult());
					break;
				}
				int itr = dd.getIterationCount()-1;
				DesignScoreBreakdown scoreBreakdown = dd.getScoreBreakdown();
				oh[0].addScore(itr,scoreBreakdown.netScore);
				oh[1].addScore(itr,scoreBreakdown.crossInteractionsOnly);
				oh[2].addScore(itr,scoreBreakdown.breathingHelixes);
				oh[3].addScore(itr,scoreBreakdown.selfFoldOnly);
			} while(dd.getIterationCount() < k);

			//Output final result:
			{
				File outFile = new File(outfile+"individualOutputs/out"+kp+".txt");
				outFile.getParentFile().mkdir();
				PrintStream out = new PrintStream(new FileOutputStream(outFile));
				out.println(dd.getResult());	
				out.close();
			}
			dd.abort();
			
			//Add this solution trajectory to the trajectories set:
			for(OutputHandler o : oh){
				File outFile = new File(outfile+o.name+".txt");
				PrintStream out = new PrintStream(new FileOutputStream(outFile));
				out.println(o);
				out.close();
			}
		}
		
	}
}
