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
package circdesigna.energy;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.Scanner;

import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;


/**
 * Implements MFE prediction and folding score functions, using the MCS algorithm as an approximation.
 */
public class NupackFolder extends CircDesigNASystemElement implements NAFolding{

	public NupackFolder(CircDesigNAConfig System) {
		super(System);
	}

	public double mfe(DomainSequence seq1, DomainSequence seq2, int[][] domain,
			int[][] problemAreas) {
		throw new RuntimeException("Not available with this folding tool.");
	}

	public double mfe(DomainSequence domainSequence, int[][] domain,
			int[][] domain_markings) {
		throw new RuntimeException("Not available with this folding tool.");
	}

	public double mfeNoDiag(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings) {
		throw new RuntimeException("Not available with this folding tool.");
	}

	public double mfeStraight(DomainSequence domainSequence,
			DomainSequence domainSequence2, int[][] domain,
			int[][] domain_markings, int markLeft, int markRight, int joffset) {
		throw new RuntimeException("Not available with this folding tool.");
	}
	
	
////////////////////////////////
//// Pair Probability functions - Warning, not maintained
////////////////////////////////
	
	public void pairPr(double[][] pairsOut, DomainSequence seq1, DomainSequence seq2, int[][] domain) {
		pairPr_viaNUPACK(pairsOut, new DomainSequence[]{seq1,seq2}, domain);
	}

	public void pairPr(double[][] pairsOut, DomainSequence seq, int[][] domain) {
		pairPr_viaNUPACK(pairsOut, new DomainSequence[]{seq}, domain);
	}
	private static int ct = 0;
	private NupackRuntime NUPACKLINK = null;
	private static class NupackRuntime {
		public Process exec;
		public Scanner in;
		public PrintWriter out;
		public void finalize() {
			exec.destroy();
			in.close();
			out.close();
		}
	}
	private void pairPr_viaNUPACK(double[][] pairsOut, DomainSequence[] seqs, int[][] domain) {
		try {
			System.out.println("Going to nupack"+ct++);
			
			if (NUPACKLINK==null){
				NUPACKLINK = new NupackRuntime();
				NUPACKLINK.exec = Runtime.getRuntime().exec("/home/Benjamin/Code/C/nupack3.0/bin/pairs -T 37 -material dna -cutoff 0.001 -multi", new String[]{"NUPACKHOME=/home/Benjamin/Code/C/nupack3.0"});
				NUPACKLINK.out = new PrintWriter(new OutputStreamWriter(new BufferedOutputStream(NUPACKLINK.exec.getOutputStream())));
				NUPACKLINK.in = new Scanner(new BufferedInputStream(NUPACKLINK.exec.getInputStream()));
			}

			NUPACKLINK.out.println("output");
			NUPACKLINK.out.println(seqs.length);
			
			int N = 0;
			for(int i = 0; i < seqs.length; i++){
				int seqLen = seqs[i].length(domain);
				for(int k = 0; k < seqLen; k++){
					NUPACKLINK.out.print(Std.monomer.displayBase(base(seqs[i], k, domain)));
				}
				N += seqLen;
				NUPACKLINK.out.println();
			}
			//Clear probability matrix
			for(int i = 0; i < N; i++){
				for(int j = 0; j < N+1; j++){
					pairsOut[i][j] = 0;
				}
			}
			for(int i = 0; i < seqs.length; i++){
				NUPACKLINK.out.print((i+1)+" ");
			}
			NUPACKLINK.out.println();
			
			if (true){
				while(NUPACKLINK.in.hasNextLine()){
					String line = NUPACKLINK.in.nextLine();
					System.out.println(line);
					if (line.equals("DONE")){
						break;
					}
				}
			}
			
			Scanner in2;
			if (seqs.length==1){
				in2 = new Scanner(new File("output.ppairs"));
			} else {
				in2 = new Scanner(new File("output.epairs"));
			}
			while(in2.hasNextLine()){
				String line[] = in2.nextLine().trim().split("\\s+");
				if (line.length==3){
					if (line[0].startsWith("%")){
						continue;
					}
					int iB = new Integer(line[0])-1;
					int jB = new Integer(line[1])-1;
					pairsOut[iB][jB] = new Double(line[2]);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		/*
		for(int k = 0; k < length1; k++){
			for(int j = k; j < length1; j++){
				pairsOut[k][j] = Math.random();
				pairsOut[j][k] = Math.random(); 
			}	
		}
		*/
	}

}
