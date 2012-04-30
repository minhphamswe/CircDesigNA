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

import static circdesigna.abstractpolymer.DnaDefinition.A;
import static circdesigna.abstractpolymer.DnaDefinition.C;
import static circdesigna.abstractpolymer.DnaDefinition.G;
import static circdesigna.abstractpolymer.DnaDefinition.T;

import java.util.ArrayList;
import java.util.Scanner;

import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.ExperimentalDuplexParams.DangleScore;
import circdesigna.energy.ExperimentalDuplexParams.InteriorLoop;
import circdesigna.energy.ExperimentalDuplexParams.NearestNeighborScore;
import circdesigna.energy.ExperimentalDuplexParams.HairpinLoop;
import circdesigna.energy.ExperimentalDuplexParams.TerminalMismatchPairScore;

public class StandardizedThermoFileLoader {
	/**
	 * Parses a dG file in the format the nupack uses.
	 * 
	 * These files are commonly:
	 * RNA_mfold2.3.dG
	 * RNA_mfold3.0.dG
	 * DNA_mfold2.3.dG
	 * 
	 * and their corresponding dH files.
	 */
	public static void makeTable(ExperimentalDuplexParams ex, String dg, String dh) {
		CircDesigNAConfig std = ex.Std;
		
		ArrayList<NearestNeighborScore> nns = new ArrayList();
		//For the very last nearest neighbor in a hairpin or an inner loop.
		ArrayList<TerminalMismatchPairScore> misHairpin = new ArrayList(); 
		ArrayList<TerminalMismatchPairScore> misInnerLoop = new ArrayList();
		ArrayList<DangleScore> topDangle = new ArrayList();
		ArrayList<DangleScore> bottomDangle = new ArrayList();
		ArrayList<HairpinLoop> specialHairpins = new ArrayList();
		ArrayList<InteriorLoop> specialInteriorLoops = new ArrayList();
		
		Scanner in = new Scanner(dg);
		
		//System.out.println(dg);
		
		int[] bases = new int[]{A,C,G,T};
		int[][] headers = new int[][]{{A,T},{C,G},{G,C},{T,A},{G,T},{T,G}};
		
		for(int k = 0; k < headers.length; k++){
			String[] nextLine = spinToNextRealLine(in).trim().split("\\s+");
			for(int u = 0; u < headers.length; u++){
				double dgv = convertFloat(nextLine[u]);
				nns.add(ex.new NearestNeighborScore(headers[k][0], headers[k][1], headers[u][0], headers[u][1], dgv));
			}
		}
		
		ex.setLoopEnergies(
						convertFloat(spinToNextRealLine(in).trim().split("\\s+")),
						convertFloat(spinToNextRealLine(in).trim().split("\\s+")),
						convertFloat(spinToNextRealLine(in).trim().split("\\s+")),
						convertFloat(spinToNextRealLine(in).trim().split("\\s+"))
						);
		
		if (!in.nextLine().startsWith(">Triloops")){
			throw new RuntimeException("Bad formatted file (>Triloops missing)");
		}
		
		for(int i = 0; i < 2; i++){
			while(in.hasNextLine()){
				String line = in.nextLine();
				if (line.startsWith(">")){
					break;
				}
				//Triloops, then tetraloops
				String[] parts = line.trim().split("\\s+");
				specialHairpins.add(ex.new HairpinLoop(convertBases(parts[0],std), convertFloat(parts[1])));
			}
		}
		
		//Mismatch hairpin
		
		for(int Y1 = 0; Y1 < bases.length; Y1++){
			for(int Y2 = 0; Y2 < bases.length; Y2++){
				String[] line = spinToNextRealLine(in).trim().split("\\s+");
				for(int lastMatch = 0; lastMatch < headers.length; lastMatch++){
					//A mismatch of Y1xY2 following the last pair of the stem, lastMatch
					double dgv = convertFloat(line[lastMatch]);
					misHairpin.add(ex.new TerminalMismatchPairScore(headers[lastMatch][0],headers[lastMatch][1],bases[Y1],bases[Y2],dgv));
				}
			}
		}
		
		//Mismatch interior
		
		for(int Y1 = 0; Y1 < bases.length; Y1++){
			for(int Y2 = 0; Y2 < bases.length; Y2++){
				String[] line = spinToNextRealLine(in).trim().split("\\s+");
				for(int lastMatch = 0; lastMatch < headers.length; lastMatch++){
					double dgv = convertFloat(line[lastMatch]);
					misInnerLoop.add(ex.new TerminalMismatchPairScore(headers[lastMatch][0],headers[lastMatch][1],bases[Y1],bases[Y2],dgv));
				}
			}
		}
		
		// 6x top dangle
		for(int k = 0; k < 6; k++){
			String[] line = spinToNextRealLine(in).trim().split("\\s+");
			for(int col = 0; col < bases.length; col++){
				double dgv = convertFloat(line[col]);
				topDangle.add(ex.new DangleScore(headers[k][0], headers[k][1], bases[col], dgv));				
			}
		}
		
		// 6x bottom dangle
		for(int k = 0; k < 6; k++){
			String[] line = spinToNextRealLine(in).trim().split("\\s+");
			for(int col = 0; col < bases.length; col++){
				double dgv = convertFloat(line[col]);
				bottomDangle.add(ex.new DangleScore(headers[k][0], headers[k][1], bases[col], dgv));				
			}
		}

		ex.setMultibranchTerms(
						convertFloat(spinToNextRealLine(in).trim().split("\\s+"))
						);
		
		ex.setATPenalty(convertFloat(spinToNextRealLine(in).trim()));
		
		// Interior loops, L1 = 1, L2 = 1
		for(int leftPair = 0; leftPair < headers.length; leftPair++){
			for(int rightPair = 0; rightPair < headers.length; rightPair++){
				spinToNextRealLine(in); // I already know this.
				for(int x1 = 0; x1 < bases.length; x1++){
					String[] line = spinToNextRealLine(in).trim().split("\\s+");
					for(int x2 = 0; x2 < bases.length; x2++){
						specialInteriorLoops.add(ex.new InteriorLoop(new int[]{
							headers[leftPair][0], bases[x1], headers[rightPair][0]
						}, new int[]{
							headers[leftPair][1], bases[x2], headers[rightPair][1]
						}, convertFloat(line[x2])));
					}
				}
			}
		}
		
		// Interior loops, L1 = 2, L2 = 2
		for(int leftPair = 0; leftPair < headers.length; leftPair++){
			for(int rightPair = 0; rightPair < headers.length; rightPair++){
				for(int x1 = 0; x1 < bases.length; x1++){
					for(int y1 = 0; y1 < bases.length; y1++){
						spinToNextRealLine(in); // I already know this.
						//Ok, now we have all combinations of the two base interior loops.
						for(int y2 = 0; y2 < bases.length; y2++){
							String[] line = spinToNextRealLine(in).trim().split("\\s+");
							for(int x2 = 0; x2 < bases.length; x2++){
								specialInteriorLoops.add(ex.new InteriorLoop(new int[]{
										headers[leftPair][0], bases[x1], bases[y1], headers[rightPair][0]
									}, new int[]{
										headers[leftPair][1], bases[x2], bases[y2], headers[rightPair][1]
									}, convertFloat(line[x2])));
							}
						}
					}
				}
			}
		}
		
		// Interior Loops, L1 = 1, L2 = 2
		for(int leftPair = 0; leftPair < headers.length; leftPair++){
			for(int rightPair = 0; rightPair < headers.length; rightPair++){
				for(int x1 = 0; x1 < bases.length; x1++){
					spinToNextRealLine(in);
					//Ok, now we have all combinations of the two bottom bases
					for(int y2 = 0; y2 < bases.length; y2++){
						String[] line = spinToNextRealLine(in).trim().split("\\s+");
						for(int x2 = 0; x2 < bases.length; x2++){
							specialInteriorLoops.add(ex.new InteriorLoop(new int[]{
									headers[leftPair][0], bases[x1], headers[rightPair][0]
								}, new int[]{
									headers[leftPair][1], bases[x2], bases[y2], headers[rightPair][1]
								}, convertFloat(line[x2])));
						}
					}
				}
			}
		}
		
		String polyC = spinToNextRealLine(in);
		String beta = spinToNextRealLine(in);
		String bimolecular = spinToNextRealLine(in);
		
		ex.setSequenceSpecificStructures(nns, misHairpin, topDangle, bottomDangle, specialHairpins, specialInteriorLoops);
	}
	private static int[] convertBases(String string, CircDesigNAConfig std) {
		int[] toRet = new int[string.length()];
		for(int i = 0; i < string.length(); i++){
			toRet[i] = std.monomer.decodeBaseChar(string.charAt(i));
		}
		return toRet;
	}
	private static double convertFloat(String string) {
		//Units of dekacals per mol
		return new Double(string)/100;
	}
	private static double[] convertFloat(String[] string) {
		//Units of dekacals per mol
		double[] toRet = new double[string.length];
		for(int i = 0; i < toRet.length; i++){
			toRet[i] = convertFloat(string[i]);
		}
		return toRet;
	}
	private static String spinToNextRealLine(Scanner in) {
		while(in.hasNextLine()){
			String line = in.nextLine();
			if (line.startsWith(">")){
				continue;
			}
			return line;
		}
		return null;
	}

}
