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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

/**
 * Used to restrict the possible sequence designs by banning certain DNA words, 
 * for example, one could ensure that the designed sequence does not contain
 * a certain restriction site.
 * 
 * This class implements a routine which checks, efficiently, whether there is a banned
 * word involving position j in a domain.
 */
public class BannedPatterns extends CircDesigNASystemElement{
	public BannedPatterns(CircDesigNAConfig system) {
		super(system);
	}
	public BannedPatterns(String text, CircDesigNAConfig system) {
		super(system);
		Scanner in = new Scanner(text);
		while(in.hasNextLine()){
			String line = in.nextLine().trim();
			if (line.equals("")){
				continue;
			}
			try {
				String[] line2 = line.split("\\s+",2);
				int[] bases = new int[line2[0].length()];
				for(int k = 0; k < bases.length; k++){
					bases[k] = Std.monomer.decodeConstraintChar(line2[0].charAt(k));
				}
				addBannedWord(bases, new Double(line2[1]));
			} catch (Throwable e){
				throw new RuntimeException("Correct format is <pattern> <penalty>. ("+line+", "+e.getMessage()+")");
			}
		}
	}
	private ArrayList<int[]> bannedWords = new ArrayList();
	private ArrayList<Double> bannedWordPenalty = new ArrayList();
	public void addBannedWord(int[] bases, Double penalty) {
		bannedWords.add(bases);
		bannedWordPenalty.add(penalty);
	}
	/**
	 * Returns true if the pattern'th pattern matches the sequence beginning at position k
	 * of seq. If a match was found, the implicated bases are marked.
	 * @param domain_markings 
	 */
	public boolean matches(int pattern, DomainSequence seq, int[][] domain, int[][] domain_markings, int k) {
		int length = seq.length(domain);
		int[] word = bannedWords.get(pattern);
		if (word.length + k > length){
			return false; //Word is too long to fit.
		}
		for(int u = 0; u < word.length; u++){
			int base = seq.base(k+u, domain, Std.monomer);
			if (!Std.monomer.allowBase(word[u], base)){
				return false;
			}
		}
		//Mark the implicated bases
		seq.mark(k, word.length, domain, domain_markings);
		
		return true;	
	}
	public String defaultBannedWords() {
		return "WWWWWW 16\nSSSSSS 16\nGGGG 16\nCCCC 16\n";
	}
	public int patternSize() {
		return bannedWords.size();
	}
	/**
	 * Gets the penalty associated with the pattern'th pattern.
	 */
	public double patternWeight(int pattern) {
		return bannedWordPenalty.get(pattern);
	}
}
