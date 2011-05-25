package edu.utexas.cssb.circdesigna;

import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static edu.utexas.cssb.circdesigna.DomainSequence.DNA_COMPLEMENT_FLAG;
import static edu.utexas.cssb.circdesigna.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.utexas.cssb.circdesigna.AbstractDomainDesignTarget.HairpinClosingTarget;
import edu.utexas.cssb.circdesigna.DomainDesigner.ScorePenalty;
import edu.utexas.cssb.circdesigna.DomainStructureBNFTree.DomainStructure;
import edu.utexas.cssb.circdesigna.DomainStructureBNFTree.HairpinStem;
import edu.utexas.cssb.circdesigna.DomainStructureBNFTree.SingleStranded;
import edu.utexas.cssb.circdesigna.DomainStructureBNFTree.ThreePFivePOpenJunc;

import DnaDesign.AbstractDesigner.PopulationDesignMember;
import DnaDesign.impl.DomainDesignPMemberImpl;
import DnaDesign.impl.DomainDesignerImpl.HairpinOpening;
import DnaDesign.impl.DomainDesignerImpl.MFEHybridScore;
import DnaDesign.impl.DomainDesignerImpl.SelfFold;
import DnaDesign.impl.DomainDesignerImpl.SelfSimilarityScore;
import DnaDesign.impl.DomainDesignerImpl.VariousSequencePenalties;
import DnaDesign.parser.CDNA2PublicParser;
import DnaDesign.parser.CDNA2Token;
import DnaDesign.parser.CDNA2PublicParser.Result;

/**
 * Utilities for 
 * 1) Walking down the structural tree of molecules looking for certain motifs (used for generating penalties)
 * 2) Cleaning up the output when reporting a solution set (not repeating strands which are equal)
 * 3) "Duplicate System" feature (for analysis purposes)
 */
public class DomainDesigner_SharedUtils {
	/*
	public static void utilJunctionSplitter(List<DomainSequence> out, String toJunctionize, DomainStructureData dsd){
		int[] bases = DomainDesigner_SharedUtils.utilReadSequence(toJunctionize,dsd);
		for(int b = 0; b < bases.length; b++){
			int nexB = b+1;
			if (nexB < bases.length){
				DomainSequence m = new DomainSequence();
				m.setDomains(bases[b],bases[nexB],dsd);
				out.add(m);
				//System.out.println(m);
			}
		}
	}
	*/

	public static int[] utilReadSequence(String rawDomains, DomainDefinitions dsd) {
		rawDomains = "TMP "+rawDomains;

		Result t = CDNA2PublicParser.parse(rawDomains);

		ArrayList<Integer> bases = new ArrayList();
		for(Object token : t.parse){
			if (token instanceof CDNA2Token.Domain){
				String name = ((CDNA2Token.Domain)token).name;
				int newBase = 0;
				if (name.endsWith("*")){
					newBase = dsd.lookupDomainName(name.substring(0,name.length()-1));
					newBase |= DNA_COMPLEMENT_FLAG;
				} else {
					newBase = dsd.lookupDomainName(name);
				}
				if (newBase >= 0){
					bases.add(newBase);
				}
			}
		}
		
		int[] bases2 = new int[bases.size()];
		for(int k = 0; k < bases.size(); k++){
			bases2[k] = bases.get(k);
		}
		return bases2;
	}


	public static void utilRemoveDuplicateSequences(
			ArrayList<HairpinClosingTarget> hairpinClosings) {
		ListIterator<HairpinClosingTarget> itr = hairpinClosings.listIterator();
		big: while(itr.hasNext()){
			HairpinClosingTarget seq = itr.next();
			for(HairpinClosingTarget q : hairpinClosings){
				if (q!=seq && q.equals(seq)){
					for(int i = 0; i < q.stemAndOpening.length; i++){
						q.stemAndOpening[i].appendMoleculeNames(seq.stemAndOpening[i]);	
					}
					itr.remove();
					continue big;
				}
			}
		}
	}
	
	public static void utilRemoveDuplicateOrComplementaryDomains(List<DomainSequence> seqs) {
		ListIterator<DomainSequence> itr = seqs.listIterator();
		big: while(itr.hasNext()){
			DomainSequence seq = itr.next();
			for(DomainSequence q : seqs){
				if (q!=seq && ((q.domainList[0] & DNA_SEQ_FLAGSINVERSE)==(seq.domainList[0] & DNA_SEQ_FLAGSINVERSE))){
					q.appendMoleculeNames(seq);	
					itr.remove();
					continue big;
				}
			}
		}
		
	}
	
	private static boolean equalsSequenceArray(DomainSequence[] q,
			DomainSequence[] seq) {
		if (q.length!=seq.length){
			return false;
		}
		for(int i = 0; i < seq.length; i++){
			if (!q[i].equals(seq[i])){
				return false;
			}
		}
		return true;
	}


	public static void utilRemoveDuplicateSequences(List<DomainSequence> seqToSynthesize) {
		ListIterator<DomainSequence> itr = seqToSynthesize.listIterator();
		big: while(itr.hasNext()){
			DomainSequence seq = itr.next();
			for(DomainSequence q : seqToSynthesize){
				if (q!=seq && q.equals(seq)){
					q.appendMoleculeNames(seq);
					itr.remove();
					continue big;
				}
			}
		}
	}

	/**
	 * Removes subsequences
	 * @param seqs
	 */
	public static void utilRemoveSubsequences(
			List<DomainSequence> seqs) {
		ListIterator<DomainSequence> itr = seqs.listIterator();
		big: while(itr.hasNext()){
			DomainSequence seq = itr.next();
			for(DomainSequence q : seqs){
				if (q!=seq && seq.isSubsequenceOf(q)){
					q.appendMoleculeNames(seq);
					itr.remove();
					continue big;
				}
			}
		}
	}

	
	public static boolean checkComplementary(DomainSequence a, DomainSequence b){
		for(int k = 0; k < a.numDomains; k++){
			for(int y = 0; y < b.numDomains; y++){
				int abase = a.domainList[k];
				int bbase = b.domainList[y];
				if ((abase&DNA_COMPLEMENT_FLAG)!=(bbase&DNA_COMPLEMENT_FLAG)){
					if ((abase&DNA_SEQ_FLAGSINVERSE)==(bbase&DNA_SEQ_FLAGSINVERSE)){
						return true;
					}
				}
			}
		}
		return false;
	}


	/**
	 * Finds domain sequences that are in danger of forming secondary structure.
	 */
	public static void utilSingleStrandedFinder(DomainStructureBNFTree dsd,
			ArrayList<DomainSequence> MustBeVanilla) {
		ArrayList<Integer> freeList = new ArrayList<Integer>();
		
		for(DomainStructure ds : dsd.structures){
			utilSingleStrandedFinder_R(dsd,ds,MustBeVanilla,freeList);
		}
		if (freeList.size()>0){
			DomainSequence freeListD = new DomainSequence();
			freeListD.setDomains(freeList,dsd);
			MustBeVanilla.add(freeListD);
		}
		
		//Check!
		for(DomainSequence q : MustBeVanilla){
			if (checkComplementary(q,q)){
				throw new RuntimeException("Sequence "+q.toString(dsd.getDomainDefs())+" contains internal complementarity.");
			}
		}
	}
	
	
	private static void utilSingleStrandedFinder_R(DomainStructureBNFTree dsd, DomainStructure ds,
			ArrayList<DomainSequence> MustBeVanilla, ArrayList<Integer> freeList) {

		//TODO: free strands need to have "hairpins" in the middle.
		//For now, we just add them all.
		if (ds instanceof SingleStranded){
			for(int k : ds.sequencePartsInvolved){
				freeList.add((dsd.domains[k]));
			}
		}
		if (ds instanceof HairpinStem){
			ArrayList<Integer> subFree = new ArrayList();
			for(DomainStructure g : ds.subStructure){
				utilSingleStrandedFinder_R(dsd,g,MustBeVanilla,subFree);
			}	
			if (subFree.size()>0){
				DomainSequence freeListD = new DomainSequence();
				freeListD.setDomains(subFree,dsd);
				MustBeVanilla.add(freeListD);
			}
		}
		
	}

	/**
	 * Returns both a new domainDefs and a Molecules text.
	 * Basically, it adds a copy of the original problem set, so it has to find the
	 * original problem set in the input, which can already have copies of same.
	 */
	public static String[] duplicateSystem(String domainDef, String molecul) {
		String[] toRet = new String[2];
		String LN = "\n";//System.getProperty("line.separator");
		//First the domainDefs
		int whichV = 0;
		StringBuffer out;
		out = new StringBuffer();
		{
			Scanner domainDefIn = new Scanner(domainDef);
			while(domainDefIn.hasNextLine()){
				String line = domainDefIn.nextLine();
				if (line.trim().length()!=0){
					out.append(line);
					out.append(LN);
				}
			}
			whichV = 1;
			for(boolean actuallyAppend : new boolean[]{false,true}){
				Scanner in = new Scanner(domainDef);
				while(in.hasNextLine()){
					String line = in.nextLine();
					if (line.trim().length()==0){
						continue;
					}
					String[] linefake = line.split("\\s+");
					//Add "vX" to first object.
					int v = linefake[0].lastIndexOf('v');
					int vInd = linefake[0].length();
					int newV = -1;
					if (v!=-1){
						vInd = v;
						try {
							newV = new Integer(linefake[0].substring(v+1));
						} catch (NumberFormatException e){
							//e.printStackTrace();
						}
					} else {
						newV = 1;
					}
					whichV = Math.max(whichV,newV);
					if (newV < whichV){
						continue;
					}
					if (actuallyAppend){
						out.append(linefake[0].substring(0,vInd)+"v"+(whichV+1));
						out.append(line.substring(linefake[0].length()));
						out.append(LN);
					}
				}
				in.close();
			}
		}
		toRet[0] = out.toString();
		out = new StringBuffer();
		{
			Scanner moleculIN = new Scanner(molecul);
			while(moleculIN.hasNextLine()){
				String line = moleculIN.nextLine();
				out.append(line);
				out.append(LN);
			}
			moleculIN.close();
			Scanner in = new Scanner(molecul);
			while(in.hasNextLine()){
				String line = in.nextLine().trim();
				if (line.length()==0){
					continue;
				}
				String[] linefake = line.split("\\s+",2);
				if (linefake[0].endsWith("v"+whichV) || whichV == 1){
					//Translate this guy!
					int subLen = ("v"+whichV).length();
					if (whichV==1){
						subLen = 0;
					}
					out.append(linefake[0].substring(0,linefake[0].length()-subLen)+"v"+(whichV+1)+" ");
					Pattern pat = Pattern.compile("\\w+|\\W+");
					Matcher m = pat.matcher(linefake[1]);
					ArrayList<String> line2S = new ArrayList<String>();
					while(m.find()){
						line2S.add(m.group(0));
					}
					String[] line2 = new String[line2S.size()]; line2S.toArray(line2);
					for(int k = 0; k < line2.length; k++){
						if (line2[k].matches("\\W+")){
							out.append(line2[k]);
						} else {
							if (whichV!=1){
								int v = line2[k].lastIndexOf('v');
								if (v==-1){
									throw new RuntimeException("Invalid duplicated molecule: "+line);
								} else {
									try {
										if (new Integer(line2[k].substring(v+1))>=0){
											if (!line2[k].endsWith("v"+whichV)){
												throw new RuntimeException("Invalid duplicated molecule: "+line);
											}
											line2[k] = line2[k].substring(0,v);
										}
									} catch (NumberFormatException e){
										throw new RuntimeException("Invalid duplicated molecule: "+line);
									}
								}
							}
							out.append(line2[k]+"v"+(whichV+1));
						}
					}
					out.append(LN);
				}
			}
			in.close();
			
		}
		toRet[1] = out.toString();
		return toRet;
	}

	//Algorithms using the Polymer Graph. 
	
	public static void utilSingleStrandedFinder(DomainPolymerGraph dsg, ArrayList<DomainSequence> mustBeVanilla) {
		LinkedList<ArrayList<Integer>> singleStrandedLoops = new LinkedList();
		singleStrandedLoops.push(new ArrayList());
		for(int k = 0; k < dsg.length(); k++){
			int domain = dsg.getDomain(k);
			int pair = dsg.getDomainPair(k);
			if (pair<0){
				singleStrandedLoops.peek().add(domain);
			} else {
				if (pair < k){
					//Loop ending. pop stack.
					ArrayList<Integer> closedCircuit = singleStrandedLoops.pop();
					//Need to assume that all molecules end in a 3prime end...
					addSingleStrandedClosedLoop(closedCircuit,mustBeVanilla, dsg);
				} else {
					//Loop beginning. push stack.
					singleStrandedLoops.push(new ArrayList());
				}
			}
		}
		addSingleStrandedClosedLoop(singleStrandedLoops.pop(),mustBeVanilla, dsg);
		if (!singleStrandedLoops.isEmpty()){
			throw new RuntimeException("Assertion error utilSingleStrandedFinder polymergraph");
		}
	}

	/**
	 * Returns any generalized single stranded regions of the complexes.
	 */
	private static void addSingleStrandedClosedLoop(ArrayList<Integer> closedCircuit,ArrayList<DomainSequence> mustBeVanilla, AbstractComplex dsd) {
		LinkedList<Integer> ThreePrimeEnds = new LinkedList();
		int ct = 0;
		for(Integer q : closedCircuit){
			if (q < 0){
				ThreePrimeEnds.add(ct);
			}
			ct++;
		}
		if (ct > 0){
			if (ThreePrimeEnds.isEmpty()){
				//Just a single, circular domainsequence
				DomainSequence made = new DomainSequence();
				made.setDomains(closedCircuit, dsd);
				made.setCircular(true);
				mustBeVanilla.add(made);
			} else {
				//Start at each of the 3 prime ends, and wrap around to the next.
				int firstIndex = ThreePrimeEnds.getFirst();
				int endIndex = ((firstIndex+1)%closedCircuit.size());
				ArrayList<Integer> subArray = new ArrayList();
				boolean wrapAroundCheck = false;
				for(int t = endIndex; ; t = (t+1)%closedCircuit.size()){
					if (t==endIndex){
						if (wrapAroundCheck){
							break;
						}
						wrapAroundCheck = true;
					}
					int domain = closedCircuit.get(t);
					if (domain < 0){
						if (!subArray.isEmpty()){
							DomainSequence newSeq = new DomainSequence();
							newSeq.setDomains(subArray, dsd);
							mustBeVanilla.add(newSeq);
							subArray.clear();
						}
					} else {
						subArray.add(domain);
					}
				}
				if (!subArray.isEmpty()){
					throw new RuntimeException("Assertion error. addSingleStrandedClosedLoop");
				}
			}
		}
	}

	/**
	 * Finds hairpin stems. Not extraordinarily useful.
	 */
	public static void utilHairpinInternalsFinder(DomainPolymerGraph dsg, ArrayList<DomainSequence> hairpinInnards) {
		ArrayList<Integer> runningHairpin = new ArrayList();
		for(int k = 0; k < dsg.length(); k++){
			int domain = dsg.getDomain(k);
			int pair = dsg.getDomainPair(k);
			if (pair >= 0){
				runningHairpin.add(domain);
			}
			//Condition to terminate running hairpin
			if (pair < 0 || !(k+1 < dsg.length())){
				if (!runningHairpin.isEmpty()){
					DomainSequence pairedSequence = new DomainSequence();
					pairedSequence.setDomains(runningHairpin,dsg);
					hairpinInnards.add(pairedSequence);
					runningHairpin.clear();
				}
			}
		}
	}

	public static void utilHairpinClosingFinder(AbstractDomainDesignTarget target, DomainPolymerGraph dsg, ArrayList<HairpinClosingTarget> hairpinLoops) {
		//Only add the hairpin closings on the closing of the hairpin (so when the later partner domain is encountered)
		for(int k = 0; k < dsg.length(); k++){
			int domain = dsg.getDomain(k);
			if (domain<0){
				continue; //Only care about hairpin domains.
			}
			int pair = dsg.getDomainPair(k);
			if (pair >= 0 && pair < k){ //only care about hairpin closing-domains
				int pairDomain = dsg.getDomain(pair);
				if (pair > 0 && (k+1)<dsg.length()){
					//Are both neighbors valid domains?
					if (dsg.getDomain(pair-1)>=0 && dsg.getDomain(k+1) >= 0){
						//Ok, can (pair-1) or (k+1) pair?
						if (canPair(dsg,pair-1) || canPair(dsg,k+1)
								&& !areComplementary(dsg.getDomain(pair-1),dsg.getDomain(k+1))){
							hairpinLoops.add(target.new HairpinClosingTarget(dsg.getDomain(pair-1),
									pairDomain,
									domain,
									dsg.getDomain(k+1),
									true,dsg));
						}
					}
				}
				if (k > 0 && (pair+1)<dsg.length()){
					//Are both neighbors valid domains?
					if (dsg.getDomain(pair+1)>=0 && dsg.getDomain(k-1) >= 0){
						//How about the "inside": (pair+1) or (k-1).
						if (canPair(dsg,pair+1) || canPair(dsg,k-1)
								&& !areComplementary(dsg.getDomain(pair+1),dsg.getDomain(k-1))){
							hairpinLoops.add(target.new HairpinClosingTarget(
									pairDomain,
									dsg.getDomain(pair+1),
									dsg.getDomain(k-1),
									domain,
									false,dsg));
						}
					}
				}
			}
		}
	}

	private static boolean canPair(DomainPolymerGraph dsg, int i) {
		//A valid domain and not paired.
		return ( dsg.getDomain(i)>=0 ) && (dsg.getDomainPair(i) < 0);
	}

	
	/**
	 * Returns a 1 if the loop (i,j) is a loop connecting a base of a domain to the corresponding base in its reverse complement, and 0 otherwise. 
	 */
	public static boolean isAlignedAndShouldPair(DomainSequence ds1,
			int i, DomainSequence ds2, int j, int[][] domains) {
		int domain1 = ds1.domainAt(i, domains);
		int domain2 = ds2.domainAt(j, domains);
		if (areComplementary(domain1,domain2)){
			int domain = domain1 & DNA_SEQ_FLAGSINVERSE;
			
			int offset1 = ds1.offsetInto(i, domains, false); //Don't uncomplement - give offset 
			int offset2 = ds2.offsetInto(j, domains, false);
			if (offset1 == (domains[domain].length - 1 - offset2)){
				return true;
			}
		}
		return false;
	}
	
	private static boolean areComplementary(int domain1, int domain2) {
		return (domain1 ^ DNA_COMPLEMENT_FLAG)==domain2;
	}

	public static void utilPairsOfDomainsFinder(DomainPolymerGraph dsg, ArrayList<DomainSequence> pairsOfDomains) {
		for(int k = 0; k < dsg.length()-1; k++){
			int domain1 = dsg.getDomain(k);
			int domain2 = dsg.getDomain(k+1);
			if (domain1 >=0 && domain2 >= 0){
				DomainSequence newSeq = new DomainSequence();
				newSeq.setDomains(domain1, domain2, dsg);
				pairsOfDomains.add(newSeq);
			}
		}
	}

	public static void utilSingleDomainsWithOverlap(DomainPolymerGraph dsg,
			ArrayList<DomainSequence> singleDomainsWithOverlap,
			int overlap_length) {
		for(int k = 0; k < dsg.length(); k++){
			int domain1 = dsg.getDomain(k);
			if (domain1<0) continue; //End of strand
			int leftOverlap = overlap_length;
			int rightOverlap = overlap_length;
			LinkedList<Integer> overlapSequence = new LinkedList<Integer>();
			for (int domainLeft = k-1; domainLeft>=0 && leftOverlap>0; domainLeft--){
				int lDomain = dsg.getDomain(domainLeft);
				if (lDomain<0) break; //End of strand
				overlapSequence.addFirst(lDomain);
				leftOverlap -= dsg.getDomainDefs().domainLengths[lDomain & DNA_SEQ_FLAGSINVERSE];
			}
			overlapSequence.add(domain1);
			for(int domainRight = k+1; domainRight<dsg.length() && rightOverlap>0; domainRight++){
				int rDomain = dsg.getDomain(domainRight);
				if (rDomain<0) break; //End of strand
				overlapSequence.add(rDomain);
				rightOverlap -= dsg.getDomainDefs().domainLengths[rDomain & DNA_SEQ_FLAGSINVERSE];
			}
			DomainSequence toAdd = new DomainSequence();
			toAdd.setDomains(overlapSequence, dsg);
			singleDomainsWithOverlap.add(toAdd);
		}
	}

	public static void utilSingleDomainsFinder(DomainPolymerGraph dsg,
			ArrayList<DomainSequence> singleDomains) {
		for(int k = 0; k < dsg.length(); k++){
			int domain1 = dsg.getDomain(k);
			if (domain1<0) continue; //End of strand
			DomainSequence toAdd = new DomainSequence();
			toAdd.setDomains(domain1, dsg);
			singleDomains.add(toAdd);
		}
	}
	
	public static DesignScoreBreakdown getScoreBreakdown(int[][] domain, AbstractDomainDesignTarget designTarget, List<ScorePenalty> penalties) {
		DesignScoreBreakdown toRet = new DesignScoreBreakdown();

		toRet.breathingHelixes = 0;
		for(int i = 0; i < domain.length; i++){
			if (domain[i].length==0){
				continue;
			}
			for(int termbase : new int[]{domain[i][0], domain[i][domain[i].length-1]}){
				switch(designTarget.Std.monomer.noFlags(termbase)){
				case G: case C:
					break;
				default:
					toRet.breathingHelixes++;
				}
			}
		}
		
		toRet.crossInteractionsOnly = 0;
		toRet.netScore = 0;
		toRet.selfFoldOnly = 0;
		toRet.selfSimilarity = 0;
		toRet.homopolymer = 0;
		for(ScorePenalty p : penalties){
			if (p instanceof MFEHybridScore || p instanceof HairpinOpening){
				toRet.crossInteractionsOnly += p.old_score;
			} else
			if (p instanceof SelfFold){
				toRet.selfFoldOnly += p.old_score;
			} else
			if (p instanceof SelfSimilarityScore){
				toRet.selfSimilarity += p.old_score;
			} else
			if (p instanceof VariousSequencePenalties){
				toRet.homopolymer += p.old_score;
			}
			toRet.netScore += p.old_score;
		}
		return toRet;
	}


	public static DomainDesignPMemberImpl getBestMember(PopulationDesignMember<DomainDesignPMemberImpl>[] population) {
		double bestScore = Double.MAX_VALUE;
		int bestIndex = 0;
		for(int i = 0; i < population.length; i++){
			double score = 0;
			DomainDesignPMemberImpl pdm = (DomainDesignPMemberImpl)population[i];
			for(ScorePenalty d : pdm.penalties){
				score += d.old_score;
			}
			if (score < bestScore){
				bestScore = score;
				bestIndex = i;
			}
		}
		return (DomainDesignPMemberImpl)population[bestIndex];
	}

	/**
	 * Gets the farthest away (by hamming distance) member in the population from the best solution, best.
	 */
	public static DomainDesignPMemberImpl getFarthestFrom(DomainDesignPMemberImpl best, PopulationDesignMember<DomainDesignPMemberImpl>[] population) {
		double bestDistance = 0;
		int bestIndex = 0;
		for(int i = 0; i < population.length; i++){
			double distance = getHammingDistance(best,(DomainDesignPMemberImpl)population[i]);
			if (distance > bestDistance){
				bestDistance = distance;
				bestIndex = i;
			}
		}
		return (DomainDesignPMemberImpl)population[bestIndex];
	}


	private static double getHammingDistance(DomainDesignPMemberImpl best, DomainDesignPMemberImpl other) {
		int H = 0;
		for(int j = 0; j < best.domain.length; j++){
			for(int i = 0; i < best.domain[j].length; i++){
				if (best.domain[j][i]!=other.domain[j][i]){
					H++;
				}
			}
		}
		return H;
	}
}
