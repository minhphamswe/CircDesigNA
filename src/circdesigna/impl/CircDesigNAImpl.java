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
package circdesigna.impl;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import circdesigna.AbstractDomainDesignTarget;
import circdesigna.AbstractDomainDesignTarget.DuplexClosingTarget;
import circdesigna.CircDesigNA;
import circdesigna.CircDesigNAOptions;
import circdesigna.CircDesigNA_SharedUtils;
import circdesigna.DesignIntermediateReporter;
import circdesigna.DomainDefinitions;
import circdesigna.DomainSequence;
import circdesigna.SequencePenalties;
import circdesigna.abstractDesigner.ParetoSort;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.ConstraintsNAFolding;
import circdesigna.energy.ConstraintsNAFoldingImpl;

/**
 * Implementation of CircDesigNA
 */
public class CircDesigNAImpl extends CircDesigNA{
	private final ConstraintsNAFolding flI;
	private final SequencePenalties sp;
	
	/**
	 * @param std 
	 * @param foldingImpl<br>
	 * The folding score functions to utilize in evaluating solution candidates.
	 * @param designSSonly<br>
	 * 4 kinds of score functions are currently implemented: Validity, SingleStrandedAssurance, Crosstalk, Dimerization.
	 * If designSSonly is true, only SingleStrandedAssurance will be used; As a result, crosstalk and dimerization may occur
	 * in solution candidates.
	 */
	public CircDesigNAImpl(ConstraintsNAFolding foldingImpl, SequencePenalties sp, CircDesigNAConfig std) {
		super(std);
		this.flI = foldingImpl;
		this.sp = sp;
	}

	public void setPhase(int phase){
		if (flI instanceof ConstraintsNAFoldingImpl){
			((ConstraintsNAFoldingImpl) flI).setScoringModel(phase);
		}
	}
	public int countPhases() {
		if (flI instanceof ConstraintsNAFoldingImpl){
			return 3;
		}
		return 1;
	}
	
	/**
	 * Returns a list of score penalties to use in evaluating population members.
	 */
	public List<ScorePenalty> listPenalties(
			AbstractDomainDesignTarget designTarget,
			DesignIntermediateReporter DIR, int[][] domain2, CircDesigNAOptions options, DomainDefinitions dsd) {

		ArrayList<DomainSequence> rawStrands = designTarget.wholeStrands;
		ArrayList<DomainSequence> makeSS = new ArrayList();
		makeSS.addAll(designTarget.generalizedSingleStranded);
		//makeSS.addAll(designTarget.singleDomains);
		ArrayList<DuplexClosingTarget> duplexClosings = designTarget.duplexClosings;

		CircDesigNA_SharedUtils.utilRemoveDuplicateSequences(rawStrands);
		CircDesigNA_SharedUtils.utilRemoveDuplicateSequences(makeSS);

		CircDesigNA_SharedUtils.utilRemoveDuplicateSequences(duplexClosings);

		List<ScorePenalty> allScores = new LinkedList<ScorePenalty>();
		//Only add this penalty if the system contains at least one base.
		int totalBases = 0;
		for(int[] row : domain2){
			totalBases += row.length;
		}
		if (totalBases>0){
			allScores.add(new VariousSequencePenalties(rawStrands,DIR));	
		}

		ArrayList<MFEHybridNonlegalScore> hybridScorings = new ArrayList<MFEHybridNonlegalScore>();
		for(int i = 0; i < makeSS.size(); i++){
			DomainSequence ds = makeSS.get(i);
			//Secondary Structure Formation
			allScores.add(new SelfFoldNonlegalScore(ds, DIR));
			hybridScorings.add(new MFEHybridNonlegalScore(ds, ds, DIR, false));

			//Hybridization
			for(int k = i+1; k < makeSS.size(); k++){ //Do only upper triangle
				DomainSequence ds2 = makeSS.get(k);
				boolean sameMol = ds.getMoleculeName().equals(ds2.getMoleculeName());
				hybridScorings.add(new MFEHybridNonlegalScore(ds, ds2, DIR, sameMol));
			}
		}

		for(int i = 0; i < hybridScorings.size(); i++){
			MFEHybridNonlegalScore target = hybridScorings.get(i);
			for(int j = i+1; j < hybridScorings.size(); j++){
				MFEHybridNonlegalScore match = hybridScorings.get(j);
				if (match.ds[0].isSubsequenceOf(target.ds[0]) && match.ds[1].isSubsequenceOf(target.ds[1])){
					hybridScorings.remove(j--);
				}
			}
		}

		allScores.addAll(hybridScorings);

		for(DuplexClosingTarget hairpin : duplexClosings){
			allScores.add(new DuplexOpening(hairpin, DIR));
		}
		
		/*
		if (options.selfSimilarityPenalty.getState() >= 0){
			//Only do each domain once.
			List<DomainSequence> domainsOnceApiece = new ArrayList();
			domainsOnceApiece.addAll(designTarget.singleDomains);
			CircDesigNA_SharedUtils.utilRemoveDuplicateOrComplementaryDomains(domainsOnceApiece);
			for(DomainSequence domain : domainsOnceApiece){
				if (domain2[domain.domainList[0] & DomainSequence.NA_COMPLEMENT_FLAGINV].length >= options.selfSimilarityPenalty.getState()){
					allScores.add(new SelfSimilarityScore(domain, DIR));
				}
			}
		}
		*/

		return allScores;
	}
	
	/**
	 * Implements the "isDominatedBy" method, required for doing pareto sorting of design population members.
	 */
	public static class DParetoSort extends ParetoSort<CircDesigNAPMemberImpl> {
		double[] scores= new double[4];
		public boolean isDominatedBy(CircDesigNAPMemberImpl t, CircDesigNAPMemberImpl t2) {
			Arrays.fill(scores,0);
			for(int i = 0; i < t.penalties.size(); i++){
				int index = 0;
				ScorePenalty sp = t.penalties.get(i);
				if (sp instanceof MFEHybridNonlegalScore || sp instanceof SelfFoldNonlegalScore){
					index = 0;
				} else {
					index = 1;
				}
				/*
				if (sp instanceof SelfSimilarityScore){
					index = 0;
				} else if (sp instanceof MFEHybridScore || sp instanceof SelfFold){
					index = 1;
				} else {
					index = 2;
				}
				*/
				scores[index] += t.penalties.get(i).cur_score - t2.penalties.get(i).cur_score;
			}
			for(int i = 0; i < scores.length; i++){
				if (scores[i] >= 0){
					return false;
				}
			}
			return true;
		}
	}

	public class MFEHybridNonlegalScore extends ScorePenalty {
		//True for intermolecular interactions
		private boolean entropicPenalty = false;
		public MFEHybridNonlegalScore(DomainSequence ds, DomainSequence ds2, DesignIntermediateReporter dir, boolean sameMolecule){
			super(dir);
			this.ds = new DomainSequence[]{ds,ds2};
			chooseScore(dir);
			entropicPenalty = !sameMolecule;
		}		
		public int getPriority(){
			return 2;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double BIMOLECULAR = options.bimolecularPenalty.getState();
			double deltaG = (flI.mfe(ds[0],ds[1],domain,domain_markings,true))+(entropicPenalty?BIMOLECULAR:0);
			//int longestHelixLength = flI.getLongestHelixLength();
			//int numBasesPaired = flI.getNumBasesPaired();
			//double normal = longestHelixLength*numBasesPaired;
			//Compute the length of the longest helix found.
			//normal = normal*Math.max(0,deltaG);
			return Math.max(-deltaG,0);
		}
		public boolean isIntramolecular(){
			return !entropicPenalty;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public ScorePenalty clone() {
			MFEHybridNonlegalScore ci = new MFEHybridNonlegalScore(ds[0],ds[1],dir,!entropicPenalty);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
	}
	
	/**
	 * Penalize complementarity at the base of a hairpin loop
	 */
	public class DuplexOpening extends ScorePenalty {
		private DuplexClosingTarget hairpin;
		private int markLeft = -1, markRight, jOffset;
		public DuplexOpening(DuplexClosingTarget hairpin, DesignIntermediateReporter dir){
			super(dir);
			this.ds = hairpin.stemAndOpening;
			this.hairpin = hairpin;
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			DuplexOpening ci = new DuplexOpening(hairpin,dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		public int getPriority(){
			return 0;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			//Prevent the score evaluator from marking bases in the "stem" region of this helix
			if (markLeft==-1){
				int start = 0;
				int end = hairpin.stemAndOpening[0].length(domain);
				
				//TODO fix this mess
				int middle = hairpin.stemOnly[0].length(domain);
				markLeft = middle;
				markRight = end;
				//Aligned at the right, so offset j by the difference
				jOffset = hairpin.stemAndOpening[0].length(domain)-hairpin.stemAndOpening[1].length(domain);
			}
			double StemAndOpeningScore =flI.mfeStraight(hairpin.stemAndOpening[0],hairpin.stemAndOpening[1],domain,domain_markings,markLeft,markRight,jOffset); 
			double OnlyStem =flI.mfeStraight(hairpin.stemOnly[0],hairpin.stemOnly[1],domain,domain_markings,0,0,0); 
			double deltaDeltaG = StemAndOpeningScore - OnlyStem;
			/*
			if (deltaDeltaG > 0){
				System.err.println("Stem+Opening increased delta g  (?)"+" "+StemAndOpeningScore+" "+OnlyStem);
				System.err.println(hairpin.toString(domain));
			}
			*/
			return Math.max(0,-deltaDeltaG);
		}
		public boolean affectedBy(int domain) {
			for(int i = 0; i < ds.length; i++){
				if (ds[i].contains(domain)){
					return true;
				}
			}
			return false;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public int getNumDomainsInvolved() {
			return 2;
		}
	}
	
	/**
	 * Calculates the Expected Minimum Free Energy for a random string of length N, and subtract that
	 * from this measure. Floor to zero to catch exceptional cases.
	 */
	public class SelfSimilarityScore extends ScorePenalty {
		public SelfSimilarityScore(DomainSequence domain, DesignIntermediateReporter dir){
			super(dir);
			DomainSequence revComp = new DomainSequence();
			domain.makeReverseComplement(revComp);
			ds = new DomainSequence[]{domain, revComp};
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			SelfSimilarityScore ci = new SelfSimilarityScore(ds[0],dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		public int getPriority(){
			return 2;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double deltaG = flI.mfeNoDiag(ds[0], ds[1], domain, domain_markings);
			return Math.max(0,-deltaG);
		}
		public boolean affectedBy(int domain) {
			return ds[0].contains(domain);
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
	}
	
	/**
	 *  A grab-bag for miscellaneous sequence constraints. Only the sequences in the provided list
	 * are analyzed. For example, penalties for unwanted sequences (GGGG and ATATAT) may be covered here.
	 * This score is not required to be in any unit, it should be used as a nonrigorous design penalty. However, it should be
	 * used, because many of the routines in this library use parameters derived from highly sequences with good G/C/A/T mixes. An
	 * inability to evaluate whether a design is good is a sign that another solution should be sought after.
	 */
	public class VariousSequencePenalties extends ScorePenalty {
		private List<DomainSequence> seqs;

		public VariousSequencePenalties(List<DomainSequence> seqToSynthesize, DesignIntermediateReporter dir) {
			super(dir);
			seqs = seqToSynthesize;
			//chooseScore(dir); Can't show up. Uses everybody.
		}
		
		public ScorePenalty clone() {
			VariousSequencePenalties ci = new VariousSequencePenalties(seqs,dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		
		public boolean affectedBy(int domain) {
			return true;
		}


		public double evalScoreSub(int[][] domain, int[][] domain_markings) {
			double sum = 0;
			for(DomainSequence seq : seqs){
				sum += sp.getSynthesizabilityScore(seq, domain,domain_markings);
			}
			return sum;
		}

		public int getNumDomainsInvolved() {
			return 1;
		}

		public int getPriority() {
			return 0;
		}

		public DomainSequence[] getSeqs() {
			return new DomainSequence[0];
		}
	}
	
	public class SelfFoldNonlegalScore extends ScorePenalty { 
		public SelfFoldNonlegalScore(DomainSequence ds, DesignIntermediateReporter dir){
			super(dir);
			this.ds = new DomainSequence[]{ds};
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			SelfFoldNonlegalScore ci = new SelfFoldNonlegalScore(ds[0],dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double deltaG = (flI.mfe(ds[0],domain,domain_markings,true));
			return Math.max(0,-deltaG);
			//int longestHelixLength = flI.getLongestHelixLength();
			//int numBasesPaired = flI.getNumBasesPaired();
			//double normal = longestHelixLength*numBasesPaired;
			//return normal*Math.max(0,-deltaG);
		}
		public int getPriority(){
			return 1;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
	}
	/*
	public class LocalDefectSSScore extends ScorePenalty {
		private DomainSequence[] ds;
		private AbstractDomainDesignTarget target;
		private double[][] probBuffer;
		public LocalDefectSSScore(DomainSequence ds,
				DesignIntermediateReporter dir, AbstractDomainDesignTarget designTarget) {
			super(dir);
			this.ds = new DomainSequence[]{ds};
			this.target = designTarget;
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			LocalDefectSSScore ci = new LocalDefectSSScore(ds[0],dir,target);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		private double[][] expandCapacity(double[][] probBuffer2, int length1,
				int length2) {
			length2++; //One extra column for unpaired
			if (probBuffer2!=null){
				if (probBuffer2.length >= length1 && probBuffer2[0].length >= length2){
					return probBuffer2;
				}
			}
			return new double[length1][length2];
		}
		public double evalScoreSub(int[][] domain, int[][] domain_markings) {
			int length1 = ds[0].length(domain);
			probBuffer = expandCapacity(probBuffer,length1,length1);
			flI.pairPr(probBuffer, ds[0], domain);
			
			double expectedMisPairedBases = 0;
			for(int i = 0; i < length1; i++){
				for(int j = i+1; j < length1; j++){
					if (CircDesigNA_SharedUtils.isAlignedAndShouldPair(ds[0], i, ds[0], j, domain)){
						expectedMisPairedBases += probBuffer[i][j];
						if (probBuffer[i][j] > .001){
							ds[0].mark(i, domain, domain_markings);
							ds[0].mark(j, domain, domain_markings);
						}
					}
				}
			}
			
			return expectedMisPairedBases;
		}
		public int getPriority() {
			return 2;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
	}
	public class PairDefectScore extends ScorePenalty {
		private DomainSequence[] ds;
		private double[][] probBuffer;
		private boolean entropicBonus;
		private AbstractDomainDesignTarget target;
		public PairDefectScore(DomainSequence ds2, DomainSequence ds3,DesignIntermediateReporter dir, boolean onSameMolecule, AbstractDomainDesignTarget designTarget) {
			super(dir);
			ds = new DomainSequence[]{ds2,ds3};
			chooseScore(dir);
			entropicBonus = onSameMolecule;
			this.target = designTarget;
		}
		public ScorePenalty clone() {
			PairDefectScore ci = new PairDefectScore(ds[0],ds[1],dir,entropicBonus,target);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		public double evalScoreSub(int[][] domain, int[][] domain_markings) {
			int length1 = ds[0].length(domain);
			int length2 = ds[1].length(domain);
			int N = length1+length2;
			probBuffer = expandCapacity(probBuffer,N,N+1);
			flI.pairPr(probBuffer, ds[0], ds[1], domain);
			
			double expectedMisPairedBases = 0;
			for(int i = 0; i < N; i++){
				for(int j = i+1; j < N; j++){
					DomainSequence sI = (i < length1?ds[0] : ds[1]);
					int iInSi = (i < length1?i : i - length1);
					DomainSequence sJ = (j < length1?ds[0] : ds[1]);
					int jInSJ = (j < length1?j : j - length1);
					if (!CircDesigNA_SharedUtils.isAlignedAndShouldPair(sI, iInSi, sJ, jInSJ, domain)){
						expectedMisPairedBases += probBuffer[i][j];
						if (probBuffer[i][j] > .001){
							sI.mark(iInSi, domain, domain_markings);
							sJ.mark(jInSJ, domain, domain_markings);
						}
					}
				}
			}
			
			return expectedMisPairedBases;
		}
		private double[][] expandCapacity(double[][] probBuffer2, int length1,
				int length2) {
			length2++; //One extra column for unpaired
			if (probBuffer2!=null){
				if (probBuffer2.length >= length1 && probBuffer2[0].length >= length2){
					return probBuffer2;
				}
			}
			return new double[length1][length2];
		}
		public int getPriority() {
			return 2;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
	}
	*/
}
