package DnaDesign.impl;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import DnaDesign.AbstractDomainDesignTarget;
import DnaDesign.DesignIntermediateReporter;
import DnaDesign.DesignerOptions;
import DnaDesign.DomainDesigner;
import DnaDesign.DomainDesigner_SharedUtils;
import DnaDesign.DomainSequence;
import DnaDesign.NAFolding;
import DnaDesign.AbstractDomainDesignTarget.HairpinClosingTarget;
import DnaDesign.Config.CircDesigNAConfig;

/**
 * Implementation of DomainDesigner
 */
public class DomainDesignerImpl extends DomainDesigner{
	
	private NAFolding flI;
	/**
	 * @param std 
	 * @param foldingImpl<br>
	 * The folding score functions to utilize in evaluating solution candidates.
	 * @param designSSonly<br>
	 * 4 kinds of score functions are currently implemented: Validity, SingleStrandedAssurance, Crosstalk, Dimerization.
	 * If designSSonly is true, only SingleStrandedAssurance will be used; As a result, crosstalk and dimerization may occur
	 * in solution candidates.
	 */
	public DomainDesignerImpl(NAFolding foldingImpl, CircDesigNAConfig std) {
		super(std);
		this.flI = foldingImpl;
	}

	public class MFEHybridScore extends ScorePenalty { 
		private boolean entropicPenalty = false;
		public MFEHybridScore(DomainSequence ds, DomainSequence ds2, DesignIntermediateReporter dir, boolean sameMolecule){
			super(dir);
			this.ds = new DomainSequence[]{ds,ds2};
			chooseScore(dir);
			entropicPenalty = !sameMolecule;
		}		
		public int getPriority(){
			return 2;
		}
		private DomainSequence[] ds;
		//This seems too low.
		private final double BIMOLECULAR = -.513;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double deltaG = (flI.mfeHybridDeltaG(ds[0],ds[1],domain,domain_markings))+(entropicPenalty?-BIMOLECULAR:0);
			//int longestHelixLength = flI.getLongestHelixLength();
			//int numBasesPaired = flI.getNumBasesPaired();
			//double normal = longestHelixLength*numBasesPaired;
			//Compute the length of the longest helix found.
			//normal = normal*Math.max(0,deltaG);
			return Math.max(-deltaG,0);
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public ScorePenalty clone() {
			MFEHybridScore ci = new MFEHybridScore(ds[0],ds[1],dir,entropicPenalty);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
	}
	
	/**
	 * Penalize complementarity at the base of a hairpin loop
	 */
	public class HairpinOpening extends ScorePenalty {
		private HairpinClosingTarget hairpin;
		private int markLeft = -1, markRight, jOffset;
		public HairpinOpening(HairpinClosingTarget hairpin, DesignIntermediateReporter dir){
			super(dir);
			this.ds = hairpin.stemAndOpening;
			this.hairpin = hairpin;
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			HairpinOpening ci = new HairpinOpening(hairpin,dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		public int getPriority(){
			return 1;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			if (markLeft==-1){
				int start = 0;
				int end = hairpin.stemAndOpening[0].length(domain);
				if (hairpin.outside){
					int middle = hairpin.stemAndOpening[0].length(domain)-hairpin.stemOnly[0].length(domain);
					markLeft = 0;
					markRight = middle;
					jOffset = 0; //Aligned at the left position.
				} else {
					int middle = hairpin.stemOnly[0].length(domain);
					markLeft = middle;
					markRight = end;
					//Aligned at the right, so offset j by the difference
					jOffset = hairpin.stemAndOpening[0].length(domain)-hairpin.stemAndOpening[1].length(domain);
				}
			}
			double StemAndOpeningScore =flI.helixDeltaG(hairpin.stemAndOpening[0],hairpin.stemAndOpening[1],domain,domain_markings,markLeft,markRight,jOffset); 
			double OnlyStem =flI.helixDeltaG(hairpin.stemOnly[0],hairpin.stemOnly[1],domain,domain_markings,0,0,0); 
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
	 * TODO:
	 * Calculate the Expected Minimum Free Energy for a random string of length N, and subtract that
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
			double deltaG = flI.mfeNoDiagonalPairing(ds[0], ds[1], domain, null);
			deltaG += ds[0].length(domain)/2; //Heuristic!?
			return Math.max(0,-deltaG);
		}
		public boolean affectedBy(int domain) {
			return ds[0].contains(domain);
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public int getNumDomainsInvolved() {
			return ds[0].numDomains;
		}
	}
	
	/**
	 * If we wish to use the validity checking as a penalty instead of an absolute condition,
	 * use this.
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
			for(int i = 0; i < domain.length; i++){
				sum += flI.affectedSequenceInvalidScore(i, seqs, domain, domain_markings);
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
	
	public class SelfFold extends ScorePenalty { 
		public SelfFold(DomainSequence ds, DesignIntermediateReporter dir){
			super(dir);
			this.ds = new DomainSequence[]{ds};
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			SelfFold ci = new SelfFold(ds[0],dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double deltaG = (flI.mfeSSDeltaG(ds[0],domain,domain_markings));
			return Math.max(0,-deltaG);
			//int longestHelixLength = flI.getLongestHelixLength();
			//int numBasesPaired = flI.getNumBasesPaired();
			//double normal = longestHelixLength*numBasesPaired;
			//return normal*Math.max(0,-deltaG);
		}
		public int getPriority(){
			return 2;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
	}
	
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
			flI.pairPrSS(probBuffer, ds[0], domain);
			
			double expectedMisPairedBases = 0;
			for(int i = 0; i < length1; i++){
				for(int j = i+1; j < length1; j++){
					if (0==DomainDesigner_SharedUtils.isAlignedAndShouldPair(ds[0], i, ds[0], j, domain)){
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
			flI.pairPrHybrid(probBuffer, ds[0], ds[1], domain);
			
			double expectedMisPairedBases = 0;
			for(int i = 0; i < N; i++){
				for(int j = i+1; j < N; j++){
					DomainSequence sI = (i < length1?ds[0] : ds[1]);
					int iInSi = (i < length1?i : i - length1);
					DomainSequence sJ = (j < length1?ds[0] : ds[1]);
					int jInSJ = (j < length1?j : j - length1);
					if (0==DomainDesigner_SharedUtils.isAlignedAndShouldPair(sI, iInSi, sJ, jInSJ, domain)){
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

	public List<ScorePenalty> listPenalties(
			AbstractDomainDesignTarget designTarget,
			DesignIntermediateReporter DIR, int[][] domain2, DesignerOptions options) {
		
		ArrayList<DomainSequence> rawStrands = designTarget.wholeStrands;
		ArrayList<DomainSequence> makeSS = new ArrayList();
		makeSS.addAll(designTarget.generalizedSingleStranded);
		makeSS.addAll(designTarget.singleDomains);
		ArrayList<HairpinClosingTarget> hairpinClosings = designTarget.hairpinClosings;
		
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(rawStrands);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(makeSS);

		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(hairpinClosings);
		
		List<ScorePenalty> allScores = new LinkedList<ScorePenalty>();
		allScores.add(new VariousSequencePenalties(rawStrands,DIR));
		
		for(int i = 0; i < makeSS.size(); i++){
			DomainSequence ds = makeSS.get(i);
			//Secondary Structure Formation
			if (DomainDesigner_SharedUtils.checkComplementary(ds, ds)){
				/*
				allScores.add(new LocalDefectSSScore(ds, DIR, designTarget));	
				//Dimerization
				allScores.add(new PairDefectScore(ds, ds, DIR, false, designTarget));
				*/
			} else {
				allScores.add(new SelfFold(ds, DIR));
				allScores.add(new MFEHybridScore(ds, ds, DIR, ds.numDomains==1));
			}

			//Hybridization
			for(int k = i+1; k < makeSS.size(); k++){ //Do only upper triangle
				DomainSequence ds2 = makeSS.get(k);
				boolean sameMol = ds.getMoleculeName().equals(ds2.getMoleculeName());
				if (DomainDesigner_SharedUtils.checkComplementary(ds, ds2)){
					//allScores.add(new PairDefectScore(ds, ds2, DIR, sameMol, designTarget));
				} else {
					allScores.add(new MFEHybridScore(ds, ds2, DIR, sameMol || (ds.numDomains==1 && ds2.numDomains==1)));
				}
			}
		}
		
		for(HairpinClosingTarget hairpin : hairpinClosings){
			allScores.add(new HairpinOpening(hairpin, DIR));
		}
		
		if (options.selfSimilarityPenalty.getState() >= 0){
			//Only do each domain once.
			List<DomainSequence> domainsOnceApiece = new ArrayList();
			domainsOnceApiece.addAll(designTarget.singleDomains);
			DomainDesigner_SharedUtils.utilRemoveDuplicateOrComplementaryDomains(domainsOnceApiece);
			for(DomainSequence domain : domainsOnceApiece){
				if (domain2[domain.domainList[0] & DomainSequence.DNA_SEQ_FLAGSINVERSE].length >= options.selfSimilarityPenalty.getState()){
					allScores.add(new SelfSimilarityScore(domain, DIR));
				}
			}
		}
		
		return allScores;
	}
}
