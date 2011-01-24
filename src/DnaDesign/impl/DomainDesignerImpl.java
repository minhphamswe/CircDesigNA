package DnaDesign.impl;

import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.LinkedList;
import java.util.List;

import DnaDesign.AbstractDomainDesignTarget;
import DnaDesign.DesignIntermediateReporter;
import DnaDesign.DnaDefinition;
import DnaDesign.DomainDesigner;
import DnaDesign.DomainDesigner_SharedUtils;
import DnaDesign.DomainSequence;
import DnaDesign.NAFolding;

/**
 * Implementation of DomainDesigner
 */
public class DomainDesignerImpl extends DomainDesigner{
	
	private NAFolding flI;
	/**
	 * @param foldingImpl<br>
	 * The folding score functions to utilize in evaluating solution candidates.
	 * @param designSSonly<br>
	 * 4 kinds of score functions are currently implemented: Validity, SingleStrandedAssurance, Crosstalk, Dimerization.
	 * If designSSonly is true, only SingleStrandedAssurance will be used; As a result, crosstalk and dimerization may occur
	 * in solution candidates.
	 */
	public DomainDesignerImpl(NAFolding foldingImpl) {
		this.flI = foldingImpl;
	}

	public class MFEHybridScore extends ScorePenalty { 
		private boolean entropicBonus = false;
		public MFEHybridScore(DomainSequence ds, DomainSequence ds2, DesignIntermediateReporter dir, boolean b){
			super(dir);
			this.ds = new DomainSequence[]{ds,ds2};
			chooseScore(dir);
			entropicBonus = b;
		}		
		public int getPriority(){
			return 0;
		}
		private DomainSequence[] ds;
		//This seems too low.
		private final double BIMOLECULAR = -.513;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			double deltaG = (flI.mfeHybridDeltaG(ds[0],ds[1],domain,null))+(entropicBonus?-BIMOLECULAR:0);
			int longestHelixLength = flI.getLongestHelixLength();
			int numBasesPaired = flI.getNumBasesPaired();
			double normal = longestHelixLength*numBasesPaired;
			//Compute the length of the longest helix found.
			normal = normal*Math.max(0,deltaG);
			return normal;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public ScorePenalty clone() {
			MFEHybridScore ci = new MFEHybridScore(ds[0],ds[1],dir,entropicBonus);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
	}
	
	/**
	 * Penalize complementarity at the base of a hairpin loop
	 */
	public class HairpinOpening extends ScorePenalty {
		public HairpinOpening(DomainSequence dsL, DomainSequence dsR, DesignIntermediateReporter dir){
			super(dir);
			this.ds = new DomainSequence[]{dsL,dsR};
			chooseScore(dir);
		}
		public ScorePenalty clone() {
			HairpinOpening ci = new HairpinOpening(ds[0],ds[1],dir);
			ci.old_score = old_score;
			ci.cur_score = cur_score;
			return ci;
		}
		public int getPriority(){
			return 0;
		}
		private DomainSequence[] ds;
		public double evalScoreSub(int[][] domain, int[][] domain_markings){
			for(int k = 0; k < 2; k++){
				if (DnaDefinition.bindScore(ds[0].base(k, domain),ds[1].base(ds[1].length(domain)-1-k, domain))>0){
					return 1;
				}
			}
			return 0;
		}
		public boolean affectedBy(int domain) {
			return (ds[0].domainList[0]& DNA_SEQ_FLAGSINVERSE)==domain ||
			(ds[1].domainList[ds[1].domainList.length-1]& DNA_SEQ_FLAGSINVERSE)==domain;
		}
		public DomainSequence[] getSeqs() {
			return ds;
		}
		public int getNumDomainsInvolved() {
			return 2;
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
			return 1;
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
			int longestHelixLength = flI.getLongestHelixLength();
			int numBasesPaired = flI.getNumBasesPaired();
			double normal = longestHelixLength*numBasesPaired;
			return normal*Math.max(0,deltaG);
		}
		public int getPriority(){
			return 0;
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
			return 1;
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
			DesignIntermediateReporter DIR) {
		
		List<DomainSequence> rawStrands = designTarget.wholeStrands;
		List<DomainSequence> eachDomainWithOverhang = designTarget.singleDomainsWithOverlap;
		
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(rawStrands);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(eachDomainWithOverhang);

		DomainDesigner_SharedUtils.utilRemoveSubsequences(eachDomainWithOverhang);
		
		List<ScorePenalty> allScores = new LinkedList<ScorePenalty>();
		allScores.add(new VariousSequencePenalties(rawStrands,DIR));
		
		for(int i = 0; i < eachDomainWithOverhang.size(); i++){
			DomainSequence ds = eachDomainWithOverhang.get(i);
			//Secondary Structure Formation
			if (DomainDesigner_SharedUtils.checkComplementary(ds, ds)){
				allScores.add(new LocalDefectSSScore(ds, DIR, designTarget));	
				//Dimerization
				allScores.add(new PairDefectScore(ds, ds, DIR, false, designTarget));
			} else {
				allScores.add(new SelfFold(ds, DIR));
				allScores.add(new MFEHybridScore(ds, ds, DIR, false));
			}

			//Hybridization
			for(int k = i+1; k < eachDomainWithOverhang.size(); k++){ //Do only upper triangle
				DomainSequence ds2 = eachDomainWithOverhang.get(k);
				boolean sameMol = ds.getMoleculeName().equals(ds2.getMoleculeName());
				if (DomainDesigner_SharedUtils.checkComplementary(ds, ds2)){
					allScores.add(new PairDefectScore(ds, ds2, DIR, sameMol, designTarget));
				} else {
					allScores.add(new MFEHybridScore(ds, ds2, DIR, sameMol));
				}
			}
		}
		
		return allScores;
	}
}
