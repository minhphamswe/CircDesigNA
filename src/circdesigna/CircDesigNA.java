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

import static circdesigna.DomainSequence.DNAMARKER_DONTMUTATE;
import static circdesigna.DomainSequence.NA_COMPLEMENT_FLAG;
import static circdesigna.DomainSequence.NA_COMPLEMENT_FLAGINV;
import static circdesigna.abstractpolymer.DnaDefinition.C;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import circdesigna.DesignIntermediateReporter.DesignIntermediateScore;
import circdesigna.SequenceDesigner.AlternativeResult;
import circdesigna.abstractDesigner.BlockDesigner;
import circdesigna.abstractDesigner.InfiniteResourceTournament;
import circdesigna.abstractDesigner.PopulationDesignMember;
import circdesigna.abstractDesigner.StandardTournament;
import circdesigna.abstractDesigner.TinyGADesigner;
import circdesigna.abstractpolymer.MonomerDefinition;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;
import circdesigna.energy.ConstraintsNAFoldingImpl;
import circdesigna.impl.CircDesigNAImpl;
import circdesigna.impl.CircDesigNAPMemberImpl;
import circdesigna.impl.CodonCode;
import circdesigna.impl.SequenceCode;
import circdesigna.impl.SequenceDesignBlockDesignerImpl;
import circdesigna.impl.SequencePenaltiesImpl;

/**
 * CircDesigNA designs a set of sequences (split into "domains", by the parameters in Domain Defs) 
 * to optimize a set of score functions over a sequence space which is itself a constrained subset of all possible dna
 * sequences. 
 * 
 * This class
 * - Provides a useful frontend to the GUI so that it can display real time design information
 * - Provides a general purpose main() method which defers to subclasses for the design routines.
 * 
 * Subclasses specify certain behaviors of the functions in main(), thus actually producing a full design algorithm.
 */
public abstract class CircDesigNA extends CircDesigNASystemElement{
	public CircDesigNA(CircDesigNAConfig std) {
		super(std);
	}
	/**
	 * Utility method for creating a designer with the default scoring algorithms on the specified design scheme
	 * 
	 * Some designer options are also taken from the CircDesigNAConfig object. 
	 */
	public static SequenceDesigner<CircDesigNAOptions> getDefaultDesigner(String Molecules,String domainDefs,CircDesigNAConfig Std) {
		ConstraintsNAFoldingImpl folder = new circdesigna.energy.ConstraintsNAFoldingImpl(Std);
		CircDesigNA r = new CircDesigNAImpl(folder, new SequencePenaltiesImpl(Std), Std);
		return new SequenceDesignAdapter(r, Molecules, domainDefs);
	}
	
	/**
	 * A coupled wrapper which allows the GUI to easily display design internals.
	 * 
	 * Why does this class exist? The CircDesigNA did not expose an API which made GUI
	 * writing as easy as it could be. To make client code simpler, this class has some adaptor methods (especially a constructor) 
	 * that will transform your input before feeding it to a CircDesigNA. 
	 * 
	 * Error checking is deferred to the CircDesigNA class itself.
	 * 
	 * Finally, this class serves to separate the public api which controls the behavior of the designer (run,
	 * pause abort) from the internal api of how to do so (mutate, run block iteration).
	 */
	public static class SequenceDesignAdapter implements SequenceDesigner<CircDesigNAOptions>{
		private boolean waitOnStart = true;
		private boolean finished = false;
		
		CircDesigNA r;
		private Thread runner;
		private List<String> inputMolecules;
		private String domainDefsBlock;
		private Map<Integer, String> lock;
		private Map<Integer, String> initial;
		private Map<Integer, DesignerCode> mutators;
		private Map<Integer, String> positionConstraints;
		private DesignIntermediateReporter dir;
		public SequenceDesignAdapter(CircDesigNA r, String molecules, String domainDefsBlock){
			ArrayList<String> inputStrands = new ArrayList<String>();
			String q = "";
			for(String newLine : molecules.split("\n")){
				q = newLine;
				if (q.trim().length()>0){
					inputStrands.add(q);
				}
			}

			this.inputMolecules = inputStrands;
			this.domainDefsBlock = domainDefsBlock;
			dir = new DesignIntermediateReporter();
			
			this.r = r;
			runOnStart.run();
			//new Thread(runOnStart).start();
		}
		/**
		 * Decodes the input, and then calls Main.
		 */
		private Runnable runOnStart = new Runnable(){
			public void run(){			
				
				//Read domainDefs
				int num_domain = 0;
				TreeMap<Integer, Integer> domain_length_t = new TreeMap<Integer, Integer>();
				final DomainDefinitions dsd = new DomainDefinitions(r.Std);
				DomainDefinitions.readDomainDefs(domainDefsBlock, dsd);
				num_domain = dsd.domainLengths.length;
				
				//Read molecules
				final AbstractDomainDesignTarget designTarget = new AbstractDomainDesignTarget(dsd,r.Std);
				for(int kq = 0; kq < inputMolecules.size(); kq++){
					String q = inputMolecules.get(kq);
					//DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
					String[] line = q.split("\\s+",2);
					if (line.length<2){
						throw new RuntimeException("Correct molecule format: <name> <molecule>");
					}
					designTarget.addTargetStructure(q);
				}
				for(int i = 0; i < dsd.domainLengths.length; i++){
					int val = dsd.domainLengths[i];
					if (val!=-1){
						domain_length_t.put(i, val);
					}
				}
				if (num_domain==0){
					throw new RuntimeException("No valid molecules to design.");
				}
				
				final int num_domain_2 = num_domain;
				final int[] domain_length = new int[num_domain_2];

				if (initial==null){
					initial = new TreeMap<Integer, String>();
				}
				if (mutators == null){
					mutators = new TreeMap<Integer, DesignerCode>();
				}
				//main() requires a nonnull set of initial and mutators (because they are maps).
				//That being said, they may very well both be empty maps.
				
				positionConstraints = new TreeMap<Integer, String>();
				
				for(int k = 0; k < domain_length.length; k++){
					Integer got = domain_length_t.get(k);
					if (got==null){
						throw new RuntimeException("No length data for domain \'"+dsd.getDomainName(k)+"\'");
					}
					domain_length[k] = got;
					
					//Constraints
					positionConstraints.put(k,dsd.getConstraint(k));
					initial.put(k,dsd.getInitialSequence(k));
					if (dsd.maintainAminos(k)){
						mutators.put(k, new CodonCode(r.Std));
					}
				}
				

				runner = new Thread(){public void run(){
					while(waitOnStart && ! r.abort){
						try {
							Thread.sleep(10);
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
					}
					if (r.abort){
						return;
					}
					ArrayList<Integer> results = new ArrayList<Integer>();
					int timesToRun = 1;
					while(timesToRun>0){
						timesToRun--;
						errorResult = null;
						try {
							results.add(r.main(num_domain_2, domain_length, Integer.MAX_VALUE, lock, initial, positionConstraints, mutators, designTarget, dir, dsd));
						} catch (OutOfMemoryError e){
							errorResult = "Error in design: " +e.getMessage()+"\nJava ran out of memory! Strategies to fix:\n\t1) Reduce the population size.\n\t2)Download this applet, and run it with the argument -Xmx700 from the command line.";//ensure nonnull
						} catch (Throwable e){
							e.printStackTrace();
							errorResult = "Error in design: " +e.getMessage()+"";//ensure nonnull
							break;
						}
					}
					//System.out.println(results);
					finished = true;

					//Help the designer out in handling listeners:
					if (r.runOnIteration!=null){
						r.runOnIteration.actionPerformed(null);
					}
				}};
				
				runner.start();
			}
		};
		
		public CircDesigNAOptions getOptions() {
			return r.options;
		}

		private String errorResult;
		
		public AlternativeResult[] getAlternativeResults(){
			return r.alternatives;
		}
		
		public String getResult() {
			if (r.alternatives==null){
				return getResult(null);
			}
			return getResult(r.alternatives[0]);
		}
		/**
		 * Returns a human-readable dump of the output.
		 */
		public String getResult(AlternativeResult res2) {
			String progress = "";
			if (r.dbesign != null){
				progress = String.format("Progress to next iteration: %.0f%%\n", 100*r.dbesign.getProgress());
			}
			if (res2==null){
				return progress + "No output. Please wait, and then click the button above to get intermediate results.";
			}
			if (res2.TYPE == AlternativeResult.ERROR){
				return errorResult;
			}
			if (res2.TYPE == AlternativeResult.LOG){
				return r.iteration_history.toString();
			}
			DomainDesignerAlternativeResult res = (DomainDesignerAlternativeResult) res2;

			ArrayList<DomainSequence> alreadyPrintedSequences = new ArrayList();
			DomainDefinitions dsd = new DomainDefinitions(r.Std);
			//DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
			DomainDefinitions.readDomainDefs(domainDefsBlock, dsd);
			StringBuffer sb = new StringBuffer();
			String lR = "\n";
			sb.append(progress);
			DesignScoreBreakdown breakdown = res.getBreakdown();
			String[] outputDomains = res.getOutputDomains();
			
			sb.append(String.format("Design Phase\t%d\nIteration\t%d\nNet 'Score'\t%.2f.\n>> Cross Interactions\t%.2f.\n>> Breathing Helix Ends\t%d.\n>> Self-Folding Interactions\t%.2f.\n>> Banned Patterns\t%.2f.\n",res.getPhase(),res.getIteration(),breakdown.netScore,breakdown.crossInteractionsOnly,breakdown.breathingHelixes,breakdown.selfFoldOnly,breakdown.bannedPatterns));
			sb.append(lR);
			sb.append("Current designer state (to resume from this point, paste as 'Domain Definition'):");
			sb.append(lR);
			sb.append("----------------");
			sb.append(lR);
			for(int k = 0; k < outputDomains.length; k++){
				sb.append(dsd.getDomainName(k));
				sb.append("\t");
				sb.append(outputDomains[k].length());
				String arg = dsd.getArguments(k); //Does not include -init()
				if (arg.length()>0){
					sb.append("\t"+dsd.getArguments(k));	
				}
				sb.append("\t");
				sb.append("-init(");
				sb.append(outputDomains[k]);
				sb.append(")");
				sb.append(lR);
			}
			sb.append("----------------");
			sb.append(lR);
			sb.append("Molecular Substrands:");
			sb.append(lR);	
			for(String q : inputMolecules){
				String a[] = q.split("\\s+",2);    
				boolean alreadyPrintedInMolecule = false;
				splitLoop: for(String subStrand : a[1].split("}")){
					if (subStrand.trim().length()==0){
						continue;
					}
					DomainSequence ds = new DomainSequence();
					ds.setDomains(subStrand,dsd,null);
					for(DomainSequence g : alreadyPrintedSequences){
						if (g.equals(ds)){
							continue splitLoop;
						}
					}
					if (!alreadyPrintedInMolecule){
						//Print out the name of the molecule:
						sb.append("Strands in molecule "+a[0]+":");
						sb.append(lR);
					}
					alreadyPrintedInMolecule = true;
					alreadyPrintedSequences.add(ds);
					sb.append("\t"+"["+subStrand.replace("[","")+"}");
					sb.append(lR+"\t");
					printSequence(sb,r,ds,true,outputDomains);
					sb.append(lR+"\t");
					printSequence(sb,r,ds,false,outputDomains);
					sb.append(lR);
				}
			}
			return sb.toString();
		}
		private void printSequence(StringBuffer sb, CircDesigNA r2, DomainSequence ds, boolean withSeperator, String[] domains) {
			if (withSeperator)
				sb.append("[");
			for(int k = 0; k < ds.numDomains; k++){
				String domain = domains[ds.domainList[k] & NA_COMPLEMENT_FLAGINV];
				if ((ds.domainList[k] & NA_COMPLEMENT_FLAG)!=0){
					domain = revComp(domain);
				}
				sb.append(domain);
				if (withSeperator && k + 1 < ds.numDomains){
					sb.append(" ");
				}
			}
			if (withSeperator)
				sb.append("}");
		}
		private String revComp(String q){
			q = q.replaceAll("\\s+","");
			char[] arr = new char[q.length()];
			for(int i = 0; i < arr.length; i++){
				arr[i] = Character.toUpperCase(q.charAt(arr.length-1-i));
				if (arr[i]=='G') arr[i] = 'C';
				else if (arr[i]=='C') arr[i] = 'G';
				else if (arr[i]=='A') arr[i] = 'T';
				else if (arr[i]=='T') arr[i] = 'A';
			}
			return new String(arr);
		}

		
		/**
		 * State functions
		 */
		public boolean isFinished() {
			return finished;
		}
		
		public boolean isEndConditionError(){
			return errorResult!=null;
		}

		public boolean isRunning() {
			return !waitOnStart && !r.waitForResume;
		}


		public void pause() {
			r.waitForResume = true;
		}

		public void resume() {
			waitOnStart = false;
			r.waitForResume = false;
		}

		public void abort() {
			r.abort = true;
			runner.interrupt();
		}
		public DesignIntermediateReporter getDir() {
			return dir;
		}
		public double getBestScore() {
			if (r.bestScore < 0){
				return 0;
			}
			return r.bestScore;
		}
		public int getCurrentIteration() {
			return r.design_iteration;
		}
		public void runIteration() {
			final boolean[] return_shared = new boolean[]{false};
			r.runOnIteration = new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					pause();
					return_shared[0] = true;
				}
			};
			resume();
			while(!return_shared[0]){
				if (isFinished()){
					break;
				}
				try {
					Thread.sleep(1);
				} catch (InterruptedException e1) {
					e1.printStackTrace();
				}
			}
			r.runOnIteration = null;
		}
		
		public DesignScoreBreakdown getScoreBreakdown(){
			if (r.alternatives==null){
				return null;
			}
			return getScoreBreakdown(r.alternatives[0]);
		}
		public DesignScoreBreakdown getScoreBreakdown(AlternativeResult ar){
			return ((DomainDesignerAlternativeResult)ar).getBreakdown();
		}

	}
	
	private class SpecialAlternativeResult extends AlternativeResult {
		private String description;
		public SpecialAlternativeResult(String description, int TYPE, int ID) {
			this.description = description;
			this.TYPE = TYPE;
			this.ID = ID;
		}
		public String getDescription() {
			return description;
		}	
	}
	
	private class DomainDesignerAlternativeResult extends AlternativeResult{
		private DesignScoreBreakdown breakdown;
		private String[] outputDomains;
		private String description;
		private int phase;
		private int iteration;
		public DomainDesignerAlternativeResult(String description, CircDesigNAPMemberImpl best, AbstractDomainDesignTarget designTarget, DomainDefinitions dsd) {
			breakdown = CircDesigNA_SharedUtils.getScoreBreakdown(best.domain,designTarget,best.penalties);
			outputDomains = displayDomains(best.domain, false, dsd);
			this.description = description;
			this.phase = CircDesigNA.this.design_phase;
			this.iteration = dbesign.getIterationCount();
		}
		public int getPhase(){
			return phase;
		}
		public int getIteration(){
			return iteration;
		}
		public DesignScoreBreakdown getBreakdown(){
			return breakdown;
		}
		public String[] getOutputDomains() {
			return outputDomains;
		}
		public String getDescription(){
			return description;
		}
	}
	
	public CircDesigNAOptions options = CircDesigNAOptions.getDefaultOptions();
	
	private int design_iteration = 0, design_phase = 0;
	private StringBuffer iteration_history = new StringBuffer();
	private ActionListener runOnIteration = null;
	private double bestScore = -1;
	private AlternativeResult[] alternatives;
	private BlockDesigner<CircDesigNAPMemberImpl> dbesign;
	private boolean waitForResume = false;

	public boolean abort = false;

	//End accessible to client.
	//Constructor
	
	/**
	 * Debug output flags
	 */
	private boolean DEBUG_PENALTIES = true;
	
	/**
	 * Use the evaluator to determine where the problem spots are, 
	 * and mutate them. Recommended, otherwise you're merely randomly mutating.
	 */
	private boolean ENABLE_MARKINGS = true;
	
	/**
	 * Use to add probability of mutation to bases involved in many penalties. 
	 */
	private boolean SORT_MARKINGS = true;
	
	/**
	 * Take the evaluators word as absolute, and never mutate bases that
	 * aren't directly involved in the "current worst penalty". (Otherwise, occasionally
	 * mutate random bases at a low frequency)
	 */
	private boolean MARKINGS_ENABLESTRICT = true;
	
	public int MAX_MUTATION_DOMAINS = 999;
	public int MAX_MUTATIONS_LIMIT = 9999; // maximum number of simultaneous mutations

	//It is sometimes advantageous to design sequences with some immutable complementarity.
	public boolean ALLOW_COMPLEMENTARY_SCORES = false;

	public abstract static class ScorePenalty {
		public static final double MAX_SCORE = 1e18;
		public double old_score, cur_score;
		//This is only true if we have evaluated this penalty, but not committed or reverted it yet.
		//Set by the USER of this class.
		public boolean in_intermediate_state;
		public DesignIntermediateScore dis;
		public DesignIntermediateReporter dir;
		public ScorePenalty(DesignIntermediateReporter dir){
			old_score = cur_score = MAX_SCORE; //A suitably large number
			this.dir = dir;
			in_intermediate_state = false;
		}
		public double getScore(int[][] domain, int[][] domain_markings){
			evalScore(domain, domain_markings);
			dedicate();
			return old_score;
		}
		public void revert(){
			cur_score = old_score;
			in_intermediate_state = false;
		}
		public void dedicate(){
			old_score = cur_score;
			in_intermediate_state = false;
		}
		public final double check(double possScore){
			return Math.min(possScore,MAX_SCORE);
		}
		public final double evalScore(int[][] domain, int[][] domain_markings){
			double new_score = check(evalScoreSub(domain, domain_markings));
			if (new_score < 0){
				throw new RuntimeException("Negative subscore.");
			}
			cur_score = new_score;
			return getCDelta();
		}
		public final double getCDelta(){
			return cur_score - old_score;
		}
		protected void chooseScore(DesignIntermediateReporter dir) {
			if (dir==null){
				return;
			}
			DomainSequence[] ds = getSeqs();
			if (ds.length==0){
				return;
			} else
			if (ds.length==1){
				dis = dir.chooseDesignIntermediateScore(ds[0].getMoleculeName(), ds[0].getMoleculeName());
			} else
			if (ds.length==2){
				dis = dir.chooseDesignIntermediateScore(ds[0].getMoleculeName(), ds[1].getMoleculeName());
			}
		}
		public abstract ScorePenalty clone();
		public abstract double evalScoreSub(int[][] domain, int[][] domain_markings);
		public boolean affectedBy(int domain){
			for(DomainSequence q : getSeqs()){
				if (q.contains(domain)){
					return true;
				}
			}
			return false;
		}
		public abstract DomainSequence[] getSeqs();
		public int getNumDomainsInvolved(){
			int sum = 0;
			for(DomainSequence q : getSeqs()){
				sum += q.numDomains;
			}
			return sum;
		}
		public abstract int getPriority();
		public String toString(DomainDefinitions dsd) {
			String myString = getClass().getSimpleName();
			DomainSequence[] seqs = getSeqs();
			if (seqs.length>0){
				myString += " : ";
				for(int k = 0; k < seqs.length; k++){
					myString += seqs[k].toString(dsd);
					if (k + 1 < seqs.length){
						myString += " vs ";
					}
				}
			} else {
				myString += " (On all molecules)";
			}
			return myString;
		}
	}
	
	private void fullScoreReport(PopulationDesignMember<CircDesigNAPMemberImpl>[] population, DesignIntermediateReporter DIR, AbstractDomainDesignTarget designTarget, DomainDefinitions dsd){
		iteration_history.append(String.format("%d %.3f",design_iteration, bestScore)+"\n");
		
		CircDesigNAPMemberImpl best = CircDesigNA_SharedUtils.getBestMember(population);
		
		//Report best member, primarily
		DIR.beginScoreReport();
		for(ScorePenalty s : best.penalties){
			if (s.dis!=null){
				s.dis.addScore(s.old_score);
			}
		}
		DIR.endScoreReport();
		DesignScoreBreakdown outputScore = CircDesigNA_SharedUtils.getScoreBreakdown(best.domain,designTarget,best.penalties);
		displayDomains(best.domain, true, dsd);
		
		//Other miscellaneous data kept include the farthest member from the best (by hamming distance)
		//and every population member.
		CircDesigNAPMemberImpl farthest = CircDesigNA_SharedUtils.getFarthestFrom(best,population);
		
		AlternativeResult[] tmp_alternatives = new AlternativeResult[AlternativeResult.OTHER + population.length];
		tmp_alternatives[AlternativeResult.BEST] = new DomainDesignerAlternativeResult("Best Design", best, designTarget, dsd);
		tmp_alternatives[AlternativeResult.FARTHEST1] = new DomainDesignerAlternativeResult("Farthest from Best", farthest, designTarget, dsd);
		tmp_alternatives[AlternativeResult.ERROR] = new SpecialAlternativeResult("Design Status", AlternativeResult.ERROR, 0);
		tmp_alternatives[AlternativeResult.LOG] = new SpecialAlternativeResult("Design Log", AlternativeResult.LOG, 0);
		for(int i = 0; i < population.length; i++){
			tmp_alternatives[AlternativeResult.OTHER+i] = new DomainDesignerAlternativeResult("Population ID "+i, (CircDesigNAPMemberImpl)population[i], designTarget, dsd);
		}
		alternatives = tmp_alternatives;
	}
	
	/**
	 * The program will attempt to prevent any interaction between the DNASequences in toSynthesize, though interactions
	 * containing complementary sequences will be ignored (i.e, [1 4* 5} will not be tested against [1 4 5}).
	 * @param initial 
	 * @param hairpinInnards 
	 * @param dsd 
	 */
	int main(int num_domain, int[] domain_length, int MAX_ITERATIONS, Map<Integer, String> lockedDomains, Map<Integer, String> initial, Map<Integer, String> positionConstraints, Map<Integer, DesignerCode> mutators, AbstractDomainDesignTarget designTarget, DesignIntermediateReporter DIR, DomainDefinitions dsd) {
		while (waitForResume && !abort){
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("----------------------------------");
		System.out.println("           CircDesigNA");
		System.out.println("----------------------------------");
		System.out.println("Designer started on "+new Date());
		
		//Initialization
		design_phase = 1;
		setPhase(design_phase);
		iteration_history.append("Phase "+design_phase+"\n");
		
		//Domains to be designed. The integer type is being used to hold DNA bases,
		//with additional metadata (multiples of 10 are added to mean certain things)
		int[][] domain = new int[num_domain][];
		int[][] domain_markings = new int[num_domain][];
		DesignerCode[] mutate = new DesignerCode[num_domain];
		//initialize domains to zeroes
		for(int i = 0; i < num_domain; i++){
			int[] cDomain = new int[domain_length[i]];
			domain[i] = cDomain; //initialize to 0s
			domain_markings[i] = new int[domain_length[i]];
		}
		//Don't need this any more.
		domain_length=null;

		//Enumerate penalty scores (see FoldingImplTestGUI for a visual of this process) via "listPenalties"
		List<ScorePenalty> allScores = listPenalties(designTarget,DIR,domain,options,dsd);
		if (allScores.isEmpty()){
			throw new RuntimeException("No scores to optimize : Algorithm has nothing to do!");
		}		
		System.out.println("Discovered "+allScores.size()+" score elements. Listing them:");
		for(ScorePenalty q : allScores){
			System.out.println(q.toString(dsd));
		}

		
		//Load the constraints
		if (positionConstraints!=null){
			for(Entry<Integer, String> init : positionConstraints.entrySet()){
				int i = init.getKey();
				for(int j = 0; j < domain[i].length; j++){
					try {
						domain[i][j] = Std.monomer.decodeConstraintChar(init.getValue().charAt(j));
					} catch (Throwable e){
						throw new RuntimeException(e.getMessage()+" (Domain "+dsd.getDomainName(i)+")");
					}
				}
			}
		}
		
		//Load the initial domains, and (possibly) check them against the constraints.
		if (initial!=null){
			for(Entry<Integer, String> init : initial.entrySet()){
				int i = init.getKey();
				for (int j = 0; j < domain[i].length; j++) {
					try {
						final char charAt = init.getValue().charAt(j);
						int base = Std.monomer.decodeInitializationChar(charAt);
						
						if (base==Std.monomer.NOBASE){
							continue;
						}
						
						if (Std.monomer.allowBase(domain[i][j], base)){
							//Don't allow initial bases unless they match!
							domain[i][j] = domain[i][j] - Std.monomer.noFlags(domain[i][j]) + base;
						} else {
							throw new RuntimeException("Base "+charAt+" conflicts with constraint at position "+(j+1));
						}
					} catch (Throwable e){
						throw new RuntimeException(e.getMessage()+" (Domain "+dsd.getDomainName(i)+")");
					}
				}
			}
		}

		//Generate a list of mutable domains (optimization)
		ArrayList<Integer> mutableDomainsL = new ArrayList();
		for(int i = 0; i < num_domain; i++){
			boolean mutable = false;
			//Optimize for completely locked domains.
			loop:for (int j = 0; j < domain[i].length; j++) {
				if (domain[i][j] - Std.monomer.noFlags(domain[i][j]) != Std.monomer.LOCK_FLAG()){
					mutable = true;
					break loop;
				}
			}
			if (mutable){ //has at least one mutable base
				mutableDomainsL.add(i);		
			}
		}

		//Derive a mutation strategy for all domains (by parsing domain constraints with the options selected)
		for(int i = 0; i < num_domain; i++){
			if (mutators!=null){
				mutate[i] = mutators.get(i);
				//If one is specified, use that one. 
				if (mutate[i]==null){//Otherwise..
					SequenceCode newCode = Std.getDefaultSequenceCode();
					DesignSequenceConstraints dsc = Std.getDefaultConstraints();
					dsd.loadConstraints(i,dsc,false);
					newCode.setConstraints(dsc);
					mutate[i] = newCode;
				}
			} else {
				throw new RuntimeException("Mutators cannot be null");
			}
		}

		//Initialize sequences, using initial seed (if initial sequence provided, no randomization occurs.)
		randomizeSequence(domain, mutate, dsd, false);

		//Locked domains. Legacy feature.
		if (lockedDomains!=null){
			for(Entry<Integer, String> lock : lockedDomains.entrySet()){
				int i;
				mutableDomainsL.remove(i = lock.getKey());
				for (int j = 0; j < domain_length[i]; j++) {
					domain[i][j] = getLockedBase(lock.getValue().charAt(j));
				}
			}
		}
		
		//Make the list of mutableDomains again (optimization)
		int[] mutableDomains = new int[mutableDomainsL.size()];
		for(int i = 0; i < mutableDomains.length; i++){
			mutableDomains[i] = mutableDomainsL.get(i);
		}
		
		//Assertion: all domains are valid according to their respective mutation strategies. 

		//Select the penalties which depend on certain domains. (optimization).
		int[][] scoredElements = new int[num_domain][];
		for(int k = 0; k < num_domain; k++){
			ArrayList<Integer> affectedBySeq = new ArrayList();
			for(int i = 0; i < allScores.size(); i++){
				ScorePenalty q = allScores.get(i);
				if (q.affectedBy(k)){
					affectedBySeq.add(i);
				}
			}
			scoredElements[k] = new int[affectedBySeq.size()];
			for(int i = 0; i  < affectedBySeq.size(); i++){
				scoredElements[k][i] = affectedBySeq.get(i);
			}
		}

		try {
			//Begin "Population member" designer of abstraction.
			//We will duplicate our initial population member POPSIZE times.
			//First though, calculate its score:
			CircDesigNAPMemberImpl initialSeed = new CircDesigNAPMemberImpl(allScores,scoredElements,domain,domain_markings);
			
			//Briefly describe the initial sequence member (because it's interesting)
			System.out.println("Randomized initial sequence (for population member 0):");
			displayDomains(domain,true,dsd);
			
			//Score the original
			scorePopulation(new CircDesigNAPMemberImpl[]{initialSeed}, 0, 0);
			selectPhaseAndUpdateScore(new CircDesigNAPMemberImpl[]{initialSeed});
			
			//Create block designer, which will produce a certain initial population from the initial sequences we chose.
			CircDesigNAPMemberImpl tempMember = initialSeed.designerCopyConstructor(-1); //needed for "reverting" mutations
			SequenceDesignBlockDesignerImpl dbesignSingle = new SequenceDesignBlockDesignerImpl(mutableDomains,mutate,dsd,this,tempMember);
			dbesign = null;
			
			if (options.standardUseGA.getState()){
				dbesign = new TinyGADesigner(dbesignSingle, new CircDesigNAImpl.DParetoSort(), true);
			} else {
				if (options.resourcePerMember.getState() < 0){
					dbesign = new InfiniteResourceTournament(dbesignSingle);
				} else {
					dbesign = new StandardTournament(dbesignSingle, options.resourcePerMember.getState());
				}
			}
			
			{
				dbesign.initialize(initialSeed, options.population_size.getState());	
			}

			//Randomize the clones?
			boolean randomizePopulation = true;
			if (randomizePopulation){  //If false, then all population members start at the same location.
				System.out.println("Randomizing "+dbesign.getPopulation().length+" - 1 population members.");
				//If true, seed each of the population members differently. 
				PopulationDesignMember<CircDesigNAPMemberImpl>[] a = dbesign.getPopulation();
				for(int i = 1; i < a.length; i++){
					CircDesigNAPMemberImpl r = (CircDesigNAPMemberImpl)a[i];
					
					for(int k = 0; k < r.domain.length; k++){
						mutate[k].reset();
						for(int j = 0; j < r.domain[k].length; j++){
							mutate[k].mutateToOther(r.domain,k,j);
						}
						if (true){
							for(int j = 0; j < r.domain[k].length; j++){
								System.out.print(Std.monomer.displayBase(r.domain[k][j]));
							}
							System.out.println();
						}
					}
					//Initial score.

					if (abort){
						return 0;
					}

					//}DIR.endScoreReport();
					dbesign.setProgress((i+1), a.length);
					
					if (abort){
						return 0;
					}
				}
			}
			
			//Compute initial scores using current phase
			scorePopulation(dbesign.getPopulation(), 1, dbesign.getPopulation().length-1);
			selectPhaseAndUpdateScore(dbesign.getPopulation());
			//Report initial scores
			fullScoreReport(dbesign.getPopulation(),DIR,designTarget,dsd);			
		
			design_iteration = 0;

			double endingScore = options.end_score_threshold.getState();
			
			int finalPhase = countPhases();
			while (!((design_phase == finalPhase && bestScore <= endingScore) || (design_iteration >= MAX_ITERATIONS) || abort)) {
				while (waitForResume && !abort){
					try {
						Thread.sleep(100);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				if (abort){
					break;
				}
				//We are currently in this iteration:
				design_iteration++;

				//Update phase and rescore, if necessary.
				selectPhaseAndUpdateScore(dbesign.getPopulation());
				if (options.random_design.getState()){
					for(CircDesigNAPMemberImpl member : dbesign.getPopulation()){
						randomizeSequence(member.domain, mutate, dsd, true); //Completely random sequence, ignoring initial seed
					}
					scorePopulation(dbesign.getPopulation(), 0, dbesign.getPopulation().length-1);
					System.out.println("Iteration "+design_iteration+" Score "+bestScore);
				} else {
					//Mutate and update scores, under the assumption that phase does not change
					dbesign.runBlockIteration(this,options.end_score_threshold.getState());
					System.out.flush();
					CircDesigNAPMemberImpl q = dbesign.getBestPerformingChild();
					if (q==null){
						System.out.println("Iteration did not complete.");
						break;
					}
					updateBestScore(dbesign.getPopulation());
				}
				System.out.flush();

				//Iteration complete.
				resetDebugPenalty();
				//Report the status of the best member.
				fullScoreReport(dbesign.getPopulation(),DIR,designTarget,dsd);

				//debugPenalty(s, q.domain, dsd);
				//System.out.println(num_domains_mut);
				
				if (runOnIteration!=null){
					runOnIteration.actionPerformed(null);
				}
			}
		} finally {
			System.out.print("Designer ended after "+design_iteration+" iterations");
			if (bestScore>=0){
				System.out.print("with a score of "+bestScore);
			}
			System.out.println();
		}
		return design_iteration;
	}
	private void randomizeSequence(int[][] domain, DesignerCode[] mutate, DomainDefinitions dsd, boolean forceNewSequence) {
		for(int i = 0; i < domain.length; i++){
			if (forceNewSequence){
				Arrays.fill(domain[i], MonomerDefinition.NOBASE);
			}
			if (options.rule_ccend_option.getState() && domain[i].length>0){
				if (domain[i][0]==0)domain[i][0] = C + Std.monomer.GCL_FLAG();
				if (domain[i][domain[i].length-1]==0)domain[i][domain[i].length-1] = C + Std.monomer.GCL_FLAG();
			}
			//Will initialize unconstrained portion of domain
			try {
				pickInitialSequence(domain,i,mutate[i]);
			} catch (Throwable e){
				throw new RuntimeException(e.getMessage()+" (Domain "+dsd.getDomainName(i)+")");
			}
		}
	}
	/**
	 * For all indices in [i, j], performs a complete rescoring of those member's scores.
	 * Then, the value of bestScore is updated to contain the best score of any member in the entire
	 * population (i.e. indices [0, pop.length-1].)
	 */
	public void scorePopulation(CircDesigNAPMemberImpl[] pop, int i, int j){
		for(int k = i; k <= j; k++){
			CircDesigNAPMemberImpl r = pop[k];
			for(ScorePenalty s : r.penalties){
				s.getScore(r.domain, r.domain_markings);
			}
		}
		updateBestScore(pop);
	}
	public void updateBestScore(CircDesigNAPMemberImpl[] pop){
		double bestScore_tmp = Double.POSITIVE_INFINITY;
		for(int k = 0; k < pop.length; k++){
			CircDesigNAPMemberImpl r = pop[k];
			double score = 0;
			for(ScorePenalty s : r.penalties){
				score += s.cur_score;
			}
			bestScore_tmp = Math.min(bestScore_tmp, score);
		}
		bestScore = bestScore_tmp;
	}
	public void selectPhaseAndUpdateScore(CircDesigNAPMemberImpl[] pop){
		double endingScore = options.end_score_threshold.getState();
		while (bestScore <= endingScore && design_phase < countPhases()){
			design_phase++;
			iteration_history.append("Phase "+design_phase+"\n");
			System.out.println("Entering Phase "+design_phase);
			setPhase(design_phase);
			//Need complete rescore because score function probably changed
			scorePopulation(pop, 0, pop.length-1);
		}
	}
	
	public abstract int countPhases();
	public abstract void setPhase(int phase);
	
	public abstract List<ScorePenalty> listPenalties(
			AbstractDomainDesignTarget designTarget, DesignIntermediateReporter DIR, int[][] domain, CircDesigNAOptions options2, DomainDefinitions dsd) ;

	public static final void deepFill(int[][] domain_markings, int i) {
		for(int[] row : domain_markings){
			Arrays.fill(row,i);
		}
	}
	/**
	 * 
	 * @param designerCode 
	 */
	private void pickInitialSequence(int[][] domain, int mut_domain, DesignerCode mutator) {
		mutator.reset();
		
		boolean[] preserveInitial = new boolean[domain[mut_domain].length];
		boolean allowModifyInitialSpecifiedSequence = false;
		initialLoop: for(int k = 0; !abort; k++){
			if (k > 10 && !allowModifyInitialSpecifiedSequence){ //10 times before allowing modifications to initial sequence.
				System.err.println("Attempting to modify initial sequence to get inside constraints...");
				allowModifyInitialSpecifiedSequence = true;
			}
			if ( k > 1000 ){
				throw new RuntimeException("Initial constraints were too strict.");
			}
			for(boolean specialFlagsOnly : new boolean[]{true,false}){
				for(int j = 0; j < domain[mut_domain].length; j++){
					int flag = domain[mut_domain][j]-Std.monomer.noFlags(domain[mut_domain][j]);
					if ((flag!=0)!=specialFlagsOnly){
						continue;
					}
					if (Std.monomer.noFlags(domain[mut_domain][j])==Std.monomer.NOBASE){
						//We must fill these in (they are NOBASE)
					} else {
						//We have a flag.
						if (flag==0){
							//Ok, initial base specified. Leave it alone.
							preserveInitial[j] = true;
							if (!allowModifyInitialSpecifiedSequence){
								continue;
							}
							//Ok, we're allowed to mutate. But still, don't completely get rid of the sequence...
							if (Math.random()>4./Math.max(domain[mut_domain].length,10)){
								continue; //Modify 10 bases at a time.
							}
							System.err.println("Allowing initial sequence override at position "+j+" in domain "+mut_domain);
						}
						//Otherwise, it's free to mutate.
					}
					if (!(mutator instanceof CodonCode)){ //no codon table: single base modifications
						mutator.mutateToOther(domain,mut_domain,j);
					} else {
						for(int i = 0; i < 3; i++){
							if (Std.monomer.noFlags(domain[mut_domain][j-(j%3)+i])==0){
								throw new RuntimeException("Protein Conserving (-p) Flag requires an initial sequence");
							}
						}
						mutator.mutateToOther(domain,mut_domain,j-(j%3));
					}
					if (Std.monomer.noFlags(domain[mut_domain][j])==0){

						for(int i : domain[mut_domain]){
							if (Std.monomer.noFlags(i)==0){
								System.err.print("?");
							} else {
								System.err.print(Std.monomer.displayBase(i));
							}
						}
						System.err.println();
						
						//Mutation failed to replace a 0. rules must have been too strict; try again from the start.
						continue initialLoop;
					}
				}
			}
			if (!mutator.isValid(domain,mut_domain)){
				/*
				for(int i : domain[mut_domain]){
					if (Std.monomer.noFlags(i)==0){
						System.err.print("?");
					} else {
						System.err.print(Std.monomer.displayBase(i));
					}
				}
				System.err.println();
				*/
				System.err.println("Sequence constraints invalidated; Trying again: "+k);
				
				continue initialLoop;
			}
			break;
		}
		
		//Ok, all quotas are satisfied. Make a randomizing pass to spread out artificial stretches
		//of single base:

		for(int j = domain[mut_domain].length-1; j >=0; j--){
			if (!preserveInitial[j]){ //Still don't overwrite initial sequence
				if (!(mutator instanceof CodonCode)){ //no codon table: single base modifications
					mutator.mutateToOther(domain,mut_domain,j);
				} else { 
					mutator.mutateToOther(domain,mut_domain,j-(j%3));
				}
			}
		}
	}


	/**
	 * Used for sorting bases we want to mutate in a max heap by priority
	 */
	private Priority_int_int[] inplacePrioritySort_shared;
	private static final class Priority_int_int implements Comparable<Priority_int_int>{
		public int a;
		public int b;
		public int compareTo(Priority_int_int o) {
			int val = o.b-b;
			if (val==0){
				return Math.random()<.5 ? -1:1;
			}
			return val;
		}
		public void set(int k, int mark) {
			a = k;
			b = mark;
		}
	}
	/**
	 * Mutates mut_domain such that it is changed at the nucleic acid sequence level. If cc is nonnull, the changes will conserve amino acid sequence. 
	 * @param seqToSynthesize 
	 * @param domain_markings 
	 * @param mutators 
	 * @param max_mutations 
	 */
	public int mutateUntilValid(int mut_domain, int[][] domain,
			int[][] domain_markings, DesignerCode mutator, int min_mutations, int max_mutations) {
		int len = domain[mut_domain].length;

		//Count the bases to mutate (1's in the "markings" array)
		int oneC = 0;
		
		//How many mutations will we try now?
		int num_mut;
		//Count the reccomended bases to mutate
		
		if (ENABLE_MARKINGS){
			if (inplacePrioritySort_shared==null || inplacePrioritySort_shared.length<len){
				inplacePrioritySort_shared = new Priority_int_int[len]; 
				for(int k = 0; k < inplacePrioritySort_shared.length; k++){
					inplacePrioritySort_shared[k] = new Priority_int_int();
				}
			}
			for(int k = 0; k < len; k++){
				inplacePrioritySort_shared[oneC].set(k,DNAMARKER_DONTMUTATE);
				int q = domain_markings[mut_domain][k];
				if (q!=DNAMARKER_DONTMUTATE || SORT_MARKINGS){
					inplacePrioritySort_shared[oneC].set(k,q);
					oneC++;
				}
			}
			if (SORT_MARKINGS) {
				num_mut = oneC;
				Arrays.sort(inplacePrioritySort_shared,0,oneC);
			} else {
				num_mut = oneC;
			}
			//Post condition: num_mut <= oneC
		} else {
			 num_mut = len*2;
		}
		num_mut = Math.min(num_mut,max_mutations);
		
		//Randomize number to mutate:
		if (num_mut > 0){
			num_mut = int_urn(min_mutations,num_mut);
		} else {
			//Was not able to mutate.
			return 0;
		}
		
		//System.out.println(num_mut+" "+min_mutations);

		int successful = 0;
		//int timesInLoop = 0;
		ArrayList<Integer> positions = new ArrayList();
		
		mutator.reset();
		mutloop: for (int k = 0; k < num_mut; k++) {
			// Attempt a mutation
			int j;
			if (ENABLE_MARKINGS){
				if (SORT_MARKINGS){
					if (k >= oneC){
						break mutloop;
					}
					j = inplacePrioritySort_shared[k].a;
				} else {
					j = inplacePrioritySort_shared[int_urn(0, oneC-1)].a;
					/*
					if (domain_markings[mut_domain][j]==DNAMARKER_DONTMUTATE){
						throw new RuntimeException("Assertion error");
					}
					*/
				}
			} else {
				j = int_urn(0, len-1);
			}
			/*
			timesInLoop++;
			if (timesInLoop > 1000){
				System.err.println("Warning. Could not mutate in 1000 guesses. Breaking");
				break;
			}
			*/
			if (!(mutator instanceof CodonCode)){ //no codon table: single base modifications
				if (!mutator.mutateToOther(domain, mut_domain, j)){
					num_mut++;
					continue;
				}
				successful ++;
			} else { //Codon table: 3 mutation points = 1 amino swap.
				//Round to codon
				j -= j%3;

				//Protein domains are multiples of 3 long, so this is fine.
				//Mutate codon!
				if (!mutator.mutateToOther(domain, mut_domain, j)){
					num_mut++;
					continue;
				}
				
				//+3.
				k+=2;
				successful += 3;
			}
		}
		return successful;
	}

	private void printDebugPenalties(){
		System.out.println("Penalty summary:");
		for(int k = 0; k < debug_keepMessages.length; k++){
			if (debug_keepMessages[k]==0.0) continue;
			System.out.println(debug_keepMessages[k]+" "+debug_messages[k]);
		}
	}
	public void resetDebugPenalty(){
		for(int k = 0; k < debug_keepMessages.length; k++){
			debug_keepMessages[k] = 0;
			debug_messages[k] = "";
		}
	}
	private static final int debugOutTopScores = 50;
	private double[] debug_keepMessages = new double[debugOutTopScores];
	private String[] debug_messages = new String[debugOutTopScores];


	public void debugPenalty(ScorePenalty s, int[][] domain, DomainDefinitions dsd) {
		double val = s.cur_score;
		for(int k = 0; k < debug_keepMessages.length; k++){
			if (k == debug_keepMessages.length-1 || val < debug_keepMessages[k+1]){
				for(int y = 0; y < k; y++){
					debug_keepMessages[y] = debug_keepMessages[y+1];
					debug_messages[y] = debug_messages[y+1];
				}
				debug_messages[k] = s.toString();
				for(DomainSequence q : s.getSeqs()){
					debug_messages[k] += " "+displaySequence(q,domain,dsd);
				}
				debug_keepMessages[k] = val;
				break;
			}
		}
	}

	
	//x1: locked. x2: GCLG / GCLC flag.
	
	public int getLockedBase(char charAt) {
		return Std.monomer.noFlags(Std.monomer.decodeConstraintChar(charAt))+Std.monomer.LOCK_FLAG();
	}
	
	private String[] displayDomains(int[][] domain, boolean toOutput, DomainDefinitions dsd) {
		String[] outputDomains = new String[domain.length];
		for(int i = 0; i < domain.length; i++){
			int[] row = domain[i];
			outputDomains[i] = "";
			for(int q : row){
				outputDomains[i]+=Std.monomer.displayBase(q);
			}
			if (toOutput){
				System.out.println(dsd.getDomainName(i)+">");
				System.out.println(wrap(outputDomains[i],200));
			}
		}
		if (toOutput){
			System.out.println();
		}
		return outputDomains;
	}
	
	/**
	 * String wrapping. Used for outputting large strings.
	 */
	private static String wrap(String q, int i){
		StringBuffer out = new StringBuffer();
		for(int d = 0; d < q.length(); d+=i){
			int end = Math.min(d+i,q.length());
			out.append(q,d,end);
			if (d+i < q.length()){
				out.append("\n");
			}
		}
		return out.toString();
	}
	private String displaySequence(DomainSequence sequence, int[][] domain, DomainDefinitions dsd) {
		StringBuffer sb = new StringBuffer();
		sb.append(sequence.toString(dsd));
		
		int len = sequence.length(domain);
		sb.append(" = ");
		for(int i = 0; i < len; i++){
			sb.append(Std.monomer.displayBase(sequence.base(i, domain, Std.monomer)));
		}
		
		return sb.toString();
	}
	public static int int_urn(int from, int to) {
		return (int)(Math.random()*(to-from+1)+from);
	}
}