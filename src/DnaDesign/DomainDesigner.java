package DnaDesign;

import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.DomainSequence.DNAMARKER_DONTMUTATE;
import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

import DnaDesign.AbstractDesigner.BlockDesigner;
import DnaDesign.AbstractDesigner.InfiniteResourceTournament;
import DnaDesign.AbstractDesigner.PopulationDesignMember;
import DnaDesign.AbstractDesigner.StandardTournament;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.Config.CircDesigNASystemElement;
import DnaDesign.DesignIntermediateReporter.DesignIntermediateScore;
import DnaDesign.impl.CodonCode;
import DnaDesign.impl.DomainDesignBlockDesignerImpl;
import DnaDesign.impl.DomainDesignPMemberImpl;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
import DnaDesign.impl.SequenceCode;

/**
 * The DomainDesigner designs a set of sequences (split into "domains", by the parameters in Domain Defs) 
 * to optimize a set of score functions over a sequence space which is itself a constrained subset of all possible dna
 * sequences. 
 * 
 * This class
 * - Provides a useful frontend to the GUI so that it can display real time design information
 * - Provides a general purpose main() method which defers to subclasses for the design routines.
 * 
 * Subclasses specify certain behaviors of the functions in main(), thus actually producing a full design algorithm.
 */
public abstract class DomainDesigner extends CircDesigNASystemElement{
	public DomainDesigner(CircDesigNAConfig std) {
		super(std);
	}
	/**
	 * Utility method for creating a designer with the default parameters and input from the GUI layer.
	 */
	public static DDSeqDesigner<DesignerOptions> getDefaultDesigner(String Molecules,String domainDefs,CircDesigNAConfig Std) {
		ArrayList<String> inputStrands = new ArrayList<String>();
		for(String q : Molecules.split("\n")){
			String[] line = q.split("\\s+");
			if (line.length == 0){
				continue;
			}
			if (line.length != 2){
				throw new RuntimeException("Correct Molecule format: <name> <molecule> @ "+line[0]);
			}
			inputStrands.add(q);
		}
		
		return new DomainDesigner_Public(inputStrands,domainDefs,null,null,null,Std);
	}
	/**
	 * A coupled, exposed wrapper which allows the GUI to easily display design internals.
	 * 
	 * Why does this class exist? The domainDesigner is WAY too flexible, and has a lot of options.
	 * To make client code simpler, this class has some adaptor methods (a constructor) that will 
	 * prepare your input to be valid to the DomainDesigner. 
	 * 
	 * Error checking is done, in the sense that some combinations of inputs are not valid. For the most
	 * part, however, input sanity is done in the DomainDesigner class itself.
	 * 
	 * Finally, this class serves to separate the public api which controls the behavior of the designer (run,
	 * pause abort) from the internal api of how to do so (mutate, run block iteration).
	 */
	public static class DomainDesigner_Public implements DDSeqDesigner<DesignerOptions>{
		private boolean waitOnStart = true;
		private boolean finished = false;
		
		DomainDesigner r;
		private List<String> inputMolecules;
		private String domainDefsBlock;
		private Map<Integer, String> lock;
		private Map<Integer, String> initial;
		private Map<Integer, DesignerCode> mutators;
		private DesignIntermediateReporter dir;
		public DomainDesigner_Public(List<String> inputStrands, String domainDefsBlock, Map<Integer, String> lock, Map<Integer, String> initial, Map<Integer, DesignerCode> mutators, CircDesigNAConfig Std){
			this.inputMolecules = inputStrands;
			this.domainDefsBlock = domainDefsBlock;
			this.lock = lock;
			this.initial = initial;
			this.mutators = mutators;
			dir = new DesignIntermediateReporter();
			
			r = new DomainDesignerImpl(new FoldingImpl(Std), Std);
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
				final DomainStructureData dsd = new DomainStructureData(r.Std);
				DomainStructureData.readDomainDefs(domainDefsBlock, dsd);
				num_domain = dsd.domainLengths.length;
				
				//Read molecules
				final AbstractDomainDesignTarget designTarget = new AbstractDomainDesignTarget(dsd,r.Std);
				for(int kq = 0; kq < inputMolecules.size(); kq++){
					String q = inputMolecules.get(kq);
					//DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
					String moleculeName = q.split("\\s+")[0];
					String inputStrand = q.split("\\s+")[1];
					designTarget.addTargetStructure(moleculeName,inputStrand);
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
				
				for(int k = 0; k < domain_length.length; k++){
					Integer got = domain_length_t.get(k);
					if (got==null){
						throw new RuntimeException("No length data for domain \'"+dsd.getDomainName(k)+"\'");
					}
					domain_length[k] = got;
					
					//Constraints are "initial sequences".
					initial.put(k,dsd.getConstraint(k));
					if (dsd.maintainAminos(k)){
						mutators.put(k, new CodonCode(r.Std));
					}
				}
				

				new Thread(){public void run(){
					while(waitOnStart && ! r.abort){
						try {
							Thread.sleep(100);
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
							results.add(r.main(num_domain_2, domain_length, Integer.MAX_VALUE, lock, initial, mutators, designTarget, dir, dsd));
						} catch (OutOfMemoryError e){
							errorResult = "Could not start designer: " +e.getMessage()+"\nJava ran out of memory! Strategies to fix:\n\t1) Reduce the population size.\n\t2)Download this applet, and run it with the argument -Xmx700 from the command line.";//ensure nonnull
						} catch (Throwable e){
							e.printStackTrace();
							errorResult = "Could not start designer: " +e.getMessage()+"";//ensure nonnull
							break;
						}
					}
					//System.out.println(results);
					finished = true;
				}}.start();
			}
		};
		
		public DesignerOptions getOptions() {
			return r.options;
		}

		private String errorResult;

		/**
		 * Returns a human-readable dump of the output.
		 */
		public String getResult() {
			if (errorResult!=null){
				return errorResult;
			}
			if (r.outputDomains==null){
				return "Output incomplete. Run designer longer";
			}

			ArrayList<DomainSequence> alreadyPrintedSequences = new ArrayList();
			DomainStructureData dsd = new DomainStructureData(r.Std);
			//DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
			DomainStructureData.readDomainDefs(domainDefsBlock, dsd);
			StringBuffer sb = new StringBuffer();
			String lR = "\n";
			sb.append("Current designer state (to resume from this point, paste as 'Domain Definition'):");
			sb.append(lR);
			sb.append("Net 'Score': "+r.best_score);
			sb.append(lR);
			sb.append("----------------");
			sb.append(lR);
			for(int k = 0; k < r.outputDomains.length; k++){
				sb.append(dsd.getDomainName(k));
				sb.append("\t");
				sb.append("[");
				sb.append(r.outputDomains[k]);
				sb.append("]");
				String arg = dsd.getArguments(k);
				if (arg.length()>0){
					sb.append("\t"+dsd.getArguments(k));	
				}
				sb.append(lR);
			}
			sb.append("----------------");
			sb.append(lR);
			sb.append("Molecular Substrands:");
			sb.append(lR);	
			for(String q : inputMolecules){
				String a[] = q.split("\\s+");    
				boolean alreadyPrintedInMolecule = false;
				splitLoop: for(String subStrand : a[1].split("}")){
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
					sb.append("\t"+"["+subStrand.replaceAll("\\s+","").replace("[","")+"}");
					sb.append(lR+"\t");
					printSequence(sb,r,ds,true);
					sb.append(lR+"\t");
					printSequence(sb,r,ds,false);
					sb.append(lR);
				}
			}

			return sb.toString();
		}
		private void printSequence(StringBuffer sb, DomainDesigner r2, DomainSequence ds, boolean withSeperator) {
			if (withSeperator)
				sb.append("[");
			for(int k = 0; k < ds.numDomains; k++){
				String domain = r.outputDomains[ds.domainList[k] & DNA_SEQ_FLAGSINVERSE];
				if ((ds.domainList[k] & DNA_COMPLEMENT_FLAG)!=0){
					domain = revComp(domain);
				}
				sb.append(domain);
				if (withSeperator && k + 1 < ds.domainList.length){
					sb.append("|");
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
		}
		public DesignIntermediateReporter getDir() {
			return dir;
		}
		public double getBestScore() {
			return r.best_score;
		}
		public int getIterationCount() {
			return r.num_mut_attempts;
		}
	}
	
	public DesignerOptions options = DesignerOptions.getDefaultOptions();

	//Output of designer. The best score, and the domains corresponding to it.
	private double best_score;
	private int num_mut_attempts = 0;
	private String[] outputDomains;
	private boolean waitForResume = false;

	public boolean abort = false;

	//End accessible to client.
	
	//Constructor
	
	
	/**
	 * Debug output flags
	 */
	public boolean DEBUG_PENALTIES = true;
	
	/**
	 * Use the evaluator to determine where the problem spots are, 
	 * and mutate them. Recommended, otherwise you're merely randomly mutating.
	 */
	private boolean ENABLE_MARKINGS = true;
	
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
		public ScorePenalty(DesignIntermediateReporter dir){
			old_score = cur_score = MAX_SCORE; //A suitably large number
			this.dir = dir;
		}
		public double getScore(int[][] domain, int[][] domain_markings){
			evalScore(domain, domain_markings);
			dedicate();
			return old_score;
		}
		public double old_score, cur_score;
		public DesignIntermediateScore dis;
		public DesignIntermediateReporter dir;
		public void revert(){
			cur_score = old_score;
		}
		public void dedicate(){
			old_score = cur_score;
		}
		public final double check(double possScore){
			return Math.min(possScore,MAX_SCORE);
		}
		public final double evalScore(int[][] domain, int[][] domain_markings){
			double new_score = check(evalScoreSub(domain, domain_markings));
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
	}
	
	/**
	 * The program will attempt to prevent any interaction between the DNASequences in toSynthesize, though interactions
	 * containing complementary sequences will be ignored (i.e, 1|4*|5 will not be tested against 1|4|5).
	 * @param initial 
	 * @param hairpinInnards 
	 * @param dsd 
	 */
	int main(int num_domain, int[] domain_length, int TOTAL_ATTEMPTS, Map<Integer, String> lockedDomains, Map<Integer, String> initial, Map<Integer, DesignerCode> mutators, AbstractDomainDesignTarget designTarget, DesignIntermediateReporter DIR, DomainStructureData dsd) {
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

		//Load the initial domains, with constraint characters.
		if (initial!=null){
			for(Entry<Integer, String> init : initial.entrySet()){
				int i = init.getKey();
				for (int j = 0; j < domain_length[i]; j++) {
					try {
						domain[i][j] = Std.monomer.decodeConstraintChar(init.getValue().charAt(j));
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
			loop:for (int j = 0; j < domain_length[i]; j++) {
				if (domain[i][j] - Std.monomer.noFlags(domain[i][j]) != Std.monomer.LOCK_FLAG()){
					mutable = true;
					break loop;
				}
			}
			if (mutable){
				mutableDomainsL.add(i);		
			}
		}

		//Derive a mutation strategy for all domains (by parsing domain constraints with the options selected)
		for(int i = 0; i < num_domain; i++){
			if (mutators!=null){
				mutate[i] = mutators.get(i);
				if (mutate[i]==null){
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

		//Using mutation strategy (which are aware of sequence space), pick valid starting sequences.
		for(int i = 0; i < num_domain; i++){
			if (options.rule_ccend_option.getState()){
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

		//Assertion: all domains are valid according to their respective mutation strategies. 

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

		/*
		//Check if this is a valid start position?
		for(i = 0; i < num_domain; i++){
			if (!isValidSequenceSetup(i,seqToSynthesize,domain)){
				throw new RuntimeException("Assertion error: Start setup invalidated rules.");
			}
		}
		 */

		//Enumerate penalty scores (see FoldingImplTestGUI for a visual of this process) via "listPenalties"
		List<ScorePenalty> allScores = listPenalties(designTarget,DIR,domain,options);		
		System.out.println("Discovered "+allScores.size()+" score elements");
		if (allScores.isEmpty()){
			throw new RuntimeException("No scores to optimize : Algorithm has nothing to do!");
		}


		/*
		for(int[] row : domain_markings){
			System.out.println(Arrays.toString(row));
		}
		 */


		//Select the penalties which depend on certain domains. Clear optimization here.
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
			//Initial score.
			double current_score = 0;
			deepFill(domain_markings,DNAMARKER_DONTMUTATE);
			DIR.beginScoreReport();{
				for(ScorePenalty q : allScores){
					current_score += q.getScore(domain, domain_markings);
					if (q.dis!=null){
						q.dis.addScore(q.old_score);
					}
					if (abort){
						return 0;
					}
				}
			}DIR.endScoreReport();
			System.out.println("Randomized initial sequence (for population member 0):");
			displayDomains(domain,true,dsd);
			System.out.println("Initial score: "+current_score);	
			best_score = current_score;

			//Create block designer, which will produce a certain initial population from the initial sequences we chose.
			DomainDesignPMemberImpl initialSeed = new DomainDesignPMemberImpl(allScores,scoredElements,domain,domain_markings);
			DomainDesignBlockDesignerImpl dbesignSingle = new DomainDesignBlockDesignerImpl(num_domain,domain_length,mutableDomains,mutate,dsd,this);
			BlockDesigner<DomainDesignPMemberImpl> dbesign;
			if (options.resourcePerMember.getState() < 0){
				dbesign = new InfiniteResourceTournament(dbesignSingle);
			} else {
				dbesign = new StandardTournament(dbesignSingle, options.resourcePerMember.getState());
			}
			dbesign.initialize(initialSeed, options.population_size.getState());

			if (true){  //If false, then all population members start at the same location.
				System.out.println("Randomizing "+dbesign.getPopulation().length+" - 1 population members.");
				//If true, seed each of the population members differently. 
				PopulationDesignMember<DomainDesignPMemberImpl>[] a = dbesign.getPopulation();
				for(int i = 1; i < a.length; i++){
					DomainDesignPMemberImpl r = (DomainDesignPMemberImpl)a[i];
					for(int k = 0; k < r.domain.length; k++){
						pickInitialSequence(r.domain,k,mutate[k]);
					}
					//Initial score.
					double current_score_pmember = 0;
					deepFill(r.domain_markings,DNAMARKER_DONTMUTATE);
					//DIR.beginScoreReport();{
					for(ScorePenalty q : r.penalties){
						current_score_pmember += q.getScore(r.domain, r.domain_markings);
					}
					//System.out.println("Initial score: "+current_score_in);
					//}DIR.endScoreReport();			
				}
			}


			//Record the time that the iteration began...
			System.out.println("Entering design loop");
			long lastDumpState = System.nanoTime();		
			num_mut_attempts = 0;
			while (!abort) {
				while (waitForResume && !abort){
					try {
						Thread.sleep(100);
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				double endingThreshold = options.end_score_threshold.getState();
				if (best_score <= endingThreshold){
					break;
				}
				dbesign.runBlockIteration(this,endingThreshold);
				System.out.flush();
				//Iteration complete.
				resetDebugPenalty();
				DomainDesignPMemberImpl q = dbesign.getBestPerformingChild();
				//Report the status of the best member.
				DIR.beginScoreReport();
				double score = 0;
				for(ScorePenalty s : q.penalties){
					score += s.old_score;
					if (s.dis!=null){
						s.dis.addScore(s.old_score);
					}
				}
				best_score = score;
				DIR.endScoreReport();
				displayDomains(q.domain,true,dsd);
				//debugPenalty(s, q.domain, dsd);
				//System.out.println(num_domains_mut);

				num_mut_attempts++;
				if (num_mut_attempts > TOTAL_ATTEMPTS){
					break;
				}
				if ((System.nanoTime()-lastDumpState)>2e9){
					lastDumpState = System.nanoTime();
					//if (DEBUG_PENALTIES){
					//System.out.println("iteration: "+num_mut_attempts+" score: "+best_score+" most significant subscore: "+(priority-1)+" "+(current_score+deltaScore));
					//}
				}
			}
		}	finally {
			System.out.println("Designed sequence output:");
			displayDomains(domain,true,dsd);
			System.out.println("Designer ended after "+num_mut_attempts+" iterations, with a score of "+best_score);
		}
		return num_mut_attempts;
	}
	public abstract List<ScorePenalty> listPenalties(
			AbstractDomainDesignTarget designTarget, DesignIntermediateReporter DIR, int[][] domain, DesignerOptions options2) ;

	public static final void deepFill(int[][] domain_markings, int i) {
		for(int[] row : domain_markings){
			Arrays.fill(row,i);
		}
	}
	/**
	 * Postcondition: any bases which were < 10 will be alternated to most ideally balance ATCG.
	 * @param designerCode 
	 */
	private void pickInitialSequence(int[][] domain, int mut_domain, DesignerCode mutator) {
		initialLoop: for(int k = 0; ; k++){
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
							continue;
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
				continue initialLoop;
			}
			break;
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
			int val = 0;
			if (val==0){
				val = o.b-b;
			}
			if (val==0){
				val = a-o.a;
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
	public boolean mutateUntilValid(int mut_domain, int[][] domain,
			int[] domain_length, int[][] domain_markings, DesignerCode mutator, int min_mutations, int max_mutations) {
		int len = domain_length[mut_domain];

		//Count the bases to mutate (1's in the "markings" array)
		int oneC = 0;
		
		//How many mutations will we try now?
		int num_mut;
		//Count the reccomended bases to mutate
		boolean SORT_MARKINGS = options.sort_markings.getState();
		
		if (ENABLE_MARKINGS){
			if (inplacePrioritySort_shared==null || inplacePrioritySort_shared.length<len){
				inplacePrioritySort_shared = new Priority_int_int[len]; 
				for(int k = 0; k < inplacePrioritySort_shared.length; k++){
					inplacePrioritySort_shared[k] = new Priority_int_int();
				}
			}
			for(int k = 0; k < len; k++){
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
			return false;
		}
		
		//System.out.println(num_mut+" "+min_mutations);

		//int timesInLoop = 0;
		for (int k = 0; k < num_mut; k++) {
			// Attempt a mutation
			int j;
			if (ENABLE_MARKINGS){
				if (SORT_MARKINGS){
					if (k >= inplacePrioritySort_shared.length){
						throw new RuntimeException("Assertion error");
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
				j = int_urn(0, oneC-1);
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
			}      
			
		}
		return true;
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


	public void debugPenalty(ScorePenalty s, int[][] domain, DomainStructureData dsd) {
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

	void displayDomains(int[][] domain, boolean toOutput, DomainStructureData dsd) {
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
		this.outputDomains = outputDomains;
		if (toOutput){
			System.out.println();
		}
	}
	
	String printf(String k){
		//System.out.print(k);
		return k;
	}
	//x1: locked. x2: GCLG / GCLC flag.
	
	public int getLockedBase(char charAt) {
		return Std.monomer.noFlags(Std.monomer.decodeConstraintChar(charAt))+Std.monomer.LOCK_FLAG();
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
	private String displaySequence(DomainSequence sequence, int[][] domain, DomainStructureData dsd) {
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