package DnaDesign;

import static DnaDesign.DesignSequenceConstraints.GCL_FLAG;
import static DnaDesign.DesignSequenceConstraints.LOCK_FLAG;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DnaDefinition.DisplayBase;
import static DnaDesign.DnaDefinition.noFlags;
import static DnaDesign.DomainSequence.DNAMARKER_DONTMUTATE;
import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

import DnaDesign.AbstractDesigner.PopulationDesignMember;
import DnaDesign.DesignIntermediateReporter.DesignIntermediateScore;
import DnaDesign.impl.CodonCode;
import DnaDesign.impl.DomainDesignBlockDesignerImpl;
import DnaDesign.impl.DomainDesignPMemberImpl;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;
import DnaDesign.impl.SequenceCode;

/**
 * The DomainDesigner engineers the input Domains (Defined in Domain Defs) to optimize a set of
 * constraint functions. Constraints on what possible dna sequences are valid is possible.
 */
public abstract class DomainDesigner {
	/**
	 * Utility method for creating a designer with the default parameters and input from the GUI layer.
	 */
	public static DDSeqDesigner<DesignerOptions> getDefaultDesigner(ArrayList<String> inputStrands,String domainDefs) {
		return new DomainDesigner_Public(inputStrands,domainDefs,null,null,null);
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
		/**
		 * Decodes the input, and then calls Main.
		 */
		private Runnable runOnStart = new Runnable(){
			public void run(){			
				final ArrayList<DomainSequence> MustBeVanilla = new ArrayList();
				final ArrayList<DomainSequence[]> hairpinLoops = new ArrayList();
				final ArrayList<DomainSequence> hairpinInnards = new ArrayList();
				int num_domain = -1;
				TreeMap<Integer, Integer> domain_length_t = new TreeMap<Integer, Integer>();
				final DomainStructureData dsd = new DomainStructureData();
				for(int kq = 0; kq < inputStrands.size(); kq++){
					String q = inputStrands.get(kq);
					//DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
					DomainStructureData.readDomainDefs(domainDefsBlock, dsd);
					String moleculeName = q.split("\\s+")[0];
					String inputStrand = q.split("\\s+")[1];
					DomainStructureData.readStructure(moleculeName, inputStrand, dsd);
					if (dsd.structures==null){
						throw new RuntimeException("Input strand invalid "+inputStrand);
					}
					DomainDesigner_SharedUtils.utilVanillaTargetFinder(dsd, MustBeVanilla);
					DomainDesigner_SharedUtils.utilHairpinInternalsFinder(dsd, hairpinInnards);
					DomainDesigner_SharedUtils.utilHairpinClosingFinder(dsd, hairpinLoops);

					num_domain = Math.max(num_domain,dsd.domainLengths.length-1);

					for(int i = 0; i < dsd.domainLengths.length; i++){
						int val = dsd.domainLengths[i];
						if (val!=-1){
							domain_length_t.put(i, val);
						}
					}
				}
				if (num_domain==-1){
					throw new RuntimeException("No valid molecules to design.");
				}
				
				final int num_domain_2 = num_domain+1;
				final int[] domain_length = new int[num_domain_2];

				if (initial==null){
					initial = new TreeMap<Integer, String>();
				}
				if (mutators == null){
					mutators = new TreeMap<Integer, DesignerCode>();
				}
				for(int k = 0; k < domain_length.length; k++){
					Integer got = domain_length_t.get(k);
					if (got==null){
						throw new RuntimeException("No length data for domain \'"+dsd.getDomainName(k)+"\'");
					}
					domain_length[k] = got;
					
					//Constraints are "initial sequences".
					initial.put(k,dsd.getConstraint(k));
					if (dsd.maintainAminos(k)){
						mutators.put(k, new CodonCode());
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
							results.add(r.main(num_domain_2, domain_length, Integer.MAX_VALUE, lock, initial, MustBeVanilla,mutators,hairpinLoops,hairpinInnards,dir, dsd));
						} catch (Throwable e){
							e.printStackTrace();
							errorResult = e.toString()+"";//ensure nonnull
							break;
						}
					}
					System.out.println(results);
					finished = true;
				}
				}.start();
			}
		};
		DomainDesigner r;
		private List<String> inputStrands;
		private String domainDefsBlock;
		private Map<Integer, String> lock;
		private Map<Integer, String> initial;
		private Map<Integer, DesignerCode> mutators;
		private DesignIntermediateReporter dir;
		public DomainDesigner_Public(List<String> inputStrands, String domainDefsBlock, Map<Integer, String> lock, Map<Integer, String> initial, Map<Integer, DesignerCode> mutators){
			this.inputStrands = inputStrands;
			this.domainDefsBlock = domainDefsBlock;
			this.lock = lock;
			this.initial = initial;
			this.mutators = mutators;
			dir = new DesignIntermediateReporter();
			
			r = new DomainDesignerImpl(new FoldingImpl());
			runOnStart.run();
		}
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

			ArrayList<DomainSequence> singleStrands = new ArrayList();
			DomainStructureData dsd = new DomainStructureData();
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
				sb.append(r.outputDomains[k]+"");
				sb.append(lR);
			}
			sb.append("----------------");
			sb.append(lR);
			sb.append("Input Strands:");
			sb.append(lR);	
			for(String q : inputStrands){
				String a = q.split("\\s+")[1]; //Erase the name
				splitLoop: for(String subStrand : a.split("}")){
					DomainSequence ds = new DomainSequence();
					ds.setDomains(subStrand,dsd);
					for(DomainSequence g : singleStrands){
						if (g.equals(ds)){
							continue splitLoop;
						}
					}
					sb.append("["+subStrand.replaceAll("\\s+","").replace("[","")+"}");
					sb.append(lR);
					singleStrands.add(ds);
					for(boolean withSeperator : new boolean[]{true,false}){
						if (withSeperator)
							sb.append("[");
						for(int k = 0; k < ds.domainList.length; k++){
							String domain = r.outputDomains[ds.domainList[k] & DNA_SEQ_FLAGSINVERSE];
							if ((ds.domainList[k] & DNA_COMPLEMENT_FLAG)!=0){
								domain = revComp(domain);
							}
							sb.append(domain);
							if (withSeperator){
								if (k + 1 < ds.domainList.length){
									sb.append("|");
								}
							}
						}
						if (withSeperator)
							sb.append("}");
						sb.append(lR);
					}
				}
			}

			return sb.toString();
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
			if (waitOnStart){
				waitOnStart = false;
			} else {
				r.waitForResume = false;
			}
		}

		public float scoreVal() {
			return (float)r.best_score;
		}

		public float statusVal() {
			return 0;
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
	}
	
	private DesignerOptions options = DesignerOptions.getDefaultOptions();

	//Output of designer. The best score, and the domains corresponding to it.
	private double best_score;
	private String[] outputDomains;
	private boolean waitForResume = false, abort = false;

	//End accessible to client.

	public static int decodeConstraintChar(char charAt){
		if (!Character.isLetter(charAt)){
			return 0;
		}
		int ret = DnaDefinition.decodeBaseChar(charAt);
		if (Character.isLowerCase(charAt)){
			ret += LOCK_FLAG;
		};
		return ret;
}
	
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
	 * How many population members to maintain if using the "Block Algorithm"
	 */
	private int populationSize = 30;
	
	/**
	 * If we are sorting markings, then markings that occur in multiple penalties are mutated first.
	 */
	private boolean SORT_MARKINGS = false;

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
		public abstract boolean affectedBy(int domain);
		public abstract DomainSequence[] getSeqs();
		public abstract int getNumDomainsInvolved();
		public abstract int getPriority();
	}
	
	/**
	 * The program will attempt to prevent any interaction between the DNASequences in toSynthesize, though interactions
	 * containing complementary sequences will be ignored (i.e, 1|4*|5 will not be tested against 1|4|5).
	 * @param initial 
	 * @param hairpinInnards 
	 * @param dsd 
	 */
	int main(int num_domain, int[] domain_length, int TOTAL_ATTEMPTS, Map<Integer, String> lockedDomains, Map<Integer, String> initial, List<DomainSequence> makeSingleStranded, Map<Integer, DesignerCode> mutators, List<DomainSequence[]> hairpinLoops, ArrayList<DomainSequence> hairpinInnards, DesignIntermediateReporter DIR, DomainStructureData dsd) {
		//SANITY OF INPUT:
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(makeSingleStranded);
		hairpinInnards.addAll(makeSingleStranded);
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(hairpinInnards);
		ArrayList<DomainSequence> preventComplementarity = hairpinInnards;
		
		/*
		System.out.println(seqToSynthesize);
		for(DomainSequence[] ds : hairpinLoops){
			System.out.print(Arrays.toString(ds));
		}
		System.out.println();
		*/
		//DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(hairpinLoops);
		
		//System.out.println(seqToSynthesize.size());

		
		//Domains to be designed. The integer type is being used to hold DNA bases,
		//with additional metadata (multiples of 10 are added to mean certain things)
		int[][] domain = new int[num_domain][];
		int[][] domain_markings = new int[num_domain][];
		DesignerCode[] mutate = new DesignerCode[num_domain];
		for(int i = 0; i < num_domain; i++){
			int[] cDomain = new int[domain_length[i]];
			domain[i] = cDomain; //initialize to 0s
			domain_markings[i] = new int[domain_length[i]];
			if (mutators!=null){
				mutate[i] = mutators.get(i);
				if (mutate[i]==null){
					SequenceCode newCode = SequenceCode.getDefaultSequenceCode();
					DesignSequenceConstraints dsc = DesignSequenceConstraints.getDefaultConstraints();
					for(int base = 0; base < DNAFLAG_ADD; base++){
						//Build the constraints by querying the DSD. Mix Default settings with DSD specification.
						final int maxComponent = dsd.getMaxComponent(i,base);
						final int minComponent = dsd.getMinComponent(i,base);
						if (maxComponent!=-2){
							dsc.setMaxConstraint(base, maxComponent);
						}
						if (minComponent!=-2){
							dsc.setMinConstraint(base, minComponent);
						}
					}
					newCode.setConstraints(dsc);
					mutate[i] = newCode;
				}
			}
		}

		//Initial domains, introducing constraints.
		if (initial!=null){
			for(Entry<Integer, String> init : initial.entrySet()){
				int i = init.getKey();
				for (int j = 0; j < domain_length[i]; j++) {
					domain[i][j] = decodeConstraintChar(init.getValue().charAt(j));
				}
			}
		}
		
		ArrayList<Integer> mutableDomainsL = new ArrayList();
		for(int i = 0; i < num_domain; i++){
			boolean mutable = false;
			//Optimize for completely locked domains.
			loop:for (int j = 0; j < domain_length[i]; j++) {
				if (domain[i][j] < DNAFLAG_ADD){
					mutable = true;
					break loop;
				}
			}
			if (mutable){
				mutableDomainsL.add(i);		
			}
		}
		
		//Unconstrained initialization.
		for(int i = 0; i < num_domain; i++){
			if (options.rule_ccend_option.getState()){
				if (domain[i][0]==0)domain[i][0] = C + GCL_FLAG;
				if (domain[i][domain[i].length-1]==0)domain[i][domain[i].length-1] = C + GCL_FLAG;
			}
			//Will initialize unconstrained portion of domain
			try {
				pickInitialSequence(domain,i,mutate[i]);
			} catch (Throwable e){
				throw new RuntimeException("Initial conditions were too strict for domain "+dsd.getDomainName(i));
			}
		}
		
		//At this point, no domain should have any bases other than ((x%10)%5)>1.
		
		//Locked domains feature.
		if (lockedDomains!=null){
			for(Entry<Integer, String> lock : lockedDomains.entrySet()){
				int i;
				mutableDomainsL.remove(i = lock.getKey());
				for (int j = 0; j < domain_length[i]; j++) {
					domain[i][j] = getLockedBase(lock.getValue().charAt(j));
				}
			}
		}
		
		//Combine locked domains...
		int[] mutableDomains = new int[mutableDomainsL.size()];
		for(int i = 0; i < mutableDomains.length; i++){
			mutableDomains[i] = mutableDomainsL.get(i);
		}

		// Initial display of starting position
		displayDomains(domain);
		
		/*
		//Is this a valid start position?
		for(i = 0; i < num_domain; i++){
			if (!isValidSequenceSetup(i,seqToSynthesize,domain)){
				throw new RuntimeException("Assertion error: Start setup invalidated rules.");
			}
		}
		*/

		//List all score calculations.
		List<ScorePenalty> allScores = listPenalties(makeSingleStranded,preventComplementarity,hairpinLoops,DIR);
		System.out.println("Checking "+allScores.size()+" score elements");

		if (allScores.isEmpty()){
			throw new RuntimeException("No scores to optimize : Algorithm has nothing to do!");
		}
				
		//Stop condition. See the scoring functions themselves for the actual threshold.
		long lastDumpState = System.nanoTime();
		
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

		//Initial score.
		float current_score = 0;
		deepFill(domain_markings,DNAMARKER_DONTMUTATE);
		DIR.beginScoreReport();{
			for(ScorePenalty q : allScores){
				current_score += q.getScore(domain, domain_markings);
			}
		}DIR.endScoreReport();
		System.out.println("Initial score: "+current_score);	
		best_score = current_score;
		
		//Create block designer
		DomainDesignPMemberImpl initialSeed = new DomainDesignPMemberImpl(allScores,scoredElements,domain,domain_markings);
		DomainDesignBlockDesignerImpl dbesign = new DomainDesignBlockDesignerImpl(num_domain,domain_length,mutableDomains,mutate,dsd,this);
		dbesign.initialize(initialSeed, populationSize);
		
		if (false){ //The children, if this block is enabled, will start at a different position than the seed. 
			PopulationDesignMember<DomainDesignPMemberImpl>[] a = dbesign.getPopulation();
			for(PopulationDesignMember r2 : a){
				DomainDesignPMemberImpl r = (DomainDesignPMemberImpl)r2;
				for(int k = 0; k < r.domain.length; k++){
					pickInitialSequence(r.domain,k,mutate[k]);
				}
				//Initial score.
				float current_score_in = 0;
				deepFill(r.domain_markings,DNAMARKER_DONTMUTATE);
				DIR.beginScoreReport();{
					for(ScorePenalty q : r.penalties){
						current_score_in += q.getScore(r.domain, r.domain_markings);
					}
				}DIR.endScoreReport();
				System.out.println("Initial score: "+current_score_in);			
			}
		}
		
		
		int num_mut_attempts = 0;
		while (!abort) {
			double endingThreshold = options.end_score_threshold.getState();
			if (best_score <= endingThreshold){
				break;
			}
			dbesign.runBlockIteration(endingThreshold);
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
			displayDomains(q.domain,false);
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
			while (waitForResume && !abort){
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		System.out.println("Designer ended after "+num_mut_attempts+" iterations, with a score of "+best_score);
		displayDomains(domain);
		return num_mut_attempts;
	}
	public abstract List<ScorePenalty> listPenalties(
			List<DomainSequence> makeSingleStranded,
			ArrayList<DomainSequence> preventComplementarity,
			List<DomainSequence[]> hairpinLoops, DesignIntermediateReporter DIR) ;

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
				throw new RuntimeException("Initial constraints for domain: "+mut_domain+ " were too strict");
			}
			for(int j = 0; j < domain[mut_domain].length; j++){
				mutator.mutateToOther(domain,mut_domain,j);
				if (noFlags(domain[mut_domain][j])==0){
					//Mutation failed to replace a 0. rules must have been too strict; try again from the start.
					continue initialLoop;
				}
			}
			break;
		}
		
		/*
		int lessATthanGC = 0;
		for(int i = 0; i < domain.length; i++){
			if (domain[i] >= DNAFLAG_ADD){
				lessATthanGC+=isGClessAT(domain[i]);
			}
		}
		for(int i = 0; i < domain.length; i++){
			int newBase = 0;
			if (i > 0){
				newBase = (((domain[i-1]%DNAFLAG_ADD))%2); //Ensure alternating
			}
			if (Math.abs(lessATthanGC)>0){
				newBase = newBase==0?T:A;
			} else { 
				newBase = newBase==0?G:C;
			}
			if (domain[i]==0){
				domain[i] = newBase;

				lessATthanGC+=isGClessAT(domain[i]);
			}
		}
		*/
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

		//Count the ones
		int oneC = 0;
		
		//How many mutations will we try now?
		int num_mut = Math.min(len,max_mutations);
		//Count the reccomended bases to mutate
		if (ENABLE_MARKINGS){
			if (SORT_MARKINGS){
				if (inplacePrioritySort_shared==null || inplacePrioritySort_shared.length<len){
					inplacePrioritySort_shared = new Priority_int_int[len]; 
					for(int k = 0; k < inplacePrioritySort_shared.length; k++){
						inplacePrioritySort_shared[k] = new Priority_int_int();
					}
				}
			}
			for(int k = 0; k < len; k++){
				int q = domain_markings[mut_domain][k];
				if (q!=DNAMARKER_DONTMUTATE){
					if (SORT_MARKINGS) inplacePrioritySort_shared[oneC].set(k,q);
					oneC++;
				}
			}
			num_mut = Math.min(num_mut,oneC);
			if (SORT_MARKINGS) Arrays.sort(inplacePrioritySort_shared,0,oneC);
			//Post condition: num_mut <= oneC
		}
		//Randomize choice:
		if (num_mut > 0){
			num_mut = int_urn(min_mutations,num_mut);
		} else {
			//Was not able to mutate.
			return false;
		}
		
		//System.out.println(num_mut+" "+min_mutations);

		int timesInLoop = 0;
		for (int k = 0; k < num_mut; k++) {
			// Attempt a mutation
			int j;
			j = int_urn(0, len-1);
			if (ENABLE_MARKINGS){
				if (SORT_MARKINGS){
					if (k >= inplacePrioritySort_shared.length){
						break;
					}
					j = inplacePrioritySort_shared[k].a;
				} else {
					if (domain_markings[mut_domain][j]==DNAMARKER_DONTMUTATE){
						num_mut++;
						continue;
					}
				}
			}
			timesInLoop++;
			if (timesInLoop > 1000){
				System.err.println("Warning. Could not mutate in 1000 guesses. Breaking");
				break;
			}
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

	void displayDomains(int[][] domain) {
		displayDomains(domain,true);
	}
	void displayDomains(int[][] domain, boolean toOutput) {
		String[] outputDomains = new String[domain.length];
		for(int i = 0; i < domain.length; i++){
			int[] row = domain[i];
			outputDomains[i] = "";
			for(int q : row){
				outputDomains[i]+=DisplayBase(q);
			}
			if (toOutput){
				System.out.print(wrap(outputDomains[i],200)+" ");
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
	
	public static int getLockedBase(char charAt) {
		return decodeConstraintChar(charAt)%DNAFLAG_ADD+DNAFLAG_ADD*1;
	}
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
	private static String displaySequence(DomainSequence sequence, int[][] domain, DomainStructureData dsd) {
		StringBuffer sb = new StringBuffer();
		sb.append(sequence.toString(dsd));
		
		int len = sequence.length(domain);
		sb.append(" = ");
		for(int i = 0; i < len; i++){
			sb.append(DisplayBase(sequence.base(i, domain)));
		}
		
		return sb.toString();
	}
	public static int int_urn(int from, int to) {
		return (int)(Math.random()*(to-from+1)+from);
	}
}