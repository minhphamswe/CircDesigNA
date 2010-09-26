package DnaDesign;

import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DnaDefinition.DisplayBase;
import static DnaDesign.DomainSequence.DNAMARKER_DONTMUTATE;
import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

import DnaDesign.DesignIntermediateReporter.DesignIntermediateScore;
import DnaDesign.impl.CodonCode;
import DnaDesign.impl.DomainDesignerImpl;
import DnaDesign.impl.FoldingImpl;

/**
 * The DomainDesigner engineers the input Domains (Defined in Domain Defs) to optimize a set of
 * constraint functions. Constraints on what possible dna sequences are valid is possible.
 */
public abstract class DomainDesigner {
	/**
	 * Utility method for creating a designer with the default parameters and input from the GUI layer.
	 */
	public static DDSeqDesigner getDefaultDesigner(ArrayList<String> inputStrands,String domainDefs) {
		return new DomainDesigner_Public(inputStrands,domainDefs,null,null,null);
	}
	/**
	 * A coupled, exposed wrapper which allows the GUI to easily display design internals.
	 */
	public static class DomainDesigner_Public implements DDSeqDesigner{
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

				/*
				for(DomainSequence q : hairpinInnards){
					System.out.println(q.toString(dsd));
				}
				*/
				
				//System.exit(0);

				final int num_domain_2 = num_domain+1;
				final int[] domain_length = new int[num_domain_2];

				if (initial==null){
					initial = new TreeMap<Integer, String>();
				}
				if (cc == null){
					cc = new TreeMap<Integer, CodonCode>();
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
						cc.put(k, new CodonCode());
					}
				}
				
				final int[] maxISO = dsd.getMaxISO();
				final int[] maxPZ = dsd.getMaxPZ();
				
				
				/**
				 * These often invalidate rules - allow this.
				if (!(initial==null || initial.isEmpty())){
					r.rule_SeqRulesAreAbsolute=0;
				}
				 */
				


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
						try {
							results.add(r.main(num_domain_2, domain_length, Integer.MAX_VALUE, lock, initial, MustBeVanilla,cc,hairpinLoops,hairpinInnards,dir, maxISO, maxPZ,dsd));
						} catch (Throwable e){
							e.printStackTrace();
							errorResult = e.toString();
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
		{
			r = new DomainDesignerImpl(new FoldingImpl());
		}
		private List<String> inputStrands;
		private String domainDefsBlock;
		private Map<Integer, String> lock;
		private Map<Integer, String> initial;
		private Map<Integer, CodonCode> cc;
		private DesignIntermediateReporter dir;
		public DomainDesigner_Public(List<String> inputStrands, String domainDefsBlock, Map<Integer, String> lock, Map<Integer, String> initial, Map<Integer,CodonCode> cc){
			this.inputStrands = inputStrands;
			this.domainDefsBlock = domainDefsBlock;
			this.lock = lock;
			this.initial = initial;
			this.cc = cc;
			options.add(handleJunctions);
			options.add(penalizeJunctions);
			options.add(performFullN4Junctions);
			options.add(rule_ccend_option);
			dir = new DesignIntermediateReporter();
			runOnStart.run();
		}
		public List<SeqDesignerOption> getOptions() {
			return options;
		}
		private ArrayList<SeqDesignerOption> options = new ArrayList();

		private SeqDesignerOption handleJunctions = new SeqDesignerOption(){
			public String getDescription() {
				return "Calculate Junction Penalties";
			}
			public boolean getState() {
				return r.CALCULATE_JUNCTION_PENALTIES;
			}
			public void toggle() {
				r.CALCULATE_JUNCTION_PENALTIES = !r.CALCULATE_JUNCTION_PENALTIES;
			}
		};

		private SeqDesignerOption penalizeJunctions = new SeqDesignerOption(){
			public String getDescription() {
				return "Scale Junction Penalties Very Low";
			}
			public boolean getState() {
				return !r.JunctionPenalty_FullPenalty;
			}
			public void toggle() {
				r.JunctionPenalty_FullPenalty = !r.JunctionPenalty_FullPenalty;
			}
		};



		private SeqDesignerOption performFullN4Junctions = new SeqDesignerOption(){
			public String getDescription() {
				return "Calculate every possible domain junctions (ignore input)";
			}
			public boolean getState() {
				return r.PERFORM_COMPLETE_JUNCTIONS;
			}
			public void toggle() {
				r.PERFORM_COMPLETE_JUNCTIONS = !r.PERFORM_COMPLETE_JUNCTIONS;
			}
		};

		private SeqDesignerOption rule_ccend_option = new SeqDesignerOption(){
			public String getDescription() {
				return "Force domains to begin / end with C";
			}
			public boolean getState() {
				return r.rule_ccend==1;
			}
			public void toggle() {
				r.rule_ccend = 1-r.rule_ccend;
			}
		};

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
	}

	//Accessible to client:
	private double best_score;
	private String[] outputDomains;
	private boolean waitForResume = false, abort = false;

	//End accessible to client.

	//Constructor
	
	
	/**
	 * Debug output flags
	 */
	private static final boolean OUTPUT_RUNTIME_STATS = false;
	private static final boolean DEBUG_PENALTIES = true;
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
	
	int MAX_MUTATION_DOMAINS = 999;
	int MAX_MUTATIONS_LIMIT = 9999; // maximunotm number of simultaneous mutations


	float INTERACTION_SCORE_MULT = .25f;
	boolean PERFORM_COMPLETE_JUNCTIONS = false;
	boolean CALCULATE_JUNCTION_PENALTIES = false;
	boolean JunctionPenalty_FullPenalty = true;
	
	public int rule_targetworst = 0; // target worst domains for mutation
	public int rule_ccend = 1; // domains MUST start and end with C
	public int rule_SeqRulesAreAbsolute = 1;
	//It is sometimes advantageous to design sequences with some immutable complementarity.
	public boolean ALLOW_COMPLEMENTARY_SCORES = false;

	protected abstract static class ScorePenalty{
		public static final double MAX_SCORE = 1e18;
		public ScorePenalty(DesignIntermediateReporter dir){
			old_score = cur_score = MAX_SCORE; //A suitably large number
		}
		public double getScore(int[][] domain, int[][] domain_markings){
			evalScore(domain, domain_markings);
			dedicate();
			return old_score;
		}
		public double old_score, cur_score;
		public DesignIntermediateScore dis;
		public void revert(){
			cur_score = old_score;
		}
		public void dedicate(){
			old_score = cur_score;
			if (dis!=null){
				dis.addScore(old_score);
			}
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
	 * @param maxPZ 
	 * @param maxISO 
	 * @param dsd 
	 */
	int main(int num_domain, int[] domain_length, int TOTAL_ATTEMPTS, Map<Integer, String> lockedDomains, Map<Integer, String> initial, List<DomainSequence> makeSingleStranded, Map<Integer, CodonCode> cc, List<DomainSequence[]> hairpinLoops, ArrayList<DomainSequence> hairpinInnards, DesignIntermediateReporter DIR, int[] maxISO, int[] maxPZ, DomainStructureData dsd) {
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
		int i,j,k;
		int[][] domain = new int[num_domain][];
		int[][] domain_markings = new int[num_domain][];
		CodonCode[] ccs = new CodonCode[num_domain];
		for(i = 0; i < num_domain; i++){
			int[] cDomain = new int[domain_length[i]];
			domain[i] = cDomain; //initialize to 0s
			domain_markings[i] = new int[domain_length[i]];
			if (cc!=null){
				ccs[i] = cc.get(i);
			}
		}

		//Initial domains, introducing constraints.
		if (initial!=null){
			for(Entry<Integer, String> init : initial.entrySet()){
				i = init.getKey();
				for (j = 0; j < domain_length[i]; j++) {
					domain[i][j] = decodeConstraintChar(init.getValue().charAt(j));
				}
			}
		}
		
		ArrayList<Integer> mutableDomainsL = new ArrayList();
		for(i = 0; i < num_domain; i++){
			boolean mutable = false;
			//Optimize for completely locked domains.
			loop:for (j = 0; j < domain_length[i]; j++) {
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
		for(i = 0; i < num_domain; i++){
			if (rule_ccend==1){
				if (domain[i][0]==0)domain[i][0] = GCLC;
				if (domain[i][domain[i].length-1]==0)domain[i][domain[i].length-1] = GCLC;
			}
			//Will initialize unconstrained portion of domain
			pickInitialSequence(domain[i]);
		}
		
		//At this point, no domain should have any bases other than ((x%10)%5)>1.
		
		//Locked domains feature.
		if (lockedDomains!=null){
			for(Entry<Integer, String> lock : lockedDomains.entrySet()){
				mutableDomainsL.remove(lock.getKey());
				for (j = 0; j < domain_length[i]; j++) {
					domain[i][j] = getLockedBase(lock.getValue().charAt(j));
				}
			}
		}
		
		//Combine locked domains...
		int[] mutableDomains = new int[mutableDomainsL.size()];
		for(i = 0; i < mutableDomains.length; i++){
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
		double EndThreshold = 0;
		if (EndThreshold < 0){
			 throw new RuntimeException("Assertion error: endingthreshold must be a positive score");
		}
		long lastDumpState = System.nanoTime();

		//Initial score.
		float current_score = 0;
		deepFill(domain_markings,DNAMARKER_DONTMUTATE);
		DIR.beginScoreReport();{
			for(ScorePenalty q : allScores){
				current_score += q.getScore(domain, domain_markings);
			}
		}DIR.endScoreReport();
		System.out.println("Initial score: "+current_score);
		
		/*
		for(int[] row : domain_markings){
			System.out.println(Arrays.toString(row));
		}
		*/

		//Select the penalties which depend on certain domains. Clear optimization here.
		ScorePenalty[][] scoredElements = new ScorePenalty[num_domain][];
		for(k = 0; k < num_domain; k++){
			ArrayList<ScorePenalty> affectedBySeq = new ArrayList();
			for(ScorePenalty q : allScores){
				if (q.affectedBy(k)){
					affectedBySeq.add(q);
				}
			}
			scoredElements[k] = new ScorePenalty[affectedBySeq.size()];
			affectedBySeq.toArray(scoredElements[k]);
		}

		//Make a back buffer for domain reverts
		int[][] old_domains = new int[num_domain][];
		int[][] old_domains_markings = new int[num_domain][];
		for(k = 0; k < num_domain; k++){
			old_domains[k] = new int[domain_length[k]];
			old_domains_markings[k] = new int[domain_length[k]];
		}
		//Blah blah blah, various temporaries. See usage below.
		int num_mut_attempts = 0, mut_domain, num_domains_mut, min_mutations, max_mutations;
		double worstPenaltyScore = -1, deltaScore = 0, deltaScoreSum = 0;
		double[] bestWorstPenaltyScore = new double[]{Double.MAX_VALUE, Double.MAX_VALUE, Double.MAX_VALUE};
		boolean revert_mutation = false, newPointReached = false;
		ScorePenalty worstPenalty = null;
		boolean[] domain_mutated = new boolean[num_domain];
		int[] mut_domains = new int[MAX_MUTATION_DOMAINS];

		double compareWithDelta = 0;
		int timesSameBeforeBump = 10;
		int timesSameCount = 0;
		double bumpAmount = 10;
		
		while (true && !abort) {

			//Mutate, in some way.

			//num_domains_mut = int_urn(1,Math.min(Math.min(worstPenalty==null?1:worstPenalty.getNumDomainsInvolved()-1,MAX_MUTATION_DOMAINS),mutableDomains.length-1));
			num_domains_mut = int_urn(1,Math.min(MAX_MUTATION_DOMAINS,mutableDomains.length));
			Arrays.fill(domain_mutated,false);
			int times = 0;
			for(k = 0; k < num_domains_mut; k++){
				min_mutations = 1;
				/*
				if (worstPenalty==null){
					mut_domain = int_urn(0, mutableDomains.length-1);
					mut_domain = mutableDomains[mut_domain];
				} else {
					
						mut_domain = int_urn(0, mutableDomains.length-1);
						mut_domain = mutableDomains[mut_domain];
					}while (!worstPenalty.affectedBy(mut_domain));
				}
				*/
				max_mutations = MAX_MUTATIONS_LIMIT;
				

				int mut_domain_search = int_urn(0, mutableDomains.length-1);
				mut_domain = 0;
				for(int o = 0; o < mutableDomains.length; o++){
					mut_domain = mutableDomains[(mut_domain_search+o)%mutableDomains.length];
					if (domain_mutated[mut_domain]){
						continue;
					}
					if (worstPenalty!=null){
						if (!worstPenalty.affectedBy(mut_domain)){
							if (Math.random()<.25f){
								continue; //We want to skew mutations towards the worst penalty.
							}
						}
					}
					break;
				}
				if(domain_mutated[mut_domain]){
					k--;
					continue;
				}
				
				/*
				if (timesSameCount > 0){
					max_mutations = Math.max(MAX_MUTATIONS_LIMIT/timesSameCount,min_mutations);
				}
				*/
				
				if (ENABLE_MARKINGS){
					boolean hasNo1 = true;
					for(int q : domain_markings[mut_domain]){
						if (q!=0){
							hasNo1 = false;
							break;
						}
					}
					if (hasNo1){
						k--;
						num_domains_mut--;
						continue;
					}
				}

				mut_domains[k] = mut_domain;
				//Backup
				System.arraycopy(domain[mut_domain],0,old_domains[mut_domain],0,domain_length[mut_domain]);
				System.arraycopy(domain_markings[mut_domain],0,old_domains_markings[mut_domain],0,domain_length[mut_domain]);
				//Mutate
				mutateUntilValid(mut_domain, domain, domain_length, makeSingleStranded, domain_markings, ccs[mut_domain], min_mutations, max_mutations, maxISO, maxPZ);
				domain_mutated[mut_domain] = true;
			}

			//Short circuiting constraint evaluation: If score was improved, keep, otherwise, break.
			//Penalties with higher priority are presumably harder to compute, thus the shortcircuiting advantage

			//Reset markers.
			for(k = 0; k < num_domains_mut; k++){
				mut_domain = mut_domains[k];
				Arrays.fill(domain_markings[mut_domain], DNAMARKER_DONTMUTATE);
			}
			deltaScoreSum = 0;
			int priority;
			
			priorityLoop: for(priority = 0; priority <= 2; priority++){
				for(k = 0; k < num_domains_mut; k++){
					mut_domain = mut_domains[k];
					for(ScorePenalty s : scoredElements[mut_domain]){
						if (s.getPriority()==priority){
							s.evalScore(domain,domain_markings); //STATE CHANGE
						}
					}
				}

				//Decide whether we improved, checking all penalties so that we can compare
				//objectively.

				worstPenaltyScore = 0;
				deltaScore = 0;
				for(ScorePenalty s : allScores){
					deltaScore += s.getCDelta();
					if (s.cur_score > worstPenaltyScore){
						worstPenalty = s;
						worstPenaltyScore = s.cur_score; //Can't have 0 penalties.
					}
				}

				revert_mutation = false;
				newPointReached = false;

				if (deltaScore > 0){
					revert_mutation = true;
				}
				/*
				if (bestWorstPenaltyScore[priority-1]<worstPenaltyScore){
					revert_mutation = true;
				} else if (bestWorstPenaltyScore[priority-1]==worstPenaltyScore){
					//Very often, there is a "wall" at a problem spot - allow the rest of the system
					//to also optimize in the background.
					if (deltaScore > 0){
						revert_mutation = true;
					}
				} else {
					newPointReached = true;
				}
				 */

				deltaScoreSum += deltaScore;
				if (revert_mutation){
					//Short circuit! Get out of there.
					priority ++;

					break priorityLoop;
				} else {
					//Keep the mutations
					bestWorstPenaltyScore[priority] = Math.min(bestWorstPenaltyScore[priority],worstPenaltyScore);
				}
			}

			if (deltaScoreSum > 0){
				revert_mutation = true;
			}
			if (deltaScoreSum < 0){
				newPointReached = true;
			}
			//Having run some or all of the affected penalty calculations,
			//Revert the states or Dedicate them as need be here.

			if (revert_mutation){
				timesSameCount++;
				//Revert
				/*
				if (ENABLE_MARKINGS){
					for(int[] row : domain_markings){
						for(int q : row){
							System.out.printf("%4d",q);
						}
						System.out.println();
					}
				}
				*/
				for(k = 0; k < num_domains_mut; k++){
					mut_domain = mut_domains[k];
					for(ScorePenalty s : scoredElements[mut_domain]){
						s.revert();
					}
					//Have to go back to old sequences..
					System.arraycopy(old_domains_markings[mut_domain],0,domain_markings[mut_domain], 0, domain[mut_domain].length);
					System.arraycopy(old_domains[mut_domain], 0, domain[mut_domain], 0, domain[mut_domain].length);
				}	

				//Which priority did we short circuit on?
				if (priority==2){
					//If we're having trouble with crosstalk (which is really crazy n^2 stuff), 
					//ignore the advice of the self fold server.
					for(k = 0; k < num_domains_mut; k++){
						mut_domain = mut_domains[k];
						//System.out.println(Arrays.toString(domain_markings[mut_domain]));
						Arrays.fill(domain_markings[mut_domain], 1);
					}
				}
				/*
				for(k = 0; k < num_domains_mut; k++){
					mut_domain = mut_domains[k];
					for(ScorePenalty s : scoredElements[mut_domain]){
						s.evalScore(domain, domain_markings);
						if (s.cur_score!=s.old_score){
							throw new RuntimeException("Assertion failure.");
						}
					}
				}
				*/
			} else {
				if (priority != 3){
					throw new RuntimeException("Assertion failure: self fold layer (layer 2) never ran.");
				}
				
				timesSameCount = 0;
				//Dedicate
				boolean DebugPenaltiesNow = DEBUG_PENALTIES && newPointReached;
				if (DebugPenaltiesNow){
					resetDebugPenalty();
				}
				current_score = 0;
				DIR.beginScoreReport();{
					for(ScorePenalty s : allScores){
						s.dedicate();
						current_score += s.cur_score;
						if (DebugPenaltiesNow){
							debugPenalty(s, domain,dsd);
						}
					}
				}DIR.endScoreReport();
				/*
				System.out.println("Current matrix:[");
				for(float[] row : DIR.currentMatrix){
					for(float val : row){
						System.out.printf("%.3f, ",val);
					}
					System.out.println();
				}
				System.out.println("]");
				*/
				if (OUTPUT_RUNTIME_STATS){
					System.out.println(num_mut_attempts+" , "+best_score);
				}
				//System.out.println(Arrays.toString(worstPenalty.getSeqs()));
				//current_score += deltaScoreSum;
				best_score = current_score;
				displayDomains(domain, false); //Updates public domains

		
				if (DebugPenaltiesNow){
					printDebugPenalties();
					System.out.println();
					displayDomains(domain, true);
					if (ENABLE_MARKINGS){
						/*
						for(int[] row : domain_markings){
							System.out.println(Arrays.toString(row));
						}
						*/
					}
				}

				double bestWorstPenaltyScoreSum = 0;
				for(double q : bestWorstPenaltyScore){
					//Disallow negative contributions. The overall EndThreshold can still be >0, however.
					bestWorstPenaltyScoreSum += Math.max(q,0);
				}
				if (bestWorstPenaltyScoreSum <= EndThreshold){
					break;
				}

				//System.out.printf("Components: CX%.3f IN%.3f JX%.3f INT%.3f\n",crosstalkSum,interactionSum,junctionPenalties,selfCrossTalk);
			}
			//Iteration complete.

			//System.out.println(num_domains_mut);

			num_mut_attempts++;
			if (num_mut_attempts > TOTAL_ATTEMPTS){
				break;
			}
			if ((System.nanoTime()-lastDumpState)>2e9){
				lastDumpState = System.nanoTime();
				//if (DEBUG_PENALTIES){
					System.out.println("iteration: "+num_mut_attempts+" score: "+best_score+" most significant subscore: "+(priority-1)+" "+(current_score+deltaScore));
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
	 */
	private void pickInitialSequence(int[] domain) {
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
	}


	private int[] mut_old_domain;
	/**
	 * Mutates mut_domain such that it is changed at the nucleic acid sequence level. If cc is nonnull, the changes will conserve amino acid sequence. 
	 * @param seqToSynthesize 
	 * @param domain_markings 
	 * @param cc 
	 * @param max_mutations 
	 */
	private boolean mutateUntilValid(int mut_domain, int[][] domain,
			int[] domain_length, List<DomainSequence> seqToSynthesize, int[][] domain_markings, CodonCode cc, int min_mutations, int max_mutations, int[] maxISO, int[] maxPZ) {
		int len = domain_length[mut_domain];
		/*
		// Save old domain
		if (!(mut_old_domain!=null && len <= mut_old_domain.length)){
			mut_old_domain = new int[len];
		}
		System.arraycopy(domain[mut_domain],0,mut_old_domain,0,len);
		int[] mut_old = mut_old_domain;
		*/
		int[] mut_new = domain[mut_domain];

		//Count the ones
		int oneC = 0;
		
		//System.out.println(num_mut+" "+min_mutations);
		
		//How many mutations will we try now?
		int num_mut = Math.min(len,max_mutations);
		//Count the reccomended bases to mutate
		if (ENABLE_MARKINGS){
			for(int q : domain_markings[mut_domain]){
				if (q!=DNAMARKER_DONTMUTATE){
					oneC++;
				}
			}
			num_mut = Math.min(num_mut,oneC);
		}
		//Randomize choice:
		if (num_mut > 0){
			num_mut = int_urn(1,num_mut);
			num_mut = Math.max(num_mut,min_mutations);
		}
		
		//System.out.println(num_mut+" "+min_mutations);
		

		for (int k = 0; k < num_mut; k++) {
			if (cc==null){ //no codon table: single base modifications
				// Attempt a mutation
				int j = int_urn(0, len-1); // select a base to mutate

				if (ENABLE_MARKINGS){
					//Is this a base we shouldn't target ? 
					if (domain_markings[mut_domain][j] == DNAMARKER_DONTMUTATE /*&& domain_markings[mut_domain][otherJ] == 0*/ && !(!MARKINGS_ENABLESTRICT && int_urn(1,1000)==1)){
						k--;
						continue;
					}
				}

				int oldJ = domain[mut_domain][j]%DNAFLAG_ADD;

				if (domain[mut_domain][j] > DNAFLAG_ADD && domain[mut_domain][j] < DNAFLAG_ADD*2) {
					// Base immutable
					continue; 
				} else if (domain[mut_domain][j] > DNAFLAG_ADD*2){
					//Alternate between G / C
					switch(domain[mut_domain][j]){
					case GCLC:
						mut_new[j] = GCLG;
						break;
					case GCLG:
						mut_new[j] = GCLC;
						break;
					}
				} else {
					//Switch to one of the other bases
					mut_new[j] = DnaDefinition.pickNewBase(mut_new,mut_new[j],maxISO[mut_domain],maxPZ[mut_domain]);
				}

				if (len > 1 && isGClessAT(mut_new[j])!=0){
					//Handle the new creation, by swapping another base to conserve the GC / AT balance.
					int otherJ;
					int candidateL3 = -1;
					int candidateL2 = -1;
					int candidateL1 = -1;
					int otherJ_startCheck = int_urn(0,len-1);
					searchOtherJ: for(int otherJi = 0; otherJi < len; otherJi++){
						otherJ = (otherJ_startCheck+otherJi)%len;
						if (otherJ != j){
							if (domain[mut_domain][otherJ] < DNAFLAG_ADD){
								//Necessary conditions over.
								if (isGC(domain[mut_domain][otherJ])!=isGC(oldJ)){
									if (domain_markings[mut_domain][otherJ]!=DNAMARKER_DONTMUTATE){
										candidateL1 = otherJ; //All conditions (L1) succeeded.
										break searchOtherJ;
									} else {
										candidateL2 = otherJ;
									}
								} else {
									candidateL3 = otherJ;
								}
							}
						}
					}
					otherJ = candidateL1;
					if (otherJ==-1){
						otherJ = candidateL2;
					}
					if (otherJ==-1){
						otherJ = candidateL3;
					}
					if (otherJ==-1){
						//It absolutely was not possible to find a partner. Oh well.
					} else {
						int newJ = domain[mut_domain][j]%DNAFLAG_ADD;
						if ((oldJ==G || oldJ==C) && (newJ==A || newJ==T)){
							domain[mut_domain][otherJ] = int_urn(0, 1)==0?G:C;
						}
						if ((oldJ==A || oldJ==T) && (newJ==G || newJ==C)){
							domain[mut_domain][otherJ] = int_urn(0, 1)==0?A:T;
						}
					}
				}
			} else { //Codon table: 3 mutation points = 1 amino swap.
				// Attempt a mutation
				int j = int_urn(0, len-1); // select a base to mutate
				//Round to codon
				j -= j%3;

				if (ENABLE_MARKINGS){
					//Is this a codon we shouldn't target ?
					boolean hasMutateMarker = false;
					for(int jb = j; jb < j+3; jb++){
						hasMutateMarker |= domain_markings[mut_domain][jb] != DNAMARKER_DONTMUTATE;
					}
					if (!hasMutateMarker){
						if (!(!MARKINGS_ENABLESTRICT && int_urn(1,1000)==1)){
							k--;
							continue;
						}
					}
				}
				
				//Mutate codon!
				cc.mutateToOther(domain[mut_domain], j);

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
	private void resetDebugPenalty(){
		for(int k = 0; k < debug_keepMessages.length; k++){
			debug_keepMessages[k] = 0;
			debug_messages[k] = "";
		}
	}
	private static final int debugOutTopScores = 50;
	private double[] debug_keepMessages = new double[debugOutTopScores];
	private String[] debug_messages = new String[debugOutTopScores];


	private void debugPenalty(ScorePenalty s, int[][] domain, DomainStructureData dsd) {
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
	
	
	//Utilities
	private static final int 
	G = DnaDefinition.G, 
	A = DnaDefinition.A, 
	T = DnaDefinition.T,
	C = DnaDefinition.C,
	D = DnaDefinition.D,
	H = DnaDefinition.H,
	GL = DnaDefinition.G+DNAFLAG_ADD,
	AL = DnaDefinition.A+DNAFLAG_ADD,
	TL = DnaDefinition.T+DNAFLAG_ADD,
	CL = DnaDefinition.C+DNAFLAG_ADD,
	DL = DnaDefinition.D+DNAFLAG_ADD,
	HL = DnaDefinition.H+DNAFLAG_ADD,
	GCLG = DnaDefinition.G+DNAFLAG_ADD*2,
	GCLC = DnaDefinition.C+DNAFLAG_ADD*2;
	
	public static int getLockedBase(char charAt) {
		return decodeConstraintChar(charAt)%DNAFLAG_ADD+DNAFLAG_ADD*1;
	}
	
	public static int decodeConstraintChar(char charAt){
		switch(charAt){
		case 'G':
			return G;
		case 'A':
			return A;
		case 'C':
			return C;
		case 'T':
			return T;
		case 'D':
			return D;
		case 'H':
			return H;
		case 'g':
			return GL;
		case 'a':
			return AL;
		case 't':
			return TL;
		case 'c':
			return CL;
		case 'd':
			return DL;
		case 'h':
			return HL;
		default:
			return 0; //In practice, pickInitialSequences will erase this anyway.
		}
	}
	private static int isGClessAT(int i) {
		return isGC(i)?1:(isAT(i)?-1:0);
	}
	private static boolean isGC(int otherJ) {
		otherJ = DnaDefinition.noFlags(otherJ);
		return (otherJ == G) || (otherJ == C);
	}
	private static boolean isAT(int otherJ) {
		otherJ = DnaDefinition.noFlags(otherJ);
		return (otherJ == A) || (otherJ == T);
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
	private static int int_urn(int from, int to) {
		return (int)(Math.random()*(to-from+1)+from);
	}
}