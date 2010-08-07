package DNASim$DB;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Map.Entry;

public class DomainDesigner_ByRandomPartialMutations {
	public static class DDSeqDesigner_ByRandomMutation_Client implements DDSeqDesigner{
		private boolean waitOnStart = true;
		private boolean finished = false;
		private Runnable runOnStart = new Runnable(){
			public void run(){
				//TODO: Alias map, so we have a contiguous list of domains.

				final ArrayList<DomainSequence> theJunctions = new ArrayList();

				int num_domain = 0;
				TreeMap<Integer, Integer> domain_length_t = new TreeMap<Integer, Integer>();
				for(String q : inputStrands){
					DomainStructureData dsd = new DomainStructureData();
					DomainDesigner_SharedUtils.utilJunctionSplitter(theJunctions, q);
					DomainStructureData.readStructure(q, dsd);

					num_domain = Math.max(num_domain,dsd.domainLengths.length-1);
					
					for(int i = 0; i < dsd.domainLengths.length; i++){
						int val = dsd.domainLengths[i];
						if (val!=-1){
							domain_length_t.put(i, val);
						}
					}
					
				}
				
				final int[] domain_length = new int[num_domain];
				
				for(int k = 0; k < domain_length.length; k++){
					Integer got = domain_length_t.get(k+1);
					if (got==null){
						throw new RuntimeException("No length data for domain "+(k+1));
					}
					domain_length[k] = got;
				}

				final Map<Integer, String> lockedDomains = new TreeMap<Integer,String>();

				
				final int num_domain_2 = num_domain;
				
				new Thread(){public void run(){
					while(waitOnStart && ! r.abort){
						try {
							Thread.sleep(100);
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
					}

					r.main(num_domain_2, domain_length, Integer.MAX_VALUE, lockedDomains, theJunctions);

					finished = true;
				}
				}.start();
			}
		};
		DomainDesigner_ByRandomPartialMutations r;
		{
			r = new DomainDesigner_ByRandomPartialMutations();
		}
		private List<String> inputStrands;
		public DDSeqDesigner_ByRandomMutation_Client(List<String> inputStrands){
			this.inputStrands = inputStrands;
			options.add(handleJunctions);
			options.add(penalizeJunctions);
			options.add(performFullN4Junctions);
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


		public String getResult() {
			if (r.outputDomains==null){
				return "Output incomplete. Run designer longer";
			}
			StringBuffer sb = new StringBuffer();
			String lR = "\n";
			for(int k = 0; k < r.outputDomains.length; k++){
				sb.append("Domain ");
				sb.append(k+1);
				sb.append(": ");
				sb.append(r.outputDomains[k]);
				sb.append(lR);
				sb.append("Complem.: ");
				sb.append(revComp(r.outputDomains[k]));
				sb.append(lR);
			}
			sb.append(lR);
			sb.append("Input Strands:");
			sb.append(lR);
			DomainSequence ds = new DomainSequence();
			for(String q : inputStrands){
				for(String subStrand : q.split("}")){
					sb.append(subStrand.replaceAll("\\s+",""));
					sb.append(lR);
					ds.setDomains(subStrand);
					sb.append("[");
					for(int k = 0; k < ds.domainList.length; k++){
						String domain = r.outputDomains[ds.domainList[k] & DNA_SEQ_FLAGSINVERSE];
						if ((ds.domainList[k] & DNA_COMPLEMENT_FLAG)!=0){
							domain = revComp(domain);
						}
						sb.append(domain);
						if (k + 1 < ds.domainList.length){
							sb.append("|");
						}
					}
					sb.append("}");
					sb.append(lR);
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
			waitOnStart = false;
			r.waitForResume = false;
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
	}
	
	public static DDSeqDesigner makeDesigner(List<String> inputStrands){
		return new DDSeqDesigner_ByRandomMutation_Client (inputStrands);
	}
	
	//Accessible to client:
	private double best_score;
	private String[] outputDomains;
	private boolean waitForResume = false, abort = false;
	
	//End accessible to client.

	int MAX_DOMAIN_LENGTH = 60 ;
	int NB_ENABLE = 1;
	int NB_DISABLE = 0;
	int MAX_MUTATIONS = 10; // maximum number of simultaneous mutations

	/**
	 * Matching (complementarity) Strengths. Negative means repulsion.
	 */
	double GCstr = 2;
	double ATstr = 1;
	double GTstr = 0;
	double MBstr = -3;
	/**
	 * "Large loop" - what?
	 */
	double LLstr = -0.5; 
	/**
	 * "Score for domain ending in a base pair" - what?
	 */
	double DHstr = 3; 
	int MAX_IMPORTANCE = 100;
	int LHbases = 4;
	double LHstart = 2;
	double LHpower = 2;
	double INTRA_SCORE = 5; // score bonus for intrastrand/dimerization interactions
	double CROSSTALK_SCORE = -5; // score bonus for crosstalk (as compared to interaction)
	double CROSSTALK_DIV = 2; // crosstalk score is divided by this much (and then score is subtracted)
	double GGGG_PENALTY = 50;
	double ATATAT_PENALTY = 20;
	boolean PERFORM_COMPLETE_JUNCTIONS = false;
	boolean CALCULATE_JUNCTION_PENALTIES = true;
	boolean JunctionPenalty_FullPenalty = true;

	int int_urn(int from, int to) {
		return (int)(Math.random()*(to-from+1)+from);
	}
	String printf(String k){
		System.out.print(k);
		return k;
	}
	private static final int G = 1, A = 2, T = 3, C = 4, GL = 11, AL = 12, TL = 13, CL = 14;
	private int getLockedBase(char charAt) {
		switch(charAt){
		case 'G':
			return GL;
		case 'A':
			return AL;
		case 'C':
			return CL;
		case 'T':
			return TL;
		}
		throw new IllegalArgumentException("Unrecognized base "+charAt);
	}
	public String DisplayBase(int base) {
		// 1 = G, 2 = A, 3 = T, 4 = C; 11 = G (locked), etc
		if (base == G)
			return printf("G");
		else if (base == A)
			return printf("A");
		else if (base == T)
			return printf("T");
		else if (base == C)
			return printf("C");
		else if (base == GL)
			return printf("G");
		else if (base == AL)
			return printf("A");
		else if (base == TL)
			return printf("T");
		else if (base == CL)
			return printf("C");
		else {
			throw new IllegalArgumentException("Unrecognized base "+base);
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
			System.out.print(" ");
		}
		this.outputDomains = outputDomains;
		System.out.println();
	}
	void displaySequence(DomainSequence sequence, int[][] domain) {
		int len = sequence.length(domain);
		System.out.print(sequence.toString()+" = ");
		for(int i = 0; i < len; i++){
			DisplayBase(sequence.base(i, domain));
		}
		System.out.println();
	}
	
	/**
	 * Interaction score.
	 */
	private double[][] Cmatrix_pairscore;
	private double[][] Smatrix_pairscore;
	private int[][] SDmatrix_pairscore;
	public static final int DNA_COMPLEMENT_FLAG = 0x8000;
	public static final int DNA_SEQ_FLAGSINVERSE = ~(DNA_COMPLEMENT_FLAG);
	private static boolean checkComplementary(DomainSequence a, DomainSequence b){
		for(int k = 0; k < a.numDomains; k++){
			for(int y = 0; y < b.numDomains; y++){
				int abase = a.domainList[k];
				int bbase = b.domainList[k];
				if ((abase&DNA_COMPLEMENT_FLAG)!=(bbase&DNA_COMPLEMENT_FLAG)){
					if ((abase&DNA_SEQ_FLAGSINVERSE)==(bbase&DNA_SEQ_FLAGSINVERSE)){
						return true;
					}
				}
			}
		}
		return false;
	}
	double pairscore(DomainSequence seq1, DomainSequence seq2, int[][] domain) {
		// Gives the score of the two sequences's crosstalk
		double score, temp;
		int i, j, k;
		int len1 = seq1.length(domain);
		int len2 = seq2.length(domain);
		if (!(Cmatrix_pairscore!=null && len1 <= Cmatrix_pairscore.length && len2 <= Cmatrix_pairscore[0].length)){
			Cmatrix_pairscore = new double[len1][len2];
			Smatrix_pairscore = new double[len1][len2];
			SDmatrix_pairscore = new int[len1][len2];
		}
		double[][] Cmatrix = Cmatrix_pairscore; // complementarity matrix
		double[][] Smatrix = Smatrix_pairscore; // score matrix
		int[][] SDmatrix = SDmatrix_pairscore; // running total of helix size, 0 if current base didn't contribute.

		// LEN1 x LEN2: Complementarity calculation .
		// Note that the binding is going 5'-3' on both strands, and assumes seq2 is in the array 3'-5'! 
		for (i = 0; i < len1; i++) {
			for (j = 0; j < len2; j++) {
				int base1 = seq1.base(i,domain);
				int base2 = seq2.base(len2-1-j,domain);
				if (((base1 + base2)%10 == 5)&&((base1 * base2)%10 == 4)) // G/C Match
					Cmatrix[i][j] = GCstr;
				else if (((base1 + base2)%10 == 5)&&((base1 * base2)%10 == 6)) // A/T Match
					Cmatrix[i][j] = ATstr;
				else if (((base1 + base2)%10 == 4)&&((base1 * base2)%10 == 3)) // G/T Wobble
					Cmatrix[i][j] = GTstr;
				else
					Cmatrix[i][j] = MBstr; // mismatch
			}
		}

		// Calculate score
		score = 0;

		//Maximize score among "starting positions" (Where to begin the window on seq1, which is 5'-3'. The window goes 3'-5' on seq2, of course.)
		for (j = 0; j < len2; j++) {
			Smatrix[0][j] = Cmatrix[0][j];
			if (Smatrix[0][j] < 0) {
				Smatrix[0][j] = 0;
				SDmatrix[0][j] = 0;
			} else {
				Smatrix[0][j] = Smatrix[0][j] + DHstr;
				SDmatrix[0][j] = 1;
			}
			if (Smatrix[0][j] > score)
				score = Smatrix[0][j];
		}


		//Slide the window, maximizing score as it goes, and giving bonus points for helices.
		//It ALSO allows for alignments that have "bulges" on one side or the other.
		for (i = 1; i < len1; i++) {
			Smatrix[i][0] = Cmatrix[i][0];
			if (Smatrix[i][0] < 0) {
				Smatrix[i][0] = 0;
				SDmatrix[i][0] = 0;
			} else {
				Smatrix[i][0] = Smatrix[i][0] + DHstr;
				SDmatrix[i][0] = 1;
			}
			if (Smatrix[i][0] > score)
				score = Smatrix[i][0];

			for (j = 1; j < len2; j++) {
				if (Cmatrix[i][j] < 0) { // destabilizing base
					SDmatrix[i][j] = 0;
					Smatrix[i][j] = 0;

					if ((SDmatrix[i-1][j-1] > 0)&&(Smatrix[i-1][j-1] + MBstr > 0)) // starting a mismatch loop
						Smatrix[i][j] = Smatrix[i-1][j-1] + MBstr;
					if ((SDmatrix[i-1][j-1] == 0)&&(Smatrix[i-1][j-1] + LLstr > 0)) // expanding a mismatch loop
						Smatrix[i][j] = Smatrix[i-1][j-1] + LLstr;

					if ((SDmatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + MBstr > 0)&&(Smatrix[i][j-1] + MBstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i][j-1] + MBstr;
					if ((SDmatrix[i][j-1] == 0)&&(Smatrix[i][j-1] + LLstr > 0)&&(Smatrix[i][j-1] + LLstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i][j-1] + LLstr;

					if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + MBstr > 0)&&(Smatrix[i-1][j] + MBstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + MBstr;
					if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] + LLstr > 0)&&(Smatrix[i-1][j] + LLstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + LLstr;

					if (Smatrix[i][j] < 0)
						Smatrix[i][j] = 0;
				} else { // stabilizing base
					Smatrix[i][j] = Cmatrix[i][j];
					SDmatrix[i][j] = 1;

					if ((SDmatrix[i-1][j-1] > 0)&&(Smatrix[i-1][j-1] > 0)) { // continuing a helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];
						SDmatrix[i][j] = SDmatrix[i-1][j-1] + 1;
					} else if ((SDmatrix[i-1][j-1] == 0)&&(Smatrix[i-1][j-1] > 0)) { // starting a new helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];
						SDmatrix[i][j] = 1;
					}	  

					if ((SDmatrix[i][j-1] > 0)&&(Smatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + Cmatrix[i][j] - Cmatrix[i][j-1] + MBstr > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i][j-1] + Cmatrix[i][j] - Cmatrix[i][j-1] + MBstr; // introducing a 1-bulge, destroying previous bond
						SDmatrix[i][j] = 1;
					} else if ((SDmatrix[i][j-1] == 0)&&(Smatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + Cmatrix[i][j] > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i][j-1] + Cmatrix[i][j]; // closing a bulge
						SDmatrix[i][j] = 1;
					}

					if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr;
						SDmatrix[i][j] = 1;
					} else if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j];
						SDmatrix[i][j] = 1;
					}

					if (SDmatrix[i][j] > LHbases) {
						// Extra points for long helices
						temp = LHstart;
						for (k = LHbases; k < SDmatrix[i][j]; k++)
							temp = temp * LHpower;
						Smatrix[i][j] = Smatrix[i][j] + temp;

						/*
						if (DEBUGNOW-->0){
							System.out.println(seq1.toString()+"|"+seq2.toString());
							System.out.println(SDmatrix[i][j]);
						}
						*/
					}
				}

				if ((SDmatrix[i][j] > 0)&&((i == (len1-1))||(j == (len2-1))))
					Smatrix[i][j] = Smatrix[i][j] + DHstr;

				//Is this the best window?
				if (Smatrix[i][j] > score)
					score = Smatrix[i][j];
			} 
		}


		//Wow. this is a complicated little thing, isn't it?
		return score;
	}


	double selfcrosstalk(DomainSequence seq, int[][] domain) {
		double score, temp;
		int i, j, k;
		int len1 = seq.length(domain);
		if (!(Cmatrix_pairscore!=null && len1 <= Cmatrix_pairscore.length && len1 <= Cmatrix_pairscore[0].length)){
			Cmatrix_pairscore = new double[len1][len1];
			Smatrix_pairscore = new double[len1][len1];
			SDmatrix_pairscore = new int[len1][len1];
		}
		double[][] Cmatrix = Cmatrix_pairscore; // complementarity matrix
		double[][] Smatrix = Smatrix_pairscore; // score matrix
		int[][] SDmatrix = SDmatrix_pairscore; // running total of helix size, 0 if current base didn't contribute.

		// NxN complementarities. 
		for (i = 0; i < len1; i++) {
			for (j = 0; j < len1; j++) {
				int base1 = seq.base(i,domain);
				int base2 = 15-seq.base(len1-1-j,domain);
				if (((base1 + base2)%10 == 5)&&((base1 * base2)%10 == 4)) // G/C Match
					Cmatrix[i][j] = GCstr;
				else if (((base1 + base2)%10 == 5)&&((base1 * (base2))%10 == 6)) // A/T Match
					Cmatrix[i][j] = ATstr;
				else if (((base1 + (base2))%10 == 4)&&((base1 * (base2))%10 == 3)) // G/T Wobble
					Cmatrix[i][j] = GTstr;
				else
					Cmatrix[i][j] = MBstr; // mismatch
			}
		}

		score = 0;

		//Unlike the pairwise scoring, which uses a sliding window,
		//this is a linear search, which only looks for places where secondary
		//structure could form (backfolding)
		Smatrix[0][0] = 0;

		for (j = 1; j < len1; j++) {
			Smatrix[0][j] = Cmatrix[0][j];
			if (Smatrix[0][j] < 0) {
				Smatrix[0][j] = 0;
				SDmatrix[0][j] = 0;
			} else {
				Smatrix[0][j] = Smatrix[0][j] + DHstr;
				SDmatrix[0][j] = 1;
			}
			if (Smatrix[0][j] > score)
				score = Smatrix[0][j];
		}


		for (i = 1; i < len1; i++) {
			Smatrix[i][0] = Cmatrix[i][0];
			if (Smatrix[i][0] < 0) {
				Smatrix[i][0] = 0;
				SDmatrix[i][0] = 0;
			} else {
				Smatrix[i][0] = Smatrix[i][0] + DHstr;
				SDmatrix[i][0] = 1;
			}
			if (Smatrix[i][0] > score)
				score = Smatrix[i][0];

			for (j = 1; j < len1; j++) {

				if (i == j) { 
					// "Main line" match, do not score
					SDmatrix[i][j] = 0;
					Smatrix[i][j] = 0;
				} else if (Cmatrix[i][j] < 0) { // destabilizing base
					SDmatrix[i][j] = 0;
					Smatrix[i][j] = 0;

					if ((SDmatrix[i-1][j-1] > 0)&&(Smatrix[i-1][j-1] + MBstr > 0)) // starting a mismatch loop
						Smatrix[i][j] = Smatrix[i-1][j-1] + MBstr;
					if ((SDmatrix[i-1][j-1] == 0)&&(Smatrix[i-1][j-1] + LLstr > 0)) // expanding a mismatch loop
						Smatrix[i][j] = Smatrix[i-1][j-1] + LLstr;

					if ((SDmatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + MBstr > 0)&&(Smatrix[i][j-1] + MBstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i][j-1] + MBstr;
					if ((SDmatrix[i][j-1] == 0)&&(Smatrix[i][j-1] + LLstr > 0)&&(Smatrix[i][j-1] + LLstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i][j-1] + LLstr;

					if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + MBstr > 0)&&(Smatrix[i-1][j] + MBstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + MBstr;
					if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] + LLstr > 0)&&(Smatrix[i-1][j] + LLstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + LLstr;

					if (Smatrix[i][j] < 0)
						Smatrix[i][j] = 0;

				} else { // stabilizing base
					Smatrix[i][j] = Cmatrix[i][j];
					SDmatrix[i][j] = 1;

					if ((SDmatrix[i-1][j-1] > 0)&&(Smatrix[i-1][j-1] > 0)) { // continuing a helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];
						SDmatrix[i][j] = SDmatrix[i-1][j-1] + 1;
					} else if ((SDmatrix[i-1][j-1] == 0)&&(Smatrix[i-1][j-1] > 0)) { // starting a new helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];
						SDmatrix[i][j] = 1;
					}	  

					if ((SDmatrix[i][j-1] > 0)&&(Smatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + Cmatrix[i][j] - Cmatrix[i][j-1] + MBstr > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i][j-1] + Cmatrix[i][j] - Cmatrix[i][j-1] + MBstr; // introducing a 1-bulge, destroying previous bond
						SDmatrix[i][j] = 1;
					} else if ((SDmatrix[i][j-1] == 0)&&(Smatrix[i][j-1] > 0)&&(Smatrix[i][j-1] + Cmatrix[i][j] > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i][j-1] + Cmatrix[i][j]; // closing a bulge
						SDmatrix[i][j] = 1;
					}

					if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr;
						SDmatrix[i][j] = 1;
					} else if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] > Smatrix[i][j])) {
						Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j];
						SDmatrix[i][j] = 1;
					}

					if (SDmatrix[i][j] > LHbases) {
						// Extra points for long helices
						temp = LHstart;
						for (k = LHbases; k < SDmatrix[i][j]; k++)
							temp = temp * LHpower;
						Smatrix[i][j] = Smatrix[i][j] + temp;
					}
				}

				if ((SDmatrix[i][j] > 0)&&((i == (len1-1))||(j == (len1-1))))
					Smatrix[i][j] = Smatrix[i][j] + DHstr;

				if (Smatrix[i][j] > score)
					score = Smatrix[i][j];

			} 
		}

		return score;
	}
	
	public static void main(String[] args){
		TreeMap<Integer, String> lock = new TreeMap<Integer, String>();
		lock.put(6,"GTTC");
		ArrayList<DomainSequence> junctions = new ArrayList();
		//utilJunctionSplitter(junctions,"[1*|4*|7|5*|6*}");
		//utilJunctionSplitter(junctions,"[3*|2*|7*|1*}");
		DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[1|7|2|3|7*|1*|4*|7|3*|2*|7*}");
		DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[3|7*|4|1|7|3*|2*|7*|1*|4*|7|5*|6*}");
		DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[4|7*|5|6|7|4*|1*|7*|6*|5*|7}");
		DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[6|7|1|4|7*|6*|5*|7|4*|1*|7*|2*|3*}");
		new DomainDesigner_ByRandomPartialMutations().main(7,new int[]{8,8,8,8,8,8,4},Integer.MAX_VALUE,lock,junctions);
	}
	/**
	 * The program will attempt to prevent any interaction between the DNASequences in toSynthesize, though interactions
	 * containing complementary sequences will be ignored (i.e, 1|4*|5 will not be tested against 1|4|5).
	 */
	void main(int num_domain, int[] domain_length, int TOTAL_ATTEMPTS, Map<Integer, String> lockedDomains, List<DomainSequence> seqToSynthesize) {
		//SANITY OF INPUT:
		DomainDesigner_SharedUtils.utilRemoveDuplicateSequences(seqToSynthesize);
		//System.out.println(seqToSynthesize.size());
		
		//CONTINUE.
		int i, j, k;
		int num_mut;
		int mut_domain; // Domain, base, old, and new values
		int[] mut_base;
		int[] mut_new;
		int[] mut_old;
		
		int temp_domain[] = new int[MAX_DOMAIN_LENGTH];
		int[] domain_importance;
		int[] domain_gatc_avail;
		int[][] domain; // 1 = G, 2 = A, 3 = T, 4 = C; 11 = G (locked), etc

		long num_mut_attempts, total_mutations;
		int rule_4g, rule_6at, rule_ccend, rule_ming, rule_init, rule_lockold, rule_targetworst, rule_gatc_avail;

		// Set up memory for mutation computations

		mut_base = new int[MAX_MUTATIONS];
		mut_old = new int[MAX_MUTATIONS];
		mut_new = new int[MAX_MUTATIONS];

		printf("\n           Domain-based sequence design\n");
		printf("                     v. 0.2\n");
		printf("                 by Dave Zhang\n\n");
		

		domain = new int[num_domain][];

		domain_gatc_avail = new int[num_domain];

		domain_importance = new int[num_domain];

		//Domain length and importance.
		for (i = 0; i < num_domain; i++) {
			domain[i] = new int[domain_length[i]];
			domain_importance[i] = 1;
		}

		rule_4g = 0; // cannot have 4 G's or 4 C's in a row
		rule_6at = 1; // cannot have 6 A/T bases in a row
		rule_ccend = 1; // domains MUST start and end with C
		rule_ming = 1; // design tries to minimize usage of G
		rule_init = 7; // 1 = polyN, 2 = poly-H, 3 = poly-Y, 4 = poly-T
		rule_targetworst = 1; // target worst domains for mutation
		rule_gatc_avail = 15; // all flags set (bases available)

		// Generate starting domain sequences
		// 1 = G, 2 = A, 3 = T, 4 = C; 11 = G (locked), etc

		for (i = 0; i < num_domain; i++) {
			domain_gatc_avail[i] = rule_gatc_avail;
			for (j = 0; j < domain_length[i]; j++) {
				domain[i][j] = 0;
				while (domain[i][j] == 0) {
					k = int_urn(1,4);
					if ((k == 4)&&(rule_init/8 == 1))
						domain[i][j] = 1;
					if ((k == 3)&&((rule_init / 4) % 2 == 1))
						domain[i][j] = 2;
					if ((k == 2)&&((rule_init / 2) % 2 == 1))
						domain[i][j] = 3;
					if ((k == 1)&&(rule_init % 2 == 1))
						domain[i][j] = 4;
				}
			}

			if (rule_ccend == 1) {
				if (rule_gatc_avail % 2 == 1) 
					domain[i][0] = 14;
				else if (rule_gatc_avail / 8 == 1)
					domain[i][0] = 11;
				else if ((rule_gatc_avail / 2) % 2 == 1)
					domain[i][0] = 13;
				else
					domain[i][0] = 12;

				if (rule_gatc_avail % 2  == 1) 
					domain[i][domain_length[i]-1] = 14;
				else if (rule_gatc_avail / 8 == 1)
					domain[i][domain_length[i]-1] = 11;
				else if ((rule_gatc_avail / 4) % 2 == 1)
					domain[i][domain_length[i]-1] = 12;
				else
					domain[i][domain_length[i]-1] = 13;
			}
		}
		
		//Locked domains?
		for(Entry<Integer, String> lock : lockedDomains.entrySet()){
			i = lock.getKey();
			for (j = 0; j < domain_length[i]; j++) {
				domain[i][j] = getLockedBase(lock.getValue().charAt(j));
			}
		}

		// Scoring temporaries
		double[][] crosstalk; // Crosstalk is for domain similarity
		double[][] interaction; // Interaction is for domain complementarity
		double[] domain_intrinsic; // intrinsic score to domains from various rules
		double[] domain_score; // domain score
		double score = Double.MAX_VALUE, old_d_intrinsic = Double.MAX_VALUE; // Score of system
		best_score = Double.MAX_VALUE;
		
		domain_score = new double[num_domain];
		Arrays.fill(domain_score, Double.MAX_VALUE);
		domain_intrinsic = new double[num_domain];
		Arrays.fill(domain_intrinsic, Double.MAX_VALUE);
		crosstalk = new double[num_domain][num_domain];
		interaction = new double[num_domain][num_domain];
		for(j = 0; j < num_domain; j++){
			for (i = 0; i < num_domain; i++) {
				updateInteractionsToMutated(i,j,domain,domain_length,interaction,crosstalk);
			}
		}
		int worst_domain = 0; // domain that causes the worst score
		
		// Initial display (Randomized start)
		displayDomains(domain);
		
		num_mut_attempts = 0;
		total_mutations = 0;

		// Main loop

		while (true && !abort) {			
			num_mut_attempts++;

			//How many mutations will we try now?
			num_mut = 0;
			while ((num_mut < MAX_MUTATIONS-1)&&(int_urn(0,1) == 1))
				num_mut++;

			num_mut++;

			// One third of all mutations are in "worst" domain if rule_targetworst active
			if (rule_targetworst==1) {
				if (int_urn(1,3) == 1)
					mut_domain = worst_domain;
				else 
					mut_domain = int_urn(0, num_domain-1); // select a domain to mutate
			} else {
				mut_domain = int_urn(0, num_domain-1); // select a domain to mutate
			}

			for (k = 0; k < num_mut; k++) {
				// Attempt a mutation
				j = int_urn(0, (domain_length[mut_domain])-1); // select a base to mutate

				if (domain[mut_domain][j] > 10) {
					// Base immutable
					mut_base[k] = j;
					mut_old[k] = domain[mut_domain][j];
					mut_new[k] = mut_old[k];
				}
				else {
					// Base is mutable

					mut_base[k] = j;
					mut_old[k] = domain[mut_domain][j];

					// hack to not mutate bases constrained to only 1 base
					if ((domain_gatc_avail[mut_domain] & 1)+(domain_gatc_avail[mut_domain] & 2)/2+(domain_gatc_avail[mut_domain] & 4)/4+(domain_gatc_avail[mut_domain] & 8)/8 > 1) {
						if (rule_ming==1) {
							// Minimize G rule active
							do {
								mut_new[k] = int_urn(1,100);
								if (mut_new[k] < 5)
									mut_new[k] = 1;
								else if (mut_new[k] < 37)
									mut_new[k] = 2;
								else if (mut_new[k] < 69)
									mut_new[k] = 3;
								else
									mut_new[k] = 4;
								// Undo mutation if new base is not allowed
								if ((mut_new[k] == 1)&&(domain_gatc_avail[mut_domain] / 8 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 2)&&((domain_gatc_avail[mut_domain] / 4)%2 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 3)&&((domain_gatc_avail[mut_domain] / 2)%2 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 4)&&(domain_gatc_avail[mut_domain] % 2== 0))
									mut_new[k] = mut_old[k];

							} while (mut_new[k] == mut_old[k]);
						} else {
							// Uniform GATC mix
							do {
								mut_new[k] = int_urn(1,100);
								if (mut_new[k] < 26)
									mut_new[k] = 1;
								else if (mut_new[k] < 51)
									mut_new[k] = 2;
								else if (mut_new[k] < 76)
									mut_new[k] = 3;
								else
									mut_new[k] = 4;
								// Undo mutation if new base is not allowed
								if ((mut_new[k] == 1)&&(domain_gatc_avail[mut_domain] / 8 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 2)&&((domain_gatc_avail[mut_domain] / 4)%2 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 3)&&((domain_gatc_avail[mut_domain] / 2)%2 == 0))
									mut_new[k] = mut_old[k];
								if ((mut_new[k] == 4)&&(domain_gatc_avail[mut_domain] % 2== 0))
									mut_new[k] = mut_old[k];
							} while (mut_new[k] == mut_old[k]);	  
						}
					}
				}
			}

			// Apply mutations, calculate score after mutations 
			for (k = 0; k < num_mut; k++)
				domain[mut_domain][mut_base[k]] = mut_new[k];

			old_d_intrinsic = domain_intrinsic[mut_domain];
			
			for (i = 0; i < num_domain; i++) {
				updateInteractionsToMutated(i,mut_domain,domain,domain_length,interaction,crosstalk);
			}  

			//Intrinsic penalties to mutated region
			domain_intrinsic[mut_domain] = 0;

			DomainSequence dnaSeq1 = new DomainSequence();
			DomainSequence dnaSeq2 = new DomainSequence();
			
			dnaSeq1.setDomains(mut_domain);
			addIntrinsicScore(dnaSeq1,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);

			// Add intrinsic scores of junctions!
			//For all junctions ... containing mut_domain ... 
			
			if (PERFORM_COMPLETE_JUNCTIONS){
				for (i = 0; i < num_domain; i++) {
					dnaSeq1.setDomains(mut_domain,i);
					addIntrinsicScore(dnaSeq1,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);
					dnaSeq1.setDomains(i,mut_domain);
					addIntrinsicScore(dnaSeq1,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);

					if (i!=mut_domain){
						dnaSeq1.setDomains(mut_domain,i | DNA_COMPLEMENT_FLAG);
						addIntrinsicScore(dnaSeq1,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);
						dnaSeq1.setDomains(mut_domain | DNA_COMPLEMENT_FLAG,i);
						addIntrinsicScore(dnaSeq1,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);

						//Add score of the interaction between all junction pairs
						dnaSeq1.setDomains(mut_domain,i);

						for(int q = 0; q < num_domain; q++){
							for(int z = 0; z < num_domain; z++){
								if (z==q) continue;

								dnaSeq2.setDomains(q,z);

								for(int swap = 0; swap < 4; swap++){
									if(swap%2==0){
										dnaSeq1.domainList[1] |= DNA_COMPLEMENT_FLAG;
										dnaSeq1.domainList[0] &= ~DNA_COMPLEMENT_FLAG;
									} else {
										dnaSeq1.domainList[1] &= ~DNA_COMPLEMENT_FLAG;
										dnaSeq1.domainList[0] &= ~DNA_COMPLEMENT_FLAG;
									}
									if (swap/2==0){
										dnaSeq2.domainList[1] |= DNA_COMPLEMENT_FLAG;
										dnaSeq2.domainList[0] &= ~DNA_COMPLEMENT_FLAG;
									} else {
										dnaSeq2.domainList[1] &= ~DNA_COMPLEMENT_FLAG;
										dnaSeq2.domainList[0] &= ~DNA_COMPLEMENT_FLAG;
									}

									if (!checkComplementary(dnaSeq1,dnaSeq2)){
										//Seems that it's not possible to do everything.
										domain_intrinsic[mut_domain] += pairscore(dnaSeq1, dnaSeq2,domain);
									}				
								}
							}
						}
					}
				}
			} else if (CALCULATE_JUNCTION_PENALTIES){
				//Only test the one's we'll synthesize
				for(DomainSequence seq : seqToSynthesize){
					if (seq.contains(mut_domain)){
						
						addIntrinsicScore(seq,mut_domain,domain,domain_intrinsic,rule_4g,rule_6at);
						
						for(DomainSequence seq2 : seqToSynthesize){
							if (seq2==seq){
								continue;
							}
							if (!checkComplementary(seq,seq2)){
								domain_intrinsic[mut_domain] += pairscore(seq,seq2,domain) * (JunctionPenalty_FullPenalty? 1 : 1f/(seqToSynthesize.size()*seqToSynthesize.size()));
							}
						}
					}
				}
			}
			
			// Domain score is max of interaction and crosstalk scores
			for(i = 0; i < num_domain; i++){
				updateDomainScore(domain_score,i,num_domain,interaction,crosstalk,domain_importance,domain_intrinsic);
			}

			score = 0;
			double worst_score = 0;
			for (i = 0; i < num_domain; i++) {
				if (domain_score[i] > worst_score) {
					worst_score = domain_score[i];
					worst_domain = i;
				}
				score += domain_score[i];
			}
			

			double compareToScore = best_score;
			if (int_urn(1,100) <= 20){
				//20% chance of "loosening" old_score
				compareToScore ++;
			}
			
			// Keep mutations if score improved; 0.2 chance of keeping mutations if score is same, otherwise revert
			boolean revert_mutation = false;
			if (score < compareToScore || compareToScore > 1e9 ) {
				//Keep the mutations.
			} else if (score == compareToScore) {
				if (int_urn(1,100) <= 80) {
					//Revert.
					revert_mutation = true;
				}
				//20% chance of keeping if no worse.
			} else if (score > compareToScore) {
				//Revert.
				revert_mutation = true;
			}
			
			if (revert_mutation){
				//Reverting the mutation, unfortunately, takes a bit of work.
				
				for (k = num_mut - 1; k >= 0; k--)
					domain[mut_domain][mut_base[k]] = mut_old[k];
				domain_intrinsic[mut_domain] = old_d_intrinsic;

				//Have to recalculate scores:
				for (i = 0; i < num_domain; i++) {
					updateInteractionsToMutated(i,mut_domain,domain,domain_length,interaction,crosstalk);				
				}  
			} else {
				//Keep the mutations

				if(false){
					// Display new bases
					for (k = 0; k < num_mut; k++) {
						DisplayBase(domain[mut_domain][mut_base[k]]);
					}
					// Display domain scores:
					for (i = 0; i < num_domain; i++) {
						System.out.printf("%f", domain_score[i]);
					}	
				}
				best_score = score;
				total_mutations = total_mutations + num_mut;
			}
			// Check for keyboard hit
			if (num_mut_attempts%1000==0){
				DEBUGNOW = 5;
				if (num_mut_attempts > TOTAL_ATTEMPTS){
					break;
				}
				displayDomains(domain);
				System.out.println(best_score);
				//displaySequence(seqToSynthesize.get(0), domain);
				/*
				for(i = 0; i < num_domain; i++){
					System.out.print(interaction[i][0]+" ");
				}
				System.out.println();
				*/
			}
			while (waitForResume && !abort){
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		displayDomains(domain);
	}
	private int DEBUGNOW = -1;
	private void updateDomainScore(double[] domain_score, int i, int num_domain,
			double[][] interaction, double[][] crosstalk,
			int[] domain_importance, double[] domain_intrinsic) {

		domain_score[i] = 0;
		for (int j = 0; j < num_domain; j++) {
			domain_score[i] = Math.max(domain_score[i],interaction[i][j] + domain_importance[i] + domain_importance[j]);
			domain_score[i] = Math.max(domain_score[i],crosstalk[i][j] + domain_importance[i] + domain_importance[j] + CROSSTALK_SCORE);
		}
		domain_score[i] += domain_intrinsic[i] + (double) (i+1) * 0.000001;
	}
	
	/**
	 * Calculates the intrinsic score of a single domain, or a junction of two domains.
	 */
	private void addIntrinsicScore(DomainSequence seq, int mut_domain, int[][] domain, double[] domain_intrinsic, int rule_4g, int rule_6at) {
		int k, q, base;

		int len = seq.length(domain);
		
		// Search for 4g, if rule applied
		if (rule_4g == 1) {
			k = 0; // G-C counter
			for (q = 0; q < len; q++) {
				base = seq.base(q, domain);
				//Look for EITHER GGGG or CCCC. GCGC or any other variant is correctly ignored.
				if ((base == G)&&(k < 100)) 
					k++;
				else if (base == G)
					k = 1;
				if ((base == C)&&(k > 100))
					k++;
				else if (base == C)
					k = 101;
				//Add the penalty to mut_domain.
				if ((k < 100)&&(k > 3))
					domain_intrinsic[mut_domain] +=  GGGG_PENALTY;
				else if (k > 103)
					domain_intrinsic[mut_domain] +=  GGGG_PENALTY;
			}
		}

		// Search for 6at, if rule applied
		if (rule_6at == 1) {
			//look for 6+ in row of A/T
			k = 0; // AT counter
			for (q = 0; q < len; q++) {
				base = seq.base(q, domain);
				if ((base == A)||(base == T))
					k++;
				else
					k = 0;
				if (k > 5)
					domain_intrinsic[mut_domain] += ATATAT_PENALTY;	    
			} // end "junction" loop
			k = 0; // GC counter
			//Also look for 6+ in row of G/C
			for (q = 0; q < len; q++) {
				base = seq.base(q, domain);
				if ((base == G)||(base == C))
					k++;
				else
					k = 0;
				if (k > 5)
					domain_intrinsic[mut_domain] += ATATAT_PENALTY;	    
			}
		}		
	}
	/**
	 * Two DNASequence used for calculating pairscores.
	 */
	private void updateInteractionsToMutated(int i, int mut_domain, int[][] domain, int[] domain_length, double[][] interaction, double[][] crosstalk) {
		DomainSequence dnaSeq1 = new DomainSequence();
		DomainSequence dnaSeq2 = new DomainSequence();
		
		if (i == mut_domain) {
			dnaSeq1.setDomains(i);
			interaction[i][i] = pairscore(dnaSeq1, dnaSeq1, domain) + INTRA_SCORE;
			//Self crosstalk is me against reverse complemented self.
			crosstalk[i][i] = selfcrosstalk(dnaSeq1,domain) / CROSSTALK_DIV + INTRA_SCORE;
			
			//dnaSeq2.setDomains(i | DNA_COMPLEMENT_FLAG);
			//crosstalk[i][i] = pairscore(dnaSeq1,dnaSeq2,domain)/ CROSSTALK_DIV + INTRA_SCORE;
		} else {
			dnaSeq1.setDomains(i);
			dnaSeq2.setDomains(mut_domain);
			interaction[i][mut_domain] = pairscore(dnaSeq1, dnaSeq2, domain);
			interaction[mut_domain][i] = pairscore(dnaSeq2, dnaSeq1, domain);
			
			dnaSeq1.setDomains(i);
			dnaSeq2.setDomains(mut_domain | DNA_COMPLEMENT_FLAG);
			crosstalk[i][mut_domain] = pairscore(dnaSeq1, dnaSeq2, domain) / CROSSTALK_DIV;
			
			dnaSeq1.setDomains(mut_domain);
			dnaSeq2.setDomains(i | DNA_COMPLEMENT_FLAG);
			crosstalk[mut_domain][i] = pairscore(dnaSeq1, dnaSeq2, domain) / CROSSTALK_DIV;
		}
	}
}