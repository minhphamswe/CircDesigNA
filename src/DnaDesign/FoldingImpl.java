package DnaDesign;

import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.DisplayBase;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.T;

import java.util.List;
import java.util.Scanner;

/**
 * Implements MFE prediction and folding score functions
 */
public class FoldingImpl implements NAFolding{


	/**
	 * "Large loop" - what?
	 */
	double LLstr = -0.5; 

	/**
	 * "Score for domain ending in a base pair" - what?
	 */
	double DHstr = 0; 

	int LHbases = 8;
	double LHstart = 2;
	double LHpower = 2;
	int MinHairpinLoopSize = 2; //Smallest a hairpin loop can be. Used in selfcrosstalk
	/**
	 * Use Unafold to evaluate selffolding score?
	 */
	private boolean FOLD_VIA_UNAFOLD = false;
	int rule_4g, rule_6at;
	{
		//DEFAULTS:
		rule_4g = 1; // cannot have 4 G's or 4 C's in a row
		rule_6at = 1; // cannot have 6 A/T bases in a row
	}
	/**
	 * Matching (complementarity) Strengths. Negative means repulsion.
	 */
	double GCstr = 3;
	double ATstr = 2;
	double GTstr = 0.1;
	double MBstr = -3;
	private static final boolean DEBUG_selfCrosstalkMethod = false;
	public FoldingImpl(){
		
	}
	
	private DomainSequence mutateSequence = new DomainSequence();

	public double affectedSequenceInvalidScore(int mut_domain, List<DomainSequence> seqToSynthesize, int[][] domain, int[][] domain_markings) {
		mutateSequence.setDomains(mut_domain,null);
		//Break if we succeed all VALIDITY checks:
		double res = 0;
		res += isValidSequence(mutateSequence,domain,domain_markings);
		//Check all junctions
		for(DomainSequence seq : seqToSynthesize){
			if (seq.contains(mut_domain)){
				//Need to make sure that changing this domain has possibility to fix.
				res += isValidSequence(seq, domain,domain_markings);
			}
		}
		return res;
	}

	/**
	 * If certain rules are applied, this routine will check whether those rules invalidate
	 * the input sequence.
	 * @param domain_markings 
	 */
	private double isValidSequence(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		int k,q,base;
		int len = seq.length(domain);

		double sumResult = 0;
		
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
				else if ((base == C)&&(k > 100))
					k++;
				else if (base == C)
					k = 101;
				else
					k = 0;
				if ((k < 100)&&(k > 3)){
					//System.out.println("4G");
					seq.mark(q, -4, domain, domain_markings);
					sumResult += 1;
					k = 0;
				} else if (k > 103) { 
					seq.mark(q, -4, domain, domain_markings);
					//System.out.println("4C");
					sumResult += 1;
					k = 0;
				}
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
				if (k > 5){
					//System.out.println("6AT");
					seq.mark(q, -6, domain, domain_markings);
					sumResult += 1;
					k = 0;
				}
			} // end "junction" loop
			k = 0; // GC counter
			//Also look for 6+ in row of G/C
			for (q = 0; q < len; q++) {
				base = seq.base(q, domain);
				if ((base == G)||(base == C))
					k++;
				else
					k = 0;
				if (k > 5){
					//System.out.println("6GC");
					seq.mark(q, -6, domain, domain_markings);
					sumResult += 1;
					k = 0;
				}
			}
		}		

		if (false){
			//Ensure A/T and G/C ratio are ~ 50
			int at = 0;
			int gc = 0;
			for (q = 0; q < len; q++) {
				base = seq.base(q, domain);
				if (base==G){
					gc++;
				}
				if (base==C){
					gc--;
				}
				if (base==A){
					at++;
				}
				if (base==T){
					at--;
				}
			}
			if (Math.abs(gc) > (len/2-2) || Math.abs(at) > (len/2-2)){
				sumResult += 1;
				at = 0;
				gc = 0;
			}
		}

		return sumResult;
	}

	public double pairscore(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] problemAreas) {
		if (FOLD_VIA_UNAFOLD){
			//0 is the target (so shift the score to make it 0) for unafold delta G output 
			return Math.max(pairscore_viaUnafold(seq1, seq2, domain, problemAreas) - (0),-1);
		} else {
			//Target is 6.25 
			return Math.max(pairscore_viaMatrix(seq1,seq2,domain,problemAreas)-6.25,-1) ;
		}
	}

	
	public double foldSingleStranded(DomainSequence seq, int[][] domain, int[][] domain_markings){
		if (FOLD_VIA_UNAFOLD){
			//0 is the target (so shift the score to make it 0) for unafold delta G output 
			return Math.max(foldSingleStranded_viaUnafold(seq, domain, domain_markings) - (0),-1);
		} else {
			//2.5 is the target for our little matrix algorithm.
			return Math.max(foldSingleStranded_viaMatrix(seq, domain, domain_markings) - (2.5), -1);
		}
	}
	
	
	/**
	 * Interaction score.
	 */
	private double[][] Cmatrix_pairscore;
	private double[][] Smatrix_pairscore;
	private byte[][] SDmatrix_pairscore;
	double pairscore_viaMatrix(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] problemAreas) {
		// Gives the score of the two sequences's crosstalk
		double score, temp;
		int i, j, k;
		int len1 = seq1.length(domain);
		int len2 = seq2.length(domain);
		if (!(Cmatrix_pairscore!=null && len1 <= Cmatrix_pairscore.length && len2 <= Cmatrix_pairscore[0].length)){
			Cmatrix_pairscore = new double[len1][len2];
			Smatrix_pairscore = new double[len1][len2];
			SDmatrix_pairscore = new byte[len1][len2];
		}
		double[][] Cmatrix = Cmatrix_pairscore; // complementarity matrix
		double[][] Smatrix = Smatrix_pairscore; // score matrix
		byte[][] SDmatrix = SDmatrix_pairscore; // running total of helix size, 0 if current base didn't contribute.

		// LEN1 x LEN2: Complementarity calculation .
		// Note that the binding is going 5'-3' on both strands, and assumes seq2 is in the array 3'-5'! 
		for (i = 0; i < len1; i++) {
			for (j = 0; j < len2; j++) {
				int base1 = seq1.base(i,domain);
				int base2 = seq2.base(len2-1-j,domain);
				Cmatrix[i][j] = DnaDefinition.bindScore(base1, base2);
				if (Cmatrix[i][j]==0.0){
					Cmatrix[i][j] = MBstr; // mismatch
				}
			}
		}

		// Calculate score
		score = 0;


		// Seems to seed the algorithm somehow?
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
						SDmatrix[i][j] = (byte)((SDmatrix[i-1][j-1] + 1)&127);
					} else if ((SDmatrix[i-1][j-1] == 0)&&(Smatrix[i-1][j-1] > 0)) { // starting a new helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];
						SDmatrix[i][j] = 1;
					} else {	  
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
				if (Smatrix[i][j] > score){
					score = Smatrix[i][j];
				}
			}
		}

		//Wow. this is a complicated little thing, isn't it?
		return score;
	}
	double pairscore_viaUnafold(DomainSequence ds, DomainSequence ds2, int[][] domain, int[][] domain_markings) {
		StringBuffer create = new StringBuffer();
		int len = ds.length(domain);
		for(int k = 0; k < len; k++){
			create.append(DisplayBase(ds.base(k, domain)));
		}
		create.append(" ");
		int len2 = ds2.length(domain);
		for(int k = 0; k < len2; k++){
			create.append(DisplayBase(ds2.base(k, domain)));
		}
		
		String str = create.toString();
		try {
			Process p = Runtime.getRuntime().exec(absPathToHybridMinMod+str);
			Scanner in = new Scanner(p.getInputStream());
			double val = 0;
			double PERFECTscore = -20;
			try {
				while(in.hasNextLine()){
					String line = in.nextLine();
					val = new Double(line.split("\\s+")[0]);
					val = Math.max(-val,PERFECTscore);
					return val;
				}
		 	} finally {
				if (val == PERFECTscore){ //"Infinite" Free Energy (?)
					return val;
				}
				for(int k = 0; k < len+len2; k++){
					char[] arr = in.nextLine().toCharArray();
					int regions = 0;
					int z, end;
					for(z = 0; z < arr.length && regions < 4; z++){
						if (arr[z]=='\t'){
							regions++;
						}
					}
					for(end = z+1; end < arr.length; end++){
						if (arr[end]=='\t'){
							break;
						}
					}
					//System.out.println(new String(arr));
					int num = new Integer(new String(arr,z,end-z));
					//System.out.println(num);
					if (num > 0){
						if (k < len){
							//Continuous numbering, for some reason, in hybrid.exe
							ds2.mark(num-1-len, domain, domain_markings);
						} else {
							ds.mark(num-1, domain, domain_markings);
						}
					}
					//Thread.sleep(100);
				}
				in.close();
				p.waitFor();
			}
		} catch (Throwable e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	private static final String absPathToHybridSSMinMod =  "\"C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\unafold\\hybrid-ss-min.exe\" -q ";
	private static final String absPathToHybridMinMod = "\"C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\unafold\\hybrid-min.exe\" -q ";
	
	double foldSingleStranded_viaUnafold(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		StringBuffer create = new StringBuffer();
		int len = seq.length(domain);
		for(int k = 0; k < len; k++){
			create.append(DisplayBase(seq.base(k, domain)));
		}
		String str = create.toString();
		try {
			Process p = Runtime.getRuntime().exec(absPathToHybridSSMinMod+str);
			Scanner in = new Scanner(p.getInputStream());
			double val = 0;
			double PERFECTscore = -20;
			try {
				while(in.hasNextLine()){
					val = new Double(in.nextLine());
					val = Math.max(-val,PERFECTscore);
					return val;
				}
			} finally {
				if (val == PERFECTscore){ //"Infinite" Free Energy (?)
					return val;
				}
				in.nextLine(); //Read off "dg" line
				for(int k = 0; k < len; k++){
					char[] arr = in.nextLine().toCharArray();
					//System.out.println(new String(arr));
					int regions = 0;
					int z, end;
					for(z = 0; z < arr.length && regions < 4; z++){
						if (arr[z]=='\t'){
							regions++;
						}
					}
					for(end = z+1; end < arr.length; end++){
						if (arr[end]=='\t'){
							break;
						}
					}
					//System.out.println(new String(arr));
					int num = new Integer(new String(arr,z,end-z));
					//System.out.println(num);
					if (num > 0){
						seq.mark(num-1, domain, domain_markings);
					}
					//Thread.sleep(100);
				}
				in.close();
				p.waitFor();
			}
		} catch (Throwable e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException();
	}
	
	double foldSingleStranded_viaMatrix(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		double score, temp;
		int i, j, k;
		int len1 = seq.length(domain);
		if (!(Cmatrix_pairscore!=null && len1 <= Cmatrix_pairscore.length && len1 <= Cmatrix_pairscore[0].length)){
			Cmatrix_pairscore = new double[len1][len1];
			Smatrix_pairscore = new double[len1][len1];
			SDmatrix_pairscore = new byte[len1][len1];
		}
		double[][] Cmatrix = Cmatrix_pairscore; // complementarity matrix
		double[][] Smatrix = Smatrix_pairscore; // score matrix
		byte[][] SDmatrix = SDmatrix_pairscore; // running total of helix size, 0 if current base didn't contribute.
		
		// NxN complementarities. 
		for (i = 0; i < len1; i++) {
			for (j = 0; j < len1; j++) {
				int base1 = seq.base(i,domain);
				int base2 = seq.base(len1-1-j,domain);
				Cmatrix[i][j] = DnaDefinition.bindScore(base1, base2);
				if (Cmatrix[i][j]==0.0){
					Cmatrix[i][j] = MBstr; // mismatch
				}
				//Make sure matrices are clear.
				SDmatrix[i][j] = 0;
				Smatrix[i][j] = 0;
			}
		}

		score = 0;
		int bestScoreX = 0, bestScoreY = 0;

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
			if (Smatrix[0][j] > score){
				score = Smatrix[0][j];
				bestScoreX = 0;
				bestScoreY = j;
			}
		}

		iloop: for (i = 1; i < len1; i++) {
			Smatrix[i][0] = Cmatrix[i][0];
			if (Smatrix[i][0] < 0) {
				Smatrix[i][0] = 0;
				SDmatrix[i][0] = 0;
			} else {
				Smatrix[i][0] = Smatrix[i][0] + DHstr;
				SDmatrix[i][0] = 1;
			}
			if (Smatrix[i][0] > score){
				score = Smatrix[i][0];
				bestScoreX = i;
				bestScoreY = 0;
			}

			//How far to go into the "right" of the matrix.
			int jloopMax = len1-i-1 - MinHairpinLoopSize;
			jloop: for (j = 1; j < jloopMax; j++) {

				if (i == j && false) { 
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

					if (true){
					if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + MBstr > 0)&&(Smatrix[i-1][j] + MBstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + MBstr;
					if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] + LLstr > 0)&&(Smatrix[i-1][j] + LLstr > Smatrix[i][j]))
						Smatrix[i][j] = Smatrix[i-1][j] + LLstr;
					}
					
					if (Smatrix[i][j] < 0)
						Smatrix[i][j] = 0;

				} else { // stabilizing base
					Smatrix[i][j] = Cmatrix[i][j];
					SDmatrix[i][j] = 1;

					if ((SDmatrix[i-1][j-1] > 0)&&(Smatrix[i-1][j-1] > 0)) { // continuing a helix
						Smatrix[i][j] = Smatrix[i-1][j-1] + Cmatrix[i][j];// * (SDmatrix[i-1][j-1] - .5);
						SDmatrix[i][j] = (byte)Math.min(SDmatrix[i-1][j-1] + 1,127);
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

					if(true){
						if ((SDmatrix[i-1][j] > 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr > Smatrix[i][j])) {
							Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j] - Cmatrix[i-1][j] + MBstr;
							SDmatrix[i][j] = 1;
						} else if ((SDmatrix[i-1][j] == 0)&&(Smatrix[i-1][j] > 0)&&(Smatrix[i-1][j] + Cmatrix[i][j] > Smatrix[i][j])) {
							Smatrix[i][j] = Smatrix[i-1][j] + Cmatrix[i][j];
							SDmatrix[i][j] = 1;
						}
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
				
				if (Smatrix[i][j] > score){
					score = Smatrix[i][j];
					bestScoreX = i;
					bestScoreY = j;
				}
				
			} 
		}

		/*
		if (DEBUG_selfCrosstalkMethod){
			for(i = 0; i < len1; i++){
				for(j = 0; j < len1; j++){
					//System.out.printf("%8d",SDmatrix[i][j]);
					System.out.printf("%8.3f",Smatrix[i][j]);
					//System.out.print((Cmatrix[i][j]>0?1:0)+" ");
				}
				System.out.println();
			}
		}
		*/
		int x = bestScoreX;
		int y = bestScoreY;
		do {
			if (x < 0 || y < 0 || x >= len1 || y >= len1){
				break;
			}
			//System.out.println(x+" "+y);
			if (SDmatrix[x][y]>0){
				//This is one of the helixes that contributes to the MFE structure.
				//Mark the bases it requires as mutable.
				
				int ty = len1-1-y;
				if (DEBUG_selfCrosstalkMethod)System.out.println(x+"-"+(x-(SDmatrix[x][y]-1))+" "+(ty)+"-"+(ty+SDmatrix[x][y]-1));
				for(int z = -SDmatrix[x][y]+1; z <= 0; z++){
					int nx = x + z;
					int ny = ty - z;
					//Store nx
					seq.mark(nx,domain,domain_markings);
					if (DEBUG_selfCrosstalkMethod)System.out.print(DisplayBase(seq.base(nx, domain)));
					seq.mark(ny,domain,domain_markings);
					if (DEBUG_selfCrosstalkMethod)System.out.print("="+DisplayBase(seq.base(ny, domain))+",");
				}
				if (DEBUG_selfCrosstalkMethod)System.out.println();

				int HelixLength = SDmatrix[x][y];

				x -= HelixLength;
				y -= HelixLength;
			} else {
				//Backtrack the mismatch region... (Right to left, if no helix found, go up the matrix)
				
				/*
				double bestScore = Smatrix[x][y];
				int bestI = y;
				for(i = y-1; i >=0; i--){
					if (Smatrix[x][i]>bestScore){
						bestI = i;
						bestScore=Smatrix[x][i];
					}
				}
				y = bestI;

				if (SDmatrix[x][y]==0){
					x--;
					y--;
				}
				*/
				x--;
				y--;
			}
		} while(true);
		if (DEBUG_selfCrosstalkMethod)System.out.println();
		if (DEBUG_selfCrosstalkMethod)System.out.println(bestScoreX+" "+bestScoreY+" "+score);


		return score;
	}

}
