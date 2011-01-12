package DnaDesign.impl;

import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.T;
import static DnaDesign.DnaDefinition.displayBase;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import DnaDesign.DnaDefinition;
import DnaDesign.DomainSequence;
import DnaDesign.ExperimentDatabase;
import DnaDesign.NAFolding;

/**
 * Implements MFE prediction and folding score functions
 */
public class FoldingImpl implements NAFolding{

	/**
	 * "Large loop" - what?
	 */
	double LLstr = -0.5; 

	/**
	 * "Score for domain ending in a base pair"
	 */
	double DHstr = 3; 

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
	public static boolean DEBUG_selfCrosstalkMethod = false;
	private ExperimentDatabase eParams;
	public FoldingImpl(){
		this(new ExperimentalDuplexParamsImpl());
	}
	public FoldingImpl(ExperimentDatabase params){
		eParams = params;
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
			return Math.max(pairscore_viaUnafold(seq1, seq2, domain, problemAreas) - (0), -1);
		} else {
			return Math.max(pairscore_viaMatrix(seq1,seq2,domain,problemAreas) - (0), 0) ;
		}
	}

	
	public double foldSingleStranded(DomainSequence seq, int[][] domain, int[][] domain_markings){
		if (FOLD_VIA_UNAFOLD){
			//0 is the target (so shift the score to make it 0) for unafold delta G output 
			return Math.max(foldSingleStranded_viaUnafold(seq, domain, domain_markings) - (0),-1);
		} else {
			return Math.max(foldSingleStranded_viaMatrix(seq, domain, domain_markings) - (0), 0);
		}
	}
	
	
	/**
	 * Interaction score.
	 */
	private DomainSequence reverseStrand_shared = new DomainSequence();
	private double[][] compMatrix_shared;
	private double[][] sMatrix_shared;
	private double[][] gamma3mat_shared;
	/**
	 * Each leaf (a pair of ints a,b) is either:
	 * a>0: a is  the length of a helix
	 * a<0: -a is the length of the left loop, -b is the length of the right loop.
	 */
	private int[][][] sdMatrix_shared;
	private int[][] sdMatrix_old_shared;
	public void ensureSharedMatrices(int len1, int len2){
		if (!(compMatrix_shared!=null && len1 <= compMatrix_shared.length && len2 <= compMatrix_shared[0].length)){
			compMatrix_shared = new double[len1][len2];
			sMatrix_shared = new double[len1][len2];
			sdMatrix_shared = new int[len1][len2][2];
			sdMatrix_old_shared = new int[len1][len2];
			gamma3mat_shared = new double[len1][len2];
		}
	}
	public double pairscore_viaUnafold(DomainSequence ds, DomainSequence ds2, int[][] domain, int[][] domain_markings) {
		StringBuffer create = new StringBuffer();
		int len = ds.length(domain);
		for(int k = 0; k < len; k++){
			create.append(displayBase(ds.base(k, domain)));
		}
		create.append(" ");
		int len2 = ds2.length(domain);
		for(int k = 0; k < len2; k++){
			create.append(displayBase(ds2.base(k, domain)));
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
	
	public double foldSingleStranded_viaUnafold(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		StringBuffer create = new StringBuffer();
		int len = seq.length(domain);
		for(int k = 0; k < len; k++){
			create.append(displayBase(seq.base(k, domain)));
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

	/**
	 * Not multithreaded.
	 * @param len2 
	 * @param seq2 
	 */
	private void foldSingleStranded_makeCMatrix(DomainSequence seq, DomainSequence seq2, int len1,int len2, int[][] domain){
		double[][] Cmatrix = compMatrix_shared; // complementarity matrix
		double[][] Smatrix = sMatrix_shared; // score matrix
		int[][][] SDmatrix = sdMatrix_shared; // running total of helix size, 0 if current base didn't contribute.
		// NxN complementarities. 
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2; j++) {
				int base1 = seq.base(i,domain);
				int base2 = seq2.base(j,domain);
				Cmatrix[i][j] = DnaDefinition.bindScore(base1, base2);
				if (Cmatrix[i][j]==0.0){
					Cmatrix[i][j] = MBstr; // mismatch
				}
				Smatrix[i][j] = 0;
				SDmatrix[i][j][0] = 0;
				SDmatrix[i][j][1] = 0;
			}
		}		
	}
	public double pairscore_viaMatrix(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] domain_markings) {
		return foldNA_viaMatrix(seq1, seq2, domain, domain_markings, true);
	}
	public double foldSingleStranded_viaMatrix(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		return foldNA_viaMatrix(seq, seq, domain, domain_markings, false);
	}
	public double foldNA_viaMatrix(DomainSequence seq, DomainSequence seq2, int[][] domain, int[][] domain_markings, boolean fullMatrix) {
		int len1 = seq.length(domain);
		int len2 = seq2.length(domain);
		ensureSharedMatrices(len1,len2);
		foldSingleStranded_makeCMatrix(seq,seq2,len1,len2,domain);
		double[][] Cmatrix = compMatrix_shared; // complementarity matrix
		double[][] Smatrix = sMatrix_shared; // score matrix
		int[][][] SDmatrix = sdMatrix_shared; // running total of helix size, 0 if current base didn't contribute.
		double[][] gamma3mat = gamma3mat_shared;
		
		double score = 0;
		int bestI = -1, bestJ = -1;
		//Positive delta G.
		//LCS, use eParams to do intelligent weighting.
		int minHairpinSize = 1+3; //Only used in singlestranding.
		//Calculate looping bounds.
		int i;
		if (fullMatrix){
			i = len1-1;
		} else {
			i = (len1-1)-minHairpinSize;
		}
		for(; i >= 0; i--){
			int j;
			if (fullMatrix){
				j = 0;
			} else {
				j = i+minHairpinSize;
			}
			for(; j < len2; j++){
				//Left loop (m), + no bonus
				double gamma1 = foldSingleStranded_calcGamma1(i,j,len1,Cmatrix,Smatrix,SDmatrix);
				//Right loop (n), + no bonus
				double gamma2 = foldSingleStranded_calcGamma2(i,j,Cmatrix,Smatrix,SDmatrix);
				//Helix, + dummy score if new helix, - dummy score if 2nd base, + nn score is length >= 2.
				//If beginning new helix, have to add the score of the last loop.
				double gamma3 = foldSingleStranded_calcGamma3(len1,len2,seq,seq2,domain,i,j,Cmatrix,Smatrix,SDmatrix,true);
				gamma3mat[i][j] = gamma3; //Store the gamma3 result for tracebacking use.
				//Greedy algorithm: take the biggest one (Is this provably correct?).
				Smatrix[i][j]=max(gamma1,gamma2,gamma3);
				//Which did we do?
				if (gamma3 == Smatrix[i][j]){
					//We are doing helix. 
					//SDmatrix[i][j] = Math.max(0,SDmatrix[i+1][j-1])+1;
					//Leave the setting of the SDMatrix up to "calcGamma3". It must therefore run LAST in the above 3 seqs.
				} else if (gamma1 == Smatrix[i][j]){
					//Continuing loop, overwrite gamma3 calculation
					SDmatrix[i][j][0] = Math.min(0,SDmatrix[i+1][j][0])-1;
					SDmatrix[i][j][1] = Math.min(0,SDmatrix[i+1][j][1]);
				} else if (gamma2 == Smatrix[i][j]){
					//Continuing loop, overwrite gamma3 calculation
					SDmatrix[i][j][0] = Math.min(0,SDmatrix[i][j-1][0]);
					SDmatrix[i][j][1] = Math.min(0,SDmatrix[i][j-1][1])-1;
				} else {
					throw new RuntimeException("Assertion failure. foldSingleStranded_viaMatrix inner loop of filling.");
				}
				//Have we found the best structure?
				if (Smatrix[i][j] > score){
					score = Smatrix[i][j];
					bestI = i;
					bestJ = j;
				}
			}
		}
		
		//Traceback?
		double overCount = foldSingleStranded_traceBack(len1,len2,Smatrix,Cmatrix,gamma3mat,bestI,bestJ,seq,seq2,domain,domain_markings,!fullMatrix);
		
		if (debugLCSAlgorithm){
			for(int k = 0; k < len1; k++){
				for(int y = 0; y < len2; y++){
					System.out.printf(" (%3d,%3d)",SDmatrix[k][y][0],SDmatrix[k][y][1]);
				}
				System.out.println();
			}
			for(int k = 0; k < len1; k++){
				for(int y = 0; y < len2; y++){
					System.out.printf(" %4.8f",Smatrix[k][y]);
				}
				System.out.println();
			}
		}
		
		return score-overCount;
	}
	/**
	 * Performs the standard nussinov tracebacking, sans bifurcation tracing (thus, no stack).
	 */
	private double foldSingleStranded_traceBack(int len1, int len2, double[][] Smatrix, double[][] Cmatrix, double[][] gamma3mat, int bestI,
			int bestJ, DomainSequence seq, DomainSequence seq2,
			int[][] domain, int[][] domain_markings, boolean isSingleStrandFold) {
		int helixLength = 0;
		MFE_numBasesPaired = 0;
		MFE_longestHelixLength = 0;
		MFE_pointlist.clear();
		double overCount = 0;
		while(true){
			//Break condition:
			//System.out.println(bestI+" "+bestJ);
			if (bestI>=len1 || bestJ < 0){
				break;
			}
			boolean isOnFringeOfMap = bestI==len1-1 || bestJ==0;
			if (isSingleStrandFold){
				if (bestJ<bestI){
					break;
				}
			}
			MFE_pointlist.add(new Point(bestI,bestJ));
			boolean inHelix = false; //only if we go diagonally do we mark a duplex pair.
			if (!isOnFringeOfMap){
				double gamma1 = Smatrix[bestI+1][bestJ];
				double gamma2 = Smatrix[bestI][bestJ-1];
				double gamma3 = Smatrix[bestI+1][bestJ-1];
				double best = max(gamma1,gamma2,gamma3);
				//This is like the Zuker algorithm, we are keeping a matrix that
				//maintains the score we would get if we always went diagonals (no bulge handling).
				if (gamma3mat[bestI][bestJ]>best){
					inHelix=true;
				}
				else if (gamma1 == best){
					//Go there.
					bestI++;
				}
				else if (gamma2 == best){
					//Go there.
					bestJ--;
				}
				else if (gamma3==best){
					inHelix = true;
				} 
				else {
					throw new RuntimeException("Assertion failure. foldSingleStranded_traceback in best check");
				}
			} else {
				inHelix = true; //Go ahead and check the last guy.
			}
			if (inHelix){
				helixLength ++;
				MFE_longestHelixLength = Math.max(MFE_longestHelixLength,helixLength);
				//Mark condition:
				if (Cmatrix[bestI][bestJ]>0){
					//Complementary pair in traceback. mark.
					seq.mark(bestI, domain, domain_markings);
					seq2.mark(bestJ, domain, domain_markings);
				}
				//Go helix!
				bestI++;
				bestJ--;
				MFE_numBasesPaired++;
			} else {
				helixLength = 0;
			}
		}
		if (MFE_longestHelixLength==1){
			overCount += foldSingleStranded_oneBaseScore;
		}
		return overCount;
	}

	private int MFE_longestHelixLength = -1, MFE_numBasesPaired = -1;
	private ArrayList<Point> MFE_pointlist = new ArrayList();
	/**
	 * WARNING: allocates a new matrix.
	 */
	public double[][] getNussinovMatrixScore(int len1, int len2) {
		double[][] nussinovScores = new double[len1][len2];
		for(int y = 0; y < len1; y++){
			for(int x = 0; x < len2; x++){
				nussinovScores[y][x] = sMatrix_shared[y][x];
			}
		}
		return nussinovScores;
	}
	/**
	 * WARNING: allocates a new matrix.
	 */
	public int[][][] getMFEStructureMatrix() {
		return sdMatrix_shared;
	}
	/**
	 * WARNING: allocates a new list.
	 */
	public ArrayList<Point> getTraceback() {
		ArrayList<Point> toRet = new ArrayList<Point>();
		toRet.addAll(MFE_pointlist);
		return toRet;
	}
	public int getLongestHelixLength() {
		return MFE_longestHelixLength;
	}
	public int getNumBasesPaired() {
		// TODO Auto-generated method stub
		return MFE_numBasesPaired;
	}
	private boolean debugLCSAlgorithm = false;
	private double foldSingleStranded_oneBaseScore = .25;
	private double foldSingleStranded_calcDummyScore = .25;
	private double foldSingleStranded_endHelixPenalty = 0;
	private double foldSingleStranded_helixBaseScore = 0;
	private double foldSingleStranded_calcGamma1(int i, int j, int len1, double[][] cMatrix, double[][] sMatrix, int[][][] sdMatrix) {
		//This is the number, if we are a "bulge" and defer to the helix in sMatrix[i+1][j].
		if (i+1>=len1){
			//Off the map.
			return 0.0;
		}
		double bulgeScore = sMatrix[i+1][j];
		//Be sure to remove dummyScore
		if (sdMatrix[i+1][j][0]==1){
			bulgeScore-=foldSingleStranded_calcDummyScore;
		}
		if (sdMatrix[i+1][j][0]>0){
			bulgeScore-=foldSingleStranded_endHelixPenalty;
		}
		return bulgeScore;
	}
	private double foldSingleStranded_calcGamma2(int i, int j, double[][] cMatrix, double[][] sMatrix, int[][][] sdMatrix) {
		if (j-1<0){
			//Off the map.
			return 0.0;
		}
		double bulgeScore = sMatrix[i][j-1];
		//Be sure to remove dummyScore
		if (sdMatrix[i][j-1][0]==1){
			bulgeScore-=foldSingleStranded_calcDummyScore;
		}
		if (sdMatrix[i][j-1][0]>0){
			bulgeScore-=foldSingleStranded_endHelixPenalty;
		}
		return bulgeScore;
	}
	/**
	 * Helix, + dummy score if new helix, - dummy score if 2nd base, + nn score is length >= 2.
	 * If beginning new helix, have to add the score of the last loop.
	 * 
	 * Both seq and seq2 should be in 5'-3' order.
	 */
	private double foldSingleStranded_calcGamma3(int len1, int len2, DomainSequence seq, DomainSequence seq2, int[][] domain, int i, int j, double[][] cMatrix, double[][] sMatrix, int[][][] sdMatrix, boolean writeToSD) {
		double dummyScore = foldSingleStranded_calcDummyScore;
		boolean onFringeOfMap = i+1>=len1 || j-1<0;
		if (cMatrix[i][j]>0){
			//This is a pair. Extend helix
			if (writeToSD){
				sdMatrix[i][j][0] = Math.max(0,onFringeOfMap?0:sdMatrix[i+1][j-1][0])+1;
				sdMatrix[i][j][1] = sdMatrix[i][j][0]; //MUST set this > 0. Have to seperate loop counters from helix counters!
			}
			//New helix
			if (onFringeOfMap || sdMatrix[i+1][j-1][0]<=0){
				if(!onFringeOfMap && (sdMatrix[i+1][j-1][0]<0 || sdMatrix[i+1][j-1][1]<0)){
					//Ending a loop, of length > 0
					int leftLoopSize = -sdMatrix[i+1][j-1][0];
					int rightLoopSize = -sdMatrix[i+1][j-1][1];
				}
				return (onFringeOfMap?0:sMatrix[i+1][j-1])+dummyScore;
			}
			//Continuing old helix
			else {
				//get NN score.
				double nn = eParams.getNNdeltaG(seq.base(i+1, domain), seq2.base(j-1,domain),seq.base(i,domain), seq2.base(j,domain));
				double helixScore = sMatrix[i+1][j-1];
				if (sdMatrix[i+1][j-1][0]==1){
					//Remove dummy score
					helixScore -= dummyScore;
					helixScore += foldSingleStranded_helixBaseScore;
				}
				//Delta gs are negative, so subtract instead.
				helixScore -= nn;
				return helixScore;
			}
		} else {
			//No helix. Extend both left and right bulges by one.
			if (writeToSD){
				sdMatrix[i][j][0] = Math.min(0,onFringeOfMap?0:sdMatrix[i+1][j-1][0])-1;
				sdMatrix[i][j][1] = Math.min(0,onFringeOfMap?0:sdMatrix[i+1][j-1][1])-1;
			}
			if (onFringeOfMap){
				return 0.0;
			}
			//Ending old helix?
			if (sdMatrix[i+1][j-1][0]>0){
				if (sdMatrix[i+1][j-1][0]==1){
					//Remove dummy score, replace with 1 base score
					return sMatrix[i+1][j-1]-dummyScore+foldSingleStranded_oneBaseScore;
				} else {
					//Add terminal score.
					double terminalMismatch = eParams.getNNdeltaGterm(seq.base(i+1, domain), seq2.base(j-1,domain),seq.base(i,domain), seq2.base(j,domain));
					//Delta gs are negative, so subtract instead.
					return sMatrix[i+1][j-1]-terminalMismatch;
				}
			} 
			//Continuing loop region
			else {
				return sMatrix[i+1][j-1];
			}
		}
	}
	private double max(double gamma1, double gamma2, double gamma3) {
		return Math.max(Math.max(gamma1,gamma2),gamma3); 
	}

}
