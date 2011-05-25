package circdesigna.energy;

import java.awt.Point;
import java.util.ArrayList;

import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.Config.CircDesigNASystemElement;
import edu.utexas.cssb.circdesigna.DomainSequence;

/**
 * Implements MFE prediction and folding score functions, using the MCS algorithm as an approximation.
 */
public class CircDesigNAMCSFolder extends CircDesigNASystemElement implements NAFolding{

	/**
	 * This is an O(N^2) approximation which essentially takes the maximum weighted subsequence
	 * of one strand that is complementary to the other, where weight is assigned using the maximum
	 * neighbor model (with terminal nearest neighbors accounted for)
	 */
	/**
	 * Constructors, define parameters and / or a configuration.
	 */
	private NAExperimentDatabase eParams;
	public CircDesigNAMCSFolder(CircDesigNAConfig sys){
		super(sys);
		eParams = new ExperimentalDuplexParamsImpl(sys);
	}
	public CircDesigNAMCSFolder(NAExperimentDatabase params, CircDesigNAConfig sys){
		super(sys);
		eParams = params;
	}
	

//////////////////////// "MCS" algorithm
//////////////////////// Essentially applies Longest Common Subsequence recursion 
//////////////////////// to two DNA strands, running one of the strands backwards (antisense)
	/**
	 * These turn any folding request into an application of a simplistic O(N^2) maximum-weighted-subsequence
	 * finding problem. foldNA_viaMatrix implements the actual algorithm. 
	 */
	public double mfeNoDiag(DomainSequence domainSequence, DomainSequence domainSequence2, int[][] domain, int[][] domain_markings){
		FoldNA_viaMatrix_Options mfeNoDiagonalPairsOpt = new FoldNA_viaMatrix_Options();
		mfeNoDiagonalPairsOpt.foldFullMatrix = true;
		mfeNoDiagonalPairsOpt.suppressDiagonalScores = true;
		//Run the generic folding algorithm under these conditions
		return foldNA_viaMatrix(domainSequence, domainSequence2, domain, domain_markings, mfeNoDiagonalPairsOpt);
	}
	
	public double mfe(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] domain_markings) {
		FoldNA_viaMatrix_Options mfe_viaMatrix_opt = new FoldNA_viaMatrix_Options();
		mfe_viaMatrix_opt.foldFullMatrix = true;
		//Run the generic folding algorithm under these conditions
		return foldNA_viaMatrix(seq1, seq2, domain, domain_markings, mfe_viaMatrix_opt);
	}
	public double mfe(DomainSequence seq, int[][] domain, int[][] domain_markings) {
		FoldNA_viaMatrix_Options mfe_viaMatrix_opt = new FoldNA_viaMatrix_Options();
		mfe_viaMatrix_opt.foldFullMatrix = false;
		//Run the generic folding algorithm under these conditions
		return foldNA_viaMatrix(seq, seq, domain, domain_markings, mfe_viaMatrix_opt);
	}
	
	/**
	 * Interaction score shared memory buffers
	 */
	private float[][] sMatrix_shared;
	/**
	 * Each leaf (a pair of ints a,b) is either:
	 * a>0: a is  the length of a helix
	 * a<0: -a is the length of the left loop, -b is the length of the right loop.
	 */
	private int[][][] sdMatrix_shared;
	public void ensureSharedMatrices(int len1, int len2){
		if (!(sMatrix_shared!=null && len1 <= sMatrix_shared.length && len2 <= sMatrix_shared[0].length)){
			sMatrix_shared = new float[len1][len2];
			sdMatrix_shared = new int[2][len2][2];
		}
	}
	/**
	 * fills the entire folding matrix with 0s. 
	 */
	private void foldSingleStranded_flushMatrixes(DomainSequence seq, DomainSequence seq2, int len1,int len2, int[][] domain){
		// NxN complementarities. 
		for (int i = 0; i < len1; i++) {
			for (int j = 0; j < len2; j++) {
				//int base1 = seq.base(i,domain);
				//int base2 = seq2.base(j,domain);
				flushFoldMatrices(i,j);
			}
		}
	}
	private void flushFoldMatrices(int i, int j){
		sMatrix_shared[i][j] = 0;
		sdMatrix(sdMatrix_shared,i,j)[0] = 0;
		sdMatrix(sdMatrix_shared,i,j)[1] = 0;
	}

	
	/**
	 * Shortcut, simply computes the free energy score of interaction assuming that seq and seq2
	 * form a base-for-base helix. Note that the two input domain sequences must have the same length. 
	 */
	public double mfeStraight(DomainSequence seq, DomainSequence seq2, int[][] domain, int[][] domain_markings, int markStart, int markEnd, int jOffset) {
		int len1 = seq.length(domain);
		int len2 = seq2.length(domain);
		ensureSharedMatrices(len1,len2);
		float[][] Smatrix = sMatrix_shared; // score matrix
		int[][][] SDmatrix = sdMatrix_shared; // running total of helix size, 0 if current base didn't contribute.
		
		double best = 0;
		int bestI = -1, bestJ = -1;
		int helix = 0;
		for(int i = len1-1, j = -jOffset; i >= 0 && j < len2;i--, j++){
			if (j < 0){
				continue;
			}
			float gamma3 = (float) foldSingleStranded_calcGamma3(len1,len2,seq,seq2,domain,i,j,Smatrix,SDmatrix,true);
			Smatrix[i][j]=gamma3;
			if (Smatrix[i][j] < best){
				helix++;
				bestI = i;
				bestJ = j;
				best = Smatrix[i][j];
				if (Std.monomer.bindScore(base(seq,i,domain), base(seq2,j,domain)) < 0){
					if (i >= markStart && i < markEnd){
						seq.mark(i, domain, domain_markings);
						seq2.mark(j, domain, domain_markings);
					}
				}
			} else {
				helix = 0;
			}
		}
		return best;
	}
	public static class FoldNA_viaMatrix_Options {
		//True for hybridizations (versus self folding)
		public boolean foldFullMatrix;
		//Used for the "Self Similarity" Check.
		public boolean suppressDiagonalScores;
	}
	/**
	 * Implements an O(N^2) maximum weighted subsequence finding algorithm. The algorithm is familiar
	 * to one who has looked at the algorithm commonly used to solve the Longest Common Subsequence problem.
	 * It is inspired by the folding algorithm used in David Zhang's Domain Designer.
	 * 
	 * This algorithm applies a simple kernel to each cell to update its value, and then returns the largest value of 
	 * any cell in the entire matrix. The kernel is equivalent to considering the minimum free energy of three cases:
	 * dg3: Base i is paired with j, and 
	 * dg1: Base i is not paired with j, break the helix from (i+1,j) or just continue that bulge
	 * dg2: Base i is not paired with j, break the helix from (i,j-1) or just continue that bulge
	 * 
	 * Thus, the algorithm will find the minimum free energy structure when loop contributions are not taken into account.
	 * To prevent impossibly small loops, self-folding evaluation begins 4 spaces right of the diagonal (so the smallest
	 * hairpin considered has size 3). This algorithm does not consider structures which contain "bifurcations", that is, 
	 * multiloops are not supported. 
	 * 
	 * In defense of this algorithm, if any long helix exists in seq with seq2, it will be located. Thus, for the purposes
	 * of removing interactions from seq and seq2, this is good enough. Additionally, it appears that remedying the flaws described above
	 * requires upping the performance of this algorithm to O(N^3), which makes designing large (>6000 base) DNA origami scaffolds impractical.
	 */
	public double foldNA_viaMatrix(DomainSequence seq, DomainSequence seq2, int[][] domain, int[][] domain_markings, FoldNA_viaMatrix_Options options) {
		int len1 = seq.length(domain);
		int len2 = seq2.length(domain);
		if (options.suppressDiagonalScores){
			if (len1!=len2){
				throw new RuntimeException("Diagonal scores can only suppressed when folding strands of equal length");
			}
		}
		ensureSharedMatrices(len1,len2);
		/*
		for(int i = 0; i < len1; i++){
			for(int j = 0; j < len2; j++){
				flushFoldMatrices(i,j);
			}
		}
		*/
		//foldSingleStranded_makeCMatrix(seq,seq2,len1,len2,domain);
		float[][] Smatrix = sMatrix_shared; // score matrix
		int[][][] SDmatrix = sdMatrix_shared; // running total of helix size, 0 if current base didn't contribute.
		
		//Minimum hairpin size of 3, so distance from diagonal is 4
		int minHairpinSize = 1+3;
		
		double score = 0, gamma1, gamma2, gamma3, pick;
		int bestI = -1, bestJ = -1;
		//Only used in the single stranded version
		//Calculate looping bounds.
		for(int i = len1-1; i >= 0; i--){
			int j;
			if (options.foldFullMatrix){
				j = 0;
			} else {
				//assertion for selffolding
				if (len1!=len2){throw new RuntimeException();};
				//warning! relies on value of minhairpinsize
				j = i+minHairpinSize;
				for(int o = i; o < j && o < len2; o++){
					flushFoldMatrices(i,o);
				}
			}
			for(; j < len2; j++){
				//Left loop (m), + no bonus
				gamma1 = foldSingleStranded_calcGamma1(i,j,len1,Smatrix,SDmatrix);
				//Right loop (n), + no bonus
				gamma2 = foldSingleStranded_calcGamma2(i,j,Smatrix,SDmatrix);
				//Helix, + dummy score if new helix, - dummy score if 2nd base, + nn score is length >= 2.
				//If beginning new helix, have to add the score of the last loop.
				boolean computeHelixScore = true;
				if (options.suppressDiagonalScores){
					if (i==len2-1-j){
						computeHelixScore = false;
					}
				}
				if(computeHelixScore){
					gamma3 = foldSingleStranded_calcGamma3(len1,len2,seq,seq2,domain,i,j,Smatrix,SDmatrix,true);
				} else { 
					gamma3 = 0;
				}
				//Greedy algorithm: take the most minimal (proof: addititivity of delta G, optimization of a sum)
				pick = Math.min(gamma1,Math.min(gamma2,gamma3));
				//If there is a tie, use the following priority:
				if (gamma3 == pick){
					//Continuing helix, calcGamma3 autoincrements the helix length.
					//SDmatrix[i][j] = Math.max(0,SDmatrix[i+1][j-1])+1;
					//Leave the setting of the SDMatrix up to "calcGamma3". It must therefore run LAST in the above 3 seqs.
				} else if (gamma1 == pick){
					//Continuing loop, fix backtracking info
					sdMatrix(SDmatrix,i,j)[0] = Math.min(0,sdMatrix(SDmatrix,i+1,j)[0])-1; //Negative means longer loop
					sdMatrix(SDmatrix,i,j)[1] = Math.min(0,sdMatrix(SDmatrix,i+1,j)[1])-1; //Negative means longer loop
				} else if (gamma2 == pick){
					//Continuing loop, fix backtracking info
					sdMatrix(SDmatrix,i,j)[0] = Math.min(0,sdMatrix(SDmatrix,i,j-1)[0])-1; //Negative means longer loop
					sdMatrix(SDmatrix,i,j)[1] = Math.min(0,sdMatrix(SDmatrix,i,j-1)[1])-1; //Negative means longer loop
				} else {
					throw new RuntimeException("Assertion failure. foldSingleStranded_viaMatrix inner loop of filling.");
				}
				//Keep track of MFE.
				Smatrix[i][j]= (float) pick;
				if (Smatrix[i][j] < score){
					score = Smatrix[i][j];
					bestI = i;
					bestJ = j;
				}
			}
		}
		
		//Traceback.
		double overCount = foldSingleStranded_traceBack(len1,len2,bestI,bestJ,seq,seq2,domain,domain_markings,!options.foldFullMatrix);
		
		if (debugLCSAlgorithm){
			/*
			for(int k = 0; k < len1; k++){
				for(int y = 0; y < len2; y++){
					System.out.printf(" (%3d,%3d)",sdMatrix(SDmatrix,k,y)[0],sdMatrix(SDmatrix,k,y)[1]);
				}
				System.out.println();
			}
			*/
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
	private double foldSingleStranded_traceBack(int len1, int len2, int bestI,
			int bestJ, DomainSequence seq, DomainSequence seq2,
			int[][] domain, int[][] domain_markings, boolean isSingleStrandFold) {
		
		int helixLength = 0;
		
		MFE_numBasesPaired = 0;
		MFE_longestHelixLength = 0;
		MFE_pointlist.clear();
		boolean inHelix = true;
		
		while(true){
			//System.out.println(inHelix+" "+bestI+" "+bestJ+" "+Arrays.toString(domain_markings[0]));
			//Break condition:
			//System.out.println(bestI+" "+bestJ);
			if (bestI>=len1 || bestJ < 0){
				break;
			}
			boolean isOnFringeOfMap;
			if (isSingleStrandFold){
				isOnFringeOfMap = bestJ<=bestI;
			} else {
				isOnFringeOfMap = bestI==len1-1 || bestJ==0;
			}
			MFE_pointlist.add(new Point(bestI,bestJ));
			if (isOnFringeOfMap){
				if (inHelix){
					if (Std.monomer.bindScore(base(seq,bestI,domain), base(seq2,bestJ,domain)) < 0){
						helixLength++;
						MFE_longestHelixLength = Math.max(MFE_longestHelixLength,helixLength);
					}
				}
			}
			if (inHelix && isOnFringeOfMap){ 
				if (Std.monomer.bindScore(base(seq,bestI,domain), base(seq2,bestJ,domain)) < 0){
					seq.mark(bestI, domain, domain_markings);
					seq2.mark(bestJ, domain, domain_markings);
				}
			}
			if (isOnFringeOfMap){
				break;
			}
			//inHelix = false;
			float gamma1 = sMatrix_shared[bestI+1][bestJ];
			float gamma2 = sMatrix_shared[bestI][bestJ-1];
			float gamma3 = sMatrix_shared[bestI][bestJ];

			float best = Math.min(gamma1,Math.min(gamma2,gamma3));
			if (gamma1 == best){
				//Go there.
				inHelix = false;
				bestI++;	
			}
			else if (gamma2 == best){
				//Go there.
				inHelix = false;
				bestJ--;
			} else if (gamma3 == best){
				if (Std.monomer.bindScore(base(seq,bestI,domain), base(seq2,bestJ,domain)) < 0){
					seq.mark(bestI, domain, domain_markings);
					seq2.mark(bestJ, domain, domain_markings);
				}
				inHelix = true;
			}
			else {
				throw new RuntimeException("Assertion failure. foldSingleStranded_traceback in best check");
			}
			if (inHelix){
				//Mark condition:
				helixLength ++;
				MFE_longestHelixLength = Math.max(MFE_longestHelixLength,helixLength);
				//Go helix!
				bestI++;
				bestJ--;
				MFE_numBasesPaired++;
			} else {
				helixLength = 0;
			}
		}	
		if (MFE_longestHelixLength==1){ //Ended on a single base.
			return foldSingleStranded_calcDummyScore;
		}
		return 0;
	}

	private int MFE_longestHelixLength = -1, MFE_numBasesPaired = -1;
	private ArrayList<Point> MFE_pointlist = new ArrayList();
	/**
	 * Allocates a new matrix. Used for visually displaying the results of the folding.
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
	 * WARNING: allocates a new list.
	 */
	public ArrayList<Point> getTraceback() {
		ArrayList<Point> toRet = new ArrayList<Point>();
		toRet.addAll(MFE_pointlist);
		return toRet;
	}
	private boolean debugLCSAlgorithm = false;
	private double foldSingleStranded_calcDummyScore = -.25;
	private double foldSingleStranded_calcGamma1(int i, int j, int len1, float[][] sMatrix, int[][][] sdMatrix) {
		//This is the number, if we are a "bulge" and defer to the helix in sMatrix[i+1][j].
		if (i+1>=len1){
			//Off the map.
			return 0.0;
		}
		double bulgeScore = sMatrix[i+1][j];
		//Be sure to remove dummyScore
		if (sdMatrix(sdMatrix,i+1,j)[0]==1){
			bulgeScore-=foldSingleStranded_calcDummyScore;
		}
		return bulgeScore;
	}
	private double foldSingleStranded_calcGamma2(int i, int j, float[][] sMatrix, int[][][] sdMatrix) {
		if (j-1<0){
			//Off the map.
			return 0.0;
		}
		double bulgeScore = sMatrix[i][j-1];
		//Be sure to remove dummyScore
		if (sdMatrix(sdMatrix,i,j-1)[0]==1){
			bulgeScore-=foldSingleStranded_calcDummyScore;
		}
		return bulgeScore;
	}
	
	/**
	 * Helix, + dummy score if new helix, - dummy score if 2nd base, + nn score is length >= 2.
	 * If beginning new helix, have to add the score of the last loop.
	 * 
	 * Both seq and seq2 should be in 5'-3' order.
	 */
	private double foldSingleStranded_calcGamma3(int len1, int len2, DomainSequence seq, DomainSequence seq2, int[][] domain, int i, int j, float[][] sMatrix, int[][][] sdMatrix, boolean writeToSD) {
		double dummyScore = foldSingleStranded_calcDummyScore;
		boolean onFringeOfMap = i+1>=len1 || j-1<0;
		if (Std.monomer.bindScore(base(seq,i,domain), base(seq2,j,domain)) < 0){
			//This is a pair. Extend helix
			if (writeToSD){
				sdMatrix(sdMatrix,i,j)[0] = Math.max(0,onFringeOfMap?0:sdMatrix(sdMatrix,i+1,j-1)[0])+1;
				sdMatrix(sdMatrix,i,j)[1] = sdMatrix(sdMatrix,i,j)[0]; //MUST set this > 0. Have to seperate loop counters from helix counters!
			}
			//New helix
			if (onFringeOfMap || sdMatrix(sdMatrix,i+1,j-1)[0]<=0){
				double addLoopOpeningPenalty = 0;
				if(!onFringeOfMap && sdMatrix(sdMatrix,i+1,j-1)[0]<0){
					//Ending a loop, of length > 0
					int leftLoopSize = -sdMatrix(sdMatrix,i+1,j-1)[0];
					int rightLoopSize = -sdMatrix(sdMatrix,i+1,j-1)[1];
				}
				return (onFringeOfMap?0:sMatrix[i+1][j-1])+dummyScore+addLoopOpeningPenalty; //Add dummy deltaG for starting helix
			}
			//Continuing old helix
			else {
				//get NN score.
				double nn = eParams.getNNdeltaG(base(seq,i+1, domain), base(seq2,j-1,domain), base(seq,i,domain), base(seq2, j,domain));
				double helixScore = sMatrix[i+1][j-1];
				if (sdMatrix(sdMatrix,i+1,j-1)[0]==1){
					//Remove dummy score
					helixScore -= dummyScore;
				}
				//Add nearest neighbor delta G
				helixScore += nn;
				return helixScore;
			}
		} else {
			//No helix. Extend both left and right bulges by one.
			if (writeToSD){
				sdMatrix(sdMatrix,i,j)[0] = Math.min(0,onFringeOfMap?0:sdMatrix(sdMatrix,i+1,j-1)[0])-1; //Negative means longer loop
				sdMatrix(sdMatrix,i,j)[1] = Math.min(0,onFringeOfMap?0:sdMatrix(sdMatrix,i+1,j-1)[1])-1;
			}
			if (onFringeOfMap){
				return 0.0;
			}
			//Ending old helix?
			if (sdMatrix(sdMatrix,i+1,j-1)[0]>0){
				if (sdMatrix(sdMatrix,i+1,j-1)[0]==1){
					//Remove dummy score
					return sMatrix[i+1][j-1]-dummyScore;
				} else {
					//Add terminal score.
					double terminalMismatch = eParams.getNNdeltaGterm(base(seq,i+1, domain), base(seq2,j-1,domain), base(seq,i,domain), base(seq2, j,domain));
					return sMatrix[i+1][j-1]+terminalMismatch;
				}
			} 
			//Continuing loop region
			else {
				return sMatrix[i+1][j-1];
			}
		}
	}
	private int[] sdMatrix(int[][][] sdMatrix, int i, int j) {
		return sdMatrix[i%2][j];
	}

////////////////////////////////
//// Pair Probability functions - Warning, not maintained
////////////////////////////////
	
	public void pairPr(double[][] pairsOut, DomainSequence seq1, DomainSequence seq2, int[][] domain) {
		throw new RuntimeException("Not available with this folding tool.");
	}

	public void pairPr(double[][] pairsOut, DomainSequence seq, int[][] domain) {
		throw new RuntimeException("Not available with this folding tool.");
	}
	
}
