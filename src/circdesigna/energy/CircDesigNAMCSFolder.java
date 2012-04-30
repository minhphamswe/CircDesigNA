package circdesigna.energy;

import java.awt.Point;
import java.util.ArrayList;

import circdesigna.DomainSequence;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

public class CircDesigNAMCSFolder extends CircDesigNASystemElement implements OneMatrixNAFolder {
	public CircDesigNAMCSFolder(CircDesigNAConfig System) {
		super(System);
		eParams = new ExperimentalDuplexParams(System);
	}
	private ExperimentalDuplexParams eParams;
	
	public double mfe(DomainSequence seq1, DomainSequence seq2, int[][] domain, int[][] domain_markings) {
		BiFoldAlgorithmConfig config = new BiFoldAlgorithmConfig();
		config.domain = domain;
		config.domain_markings = domain_markings;
		return biFoldAlgorithm(seq1, seq2, domain, domain_markings, config);
	}
	public double mfe(DomainSequence domainSequence, int[][] domain, int[][] domain_markings) {
		BiFoldAlgorithmConfig config = new BiFoldAlgorithmConfig();
		config.foldSingleStranded = true;
		config.domain = domain;
		config.domain_markings = domain_markings;
		return biFoldAlgorithm(domainSequence, domainSequence, domain, domain_markings, config);
	}
	public double mfeNoDiag(DomainSequence domainSequence, DomainSequence domainSequence2, int[][] domain, int[][] domain_markings) {
		double score = mfeNoDiag_NoBaseline(domainSequence, domainSequence2, domain, domain_markings);
		
		int len1 = domainSequence.length(domain);
		//SEE data / CircDesigNAMCS for the source R file that runs this regression! Must be updated
		//whenever interaction - related behavior changes.
		if (Std.isDNAMode()){
			//DNA parameters, cut off at 0       
			return Math.min(score - (15.8896 + (-0.6739)*len1),0);
		} else {
			//RNA parameters, cut off at 0         
			return Math.min(score - (22.366 + (-1.134)*len1),0);
		}
	}
	public double mfeNoDiag_NoBaseline(DomainSequence domainSequence, DomainSequence domainSequence2, int[][] domain, int[][] domain_markings) {
		BiFoldAlgorithmConfig config = new BiFoldAlgorithmConfig();
		config.ignoreDiagonal = true;
		config.domain = domain;
		config.domain_markings = domain_markings;
		double score = biFoldAlgorithm(domainSequence, domainSequence2, domain, domain_markings, config);
		
		return score;
	}
	//Minimum hairpin size of 3, so distance from diagonal is 4
	private static final int minHairpinSize = 1+3;
	/**
	 * Negative values in traceback means a (-k + 1) bifurcation
	 */
	private static final int STRAIGHT = 1, LEFT_BULGE = STRAIGHT + 1, RIGHT_BULGE = LEFT_BULGE + 1, STRAIGHT_ONEBASE = RIGHT_BULGE + 1;
	private double biFoldAlgorithm(DomainSequence ds1, DomainSequence ds2, int[][] domain, int[][] domain_markings, BiFoldAlgorithmConfig options) {
		int len1 = ds1.length(options.domain);
		int len2 = ds2.length(options.domain);
		if (options.ignoreDiagonal || options.foldSingleStranded){
			if (len1!=len2){
				throw new RuntimeException("Options "+options+" are only valid when folding strands of equal length");
			}
		}
		ensureMatrixes(len1, len2);

		int best = 0, bestI = -1, bestJ = -1;
		
		for(int i = len1-1; i>=0; i--){
			int j = options.getFirstJOnRow(i, len2);
			
			for(; j < len2; j++){
				int score = 0;
				int tb = 0;
				
				int pairedScore = pairedScore(len1, len2, ds1, ds2, i, j, options);
				if (pairedScore < 0){
					if (TB_shared[i+1][j-1]==STRAIGHT || TB_shared[i+1][j-1]==STRAIGHT_ONEBASE){
						tb = STRAIGHT;
						score = W_shared[i+1][j-1] + pairedScore;
					} else {
						tb = STRAIGHT_ONEBASE;
						score = options.getMatrixScore(i+1,j-1,len1,len2) + pairedScore;
					}
				} else {
					score = Integer.MAX_VALUE;
				}
				
				int unpaired;
				if (options.foldSingleStranded){
					bigloop: for(int index = 0; ; index++){
						int k;
						if (options.shortcutBifurcations){
							switch(index){
							case 0: k = i; break; 
							case 1: k = j-1; break;
							default: k = (j+i)/2; break;
							case 3: break bigloop;
							}
						} else {
							k = index + i;
							if (k >= j){
								break bigloop;
							}
						}
						
						int bifurc = options.getMatrixScore(i, k, len1, len2) + options.getMatrixScore(k+1, j, len1, len2);
						
						if (bifurc < score){
							tb = -(k+1);
							score = bifurc;
						}
					}
				} else {
					int leftBulge = options.getMatrixScore(i, j-1, len1, len2);
					int rightBulge = options.getMatrixScore(i+1, j, len1, len2);
					
					if (leftBulge < score){
						tb = LEFT_BULGE;
						score = leftBulge;
					}
					if (rightBulge < score){
						tb = RIGHT_BULGE;
						score = rightBulge;
					}
				}
				
				W_shared[i][j] = score;
				TB_shared[i][j] = tb;
				
				if (score < best){
					best = score;
					bestI = i;
					bestJ = j;
				}
			}
		}

		if (bestI != -1){
			traceback(bestI, bestJ, len1, len2, ds1, ds2, options);
		}
		
		//Convert to kcal / mol for return value
		return best / 100.0;
	}

	private ArrayList<Point> MFE_pointlist = new ArrayList();
	private void traceback(int bestI, int bestJ, int len1, int len2, DomainSequence ds1, DomainSequence ds2, BiFoldAlgorithmConfig options) {
		MFE_pointlist.clear();
		final int[][] domain = options.domain;
		final int[][] domain_markings = options.domain_markings;
		
		ArrayList<Point> stack = new ArrayList();
		stack.add(new Point(bestI, bestJ));
		while(!stack.isEmpty()){
			MFE_pointlist.add(null);
			Point got = stack.remove(stack.size()-1);
			if (Std.monomer.bindScore(base(ds1,got.x,domain), base(ds2,got.y,domain)) < 0){
				ds1.mark(got.x, domain, domain_markings);
				ds2.mark(got.y, domain,domain_markings);
			}
			while(true){
				MFE_pointlist.add(got);
				final int i = got.x;
				final int j = got.y;
				if (i == len1-1 || j == 0){
					break;
				}
				int move = TB_shared[i][j];
				
				if (move==0){
					break;
				}

				if (move == LEFT_BULGE){
					got = new Point(i,j-1);
				}
				if (move == RIGHT_BULGE){
					got = new Point(i+1,j);
				}
				if (move == STRAIGHT || move == STRAIGHT_ONEBASE){
					if (Std.monomer.bindScore(base(ds1,i+1,domain), base(ds2,j-1,domain)) < 0){
						ds1.mark(i+1, domain, domain_markings);
						ds2.mark(j-1, domain,domain_markings);
					}
					got = new Point(i+1,j-1);
				}
				if (move < 0){
					int k = -(move+1);
					//i,k and k+1, j
					Point leftBifurc = new Point(i,k);
					Point rightBifurc = new Point(k+1,j);

					MFE_pointlist.add(null);
					got = leftBifurc;
					stack.add(rightBifurc);
				}
			} //end while true
		}
	}
	private int min(int a, int b){
		return a < b ? a : b;
	}
	public double mfeStraight(DomainSequence ds1, DomainSequence ds2, int[][] domain, int[][] domain_markings, int markLeft, int markRight, int jOffset) {
		//if (ds1.length(domain)!=ds2.length(domain)){
		//	throw new RuntimeException("mfeStraight invalid arguments - must be two domain sequences of equal length.");
		//}
		int len1 = ds1.length(domain);
		int len2 = ds2.length(domain);
		
		BiFoldAlgorithmConfig config = new BiFoldAlgorithmConfig();
		config.domain = domain;
		config.domain_markings = domain_markings;
		
		//Need only to keep one value - we are traversing up the bottom left to upper right diagonal.
		int score = 0;
		for(int i = len1-1, j = -jOffset; i >= 0 && j < len2; i--, j++){
			if (j < 0)
				continue;
			score = score + pairedScore(len1, len2, ds1, ds2, i, j, config);
			
			if (Std.monomer.bindScore(base(ds1,i,domain), base(ds2,j,domain)) < 0){
				if (i >= markLeft && i < markRight){
					ds1.mark(i, domain, domain_markings);
					ds2.mark(j, domain, domain_markings);
				}
			}		
		}
		
		//Convert to kcal / mol for return value
		return score / 100.0;
	}
	
	private int[][] W_shared;
	private int[][] TB_shared;
	private void flushMatrixes(int i, int j) {
		W_shared[i][j] = 0;
		TB_shared[i][j] = 0;
	}
	private class BiFoldAlgorithmConfig {
		public boolean ignoreDiagonal = false;
		public boolean foldSingleStranded = false;
		public boolean shortcutBifurcations = true;
		public int[][] domain;
		public int[][] domain_markings;
		public String toString(){
			return ignoreDiagonal ? "IGNORE_DIAGONALS" : (foldSingleStranded ? "SINGLE_STRANDED" : "???");
		}
		public int getFirstJOnRow(int i, int len2) {
			int j;
			if (foldSingleStranded){
				j = i+minHairpinSize;
				for(int o = i; o < j && o < len2; o++){
					flushMatrixes(i,o);
				}				
			} else {
				j = 0;
			}
			return j;
		}
		public int getMatrixScore(int i, int j, int len1, int len2) {
			if (i >= len1 || j < 0)
				return 0;
			if (foldSingleStranded){
				if (j < i)
					return 0;
			}
			if (TB_shared[i][j]==STRAIGHT_ONEBASE || TB_shared[i][j]==STRAIGHT){
				return W_shared[i+1][j-1];
			}
			return W_shared[i][j];
		}
	}
	private void ensureMatrixes(int h, int w){
		if (h==0 || w==0){
			throw new RuntimeException("EnsureMatrixes called with 0 length or width");
		}
		if (W_shared==null || W_shared.length < h || W_shared[0].length < w){
			W_shared = null;
			TB_shared = null;
			W_shared = new int[h][w];
			TB_shared = new int[h][w];
		}
	}
	
	private int pairedScore(int len1, int len2, DomainSequence ds1, DomainSequence ds2, int i, int j, BiFoldAlgorithmConfig options) {
		int[][] domain = options.domain;
		int[][] domain_markings = options.domain_markings;
		
		if (options.ignoreDiagonal){
			if (i==len2-1-j){
				return 0;
			}
		}
		
		final int bi = base(ds1,i,domain);
		final int bj = base(ds2,j,domain);

		//Do we have a prior base?
		final int ni = i+1;
		final int nj = j-1;
		//Is increase 5' - 3'
		//Js decrease 3' - 5' (i.e. i/j point in the same direction)
		if (ni >= len1 || nj < 0){
			return 0;
		}
		//Add a stack / mismatch score.
		final int bni = base(ds1,ni,domain);
		final int bnj = base(ds2,nj,domain);
		
		if (Std.monomer.bindScore(bi, bj) < 0){
			if (Std.monomer.bindScore(bni, bnj) < 0){
				//Stack.
				return eParams.getNN_deci(bi, bj, bni, bnj);
			} else {
				//Mismatch
				return eParams.getInteriorNNTerminal_deci(bi, bj, bni, bnj);
			}
		} else {
			if (Std.monomer.bindScore(bni, bnj) < 0){
				//Mismatch
				return eParams.getInteriorNNTerminal_deci(bni, bnj, bi, bj);
			}
		}
		return 0;
	}
	public void pairPr(double[][] pairsOut, DomainSequence seq1, DomainSequence seq2, int[][] domain) {
		throw new RuntimeException("Not implemented");
	}
	public void pairPr(double[][] pairsOut, DomainSequence seq1, int[][] domain) {
		throw new RuntimeException("Not implemented");
	}
	public double[][] getScoreMatrix(int len1, int len2) {
		double[][] nussinovScores = new double[len1][len2];
		for(int y = 0; y < len1; y++){
			for(int x = 0; x < len2; x++){
				nussinovScores[y][x] = W_shared[y][x] / 100.0;
			}
		}
		return nussinovScores;
	}
	public ArrayList<Point> getTraceback() {
		return MFE_pointlist;
	}
}
