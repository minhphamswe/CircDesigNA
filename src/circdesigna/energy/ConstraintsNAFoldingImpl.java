package circdesigna.energy;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import circdesigna.CircDesigNA_SharedUtils;
import circdesigna.Connector;
import circdesigna.GeneralizedInteractiveRegion;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;

/**
 * Uses an N^3 DP algorithm to compute the MFE among all unpseudoknotted folded structures
 * of one or two sequences.
 */
public class ConstraintsNAFoldingImpl extends CircDesigNASystemElement implements ConstraintsNAFolding, OneMatrixNAFolding{
	private ExperimentalDuplexParams eParams;
	private int[][][] memo_shared;
	private int[] seq_shared;
	private int[][] seq_origin_shared;
	private ConnectorFold[] connectors_shared;
	private FoldingConstraints constraints_shared;
	private int scoringModel;

	public ConstraintsNAFoldingImpl(CircDesigNAConfig System) {
		super(System);
		eParams = new ExperimentalDuplexParams(System);
		setScoringModel(3);
	}
	
	private boolean lockClaimed = false;
	private void claimLock(){
		if (lockClaimed){
			throw new RuntimeException("Attempt to multithread nonmultithreadsafe method MFE");
		}
		lockClaimed = true;
	}
	private void returnLock(){
		lockClaimed = false;
	}
	private int[] getSeq(int N){
		if (seq_shared == null || seq_shared.length < N){
			seq_shared = new int[N];
		}
		return seq_shared;
	}
	private int[][] getSeqOrigin(int N){
		if (seq_origin_shared == null || seq_origin_shared[0].length < N){
			seq_origin_shared = new int[N][2];
		}
		return seq_origin_shared;
	}
	private ConnectorFold[] getConnectors(int N){
		if (connectors_shared == null || connectors_shared.length < N){
			connectors_shared = new ConnectorFold[N];
			for(int k = 0; k < N; k++){
				connectors_shared[k] = new ConnectorFold();
			}
		}
		return connectors_shared;
		
	}
	private int[][][] getMemo(int N){
		if (memo_shared == null){
			memo_shared = new int[EXTERNAL_EQ1_0S+1][][];
		}
		return memo_shared;
	}
	private FoldingConstraints getConstraints(int N){
		if (constraints_shared == null || constraints_shared.preventPairing[0].length < N){
			constraints_shared = new FoldingConstraints(N);
		}
		return constraints_shared;
	}
	private class ConnectorFold {
		public void init(GeneralizedInteractiveRegion gir, int i, int[][] domain) {
			leftN = 0;
			rightN = 0;
			//Scoring a set of connectors involves entirely coaxial stacking
			//which in turn is entirely made up of dangle penalties
			{
				List<Connector> connectors = gir.getConnectorsTo5Of(i, domain);
				if (connectors!=null && !connectors.isEmpty()){
					int N = connectors.size();
					leftN = N;
					leftPairBase0 = base(connectors.get(N-1), 0, domain);
					leftPairBase1 = base(connectors.get(N-1), 1, domain);
					leftScore00 = scoreConnectors(connectors,0,N-1,domain);
					leftScore01 = scoreConnectors(connectors,0,N-2,domain);
					leftScore10 = scoreConnectors(connectors,1,N-1,domain);
					leftScore11 = scoreConnectors(connectors,1,N-2,domain);
				}
			}
			{
				List<Connector> connectors = gir.getConnectorsTo3Of(i, domain);
				if (connectors!=null && !connectors.isEmpty()){
					int N = connectors.size();
					rightN = N;
					rightPairBase0 = base(connectors.get(0), 0, domain);
					rightPairBase1 = base(connectors.get(0), 1, domain);
					rightScore00 = scoreConnectors(connectors,0,N-1,domain);
					rightScore01 = scoreConnectors(connectors,0,N-2,domain);
					rightScore10 = scoreConnectors(connectors,1,N-1,domain);
					rightScore11 = scoreConnectors(connectors,1,N-2,domain);
				}
			}
		}
		/**
		 * Dynamic programming folding algorithm which finds the minimum free energy of coaxial stacking of the sequence of connectors
		 * i through j.
		 */
		private int scoreConnectors(List<Connector> connectors, int i, int j, int[][] domain) {
			int bestKmin2 = 0;
			int bestKmin1 = 0;
			for(int k = 1; k <= j; k++){
				//Either connector k is coaxially stacked to k-1,
				//or it isn't
				int bestK = alt(combine(bestKmin2,
							getCoaxialStackBonus(
												base(connectors.get(k-1),1,domain),
												base(connectors.get(k-1),0,domain),
												base(connectors.get(k),0,domain),
												base(connectors.get(k),1,domain))
						),
						bestKmin1);
				//Move ahead
				bestKmin2 = bestKmin1;
				bestKmin1 = bestK;
			}
			return bestKmin1;
		}
		public int leftN, rightN;
		public int leftPairBase0, leftPairBase1, rightPairBase0, rightPairBase1;
		//Free energy of the connectors on the left, due to coaxial stacking. 10 means that the leftmost pair in the left connectors
		//is engaged with coaxial stacking with a pair external to the set of connectors.
		//01 means that the rightmost pair in the left connectors is engaged with coaxial stacking with a pair external to the set of connectors.
		public int leftScore00, leftScore01, leftScore10, leftScore11;
		//See description of leftScore, but apply to the set of connectors on the right of this base.
		public int rightScore00, rightScore01, rightScore10, rightScore11;
	}
	public class ConnectorSummary {
		/**
		 * Extends the interval of bases this summary represents one to the right
		 * Note that the bases i and j are excluded from the interval itself (in terms of "unpaired bases", etc)
		 */
		public void extendInterval(int basei, ConnectorFold i, int baseip1, ConnectorFold ip1, int baseip2, ConnectorFold ip2){
			if (i.rightN != ip1.leftN || ip1.rightN != ip2.leftN){
				throw new RuntimeException("Assertion error");
			}
			unpairedBases++;
			
			//Add loops between i and ip1
			if (i.rightN > 0){
				if (numConnectors == 0){
					leftPairBase0 = i.rightPairBase0;
					leftPairBase1 = i.rightPairBase1;
				}
				rightPairBase0 = ip1.leftPairBase0;
				rightPairBase1 = ip1.leftPairBase1;
				numConnectors += i.rightN;
				//add dangle on left of i+1
				//interiorDangles = combine(interiorDangles, eParams.getDangleTop_deci(ip1.leftPairBase1, ip1.leftPairBase0, baseip1));
			}
			if (i.rightN == 0 && numConnectors == 0){
				unpairedBasesLeft++;
			}
			
			//Add loops between ip1 and ip2
			if (ip1.rightN > 0){
				if (numConnectors == 0){
					leftPairBase0 = ip1.rightPairBase0;
					leftPairBase1 = ip1.rightPairBase1;
				}
				rightPairBase0 = ip2.leftPairBase0;
				rightPairBase1 = ip2.leftPairBase1;
				numConnectors += ip1.rightN;
				//add dangle on right of i+1
				//interiorDangles = combine(interiorDangles, eParams.getDangleBottom_deci(ip1.rightPairBase1, ip1.rightPairBase0, baseip1));
			}
			if (ip1.rightN == 0 && numConnectors > 0){
				unpairedBasesRight++;
			}
		}
		public int unpairedBases = 0; 
		public int unpairedBasesLeft = 0, unpairedBasesRight = 0;
		public int leftPairBase0, leftPairBase1, rightPairBase0, rightPairBase1;
		public int numConnectors = 0;
	}
	
	public void setScoringModel(int i) {
		scoringModel = i;
	}
	

	public double mfe(GeneralizedInteractiveRegion seq1, GeneralizedInteractiveRegion seq2, int[][] domain, int[][] domain_markings) {
		return mfe(seq1, seq2, domain, domain_markings, false);
	}
	/**
	 * Computes the minimum free energy of the minimum free energy structure of all structures which are
	 * nonpseudoknotted and which do not make any allowed pairs (so complementary domains will not pair.)
	 * 
	 * If order is 2, the optimization searches only over all structures which have the nested loops property.
	 * (So, if a is paired to b, c is paired to d, either a < c < d < b or c < a < b < d.)
	 */
	public double mfe(GeneralizedInteractiveRegion seq1, GeneralizedInteractiveRegion seq2, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing) {		
		if (!seq1.isCircular() && seq2.isCircular()){
			//swap the sequences so that 1 is circular and 2 is not.
			GeneralizedInteractiveRegion tmp = seq1;
			seq1 = seq2;
			seq2 = tmp;
		}
		
		int N1 = seq1.length(domain);
		int N2 = seq2.length(domain);
		int N = N1 + N2;
		
		int[] nicks;
		if (seq1.isCircular() && seq2.isCircular()){
			//Such pairing is practically impossible.
			return 0;
		} else {
			if (seq1.isCircular()){
				//1 is circular, 2 is not. TODO
				//throw new RuntimeException();
				//Treat as if neither are circular.
				nicks = new int[]{N1-1,N-1};
			} else {
				//Neither are circular
				nicks = new int[]{N1-1,N-1};
			}
		}

		claimLock();
		try {
			int[] seq = getSeq(N);
			int[][] seq_origin = getSeqOrigin(N);
			ConnectorFold[] connectors = getConnectors(N);
			int[][][] memo2 = getMemo(N);
			FoldingConstraints constraints = new FoldingConstraints(0);
			SequenceMarker marker = new SequenceMarker(N, seq_origin, domain_markings);
			
			for(int k = 0; k < N1; k++){
				seq[k] = seq1.base(k, domain, Std.monomer);
				seq_origin[k][0] = seq1.domainAt(k, domain);
				seq_origin[k][1] = seq1.offsetInto(k, domain, true);
				connectors[k].init(seq1, k, domain);
			}
			for(int k = 0; k < N2; k++){
				seq[N1+k] = seq2.base(k, domain, Std.monomer);
				seq_origin[N1+k][0] = seq2.domainAt(k, domain);
				seq_origin[N1+k][1] = seq2.offsetInto(k, domain, true);
				connectors[N1+k].init(seq2, k, domain);
			}

			if (onlyIllegalPairing){
				if (CircDesigNA_SharedUtils.checkComplementary(seq1, seq1) || 
						CircDesigNA_SharedUtils.checkComplementary(seq1, seq2) || 
						CircDesigNA_SharedUtils.checkComplementary(seq2, seq2)){
					constraints = getConstraints(N);

					for(int i = 0; i < N; i++){
						for(int j = 0; j < N; j++){
							constraints.preventPairing[i][j] = (seq_origin[i][1] == seq_origin[j][1]) && 
									CircDesigNA_SharedUtils.checkComplementary(seq_origin[i][0], seq_origin[j][0]);
						}
					}
				}
			}  

			NXFold(memo2, seq, N, nicks, connectors, true, constraints, marker);
			//.51 is bimolecular penalty (TODO for DNA? RNA?)
			int nStrands = 2;
			double toRet = getQe(memo2, 0, N-1, 0, 0) / 100.0 - .51 * (nStrands - 1);

			return toRet;
		} finally {
			returnLock();
		}
	}
	public double mfe(GeneralizedInteractiveRegion GeneralizedInteractiveRegion, int[][] domain, int[][] domain_markings) {
		return mfe(GeneralizedInteractiveRegion, domain, domain_markings, false);
	}

	public double mfe(GeneralizedInteractiveRegion gir, int[][] domain, int[][] domain_markings, boolean onlyIllegalPairing) {
		claimLock();
		try {

			int N = gir.length(domain);
			int[] nicks;
			if (gir.isCircular()){
				nicks = new int[]{};
			} else {
				nicks = new int[]{N-1};
			}	


			int[] seq = getSeq(N);
			int[][] seq_origin = getSeqOrigin(N);
			ConnectorFold[] connectors = getConnectors(N);
			int[][][] memo2 = getMemo(N);
			FoldingConstraints constraints = new FoldingConstraints(0);
			SequenceMarker marker = new SequenceMarker(N, seq_origin, domain_markings);
			for(int k = 0; k < N; k++){
				seq[k] = gir.base(k, domain, Std.monomer);
				seq_origin[k][0] = gir.domainAt(k, domain);
				seq_origin[k][1] = gir.offsetInto(k, domain, true);
				connectors[k].init(gir, k, domain);
			}

			if (onlyIllegalPairing){
				if (CircDesigNA_SharedUtils.checkComplementary(gir, gir)){
					constraints = getConstraints(N);
					for(int i = 0; i < N; i++){
						for(int j = 0; j < N; j++){
							constraints.preventPairing[i][j] = (seq_origin[i][1] == seq_origin[j][1]) && 
									CircDesigNA_SharedUtils.checkComplementary(seq_origin[i][0], seq_origin[j][0]);
						}
					}
				}
			}  

			double toRet;

			if (gir.isCircular()){
				NXFold(memo2, seq, N, nicks, connectors, false, constraints, marker);
				int[][] Qb = memo2[PAIRED];
				//If the circular DNA has no pairs
				int best = 0;
				//If the circular DNA has one pair 
				for(int i = 0; i <= N-2; i++){
					for(int j = i+1; j <= N-1; j++){
						best = alt(best, combine(Qb[i][j],Qb[j][i]));
					}
				}
				toRet = best / 100.0;
			} else {
				NXFold(memo2, seq, N, nicks, connectors, true, constraints, marker);
				toRet = getQe(memo2, 0, N-1, 0, 0) / 100.0;
			}
			return toRet;
			
		} finally {
			returnLock();
		}
	}
	
	private void NXFold(int[][][] memo2, int[] seq, int N, int[] nicks, ConnectorFold[] connectors, boolean onlyUpperTriangle, FoldingConstraints constraints, SequenceMarker marker) {
		switch(scoringModel){
		//3: Finds the MFE over all unpseudoknotted structures, implemented in O(n^3) time.
		case 3:
			N3Fold(memo2, seq, N, nicks, onlyUpperTriangle, constraints, marker);
			return;
		//2: Finds the MFE over all unpseudoknotted structures with no multiloops, scoring interior loops with a linear overapproximation
		case 2:
			N2Fold(memo2, seq, N, nicks, onlyUpperTriangle, constraints, marker);
			/* Refinement is incorrect, does not take into account coaxial stacking on external helixes
			//Refinement: Improve the values of Qe[0][x] for x going from 1 to N-1.
			int[][] Qe00 = memo2[EXTERNAL_00];
			int[][] Qe01 = memo2[EXTERNAL_01];
			int[][] Qe10 = memo2[EXTERNAL_10];
			int[][] Qe11 = memo2[EXTERNAL_11];
			for(int x = 2; x < N; x++){
				for(int m = 1; m < x; m++){
					Qe00[0][x] = alt(alt(Qe00[0][x], combine(Qe00[0][m], Qe00[m][x])), combine(Qe01[0][m], Qe10[(m+1)%N][x]));
					Qe01[0][x] = alt(alt(Qe01[0][x], combine(Qe00[0][m], Qe01[m][x])), combine(Qe01[0][m], Qe11[(m+1)%N][x]));
					Qe10[0][x] = alt(alt(Qe10[0][x], combine(Qe10[0][m], Qe00[m][x])), combine(Qe11[0][m], Qe10[(m+1)%N][x]));
					Qe11[0][x] = alt(alt(Qe11[0][x], combine(Qe10[0][m], Qe01[m][x])), combine(Qe11[0][m], Qe11[(m+1)%N][x]));
				}
			}
			*/
			return;
		//1: Finds the MFE over all unpseudoknotted structures with no interior loops or bulges
		case 1:
			N2Fold_NoLoops(memo2, seq, N, nicks, connectors, onlyUpperTriangle, constraints, marker);
			return;
		}
		throw new RuntimeException("Not a valid scoring model: "+scoringModel);
	}

	//Indeces in memo
	private static final int 
	EXTERNAL_00 = 0,
	EXTERNAL_01 = EXTERNAL_00+1,
	EXTERNAL_10 = EXTERNAL_01+1,
	EXTERNAL_11 = EXTERNAL_10+1,
	PAIRED = EXTERNAL_11 + 1,
	LOOP_EXTENSION_J = PAIRED+1,
	LOOP_EXTENSION_JMIN1 = LOOP_EXTENSION_J+1,
	LOOP_EXTENSION_JMIN2 = LOOP_EXTENSION_JMIN1+1,
	RIGHT_BULGE_EXTENSION_J = LOOP_EXTENSION_JMIN2+1,
	RIGHT_BULGE_EXTENSION_JMIN1 = RIGHT_BULGE_EXTENSION_J+1,
	LEFT_BULGE_EXTENSION_J = RIGHT_BULGE_EXTENSION_JMIN1+1,
	LEFT_BULGE_EXTENSION_JMIN1 = LEFT_BULGE_EXTENSION_J+1,
	MULTIBRANCH_EQ1_00 = LEFT_BULGE_EXTENSION_JMIN1+1,
	MULTIBRANCH_EQ1_01 = MULTIBRANCH_EQ1_00+1,
	MULTIBRANCH_EQ1_10 = MULTIBRANCH_EQ1_01+1,
	MULTIBRANCH_EQ1_0S = MULTIBRANCH_EQ1_10+1,
	MULTIBRANCH_EQ1_S0 = MULTIBRANCH_EQ1_0S+1,
	MULTIBRANCH_EQ2COAX_01 = MULTIBRANCH_EQ1_S0+1,
	MULTIBRANCH_EQ2COAX_11 = MULTIBRANCH_EQ2COAX_01+1,
	MULTIBRANCH_EQ2COAX_S1 = MULTIBRANCH_EQ2COAX_11+1,
	MULTIBRANCH_EQ2COAX_0S = MULTIBRANCH_EQ2COAX_S1+1,
	MULTIBRANCH_EQ2COAX_1S = MULTIBRANCH_EQ2COAX_0S+1,
	MULTIBRANCH_EQ2COAX_SS = MULTIBRANCH_EQ2COAX_1S+1,
	MULTIBRANCH_GE2_00 = MULTIBRANCH_EQ2COAX_SS+1,
	MULTIBRANCH_GE2_01 = MULTIBRANCH_GE2_00+1,
	MULTIBRANCH_GE2_10 = MULTIBRANCH_GE2_01+1,
	MULTIBRANCH_GE2_11 = MULTIBRANCH_GE2_10+1,
	EXTERNAL_EQ2COAX_01 = MULTIBRANCH_GE2_11+1,
	EXTERNAL_EQ2COAX_11 = EXTERNAL_EQ2COAX_01+1,
	EXTERNAL_EQ2COAX_S1 = EXTERNAL_EQ2COAX_11+1,
	EXTERNAL_EQ2COAX_0S = EXTERNAL_EQ2COAX_S1+1,
	EXTERNAL_EQ2COAX_1S = EXTERNAL_EQ2COAX_0S+1,
	EXTERNAL_EQ2COAX_SS = EXTERNAL_EQ2COAX_1S+1,
	EXTERNAL_EQ1_01 = EXTERNAL_EQ2COAX_SS+1,
	EXTERNAL_EQ1_0S = EXTERNAL_EQ1_01+1
	;

	/**
	 * Performs the recursion to compute the MFE energy over all foldings of the given strands which are
	 * connected foldings and which obey the folding constraints.
	 * 
	 * If onlyUpperTriangle is true, then only recursions on subsequences [i,j], where i <= j, are computed.
	 * This is essentially always on, with the exception of circular folding.
	 * @param domain_markings 
	 * @param domain 
	 * @param ds 
	 */
	private void N3Fold(int[][][] memo2, int[] seq, int N, int[] nicks, boolean onlyUpperTriangle, FoldingConstraints constraints, SequenceMarker marker) {
		if (true) throw new RuntimeException("Phase 3 not yet written to fold GIR's");
		
		//Initialization: all is impossible.
		deepFill3(memo2, N, N, Integer.MAX_VALUE);

		//Initialization: A single, external base has 0 free energy.
		int[][] Qe00 = memo2[EXTERNAL_00];
		int[][] Qe01 = memo2[EXTERNAL_01];
		int[][] Qe10 = memo2[EXTERNAL_10];
		int[][] Qe11 = memo2[EXTERNAL_11];
		for(int i = 0; i <= N-1; i++){
			Qe00[i][i] = 0;
		}
		
		int[][] Qb = memo2[PAIRED];
		int[][] QbNoStk = Qb; 
		int[][] QbStk = Qb; 

		int a1 = eParams.getMultibranchBase_deci();
		int a2 = eParams.getMultibranchBranch_deci();
		int a3 = eParams.getMultibranchUnpairedBase_deci();
		
		//Multibranch terms
		int[][] QmEq1_00 = memo2[MULTIBRANCH_EQ1_00];
		int[][] QmEq1_01 = memo2[MULTIBRANCH_EQ1_01];
		int[][] QmEq1_10 = memo2[MULTIBRANCH_EQ1_10];
		int[][] QmEq1_0S = memo2[MULTIBRANCH_EQ1_0S];
		int[][] QmEq1_S0 = memo2[MULTIBRANCH_EQ1_S0];
		
		int[][] QmEq2Coax_01 = memo2[MULTIBRANCH_EQ2COAX_01];
		int[][] QmEq2Coax_11 = memo2[MULTIBRANCH_EQ2COAX_11];
		int[][] QmEq2Coax_S1 = memo2[MULTIBRANCH_EQ2COAX_S1];
		int[][] QmEq2Coax_0S = memo2[MULTIBRANCH_EQ2COAX_0S];
		int[][] QmEq2Coax_1S = memo2[MULTIBRANCH_EQ2COAX_1S];
		int[][] QmEq2Coax_SS = memo2[MULTIBRANCH_EQ2COAX_SS];
		
		int[][] QmGe2_00 = memo2[MULTIBRANCH_GE2_00];
		int[][] QmGe2_01 = memo2[MULTIBRANCH_GE2_01];
		int[][] QmGe2_10 = memo2[MULTIBRANCH_GE2_10];
		int[][] QmGe2_11 = memo2[MULTIBRANCH_GE2_11];

		//External coaxial stacking
		int[][] QeEq2Coax_01 = memo2[EXTERNAL_EQ2COAX_01];
		int[][] QeEq2Coax_11 = memo2[EXTERNAL_EQ2COAX_11];
		int[][] QeEq2Coax_S1 = memo2[EXTERNAL_EQ2COAX_S1];
		int[][] QeEq2Coax_0S = memo2[EXTERNAL_EQ2COAX_0S];
		int[][] QeEq2Coax_1S = memo2[EXTERNAL_EQ2COAX_1S];
		int[][] QeEq2Coax_SS = memo2[EXTERNAL_EQ2COAX_SS];
	
		int[][] QeEq1_01 = memo2[EXTERNAL_EQ1_01];
		int[][] QeEq1_0S = memo2[EXTERNAL_EQ1_0S];

		//Build up to larger subproblems:
		for( int L = 2; L <= N; L++ ){
			//Swap the Qx buffers
			swap(memo2, LOOP_EXTENSION_JMIN2, LOOP_EXTENSION_JMIN1);
			swap(memo2, LOOP_EXTENSION_JMIN1, LOOP_EXTENSION_J);
			
			int[][] Qxj = memo2[LOOP_EXTENSION_J];
			int[][] Qxj1 = memo2[LOOP_EXTENSION_JMIN1];
			int[][] Qxj2 = memo2[LOOP_EXTENSION_JMIN2];
			for(int[] row : Qxj){
				Arrays.fill(row, Integer.MAX_VALUE);
			}
			
			for( int i = 0; i <= N-1; i++){
				int j = (i + L - 1) % N;
				if (onlyUpperTriangle && j < i){
					break;
				}

				//Structure of the recursion:
				//001) i and j are paired with eachother
					//001) They form an interior structure
					//001.5) They form an exterior structure
				//002) [i,j] is inside a multiloop (or exterior, as the calculation is highly similar)
					//002) [i,j] contains exactly one pair
					//002.25) [i,j] contains 2 pairs exactly, coaxially stacked, and j is paired or j unpaired and j-1 paired
					//002.5) External loop general case (>= 0 pairs)
					//002.75) Multibranch general case (>=2 pairs)
								
				//001)
				//if i and j are paired to one another, and they close an interior region
				if (L >= 4 && !(containsNick(i, (i+1)%N, nicks) && containsNick((j-1+N)%N,j,nicks))){
					//if the region contains no pairs, i.e. hairpin
					if (!containsNick(i, j, nicks)){
						QbNoStk[i][j] = alt(QbNoStk[i][j],eParams.getHairpinLoopDeltaG_deci(seq, N, i, L));
					}
					
					//if the region contains 1 pair
					//Step 1: Find the optimal interior loop and store it in Qxj
					//If L1 > 4 and L2 > 4, there is no sequence dependence and hence we use memo 
					for(int s = 10; s <= L - 4; s++){
						//symmetry term is unmodified. However, the length of the interior loop is increased.
						Qxj[i][s] = combine(Qxj2[(i+1)%N][s-2], eParams.getInteriorLoopSizeTerm_deci(s) - eParams.getInteriorLoopSizeTerm_deci(s-2));
					}
					//If L2 >= 4 and L1 = 4, 
					for(int L2 = 4, L1 = 4; L1 + L2 <= L - 4; L2++){
						int d = (i + L1 + 1)%N;
						int e = (j - L2 - 1 + N)%N;
						int s = L1 + L2;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int closingPair = eParams.getInteriorNNTerminal_deci(seq[e], seq[d], seq[(e+1)%N], seq[(d-1+N)%N]);
						Qxj[i][s] = alt(Qxj[i][s], 
								combine(Qb[d][e], 
										combine(eParams.getInteriorLoopSizeTerm_deci(s) + eParams.getNINIOAssymetry(L1, L2),
												closingPair
												)));
					}
					//If L1 > 4 and L2 = 4, 
					for(int L1 = 5, L2 = 4; L1 + L2 <= L - 4; L1++){
						int d = (i + L1 + 1)%N;
						int e = (j - L2 - 1 + N)%N;
						int s = L1 + L2;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int closingPair = eParams.getInteriorNNTerminal_deci(seq[e], seq[d], seq[(e+1)%N], seq[(d-1+N)%N]);
						Qxj[i][s] = alt(Qxj[i][s], 
								combine(Qb[d][e], 
										combine(eParams.getInteriorLoopSizeTerm_deci(s) + eParams.getNINIOAssymetry(L1, L2),
												closingPair
												)));
					}
					//Step 2: Translate the value in Qxj as a possibility for Qb
					int internalLoopLeftClosing = eParams.getInteriorNNTerminal_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N]);
					for(int s = 2; s <= L - 4; s++){
						QbNoStk[i][j] = alt(QbNoStk[i][j], combine(internalLoopLeftClosing, Qxj[i][s]));
					}
					
					//If L1 < 4.
					for(int L1 = 0; L1 < 4; L1++){
						for(int L2 = 0; L1 + L2 <= L - 4; L2++){
							int d = (i+L1+1)%N;
							int e = (j-L2-1+N)%N;
							if (containsNick(i, d, nicks) || containsNick(e,j, nicks)){
								break; //Only inner loop!
							}
							
							int stack = eParams.getNN_deci(seq[i], seq[j], seq[d], seq[e]);
							if (L1 == 0 && L2 == 0){
								//Continue a stack
								QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], stack));
								//Possibly terminate stack as well?
								QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], stack));
							} 
							//In the standard model, bulge loops have no sequence dependence on the unpaired bases. 
							if (L1 == 0 && L2 > 0){
								int bulge = eParams.getBulgeLoop_deci(L2); 
								//For L2 == 1 we consider that the stack is not terminated
								if (L2 == 1){
									bulge = combine(bulge, stack);
									QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], bulge));
									//Possibly terminate stack
									QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], bulge));
								} else {
									int ATPenalty = eParams.getATPenalty_deci(seq[d], seq[e]) + eParams.getATPenalty_deci(seq[i], seq[j]);
									QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], bulge + ATPenalty));
								}
							}
							if (L2 == 0 && L1 > 0){
								//For L1 == 1 we consider that the stack is not terminated
								int bulge = eParams.getBulgeLoop_deci(L1);
								if (L1 == 1){
									bulge = combine(bulge, stack);
									QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], bulge));
									//Possibly terminate stack
									QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], bulge));
								} else {
									int ATPenalty = eParams.getATPenalty_deci(seq[d], seq[e]) + eParams.getATPenalty_deci(seq[i], seq[j]);
									QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], bulge + ATPenalty));
								}
							}
							if (L1 > 0 && L2 > 0){
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], eParams.getInteriorLoop_deci(seq, N, i, j, L1, L2)));
							}
						}
					} //end L1 < 4
					
					//If L1 >= 4 but L2 < 4
					for(int L2 = 0; L2 < 4; L2++){
						for(int L1 = 4; L1 + L2 <= L - 4; L1++){
							int d = (i+L1+1)%N;
							int e = (j-L2-1+N)%N;
							if (containsNick(i, d, nicks) || containsNick(e,j, nicks)){
								break; //only inner loop!
							}
							if (L2 == 0 && !containsNick(i, d, nicks) && !containsNick(e,j, nicks)){
								int bulge = eParams.getBulgeLoop_deci(L1);
								int ATPenalty = eParams.getATPenalty_deci(seq[d], seq[e]) + eParams.getATPenalty_deci(seq[i], seq[j]);
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], bulge + ATPenalty));
							}
							if (L2 > 0){
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], eParams.getInteriorLoop_deci(seq, N, i, j, L1, L2)));
							}
						}
					}
				
					//Check: the above cases have checked
					//(L1 > 4, L2 > 4), (L1 == 4, L2 >= 4), (L2 == 4, L1 > 4), (L1 < 4), (L2 < 4, L1 >= 4).
					
					//if the region contains >= 2 pairs (multibranch). 
					//Remark: multi-branch scoring is the most shaky (experimentally) part of the model.
					if (L >= 6){
						//the pair i, j coaxially stacks with an internal loop?
						//Search for 5'-of-i coaxial stacks
						for(int d = (i+1)%N, y = 0; y <= 1; d=(d+1)%N, y++){
							for(int e = (d+1)%N; e!=(y-2+N)%N; e=(e+1)%N){
								int internalDangle = eParams.getDangleTop_deci(seq[e],seq[d],seq[(e+1)%N]);
								int externalDangle = eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]);
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(
										getQmGe1(memo2, seq, N, nicks, (e+1)%N,(j-1+N)%N, internalDangle, externalDangle),
										getCoaxialStackBonus(memo2, seq, N, nicks, i, j, d, e, y),
										Qb[d][e], a1+2*a2+a3*y
										));
							}
						}
						//Search for 3'-of-j coaxial stacks
						for(int e = (j-1+N)%N, y = 0; y <= 1; e=(e-1+N)%N, y++){
							for(int d = (i+3)%N; d!=e; d=(d+1)%N){
								int internalDangle = eParams.getDangleBottom_deci(seq[e], seq[d], seq[(d-1+N)%N]);
								int externalDangle = eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]);
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(
										getQmGe1(memo2, seq, N, nicks, (i+1)%N, (d-1+N)%N, externalDangle, internalDangle),
										getCoaxialStackBonus(memo2, seq, N, nicks, e, d, j, i, y),
										Qb[d][e], a1+2*a2+a3*y
										));
							}	
						}
						//no coaxial stacking of i,j.
						int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
						int leftDangle = eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]);
						int rightDangle = eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]);
						QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, a1 + a2, getQmGe2(memo2, (i+1)%N, (j-1+N)%N, leftDangle, rightDangle)));
					}
				} //end L >= 4
				
				//001.5)
				//if i and j are paired to one another, and their pair is exterior
				if (L == 2 && containsNick(i, j, nicks)){
					int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
					QbNoStk[i][j] = ATPenalty;  //Blunt end
				}
				//Read this line carefully: We rule out only the case where BOTH (i,i+1) and (j-1, j) are nicked.
				if (L >= 3 && !(containsNick(i, (i+1)%N, nicks) && containsNick((j-1+N)%N,j,nicks))){
					for(int k = 0; k < nicks.length; k++){
						int nick = nicks[k];
						if (containsNick(i,j,nick)){
							int leftloopNoCoax = Integer.MAX_VALUE;
							int leftloopCoax = Integer.MAX_VALUE;
							if (nick==i){
								//empty left loop
							} else { 
								int topDangle = eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]);
								leftloopNoCoax = getQe(memo2, (i+1)%N, nick, topDangle, 0);
								//coaxial stack of [i,j] with a stack on the continuous backbone of the i side (the left)
								for(int e = (i+2)%N; e!=(nick+1)%N; e=(e+1)%N){ 
									for(int d = (i+1)%N, y = 0; d!=e && y <= 1; d=(d+1)%N, y++){
										int remainder = 0;
										if (e!=nick){
											int midDangle = eParams.getDangleTop_deci(seq[e], seq[d], seq[(e+1)%N]);
											remainder = getQe(memo2, (e+1)%N, nick, midDangle, 0);
										}
										int pairedScore = Qb[d][e];
										leftloopCoax = alt(leftloopCoax, combine(remainder, pairedScore, getCoaxialStackBonus(memo2, seq, N, nicks, i, j, d, e, y)));
									}
								}
							}
							int rightloopNoCoax = Integer.MAX_VALUE;
							int rightloopCoax = Integer.MAX_VALUE;
							if (nick==(j-1+N)%N){
								//empty right loop
							} else {
								int bottomDangle = eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]);
								rightloopNoCoax = getQe(memo2, (nick+1)%N, (j-1+N)%N, 0, bottomDangle);
								//coaxial stack of [i,j] with a stack on the continuous backbone of the j side (the right)
								for(int d = (nick+1)%N; d!=j; d=(d+1)%N){ 
									for(int e = (j-1+N)%N, y = 0; e!=j && y <= 1; e=(e-1+N)%N, y++){
										int remainder = 0;
										if (d!=(nick+1)%N){
											int midDangle = eParams.getDangleBottom_deci(seq[e], seq[d], seq[(d-1+N)%N]);
											remainder = getQe(memo2, (nick+1)%N, (d-1+N)%N, 0, midDangle);
										}
										int pairedScore = Qb[d][e];
										rightloopCoax = alt(rightloopCoax, combine(remainder, pairedScore, getCoaxialStackBonus(memo2, seq, N, nicks, e, d, j, i, y)));
									}
								}
							}

							int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
							if (nick==i){
								//empty left loop
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, alt(rightloopCoax, rightloopNoCoax)));
							} else
							if (nick==(j-1+N)%N){
								//empty right loop
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, alt(leftloopCoax, leftloopNoCoax)));
							} else {
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, combine(leftloopNoCoax, rightloopNoCoax)));
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, combine(leftloopCoax, rightloopNoCoax)));
								QbNoStk[i][j] = alt(QbNoStk[i][j], combine(ATPenalty, combine(leftloopNoCoax, rightloopCoax)));
							}
						}
					}
				}
				
				//If i and j are paired to one another, and the pair is either a new stack, or a single pair
				//Qb[i][j] = alt(QbNoStk[i][j], QbStk[i][j]);
				
				//Overrides: Can i and j pair?
				if (constraints.preventPairing(i,j)){
					Qb[i][j] = Integer.MAX_VALUE;
				}
				
				if (Qb[i][j] < 0){
					marker.mark(i);
					marker.mark(j);
				}

				//002) If [i,j] is in a multiloop or is external
				//If [i,j] contains exactly 1 pair
				if (L >= 4){
					//0S
					for(int e = (j-1+N)%N, d = (i+1)%N, unpairedOnLeft = 1; d!=e; d=(d+1)%N, unpairedOnLeft++){
						if (containsNick(i,d,nicks)){
							break;
						}
						int pairedScoreWithDangles = getQbFull(memo2, seq, N, nicks, d, e, 1, 1);
						QmEq1_0S[i][j] = alt(QmEq1_0S[i][j], combine(pairedScoreWithDangles, a2+unpairedOnLeft*a3));
						QeEq1_0S[i][j] = alt(QeEq1_0S[i][j], pairedScoreWithDangles);
					}
					//00, S0
					for(int e = (i+2)%N, unpairedOnRight = L - 3; e!=j; e=(e+1)%N, unpairedOnRight--){
						if (containsNick(e,j,nicks)){
							continue;
						}
						QmEq1_00[i][j] = alt(QmEq1_00[i][j],combine(QmEq1_0S[i][(e+1)%N],unpairedOnRight*a3));
						int d = (i+1)%N;
						int pairedScoreWithDangles = getQbFull(memo2, seq, N, nicks, d, e, 1, 1);
						QmEq1_S0[i][j] = alt(QmEq1_S0[i][j],combine(pairedScoreWithDangles, a2+unpairedOnRight*a3));
					}
				}
				if (L >= 3){
					//01
					for(int e = j, d=(i+1)%N, unpairedOnLeft = 0; d!=e; d=(d+1)%N, unpairedOnLeft++){
						if (containsNick(i,d,nicks)){
							break;
						}
						int pairedScoreWithLeftDangle = getQbFull(memo2, seq, N, nicks, d, e, 1, 0);
						QmEq1_01[i][j] = alt(QmEq1_01[i][j], combine(pairedScoreWithLeftDangle, a2 + unpairedOnLeft*a3));
						QeEq1_01[i][j] = alt(QeEq1_01[i][j], pairedScoreWithLeftDangle);
					}
					//10
					for(int d = i, e =(i+1)%N, unpairedOnRight = L - 2; e!=j; e=(e+1)%N, unpairedOnRight++){
						if (containsNick(e, j,nicks)){
							break;
						}
						int pairedScoreWithRightDangle = getQbFull(memo2, seq, N, nicks, d, e, 0, 1);
						QmEq1_10[i][j] = alt(QmEq1_10[i][j],combine(pairedScoreWithRightDangle, a2 + unpairedOnRight*a3));
					}
				}
				//002.25)
				//if [i,j] contains exactly 2, coaxially stacked, pairs and j is paired or j - 1 is paired and j is unpaired
				if (L >= 4){
					//11
					QmEq2Coax_11[i][j] = getQmEq2Coax(memo2, seq, N, nicks, i, j, 0, 0);
					QeEq2Coax_11[i][j] = getQeEq2Coax(memo2, seq, N, nicks, i, j, 0, 0);
				}
				if(L >= 5){
					//1S
					QmEq2Coax_1S[i][j] = getQmEq2Coax(memo2, seq, N, nicks, i, (j-1+N)%N, 0, 1);
					QeEq2Coax_1S[i][j] = getQeEq2Coax(memo2, seq, N, nicks, i, (j-1+N)%N, 0, 1);
					//S1
					QmEq2Coax_S1[i][j] = getQmEq2Coax(memo2, seq, N, nicks, (i+1)%N, j, 1, 0);
					QeEq2Coax_S1[i][j] = getQeEq2Coax(memo2, seq, N, nicks, (i+1)%N, j, 1, 0);
					//01
					for(int d = (i+1)%N, unpairedOnLeft = 1; d!=(j-1+N)%N; d=(d+1)%N, unpairedOnLeft++){
						QmEq2Coax_01[i][j] = alt(QmEq2Coax_01[i][j], combine(QmEq2Coax_S1[(d-1+N)%N][j], unpairedOnLeft*a3));
						QeEq2Coax_01[i][j] = alt(QeEq2Coax_01[i][j], QeEq2Coax_S1[(d-1+N)%N][j]);
					}
				}
				if (L >= 6){
					//SS
					QmEq2Coax_SS[i][j] = getQmEq2Coax(memo2, seq, N, nicks, (i+1)%N, (j-1+N)%N, 1, 1);
					QeEq2Coax_SS[i][j] = getQeEq2Coax(memo2, seq, N, nicks, (i+1)%N, (j-1+N)%N, 1, 1);
					//0S 
					for(int d = (i+1)%N, e = (j-1+N)%N, unpairedOnLeft = 1; d!=e; d=(d+1)%N, unpairedOnLeft++){
						QmEq2Coax_0S[i][j] = alt(QmEq2Coax_0S[i][j], combine(QmEq2Coax_SS[(d-1+N)%N][(e+1)%N], unpairedOnLeft*a3));
						QeEq2Coax_0S[i][j] = alt(QeEq2Coax_0S[i][j], QeEq2Coax_SS[(d-1+N)%N][(e+1)%N]);
					}
				}
				
				//002.5)
				//external, general case
				//If [i,j] contains 0 pairs
				if (!containsNick(i, j, nicks)){
					Qe00[i][j] = alt(Qe00[i][j], 0);
				}
				//The leftmost pair involves j, and it is not coaxially stacked
				Qe01[i][j] = alt(Qe01[i][j], QeEq1_01[i][j]);
				Qe11[i][j] = alt(Qe11[i][j], getQbFull(memo2, seq, N, nicks, i, j, 0, 0));
				//The leftmost pair does not involve j, and it is not coaxially stacked
				if (L >= 3){
					for(int e = (i+1)%N; e!=j; e=(e+1)%N){
						//e+1 is paired
						Qe00[i][j] = alt(Qe00[i][j], combine(QeEq1_01[i][e], Qe10[(e+1)%N][j]));
						Qe01[i][j] = alt(Qe01[i][j], combine(QeEq1_01[i][e], Qe11[(e+1)%N][j]));
						Qe11[i][j] = alt(Qe11[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 0), Qe11[(e+1)%N][j]));
						Qe10[i][j] = alt(Qe10[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 0), Qe10[(e+1)%N][j]));
						//e+1 is not paired
						Qe00[i][j] = alt(Qe00[i][j], combine(QeEq1_0S[i][(e+1)%N], Qe00[(e+1)%N][j]));
						Qe01[i][j] = alt(Qe01[i][j], combine(QeEq1_0S[i][(e+1)%N], Qe01[(e+1)%N][j]));
						Qe10[i][j] = alt(Qe10[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 1), Qe00[(e+1)%N][j]));
						Qe11[i][j] = alt(Qe11[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 1), Qe01[(e+1)%N][j]));
					}
				}
				//The leftmost 2 pairs are coaxially stacked, and j is involved (i.e. paired as part of the right helix)
				if (L >= 4){
					Qe01[i][j] = alt(Qe01[i][j], QeEq2Coax_01[i][j]);
					Qe11[i][j] = alt(Qe11[i][j], QeEq2Coax_11[i][j]);
				}
				//The leftmost 2 pairs are coaxially stacked, and j is not involved in this coaxial stacking
				if (L >= 5){
					for(int e = (i+3)%N; e!=j; e=(e+1)%N){
						//e+1 is paired
						Qe00[i][j] = alt(Qe00[i][j], combine(QeEq2Coax_01[i][e], Qe10[(e+1)%N][j]));
						Qe01[i][j] = alt(Qe01[i][j], combine(QeEq2Coax_01[i][e], Qe11[(e+1)%N][j]));
						Qe11[i][j] = alt(Qe11[i][j], combine(QeEq2Coax_11[i][e], Qe11[(e+1)%N][j]));
						Qe10[i][j] = alt(Qe10[i][j], combine(QeEq2Coax_11[i][e], Qe10[(e+1)%N][j]));
						//e+1 is not paired
						Qe00[i][j] = alt(Qe00[i][j], combine(QeEq2Coax_0S[i][(e+1)%N], Qe00[(e+1)%N][j]));
						Qe01[i][j] = alt(Qe01[i][j], combine(QeEq2Coax_0S[i][(e+1)%N], Qe01[(e+1)%N][j]));
						Qe10[i][j] = alt(Qe10[i][j], combine(QeEq2Coax_1S[i][(e+1)%N], Qe00[(e+1)%N][j]));
						Qe11[i][j] = alt(Qe11[i][j], combine(QeEq2Coax_1S[i][(e+1)%N], Qe01[(e+1)%N][j]));
					}
				}
				
				//002.75)
				//multibranch, >=2 pairs
				if (L >= 4){
					//the leftmost pair is not coaxially stacked with the second
					for(int e = (i+1)%N; e!=j; e=(e+1)%N){
						//e+1 is paired
						QmGe2_00[i][j] = alt(QmGe2_00[i][j], combine(QmEq1_01[i][e], getQmGe1_10(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_01[i][j] = alt(QmGe2_01[i][j], combine(QmEq1_01[i][e], getQmGe1_11(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_10[i][j] = alt(QmGe2_10[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 0), a2, getQmGe1_10(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_11[i][j] = alt(QmGe2_11[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 0), a2, getQmGe1_11(memo2, seq, N, nicks, (e+1)%N, j)));
						//e+1 is not paired
						QmGe2_00[i][j] = alt(QmGe2_00[i][j], combine(QmEq1_0S[i][(e+1)%N], getQmGe1_00(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_01[i][j] = alt(QmGe2_01[i][j], combine(QmEq1_0S[i][(e+1)%N], getQmGe1_01(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_10[i][j] = alt(QmGe2_10[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 1), a2, getQmGe1_00(memo2, seq, N, nicks, (e+1)%N, j)));
						QmGe2_11[i][j] = alt(QmGe2_11[i][j], combine(getQbFull(memo2, seq, N, nicks, i, e, 0, 1), a2, getQmGe1_01(memo2, seq, N, nicks, (e+1)%N, j)));
					}
					//the leftmost pair is coaxially stacked with the second and j is involved
					QmGe2_01[i][j] = alt(QmGe2_01[i][j], QmEq2Coax_01[i][j]);
					QmGe2_11[i][j] = alt(QmGe2_11[i][j], QmEq2Coax_11[i][j]); 
					//or j is not involved
					if (L >= 5){
						for(int e = (i+3)%N, basesOnRight = L - 4; e!=j; e=(e+1)%N, basesOnRight--){
							//e+1 is paired
							QmGe2_00[i][j] = alt(QmGe2_00[i][j], combine(QmEq2Coax_01[i][e], getQmGe1_10(memo2, seq, N, nicks, (e+1)%N, j)));
							QmGe2_01[i][j] = alt(QmGe2_01[i][j], combine(QmEq2Coax_01[i][e], getQmGe1_11(memo2, seq, N, nicks, (e+1)%N, j)));
							QmGe2_10[i][j] = alt(QmGe2_10[i][j], combine(QmEq2Coax_11[i][e], getQmGe1_10(memo2, seq, N, nicks, (e+1)%N, j)));
							QmGe2_11[i][j] = alt(QmGe2_11[i][j], combine(QmEq2Coax_11[i][e], getQmGe1_11(memo2, seq, N, nicks, (e+1)%N, j)));
							//e+1 is not paired
							QmGe2_00[i][j] = alt(QmGe2_00[i][j], combine(QmEq2Coax_0S[i][(e+1)%N], getQmGe0_00(memo2, seq, N, nicks, (e+1)%N, j, basesOnRight)));
							QmGe2_01[i][j] = alt(QmGe2_01[i][j], combine(QmEq2Coax_0S[i][(e+1)%N], getQmGe1_01(memo2, seq, N, nicks, (e+1)%N, j)));
							QmGe2_10[i][j] = alt(QmGe2_10[i][j], combine(QmEq2Coax_1S[i][(e+1)%N], getQmGe0_00(memo2, seq, N, nicks, (e+1)%N, j, basesOnRight)));
							QmGe2_11[i][j] = alt(QmGe2_11[i][j], combine(QmEq2Coax_1S[i][(e+1)%N], getQmGe1_01(memo2, seq, N, nicks, (e+1)%N, j)));
						}
					}
				} //end multibranch, general case (>=2 pairs)
							
			}//End scoring for [i,j].
		}
	}
	
	/**
	 * Performs the recursion to compute the MFE energy over all foldings of the given strands which are
	 * connected foldings and which do not contain any of the pairs in the constraint matrix.
	 * 
	 * If onlyUpperTriangle is true, then only recursions on subsequences [i,j], where i <= j, are computed.
	 * This is essentially always on, with the exception of circular folding.
	 * @param domain_markings 
	 * @param domain 
	 * @param ds 
	 */
	private void N2Fold(int[][][] memo2, int[] seq, int N, int[] nicks, boolean onlyUpperTriangle, FoldingConstraints constraints, SequenceMarker marker) {
		if (true) throw new RuntimeException("Phase 2 not yet written to fold GIR's");
		
		//Initialization: all is impossible.
		int maxLoopLen = 30;
		deepFill3(memo2, 1, N, Integer.MAX_VALUE, 
				EXTERNAL_00,
				EXTERNAL_01,
				EXTERNAL_10,
				EXTERNAL_11);
		deepFill3(memo2, N, N, Integer.MAX_VALUE,
				PAIRED);
		deepFill3(memo2, N, Math.min(N, maxLoopLen+1), Integer.MAX_VALUE, 
				LOOP_EXTENSION_J,
				LOOP_EXTENSION_JMIN1,
				LOOP_EXTENSION_JMIN2,
				RIGHT_BULGE_EXTENSION_J,
				RIGHT_BULGE_EXTENSION_JMIN1,
				LEFT_BULGE_EXTENSION_J,
				LEFT_BULGE_EXTENSION_JMIN1
				);
		
		int[][] Qb = memo2[PAIRED];
		int[][] QbNoStk = Qb; 
		int[][] QbStk = Qb;

		int slope_loop = eParams.getInteriorLoopSizeTermDx_deci(maxLoopLen);
		int slope_bulge = eParams.getBulgeLoopDx_deci(maxLoopLen);

		//Build up to larger subproblems:
		for( int L = 2; L <= N; L++ ){
			//Swap the Qx buffers
			swap(memo2, LOOP_EXTENSION_JMIN2, LOOP_EXTENSION_JMIN1);
			swap(memo2, LOOP_EXTENSION_JMIN1, LOOP_EXTENSION_J);

			int[][] Qxj = memo2[LOOP_EXTENSION_J];
			int[][] Qxj1 = memo2[LOOP_EXTENSION_JMIN1];
			int[][] Qxj2 = memo2[LOOP_EXTENSION_JMIN2];
			
			for(int[] row : Qxj){
				for(int i = 0; i <= maxLoopLen && i < row.length; i++) row[i] = Integer.MAX_VALUE;
			}

			//Swap the QxBulge buffers
			swap(memo2, RIGHT_BULGE_EXTENSION_J, RIGHT_BULGE_EXTENSION_JMIN1);
			swap(memo2, LEFT_BULGE_EXTENSION_J, LEFT_BULGE_EXTENSION_JMIN1);
			int[][] QxjBulgeR = memo2[RIGHT_BULGE_EXTENSION_J];
			int[][] QxjBulgeR1 = memo2[RIGHT_BULGE_EXTENSION_JMIN1];
			int[][] QxjBulgeL = memo2[LEFT_BULGE_EXTENSION_J];
			int[][] QxjBulgeL1 = memo2[LEFT_BULGE_EXTENSION_JMIN1];
			for(int[] row : QxjBulgeR){
				for(int i = 0; i <= maxLoopLen && i < row.length; i++) row[i] = Integer.MAX_VALUE;
			}
			for(int[] row : QxjBulgeL){
				for(int i = 0; i <= maxLoopLen && i < row.length; i++) row[i] = Integer.MAX_VALUE;
			}
			
			for( int i = 0; i <= N-1; i++){
				int j = (i + L - 1) % N;
				if (onlyUpperTriangle && j < i){
					break;
				}
				
				//Structure of the recursion:
				//001-Prelim) There is one top level pair in the interval i+1, j-1.
				//001) i and j are paired with eachother
					//001) They form an interior structure
					//001.5) They form an exterior structure
				//002) [i,j] is inside a multiloop (or exterior)
					//002) [i,j] contains exactly one pair
					//002.25) [i,j] contains 2 pairs exactly, coaxially stacked, and j is paired or j unpaired and j-1 paired
					//002.5) External loop general case (>= 0 pairs)
					//002.75) Multibranch general case (>=2 pairs)
								
				//001-Prelim) There is one top level pair in the interval i+1, j-1.
				if (L >= 4 && !(containsNick(i, (i+1)%N, nicks) && containsNick((j-1+N)%N,j,nicks))){					
					//Find the optimal interior loop and store it in Qxj

					//s > maxLoopLen, overestimate
					if (maxLoopLen <= L - 4){
						Qxj[i][maxLoopLen] = alt(Qxj[i][maxLoopLen], 
								combine(slope_loop, alt(Qxj1[i][maxLoopLen], Qxj1[(i+1)%N][maxLoopLen])));
					}
					if (maxLoopLen <= L - 3){
						QxjBulgeL[i][maxLoopLen] = alt(QxjBulgeL[i][maxLoopLen], 
								combine(slope_bulge, QxjBulgeL1[(i+1)%N][maxLoopLen]));
						QxjBulgeR[i][maxLoopLen] = alt(QxjBulgeR[i][maxLoopLen], 
								combine(slope_bulge, QxjBulgeR1[i][maxLoopLen]));
					}
					
					//s <= maxLoopLen, L1 > 4 and L2 > 4. 
					for(int s = 10; s <= L - 4 && s <= maxLoopLen; s++){
						//symmetry term is unmodified. However, the length of the interior loop is increased.
						Qxj[i][s] = alt(Qxj[i][s], combine(Qxj2[(i+1)%N][s-2], eParams.getInteriorLoopSizeTerm_deci(s) - eParams.getInteriorLoopSizeTerm_deci(s-2)));
					}
					//s <= maxLoopLen, L2 >= 4 and L1 = 4, 
					for(int L2 = 4, L1 = 4; L1 + L2 <= L - 4 && L1 + L2 <= maxLoopLen; L2++){
						int d = (i + L1 + 1)%N;
						int e = (j - L2 - 1 + N)%N;
						int s = L1 + L2;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int closingPair = eParams.getInteriorNNTerminal_deci(seq[e], seq[d], seq[(e+1)%N], seq[(d-1+N)%N]);
						Qxj[i][s] = alt(Qxj[i][s], 
								combine(Qb[d][e], 
										combine(eParams.getInteriorLoopSizeTerm_deci(s) + eParams.getNINIOAssymetry(L1, L2),
												closingPair
												)));
					}
					//s <= maxLoopLen, L1 > 4 and L2 = 4, 
					for(int L1 = 5, L2 = 4; L1 + L2 <= L - 4 && L1 + L2 <= maxLoopLen; L1++){
						int d = (i + L1 + 1)%N;
						int e = (j - L2 - 1 + N)%N;
						int s = L1 + L2;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int closingPair = eParams.getInteriorNNTerminal_deci(seq[e], seq[d], seq[(e+1)%N], seq[(d-1+N)%N]);
						Qxj[i][s] = alt(Qxj[i][s], 
								combine(Qb[d][e], 
										combine(eParams.getInteriorLoopSizeTerm_deci(s) + eParams.getNINIOAssymetry(L1, L2),
												closingPair
												)));
					}
					
					//s <= maxLoopLen, L1 = 0 and L2 >= 2
					for(int L1 = 0, L2 = 2; L1 + L2 <= L-4 && L1 + L2 <= maxLoopLen; L2++){
						int s = L1 + L2;
						int d = (i+L1+1)%N;
						int e = (j-L2-1+N)%N;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int bulge = eParams.getBulgeLoop_deci(L2); 
						int ATPenalty = eParams.getATPenalty_deci(seq[d], seq[e]);
						int rightBulge = combine(Qb[d][e], bulge + ATPenalty);
						QxjBulgeR[i][s] = alt(QxjBulgeR[i][s], rightBulge);
					}
					//s <= maxLoopLen, L2 = 0 and L1 >= 2
					for(int L1 = 2, L2 = 0; L1 + L2 <= L-4 && L1 + L2 <= maxLoopLen; L1++){
						int s = L1 + L2;
						int d = (i+L1+1)%N;
						int e = (j-L2-1+N)%N;
						if (containsNick(i, d, nicks) || containsNick(e, j, nicks)){
							break;
						}
						int bulge = eParams.getBulgeLoop_deci(L1); 
						int ATPenalty = eParams.getATPenalty_deci(seq[d], seq[e]);
						int leftBulge = combine(Qb[d][e], bulge + ATPenalty);
						QxjBulgeL[i][s] = alt(QxjBulgeL[i][s], leftBulge);
					}
				}
				
				//001)
				boolean canPair = true;
				//Overrides: Can i and j pair?
				if (constraints.preventPairing(i,j)){
					canPair = false;
				}
				//Based on base i and j, can i and j pair?
				if (Std.monomer.bindScore(seq[i], seq[j]) == 0){
					canPair = false;
				}
				if (canPair){
					if (L >= 4 && !(containsNick(i, (i+1)%N, nicks) && containsNick((j-1+N)%N,j,nicks))){
						//if the region contains no pairs, i.e. hairpin
						if (!containsNick(i, j, nicks)){
							QbNoStk[i][j] = alt(QbNoStk[i][j],eParams.getHairpinLoopDeltaG_deci(seq, N, i, L));
						}

						//if the region contains one pair
						
						//Bulges, interior loops where L1 >= 4 and L2 >= 4
						int internalLoopLeftClosing = eParams.getInteriorNNTerminal_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N]);
						for(int s = 8; s <= L - 4 && s <= maxLoopLen; s++){
							QbNoStk[i][j] = alt(QbNoStk[i][j], combine(internalLoopLeftClosing, Qxj[i][s]));
						}

						//Bulges, interior loops, where L1 < 4.
						for(int L1 = 0; L1 < 4; L1++){
							for(int L2 = 0; L1 + L2 <= L - 4 && L1 + L2 <= maxLoopLen; L2++){
								int d = (i+L1+1)%N;
								int e = (j-L2-1+N)%N;
								int s = L1 + L2;
								if (containsNick(i, d, nicks) || containsNick(e,j, nicks)){
									break; //there must be L1 (L2) connected bases on the left (right).
								}

								int stack = eParams.getNN_deci(seq[i], seq[j], seq[d], seq[e]);
								if (L1 == 0 && L2 == 0){
									//Continue a stack
									QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], stack));
									//Possibly terminate stack as well?
									QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], stack));
								} 
								//In the standard model, bulge loops have no sequence dependence on the unpaired bases. 
								if (L1 == 0 && L2 > 0){
									int bulge = eParams.getBulgeLoop_deci(L2); 
									//For L2 == 1 we consider that the stack is not terminated
									if (L2 == 1){
										bulge = combine(bulge, stack);
										QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], bulge));
										//Possibly terminate stack
										QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], bulge));
									} else {
										int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
										QbStk[i][j] = alt(QbStk[i][j], combine(ATPenalty, QxjBulgeR[i][s]));
									}
								}
								if (L2 == 0 && L1 > 0){
									//For L1 == 1 we consider that the stack is not terminated
									int bulge = eParams.getBulgeLoop_deci(L1);
									if (L1 == 1){
										bulge = combine(bulge, stack);
										QbStk[i][j] = alt(QbStk[i][j], combine(QbStk[d][e], bulge));
										//Possibly terminate stack
										QbStk[i][j] = alt(QbStk[i][j], combine(QbNoStk[d][e], bulge));
									} else {
										int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
										QbStk[i][j] = alt(QbStk[i][j], combine(ATPenalty, QxjBulgeL[i][s]));
									}
								}

								if (L1 > 0 && L2 > 0){
									QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], eParams.getInteriorLoop_deci(seq, N, i, j, L1, L2)));
								}
							}
						} //end L1 < 4

						//Bulges, interior loops where L1 >= 4 but L2 < 4
						for(int L2 = 0; L2 < 4; L2++){
							for(int L1 = 4; L1 + L2 <= L - 4 && L1 + L2 <= maxLoopLen; L1++){
								int d = (i+L1+1)%N;
								int e = (j-L2-1+N)%N;
								int s = L1 + L2;
								if (containsNick(i, d, nicks) || containsNick(e,j, nicks)){
									break; //only inner loop!
								}
								if (L2 == 0){
									int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
									QbStk[i][j] = alt(QbStk[i][j], combine(ATPenalty, QxjBulgeL[i][s]));
								}
								if (L2 > 0){
									QbNoStk[i][j] = alt(QbNoStk[i][j], combine(Qb[d][e], eParams.getInteriorLoop_deci(seq, N, i, j, L1, L2)));
								}
							}
						}

						//Check: the above cases have checked
						//(L1 > 4, L2 > 4), (L1 == 4, L2 >= 4), (L2 == 4, L1 > 4), (L1 < 4), (L2 < 4, L1 >= 4).

					} //end i and j are paired, and their pair is internal

					//001.5)
					//if i and j are paired to one another, and their pair is exterior
					if (L == 2 && containsNick(i, j, nicks)){
						int ATPenalty = eParams.getATPenalty_deci(seq[i], seq[j]);
						QbNoStk[i][j] = ATPenalty;  //Blunt end
					}
					//Read this line carefully: We rule out only the case where BOTH (i,i+1) and (j-1, j) are nicked.
					if (L >= 3 && !(containsNick(i, (i+1)%N, nicks) && containsNick((j-1+N)%N,j,nicks))){
						for(int k = 0; k < nicks.length; k++){
							int nick = nicks[k];
							if (containsNick(i,j,nick)){
								int score = 0;
								if (nick==i){
									//empty left loop
								} else { 
									score = combine(score, eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]));
								}
								if (nick==(j-1+N)%N){
									//empty right loop
								} else {
									score = combine(score, eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]));
								}

								score = combine(score, eParams.getATPenalty_deci(seq[i], seq[j]));

								QbNoStk[i][j] = alt(QbNoStk[i][j], score);
							}
						}
					}
					
					//If i and j are paired to one another, and the pair is either a new stack, or a single pair
					//Qb[i][j] = alt(QbNoStk[i][j], QbStk[i][j]);
					
					if (Qb[i][j] < 0){
						marker.mark(i);
						marker.mark(j);
					}
				} //end i, j are paired.
			}//End scoring for [i,j].
		}
		
		if (onlyUpperTriangle){
			//Only fill in [0][N-1]
			int[][] Qe00 = memo2[EXTERNAL_00];
			int[][] Qe01 = memo2[EXTERNAL_01];
			int[][] Qe10 = memo2[EXTERNAL_10];
			int[][] Qe11 = memo2[EXTERNAL_11];
			
			if (!containsNick(0, N-1, nicks)){
				Qe00[0][N-1] = 0; //No pairs, allowed only if folding a single strand.
			}
			
			for(int i = 0; i < N; i++){
				for(int j = i+1; j < N; j++){
					if (i == 0){
						if (j == N-1){
							//11
							Qe11[0][N-1] = alt(Qe11[0][N-1], 
									getQbFull(memo2, seq, N, nicks, i, j, 0, 0));
						} else {
							//10
							Qe10[0][N-1] = alt(Qe10[0][N-1], 
									getQbFull(memo2, seq, N, nicks, i, j, 0, 1));
						}
					} else {
						if (j == N-1){
							//01
							Qe01[0][N-1] = alt(Qe01[0][N-1], 
									getQbFull(memo2, seq, N, nicks, i, j, 1, 0));
						} else {
							//00
							Qe00[0][N-1] = alt(Qe00[0][N-1], 
									getQbFull(memo2, seq, N, nicks, i, j, 1, 1));
						}
					}
				}
			}
		}
	}
	
	/**
	 * Performs the recursion to compute the MFE energy over all foldings of the given strands which are
	 * connected foldings and which do not contain any of the pairs in the constraint matrix and which do
	 * not contain any interior loops, bulge loops, or multiloops (what remains is really just stems and hairpins)
	 * Furthermore, any introduced pairs are assumed to NOT engage in coaxial stacking with other helixes.
	 * 
	 * If onlyUpperTriangle is true, then only recursions on subsequences [i,j], where i <= j, are computed.
	 * This is essentially always on, with the exception of circular folding.
	 * @param domain_markings 
	 * @param domain 
	 * @param ds 
	 */
	private void N2Fold_NoLoops(int[][][] memo2, int[] seq, int N, int[] nicks, ConnectorFold[] connectors, boolean onlyUpperTriangle, FoldingConstraints constraints, SequenceMarker marker) {
		//Initialization: All structures are impossible.
		deepFill3(memo2, 1, N, Integer.MAX_VALUE, 
				EXTERNAL_00,
				EXTERNAL_01,
				EXTERNAL_10,
				EXTERNAL_11);
		
		int QbMemory = 4; //Remember 4 rows of the Qb matrix (i through i+3)
		
		if (memo2[PAIRED] == null || memo2[PAIRED].length < N){
			memo2[PAIRED] = new int[Math.max(N, QbMemory)][];
		}
		deepFill3(memo2, onlyUpperTriangle?QbMemory:N, N, Integer.MAX_VALUE, 
				PAIRED);
		int[][] Qb = memo2[PAIRED];

		//Only fill in [0][N-1]
		int[] Qe00 = memo2[EXTERNAL_00][0];
		int[] Qe01 = memo2[EXTERNAL_01][0];
		int[] Qe10 = memo2[EXTERNAL_10][0];
		int[] Qe11 = memo2[EXTERNAL_11][0];			
		
		//Zero free energy is defined to be the free energy of the structure with no pairing of the numbered bases of the GIR
		//However, we only measure the free energy of connected structures.
		if (!containsNick(0, N-1, nicks)){
			Qe00[N-1] = 0; //No pairs defined to have 0 delta G (connected structure only if a single strand (no nicks in [0, N-1]).
		}

		int maxI = N-1; //Maximum possible value of i, trivially.
		if (nicks.length >= 2 && onlyUpperTriangle){
			//Then, the outermost pair must have as its left base a base in the first strand 
			maxI = nicks[0];
			//So that when we enter the loop, we fill Qb[maxI].
			for(int k = 0; k < QbMemory; k++){
				swap(Qb, k%N, (maxI + 1 + k)%N);
			}
		}

		for( int i = maxI; i >= 0; i--){
			if (onlyUpperTriangle){
				//Memory saving trick: Qb[i][xxx] depends only on Qb[i+1][xxx], so drop all but 2 rows of Qb.
				//Collect row i+QbMemory to use as row i
				swap(Qb, (i + QbMemory)%N, i);
			}
			int[] Qbi = Qb[i];
			
			int minL = 2; //Minimum possible value of L so that i and j are different bases, trivially.
			if (nicks.length >= 2 && onlyUpperTriangle){
				//Then, the outermost pair must have as its right base a base in the last strand 
				int minJ = nicks[nicks.length-2]+1;
				//It holds that minJ > maxI.
				minL = minJ - maxI + 1;
			}
			
			ConnectorSummary conSum_ij = new ConnectorSummary();
			for( int L = minL; L <= N; L++ ){
				int j = (i + L - 1) % N;
				if (onlyUpperTriangle && j < i){
					break;
				}
				int countNicks = countNicks(i,j,nicks);
				boolean hasTopNick = containsNick(i, (i+1)%N, nicks);
				boolean hasBottomNick = containsNick((j-1+N)%N,j,nicks);
			
				if (countNicks > 0){
					conSum_ij = null;
				} else {
					if (L >= 3){
						conSum_ij.extendInterval(seq[(j-2+N)%N], connectors[(j-2+N)%N], 
								seq[(j-1+N)%N], connectors[(j - 1 + N)%N], 
								seq[j], connectors[j]);
					}
				}
				Qbi[j] = Integer.MAX_VALUE;
				
				//i and j are paired with eachother
				if (!constraints.preventPairing(i,j) && !(Std.monomer.bindScore(seq[i], seq[j]) == 0)){
					//i and j pair, enclose an interior region, and all bases between i and j are unpaired
					if (L >= 4 && (countNicks==0)){
						Qbi[j] = alt(Qbi[j], intLoopDG_1pair(seq, N, i, L, conSum_ij));
					}
					//i and j pair, enclose a nicked region, and all bases between i and j are unpaired
					if (countNicks > 0){
						if (countNicks >= 2){
							throw new RuntimeException("Folding of >=2 strands not yet implemented");
						}
						Qbi[j] =  alt(Qbi[j], extLoopDG_1pair(seq, N, i, L, nicks, connectors));
					}
					
					//i and j pair, and [i,j] contains one non-connector pair, [d,e]. 
					//The region formed between these two pairs may be interior or exterior. 
					//Recall that no disconnected structures are allowed.
					for(int L1 = 0; L1 <= QbMemory - 2; L1++){
						for(int L2 = 0; L2 <= QbMemory - 2; L2++){
							if (L1 + L2 + 4 <= L){
								int d = (i+L1+1)%N;
								int e = (j-L2-1+N)%N;
								int score = Qb[d][e];
								if (score != Integer.MAX_VALUE){
									score = combine(score, loopDG_2pair(seq, N, nicks, i, L, L1, L2, connectors));
									Qbi[j] =  alt(Qbi[j], score);		
								}
							}
						}
					}
					
					if (Qbi[j] < 0){
						marker.mark(i);
						marker.mark(j);
					}

					if (onlyUpperTriangle){
						//Outermost pair is i,j
						int score = Qbi[j];
						//Add score of external loop j to i
						if (containsNick(0, i, nicks)){
							//Disconnected!
							score = Integer.MAX_VALUE;
						} else
						if (containsNick(j, N-1, nicks)){
							//Disconnected!
							score = Integer.MAX_VALUE;
						} else
						{
							score =  combine(score, extLoopDG_1pair(seq, N, j, N - L + 2, nicks, connectors));
						}
						
						if (i == 0){
							if (j == N-1){
								//11
								Qe11[N-1] = alt(Qe11[N-1],score);
							} else {
								//10
								Qe10[N-1] = alt(Qe10[N-1],score);
							}
						} else {
							if (j == N-1){
								//01
								Qe01[N-1] = alt(Qe01[N-1],score);
							} else {
								//00
								Qe00[N-1] = alt(Qe00[N-1],score);
							}
						}
					}
				}//End "can pair"
			}//End scoring for [i,j].
		}
	}	

	/**
	 * Method assumes the interval described contains no nicks
	 */
	private int intLoopDG_1pair(int[] seq, int N, int i, int L, ConnectorSummary interval) {
		int a1 = eParams.getMultibranchBase_deci();
		int a2 = eParams.getMultibranchBranch_deci();
		int a3 = eParams.getMultibranchUnpairedBase_deci();
		int j = (i + L - 1)%N;
		if (interval.numConnectors == 0){
			//Hairpin
			int score = eParams.getHairpinLoopDeltaG_deci(seq, N, i, L);
			return score;
		}
		
		if (interval.numConnectors == 1){
			//Interior loop scoring uses no dangles, so we lose all dangles!
			int lostDangles = 0;
			if (interval.numConnectors > 0){
				//we lose the leftmost connector's left dangle
				lostDangles -= eParams.getDangleBottom_deci(interval.leftPairBase1, interval.leftPairBase0, seq[(i + interval.unpairedBasesLeft)%N]);
			}
			if (interval.numConnectors > 0){
				//we lose the rightmost connector's right dangle
				lostDangles -= eParams.getDangleTop_deci(interval.rightPairBase1, interval.rightPairBase0, seq[(j - interval.unpairedBasesRight + N)%N]);
			}
			
			//Internal loop
			int L1 = interval.unpairedBasesLeft;
			int L2 = interval.unpairedBasesRight;

			int based = interval.leftPairBase0;
			int basee = interval.leftPairBase1;
			
			int loopScore = intLoopDG_2pair_noConnectors(seq, N, i, j, based, basee, L1, L2);
			return combine(loopScore, lostDangles);
		}

		if (interval.numConnectors >= 2){
			int lostDangles = 0;
			if (interval.numConnectors > 0 && interval.unpairedBasesLeft == 0){
				//we lose the leftmost connector's left dangle because i,j are paired
				lostDangles -= eParams.getDangleBottom_deci(interval.leftPairBase1, interval.leftPairBase0, seq[i]);
			}
			if (interval.numConnectors > 0 && interval.unpairedBasesRight == 0){
				//we lose the rightmost connector's right dangle because i,j are paired
				lostDangles -= eParams.getDangleTop_deci(interval.rightPairBase1, interval.rightPairBase0, seq[j]);
			}
			
			//Multiloop
			int score = a1;
			score += a2*(interval.numConnectors + 1);
			score += a3*interval.unpairedBases; //which are unpaired
			
			//two new dangles
			if (interval.unpairedBasesLeft > 0){
				score = combine(score, eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]));
			}
			if (interval.unpairedBasesRight > 0){
				score = combine(score, eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]));
			}
			score = combine(lostDangles, score);

			//AT penalty added for the new closing of one new stack
			score = combine(score, eParams.getATPenalty_deci(seq[i], seq[j]));
			
			return score;
		}
		throw new RuntimeException("Assertion error");
	}

	/**
	 * Method assumes the region described contains no nicks
	 */
	private int intLoopDG_2pair_noConnectors(int[] seq, int N, int i, int j, int based, int basee, int L1, int L2) {
		int stack = eParams.getNN_deci(seq[i], seq[j], based, basee);
		if (L1 == 0 && L2 == 0){
			return stack;
		} 
		//In the standard model, bulge loops have no sequence dependence on the unpaired bases. 
		if (L1 == 0 && L2 > 0){
			int bulge = eParams.getBulgeLoop_deci(L2); 
			//For L2 == 1 we consider that the stack is not terminated
			if (L2 == 1){
				return combine(bulge, stack);
			} else {
				int ATPenalty = eParams.getATPenalty_deci(based, basee) + eParams.getATPenalty_deci(seq[i], seq[j]);
				return combine(bulge, ATPenalty);
			}
		}
		if (L2 == 0 && L1 > 0){
			//For L1 == 1 we consider that the stack is not terminated
			int bulge = eParams.getBulgeLoop_deci(L1);
			if (L1 == 1){
				return combine(bulge, stack);
			} else {
				int ATPenalty = eParams.getATPenalty_deci(based, basee) + eParams.getATPenalty_deci(seq[i], seq[j]);
				return combine(bulge, ATPenalty);
			}
		}
		//true interior loops
		if (L1 > 0 && L2 > 0){
			if (L1 == 1 && L2 == 1){
				return eParams.get1x1InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], based, basee);
			} else if (L1 == 1 && L2 == 2){
				return eParams.get1x2InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(j-2+N)%N], based, basee);
			} else if (L1 == 2 && L2 == 1){
				return eParams.get2x1InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(i+2)%N], based, basee);
			} else if (L1 == 2 && L2 == 2){
				return eParams.get2x2InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(i+2)%N], seq[(j-2+N)%N], based, basee);
			} else {
				int NINIOAssym = eParams.getNINIOAssymetry(L1, L2);
				int leftTerminator = eParams.getInteriorNNTerminal_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N]);
				if (leftTerminator == Integer.MAX_VALUE){
					return Integer.MAX_VALUE;
				}
				int ep1 = (j - L2)%N;
				int dm1 = (i + L1)%N;
				int rightTerminator = eParams.getInteriorNNTerminal_deci(basee, based, seq[ep1], seq[dm1]);
				if (rightTerminator == Integer.MAX_VALUE){
					return Integer.MAX_VALUE;
				}
				return combine(
						combine(
							combine(eParams.getInteriorLoopSizeTerm_deci(L1 + L2), NINIOAssym), 
							leftTerminator), rightTerminator);
			}
		}
		throw new RuntimeException("Assertion error");
	}
	
	/**
	 * Generic (can contain nicks) 2-pair defined region scorer.
	 */
	private int loopDG_2pair(int[] seq, int N, int[] nicks, int i, int L, int L1, int L2, ConnectorFold[] connectors) {
		int j = (i + L - 1)%N;
		int d = (i + L1 + 1)%N;
		int e = (j - L2 - 1 + N)%N;
		
		int lostDangles = 0;
		int addedDangles = 0; //added to multiloop or exterior loop scores
		
		if (connectors[i].rightN > 0){
			//we lose the dangle from i
			lostDangles -= eParams.getDangleBottom_deci(connectors[i].rightPairBase1, connectors[i].rightPairBase0, seq[i]);
		} else if (L1 > 0){
			//we gain the dangle from i+1
			addedDangles += eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]);
		}
		if (connectors[d].leftN > 0){
			//we lose the dangle from d
			lostDangles -= eParams.getDangleTop_deci(connectors[d].leftPairBase1, connectors[d].leftPairBase0, seq[d]);
		} else if (L1 > 0){
			//we gain the dangle from d-1
			addedDangles += eParams.getDangleBottom_deci(seq[e], seq[d], seq[(d-1+N)%N]);
		}
		if (connectors[j].leftN > 0){
			//we lose the dangle from j
			lostDangles -= eParams.getDangleTop_deci(connectors[j].leftPairBase1, connectors[j].leftPairBase0, seq[j]);
		} else if (L2 > 0){
			//we gain the dangle from j-1
			addedDangles += eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]);
		}
		if (connectors[e].rightN > 0){
			//we lose the dangle from e
			lostDangles -= eParams.getDangleBottom_deci(connectors[e].rightPairBase1, connectors[e].rightPairBase0, seq[e]);
		} else if (L2 > 0){
			//we gain the dangle from e+1
			addedDangles += eParams.getDangleTop_deci(seq[e], seq[d], seq[(e+1)%N]);
		}

		boolean leftNicked = containsNick(i,d,nicks);
		boolean rightNicked = containsNick(e,j,nicks);
		
		if (leftNicked && rightNicked){
			//Disconnected!
			return Integer.MAX_VALUE;
		} else if (leftNicked && !rightNicked){
			//Loop is exterior, but we may be able to coaxially stack the two stems
			if (L2 <= 1){ //one intermediate base is O.K.
				return combine(lostDangles, addedDangles, eParams.getNN_deci(seq[i], seq[j], seq[d], seq[e]));
			} else {
				//Two stacks
				int score = eParams.getATPenalty_deci(seq[i], seq[j]);
				score = combine(score, eParams.getATPenalty_deci(seq[d], seq[e]));
				return combine(lostDangles, addedDangles, score);
			}
		} else if (rightNicked && !leftNicked){
			//Loop is exterior, but we may be able to coaxially stack the two stems
			if (L1 <= 1){ //one intermediate base is O.K.
				return combine(lostDangles, addedDangles, eParams.getNN_deci(seq[i], seq[j], seq[d], seq[e]));
			} else {
				//Two stacks
				int score = eParams.getATPenalty_deci(seq[i], seq[j]);
				score = combine(score, eParams.getATPenalty_deci(seq[d], seq[e]));
				return combine(lostDangles, addedDangles, score);
			}
		} else {
			//Interior region.
			int numLoops = 2; //For the two introduced pairs
			for(int k = i; k != d; k = (k+1)%N){
				numLoops += connectors[k].rightN;
			}
			for(int k = e; k != j; k = (k+1)%N){
				numLoops += connectors[k].rightN;
			}

			if (numLoops >= 3){	
				int a1 = eParams.getMultibranchBase_deci();
				int a2 = eParams.getMultibranchBranch_deci();
				int a3 = eParams.getMultibranchUnpairedBase_deci();
				//Interior multiloop
				int loop = a1;
				loop += a2*(numLoops);
				loop += a3*(L1 + L2);
				//AT penalty added for the new closing of two new stacks
				int score = combine(loop, eParams.getATPenalty_deci(seq[i], seq[j]));
				score = combine(score, eParams.getATPenalty_deci(seq[d], seq[e]));
				return combine(lostDangles, addedDangles, score);
			} else {
				if (lostDangles != 0){
					throw new RuntimeException("Assertion error");
				}
				//Interior, non-multiloop
				int score = intLoopDG_2pair_noConnectors(seq, N, i, j, seq[d], seq[e], L1, L2);
				return score;
			}
		}
	}

	/**
	 * Method assumes there is precisely one nick in [i, i+L].
	 */
	private int extLoopDG_1pair(int[] seq,  int N, int i, int L, int[] nicks, ConnectorFold[] connectors) {
		int j = (i + L - 1)%N;
		boolean hasTopNick = containsNick(i, (i+1)%N, nicks);
		boolean hasBottomNick = containsNick((j-1+N)%N,j,nicks);
		//AT penalty
		int score = eParams.getATPenalty_deci(seq[i], seq[j]);
		int lostDangles = 0;
		
		if (connectors[i].rightN > 0){
			//we lose i's dangle
			lostDangles -= eParams.getDangleBottom_deci(connectors[i].rightPairBase1, connectors[i].rightPairBase0, seq[i]);
		} else {
			//we gain a dangle on i
			if (L > 2 && !hasTopNick){
				score = combine(score, eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]));
			}
		}
		
		if (connectors[j].leftN > 0){
			//we lose j's dangle
			lostDangles -= eParams.getDangleTop_deci(connectors[j].leftPairBase1, connectors[j].leftPairBase0, seq[j]);
		} else {
			//we gain a dangle on j
			if (L > 2 && !hasBottomNick){
				score = combine(score, eParams.getDangleBottom_deci(seq[i], seq[j], seq[(j-1+N)%N]));
			}
		}
		
		return combine(lostDangles, score);
	}
	
	private int getQmGe1(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int bonusIfIUnpaired, int bonusIfJUnpaired) {
		int best = combine(getQmGe1_00(memo2, seq, N, nicks, i, j), bonusIfIUnpaired, bonusIfJUnpaired);
		best = alt(best, combine(getQmGe1_01(memo2, seq, N, nicks, i, j), bonusIfIUnpaired));
		best = alt(best, combine(getQmGe1_10(memo2, seq, N, nicks, i, j), bonusIfJUnpaired));
		best = alt(best, getQmGe1_11(memo2, seq, N, nicks, i, j));
		return best;
	}

	private int getQmGe0_00(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int L) {
		int a3 = eParams.getMultibranchUnpairedBase_deci();
		return alt(alt(memo2[MULTIBRANCH_EQ1_00][i][j], memo2[MULTIBRANCH_GE2_00][i][j]), a3*L);
	}
	private int getQmGe1_00(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j) {
		return alt(memo2[MULTIBRANCH_EQ1_00][i][j], memo2[MULTIBRANCH_GE2_00][i][j]);
	}
	private int getQmGe1_10(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j) {
		return alt(memo2[MULTIBRANCH_EQ1_10][i][j], memo2[MULTIBRANCH_GE2_10][i][j]);
	}
	private int getQmGe1_01(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j) {
		return alt(memo2[MULTIBRANCH_EQ1_01][i][j], memo2[MULTIBRANCH_GE2_01][i][j]);
	}
	private int getQmGe1_11(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j) {
		int a2 = eParams.getMultibranchBranch_deci();
		return alt(combine(getQbFull(memo2, seq, N, nicks, i, j, 0, 0), a2), memo2[MULTIBRANCH_GE2_11][i][j]);
	}


	private int getQmEq2Coax(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int dangleOnLeft, int dangleOnRight) {
		int[][] Qb = memo2[PAIRED];

		int a2 = eParams.getMultibranchBranch_deci();
		int a3 = eParams.getMultibranchUnpairedBase_deci();
		
		int d = i;
		int g = j;
		if (dangleOnLeft==1){
			if (containsNick((i-1+N)%N, i, nicks)){
				return Integer.MAX_VALUE;
			}
		}
		if (dangleOnRight==1){
			if (containsNick(j, (j+1)%N, nicks)){
				return Integer.MAX_VALUE;
			}
		}
		//d is paired with e, f is paired with g, and the two helixes coaxially stack. Apply dangles.
		int best = Integer.MAX_VALUE;
		for(int e = (i+1)%N; e!=(g-1+N)%N; e=(e+1)%N){
			for(int f = (e+1)%N, y=0; y<=1; f=(f+1)%N, y++){
				int dangles = 0;
				if (dangleOnLeft==1){
					dangles = combine(dangles, eParams.getDangleBottom_deci(seq[e], seq[i], seq[(i-1+N)%N]));
				}
				if (dangleOnRight==1){
					dangles = combine(dangles, eParams.getDangleTop_deci(seq[j], seq[f], seq[(j+1)%N]));
				}
				int coax = getCoaxialStackBonus(memo2, seq, N, nicks, e, d, f, g, y);
				best = alt(best, combine(dangles, coax, Qb[d][e], Qb[f][g], y*a3+2*a2));
			}
		}
		return best;
	}
	private int getQeEq2Coax(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int dangleOnLeft, int dangleOnRight) {
		int[][] Qb = memo2[PAIRED];
		
		int d = i;
		int g = j;
		if (dangleOnLeft==1){
			if (containsNick((i-1+N)%N, i, nicks)){
				return Integer.MAX_VALUE;
			}
		}
		if (dangleOnRight==1){
			if (containsNick(j, (j+1)%N, nicks)){
				return Integer.MAX_VALUE;
			}
		}
		//d is paired with e, f is paired with g, and the two helixes coaxially stack. Apply dangles.
		int best = Integer.MAX_VALUE;
		for(int e = (i+1)%N; e!=(g-1+N)%N; e=(e+1)%N){
			for(int f = (e+1)%N, y=0; y<=1; f=(f+1)%N, y++){
				int dangles = 0;
				if (dangleOnLeft==1){
					dangles = combine(dangles, eParams.getDangleBottom_deci(seq[e], seq[i], seq[(i-1+N)%N]));
				}
				if (dangleOnRight==1){
					dangles = combine(dangles, eParams.getDangleTop_deci(seq[j], seq[f], seq[(j+1)%N]));
				}
				int coax = getCoaxialStackBonus(memo2, seq, N, nicks, e, d, f, g, y);
				best = alt(best, combine(dangles, coax, Qb[d][e], Qb[f][g]));
			}
		}
		return best;
	}



	private int getQbFull(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int dangleOnLeft, int dangleOnRight) {
		int score = memo2[PAIRED][i][j];
		if (dangleOnLeft==1){
			if (containsNick((i-1+N)%N, i, nicks)){
				return Integer.MAX_VALUE;
			}
			score = combine(score, eParams.getDangleBottom_deci(seq[j], seq[i], seq[(i-1+N)%N]));
		}
		if (dangleOnRight==1){
			if (containsNick(j, (j+1)%N, nicks)){
				return Integer.MAX_VALUE;
			}
			score = combine(score, eParams.getDangleTop_deci(seq[j], seq[i], seq[(j+1)%N]));
		}
		score = combine(score, eParams.getATPenalty_deci(seq[i], seq[j]));
		return score;
	}



	/**
	 * Returns the MFE structure of an exterior loop [i,j], applying appropriate bonuses applied bonusIfIUnpaired (if i is unpaired)
	 * and bonusIfJUnpaired (if j is unpaired)
	 */
	private int getQe(int[][][] memo2, int i, int j, int bonusIfIUnpaired, int bonusIfJUnpaired) {
		int best = combine(memo2[EXTERNAL_00][i][j], bonusIfIUnpaired, bonusIfJUnpaired);
		best = alt(best, combine(memo2[EXTERNAL_01][i][j], bonusIfIUnpaired));
		best = alt(best, combine(memo2[EXTERNAL_10][i][j], bonusIfJUnpaired));
		best = alt(best, memo2[EXTERNAL_11][i][j]);
		return best;
	}
	/**
	 * Returns the MFE structure of an multibranch interior with at least two internal pairs spanning [i,j], applying appropriate bonuses applied bonusIfIUnpaired (if i is unpaired)
	 * and bonusIfJUnpaired (if j is unpaired)
	 */
	private int getQmGe2(int[][][] memo2, int i, int j, int bonusIfIUnpaired, int bonusIfJUnpaired) {
		int best = combine(memo2[MULTIBRANCH_GE2_00][i][j], bonusIfIUnpaired, bonusIfJUnpaired);
		best = alt(best, combine(memo2[MULTIBRANCH_GE2_01][i][j], bonusIfIUnpaired));
		best = alt(best, combine(memo2[MULTIBRANCH_GE2_10][i][j], bonusIfJUnpaired));
		best = alt(best, memo2[MULTIBRANCH_GE2_11][i][j]);
		return best;
	}
	


	/**
	 * [d,e] is coaxial to [i,j], with y bases on the continuous backbone of the coaxial stack, and the left backbone of [i,j] is continuous to d.
	 * In other words, d = (i+y+1)%N. 
	 */
	private int getCoaxialStackBonus(int[][][] memo2, int[] seq, int N, int[] nicks, int i, int j, int d, int e, int y) {
		if (containsNick(i, d, nicks)){
			return Integer.MAX_VALUE;
		}
		int dangles = 0;
		if (y > 0){
			dangles = combine(eParams.getDangleTop_deci(seq[i], seq[j], seq[(i+1)%N]), 
					eParams.getDangleBottom_deci(seq[e], seq[d], seq[(d-1+N)%N]));
		} else {
			dangles = combine(dangles, 
					eParams.getDangleBottom_deci(seq[i], seq[j], seq[e]), 
					eParams.getDangleTop_deci(seq[e], seq[d], seq[j]));
		}
		return dangles;
	}
	/**
	 * [d,e] is coaxial to [i,j], with y bases on the continuous backbone of the coaxial stack, and the left backbone of [i,j] is continuous to d.
	 * In other words, d = (i+y)%N. 
	 * 
	 * This is a special case of getCoaxialStackBonus, where the bases are provided and y = 0.
	 */
	private int getCoaxialStackBonus(int basei, int basej, int based, int basee) {
		int dangles = combine(eParams.getDangleBottom_deci(basei, basej, basee), 
					eParams.getDangleTop_deci(basee, based, basej));
		return dangles;
	}
	
	private int countNicks(int i, int j, int[] nicks){
		int count = 0;
		for(int nick = 0; nick < nicks.length; nick++){
			if (containsNick(i, j, nicks[nick])){
				count++;
			}
		}
		return count;
	}
	
	private boolean containsNick(int i, int j, int[] nicks) {
		for(int k = 0; k < nicks.length; k++){
			if (containsNick(i,j,nicks[k])){
				return true;
			}
		}
		return false;
	}
	private boolean containsNick(int i, int j, int nick){
		if (i <= j){
			if (nick >= i && nick < j){
				return true;
			}
		} else {
			if (nick >= i || nick < j){
				return true;
			}
		}
		return false;
	}

	private static void deepFill3(int[][][] memo, int rows, int cols, int value) {
		for(int index = 0; index < memo.length; index++){
			if (memo[index]==null || memo[index].length < rows || memo[index][0].length < cols){
				memo[index] = new int[rows][cols];
			}
			for(int row = 0; row < rows; row++){
				for(int col = 0; col < cols; col++){
					memo[index][row][col] = value;
				}
			}
		}	
	}
	private static void deepFill3(int[][][] memo, int rows, int cols, int value, int ... indeces) {
		for(int index : indeces){
			if (memo[index]==null || memo[index].length < rows){
				memo[index] = new int[rows][];
			}
			for(int row = 0; row < rows; row++){
				if (memo[index][row] == null || memo[index][row].length < cols){
					memo[index][row] = new int[cols];
				}
				for(int col = 0; col < cols; col++){
					memo[index][row][col] = value;
				}
			}
		}	
	}
	private static void swap(Object[] ar, int i, int j) {
		Object tmp = ar[i];
		ar[i] = ar[j];
		ar[j] = tmp;
	}

	/**
	 * Combine two choices, which must be taken together
	 */
	private static int combine(int i, int j) {
		if (i == Integer.MAX_VALUE || j == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return i + j;
	}
	/**
	 * Combine three choices, which must be taken together
	 */
	private static int combine(int i, int j, int k) {
		if (i == Integer.MAX_VALUE || j == Integer.MAX_VALUE || k == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return i + j + k;
	}
	/**
	 * Combine four choices, which must be taken together
	 */
	private static int combine(int i, int j, int k, int m) {
		if (i == Integer.MAX_VALUE || j == Integer.MAX_VALUE || k == Integer.MAX_VALUE || m == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return i + j + k + m;
	}
	/**
	 * Combine five choices, which must be taken together
	 */
	private static int combine(int i, int j, int k, int m, int n) {
		if (i == Integer.MAX_VALUE || j == Integer.MAX_VALUE || k == Integer.MAX_VALUE || m == Integer.MAX_VALUE || n == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return i + j + k + m + n;
	}

	/**
	 * Choose between two alternatives
	 */
	private static int alt(int i, int j) {
		/*
		if (i < -1e7 || j < -1e7){
			throw new RuntimeException("UNDERFLOW");
		}
		*/
		if (i < j){
			return i;
		}
		return j;
	}
	private static int altALERT(String message, int i, int j) {
		if (i < j){
			System.out.println(message+" : "+i);
			return i;
		}
		return j;
	}

	public double mfeNoDiag(GeneralizedInteractiveRegion GeneralizedInteractiveRegion,
			GeneralizedInteractiveRegion GeneralizedInteractiveRegion2, int[][] domain,
			int[][] domain_markings) {
		throw new RuntimeException("TODO");
	}

	public double mfeStraight(GeneralizedInteractiveRegion seq1, GeneralizedInteractiveRegion seq2, int[][] domain, int[][] domain_markings, int markLeft, int markRight, int joffset) {
		int score = 0;
		int N1 = seq1.length(domain);
		int N2 = seq2.length(domain);
		int L = 2;
		for(int i = 0; i + L <= N1; i++){
			int j = N2 - 1 - i;
			if (j - L + 1 >= 0){
				score = alt(score, combine(score, 
						eParams.getNN_deci(base(seq1, i, domain), base(seq2, j, domain), base(seq1, i+1, domain), base(seq2, j-1, domain))
								));
			} else {
				break;
			}
		}
		return score / 100.0;
	}
	public ArrayList<Point> getTraceback() {
		throw new RuntimeException("TODO");
	}
	public double[][] getScoreMatrix(int len1, int len2) {
		throw new RuntimeException("TODO");
	}
}
