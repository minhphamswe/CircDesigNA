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
package circdesigna.energy;

import static circdesigna.abstractpolymer.DnaDefinition.A;
import static circdesigna.abstractpolymer.DnaDefinition.C;
import static circdesigna.abstractpolymer.DnaDefinition.G;
import static circdesigna.abstractpolymer.DnaDefinition.T;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import circdesigna.ZipExtractor;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;



/**
 * The actual Experimental Parameters Database.
 */
public class ExperimentalDuplexParams extends CircDesigNASystemElement {
	public ExperimentalDuplexParams(CircDesigNAConfig config){
		super(config);
		//System.out.print("Unpacking Thermo Parameters ... ");
		
		ZipInputStream paramZip = ZipExtractor.getFile("parameters.zip");
		ZipEntry nextEntry;
		String dG = null, dH = null;
		try {
			while((nextEntry= paramZip.getNextEntry())!=null){
				if (nextEntry.getName().startsWith(config.getParameterName())){
					ByteArrayOutputStream baos = ZipExtractor.readFully(paramZip);
					if (nextEntry.getName().endsWith(".dG")){
						dG = baos.toString();	
					}
					if (nextEntry.getName().endsWith(".dH")){
						dH = baos.toString();
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		StandardizedThermoFileLoader.makeTable(this,dG,dH);
	}
	private static final int getNormalBase(int nonnormalBase){
		return nonnormalBase - 1;
		//return Std.monomer.getNormalBaseFromZero(nonnormalBase);
	}
	private static final int D2DECI(double value){
		//Round.
		return (int)Math.round(value * 100);
	}
	private static final ArrayList<Integer> D2DECIList(double[] value){
		ArrayList<Integer> toRet = new ArrayList<Integer>();
		for(int i = 0; i < value.length; i++){
			toRet.add(D2DECI(value[i]));
		}
		return toRet;
	}
	private static final int[] D2DECI(double[] value){
		int[] toRet = new int[value.length];
		for(int i = 0; i < toRet.length; i++){
			toRet[i] = D2DECI(value[i]);
		}
		return toRet;
	}
	
	private int[][][][] getNN_deci;
	public double getNN(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNN_deci[W][X][Y][Z]/100.0;
	}
	/**
	 * Returns (in deci-kcal/mol) the stack energy of the nearest neighbor 5'-WY-3' over 3'-XZ-5'
	 */
	public int getNN_deci(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNN_deci[W][X][Y][Z];
	}
	
	private int[][][] getDangleTop_deci;
	public int getDangleTop_deci(int X1, int X2, int Y) {
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		Y = getNormalBase(Y);
		return getDangleTop_deci[X1][X2][Y];
	}
	private int[][][] getDangleBottom_deci;
	public int getDangleBottom_deci(int X1, int X2, int Y) {
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		Y = getNormalBase(Y);
		return getDangleBottom_deci[X1][X2][Y];
	}
	private int[][][][] getInteriorNNTerminal_deci;
	public double getInteriorNNTerminal(int W, int X, int Y, int Z) {
		return getInteriorNNTerminal_deci(W,X,Y,Z) / 100.0;
	}
	/**
	 * Returns (in deci-kcal/mol) the energy contribution of a a stack whose terminal bases are inside an interior loop
	 * (or hairpin). The energy is of the pair 5'-WY-3' over 3'-XZ-5', where W and X are paired but Y and Z are not paired.
	 */
	public int getInteriorNNTerminal_deci(int W, int X, int Y, int Z) {
		return getHairpinNNTerminal_deci(W, X, Y, Z);
		/*
		int ATPenalty = getATPenalty_deci(W, X);
		int loop = getHairpinNNTerminal_deci(W, X, Y, Z);
		if (loop == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return loop + ATPenalty;
		*/
	}
	public int getHairpinNNTerminal_deci(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		int loop = getInteriorNNTerminal_deci[W][X][Y][Z]; 
		if (loop == Integer.MAX_VALUE){
			return Integer.MAX_VALUE;
		}
		return loop;
	}
 
	/**
	 * The loop consisting of bases seq[i], seq[(i+1)%N], ... seq[(i+L)%N], where the first and last bases are 
	 * the closing pair of the hairpin.
	 */
	private int[][][][][][] getTetraLoop_deci;
	private int[][][][][] getTriLoop_deci;
	private ArrayList<Integer> getHairpinLoopGeneral_deci;
	public int getHairpinLoopDeltaG_deci(int[] seq, int N, int i, int L) {
		if (L < 5){
			return Integer.MAX_VALUE;
		}
		//Loop contribution

		//See Matthews, Mol Biol, 2009
		int s = L - 2 - 1; //0 index corresponds to hairpin loop of size 1 (which is impossible, incidentally)
		while (s >= getHairpinLoopGeneral_deci.size() && s > 8){
			int n = getHairpinLoopGeneral_deci.size() + 1;
			double T = 310.5;
			getHairpinLoopGeneral_deci.add(getHairpinLoopGeneral_deci.get(8) + D2DECI(1.75 * getR_kcalmol() * T * Math.log(n/9.0)));
		}
		
		int energy = getHairpinLoopGeneral_deci.get(s);
		
		//Bonuses.
		int j = (i+L-1)%N;
		if (L == 5){
			energy += getTriLoop_deci[base(seq,i,N)][base(seq,i+1,N)][base(seq,i+2,N)][base(seq,i+3,N)][base(seq,i+4,N)];
		}
		if (L == 6){
			energy += getTetraLoop_deci[base(seq,i,N)][base(seq,i+1,N)][base(seq,i+2,N)][base(seq,i+3,N)][base(seq,i+4,N)][base(seq,i+5,N)];
		}
		if (L > 5){
			int terminalBonus = getHairpinNNTerminal_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N]);
			if (terminalBonus == Integer.MAX_VALUE){
				return Integer.MAX_VALUE;
			}
			energy += terminalBonus;
		} else {
			int ATPenalty = getATPenalty_deci(seq[i], seq[j]);
			energy += ATPenalty;
		}
		return energy;
	}

	private ArrayList<Integer> getBulgeLoop_deci;
	public int getBulgeLoop_deci(int s) {
		s--; //0-index corresponds to s==1
		//See Matthews, Mol Biol, 2009
		while (s >= getBulgeLoop_deci.size() && s > 5){
			int n = getBulgeLoop_deci.size() + 1;
			double T = 310.5;
			getBulgeLoop_deci.add(getBulgeLoop_deci.get(5) + D2DECI(1.75 * getR_kcalmol() * T * Math.log(n/6.0)));
		}
		return getBulgeLoop_deci.get(s);
	}
	private ArrayList<Integer> getInteriorLoopSizeTerm_deci;
	public int getInteriorLoopSizeTerm_deci(int s) {
		s--; //0 index corresponds to interior loop of size 1 (which is impossible, incidentally)
		//See Matthews, Mol Biol, 2009. 
		while (s >= getInteriorLoopSizeTerm_deci.size() && s > 5){
			int n = getInteriorLoopSizeTerm_deci.size() + 1;
			double T = 310.5;
			getInteriorLoopSizeTerm_deci.add(getInteriorLoopSizeTerm_deci.get(5) + D2DECI(1.75 * getR_kcalmol() * T * Math.log(n/6.0)));
		}
		return getInteriorLoopSizeTerm_deci.get(s);
	}
	
	private int[][][][][][] get1x1InteriorLoop_deci;
	/**
	 * Loop of 5'-A X1 C-3' over 3'-B X2 D-5'
	 */
	private int get1x1InteriorLoop_deci(int A, int B, int X1, int X2, int C, int D) {
		A = getNormalBase(A);
		B = getNormalBase(B);
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		C = getNormalBase(C);
		D = getNormalBase(D);
		return get1x1InteriorLoop_deci[A][B][X1][X2][C][D];
	}
	private int[][][][][][][] get1x2InteriorLoop_deci;
	/**
	 * Loop of 5'-A X1 C-3' over 3'-B X2 Y2 D-5'
	 */
	private int get1x2InteriorLoop_deci(int A, int B, int X1, int X2, int Y2, int C, int D) {
		A = getNormalBase(A);
		B = getNormalBase(B);
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		Y2 = getNormalBase(Y2);
		C = getNormalBase(C);
		D = getNormalBase(D);
		return get1x2InteriorLoop_deci[A][B][X1][X2][Y2][C][D];
	}
	/**
	 * Loop of 5'-A X1 Y1 C-3' over 3'-B X2 D-5'
	 */
	private int get2x1InteriorLoop_deci(int A, int B, int X1, int X2, int Y1, int C, int D) {
		A = getNormalBase(A);
		B = getNormalBase(B);
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		Y1 = getNormalBase(Y1);
		C = getNormalBase(C);
		D = getNormalBase(D);
		//180* rotation to get score
		return get1x2InteriorLoop_deci[D][C][X2][X1][Y1][B][A];
	}
	private int[][][][][][][][] get2x2InteriorLoop_deci;
	/**
	 * Loop of 5'-A X1 Y1 C-3' over 3'-B X2 Y2 D-5'
	 */
	private int get2x2InteriorLoop_deci(int A, int B, int X1, int X2, int Y1, int Y2, int C, int D) {
		A = getNormalBase(A);
		B = getNormalBase(B);
		X1 = getNormalBase(X1);
		X2 = getNormalBase(X2);
		Y1 = getNormalBase(Y1);
		Y2 = getNormalBase(Y2);
		C = getNormalBase(C);
		D = getNormalBase(D);
		return get2x2InteriorLoop_deci[A][B][X1][X2][Y1][Y2][C][D];
	}
	
	public int getInteriorLoop_deci(int[] seq, int N, int i, int j, int L1, int L2) {
		i %= N;
		j %= N;
		int d = (i+L1+1)%N;
		int e = (j-L2-1+N)%N;
		if (d==i || e == j || i < 0 || j < 0){
			return Integer.MAX_VALUE;
		}
		
		int loop;
		if (L1 == 1 && L2 == 1){
			loop = get1x1InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[d], seq[e]);
		} else if (L1 == 1 && L2 == 2){
			loop = get1x2InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(j-2+N)%N], seq[d], seq[e]);
		} else if (L1 == 2 && L2 == 1){
			loop = get2x1InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(i+2)%N], seq[d], seq[e]);
		} else if (L1 == 2 && L2 == 2){
			loop = get2x2InteriorLoop_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N], seq[(i+2)%N], seq[(j-2+N)%N], seq[d], seq[e]);
		} else {
			int NINIOAssym = getNINIOAssymetry(L1, L2);
			int leftTerminator = getInteriorNNTerminal_deci(seq[i], seq[j], seq[(i+1)%N], seq[(j-1+N)%N]);
			if (leftTerminator == Integer.MAX_VALUE){
				return Integer.MAX_VALUE;
			}
			int rightTerminator = getInteriorNNTerminal_deci(seq[e], seq[d], seq[(e+1)%N], seq[(d-1+N)%N]);
			if (rightTerminator == Integer.MAX_VALUE){
				return Integer.MAX_VALUE;
			}
			loop = getInteriorLoopSizeTerm_deci(L1 + L2) + NINIOAssym + leftTerminator + rightTerminator;
		}
		return loop;
	}

	private int[] NINIOAssym_deci;
	//Tested: I Verified that nupack does it this way.
	public int getNINIOAssymetry(int L1, int L2) {
		int m = Math.min(Math.min(4,L1),L2);
		return Math.min(NINIOAssym_deci[4], NINIOAssym_deci[m-1]*Math.abs(L1 - L2));
	}
	private int[] MBterms_deci;
	public int getMultibranchBase_deci() {
		return MBterms_deci[0];
	}	
	public int getMultibranchBranch_deci() {
		return MBterms_deci[1];
	}
	public int getMultibranchUnpairedBase_deci() {
		return MBterms_deci[2];
	}	
	private int ATPenalty_deci;
	/**
	 * Returns either 0 or the penalty incurred by a stack terminating with the pair X1 X2.
	 */
	public int getATPenalty_deci(int X1, int X2){
		if ((X1 == C && X2 == G) || (X1 == G && X2 == C)){
			return 0;	
		}
		return ATPenalty_deci;
	}

	
	private final int base(int[] seq, int i, int N) {
		return getNormalBase(seq[i%N]);
	}
	
	public static void main(String[] args) throws Throwable{
		CircDesigNAConfig config = new CircDesigNAConfig();
		//config.setMode(CircDesigNAConfig.RNA_MODE);
		final ExperimentalDuplexParams x = new ExperimentalDuplexParams(config);

		//parseMFoldParamsFile("dangles.add(new DangleScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\dangle.dgd",true,"));");
		//parseMFoldParamsFile("nns.add(new NearestNeighborScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\stack.dgd",false,"));");
		//parseMFoldParamsFile("tns.add(new TerminalMismatchPairScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\tstackh.dgd",false,"));");

		//x.parseMFoldParamsFile("dangles.add(new DangleScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\dangle.dg",true,"));");
		//x.parseMFoldParamsFile("nns.add(new NearestNeighborScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\stack.dg",false,"));");
		//x.parseMFoldParamsFile("tns.add(new TerminalMismatchPairScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\tstackh.dg",false,"));");
		
		//System.out.println(x.getDeltaGAssoc(2, 310.15));
		System.out.println(x.getNN(A,T,G,T));
		System.out.println(x.getInteriorNNTerminal(A,T,G,G));
		
		//OK!
	}

	/**
	 * Called after makeTable*
	 */
	public void setSequenceSpecificStructures(Collection<NearestNeighborScore> nns,
			Collection<TerminalMismatchPairScore> tns,
			Collection<DangleScore> topDangles,
			Collection<DangleScore> bottomDangles,
			Collection<HairpinLoop> specialHairpins, 
			Collection<InteriorLoop> specialInteriorLoops
			) {
		//Ok! parse the arraylists to actual tables.
		getNN_deci = new int[4][4][4][4];
		deepFill4(getNN_deci, Integer.MAX_VALUE);
		getInteriorNNTerminal_deci = new int[4][4][4][4];
		deepFill4(getInteriorNNTerminal_deci, Integer.MAX_VALUE);
		getDangleTop_deci = new int[4][4][4];
		deepFill3(getDangleTop_deci, Integer.MAX_VALUE);
		getDangleBottom_deci = new int[4][4][4];
		deepFill3(getDangleBottom_deci, Integer.MAX_VALUE);
		getTetraLoop_deci = new int[4][4][4][4][4][4];
		//deepFill6(getTetraLoop_deci, 0);
		getTriLoop_deci = new int[4][4][4][4][4];
		//deepFill5(getTriLoop_deci, 0);
		get1x1InteriorLoop_deci = new int[4][4][4][4][4][4];
		deepFill6(get1x1InteriorLoop_deci, Integer.MAX_VALUE);
		get1x2InteriorLoop_deci = new int[4][4][4][4][4][4][4];
		deepFill7(get1x2InteriorLoop_deci, Integer.MAX_VALUE);
		get2x2InteriorLoop_deci = new int[4][4][4][4][4][4][4][4];
		deepFill8(get2x2InteriorLoop_deci, Integer.MAX_VALUE);
		
		for(NearestNeighborScore nnsi : nns){
			getNN_deci[nnsi.W]
			            [nnsi.X]
			             [nnsi.Y]
			              [nnsi.Z] = D2DECI(nnsi.score);
		}
		for(TerminalMismatchPairScore tnsi : tns){
			getInteriorNNTerminal_deci[tnsi.W]
			                [tnsi.X]
			                 [tnsi.Y]
			                  [tnsi.Z] = D2DECI(tnsi.score);
		}
		for(DangleScore dang : topDangles){
			getDangleTop_deci[dang.X1]
			          [dang.X2]
			           [dang.Y] 
			        		= D2DECI(dang.score);
		}
		for(DangleScore dang : bottomDangles){
			getDangleBottom_deci[dang.X1]
			          [dang.X2]
			           [dang.Y] 
			        		= D2DECI(dang.score);
		}
		for(HairpinLoop spec : specialHairpins){
			if (spec.bases.length==5){
				getTriLoop_deci[spec.bases[0]]
						[spec.bases[1]]
							[spec.bases[2]]
								[spec.bases[3]]
									[spec.bases[4]] = D2DECI(spec.score);
			} else
			if (spec.bases.length==6){
				getTetraLoop_deci[spec.bases[0]]
						[spec.bases[1]]
							[spec.bases[2]]
								[spec.bases[3]]
									[spec.bases[4]]
										[spec.bases[5]] = D2DECI(spec.score);
			} else {
				throw new RuntimeException("Only special triloops and tetraloops supported");
			}
		}
		for(InteriorLoop in : specialInteriorLoops){
			int L1 = in.basesTop.length - 2;
			int[] A1 = in.basesTop;
			int L2 = in.basesBottom.length - 2;
			int[] A2 = in.basesBottom;
			if (L1 == 1 && L2 == 1){
				get1x1InteriorLoop_deci[A1[0]][A2[0]][A1[1]][A2[1]][A1[2]][A2[2]] = D2DECI(in.score);
			} else
			if (L1 == 1 && L2 == 2){
				get1x2InteriorLoop_deci[A1[0]][A2[0]]
							[A1[1]][A2[1]]
									[A2[2]]
										[A1[2]][A2[3]] = D2DECI(in.score);
			} else
			if (L1 == 2 && L2 == 2){
				get2x2InteriorLoop_deci[A1[0]][A2[0]]
						[A1[1]][A2[1]]
								[A1[2]][A2[2]]
									[A1[3]][A2[3]] = D2DECI(in.score);
			} else {
				throw new RuntimeException("Only 1x1, 1x2, and 2x2 loops supported in input (2x1 are derived from 1x2)");
			}
		}
		
		//Print them out:
		//System.out.println("NNDelta: "+Arrays.deepToString(getNNdeltaG));
		//System.out.println("Terminating: "+Arrays.deepToString(getNNdeltaGterm));
		//System.out.println("Dangles: "+Arrays.deepToString(getDangle));
	}

	public void setLoopEnergies(double[] hairpinLoop, double[] bulgeLoop,
			double[] interiorLoop, double[] assymetry) {
		getHairpinLoopGeneral_deci = D2DECIList(hairpinLoop);
		getBulgeLoop_deci = D2DECIList(bulgeLoop);
		getInteriorLoopSizeTerm_deci = D2DECIList(interiorLoop);
		NINIOAssym_deci = D2DECI(assymetry);
	}
	public void setMultibranchTerms(double[] MLterms) {
		MBterms_deci = D2DECI(MLterms);
	}

	public void setATPenalty(double atPenalize) {
		ATPenalty_deci = D2DECI(atPenalize);
	}
	private void deepFill4(int[][][][] mat, int val) {
		for(int[][][] sblock : mat){
			for(int[][] block : sblock){
				for(int[] row : block){
					Arrays.fill(row,val);
				}
			}
		}
	}
	private void deepFill3(int[][][] sblock, int val) {
		for(int[][] block : sblock){
			for(int[] row : block){
				Arrays.fill(row,val);
			}
		}
	}
	private void deepFill6(int[][][][][][] sblock, int val) {
		for(int[][][][][] block : sblock){
			for(int[][][][] subblock : block){
				for(int[][][] hypermatrix : subblock){
					deepFill3(hypermatrix, val);
				}
			}
		}
	}
	private void deepFill7(int[][][][][][][] sblock, int val) {
		for(int[][][][][][] block : sblock){
			deepFill6(block, val);
		}
	}
	private void deepFill8(int[][][][][][][][] block, int val) {
		for(int[][][][][][][] subblock : block){
			deepFill7(subblock, val);
		}
	}
	
	public class EnergyList {
		public EnergyList(double[] energies){
			this.energies = energies;
		}
		public final double[] energies; 
	}
	public class HairpinLoop {
		public HairpinLoop(int[] bases_, double d){
			this.bases = new int[bases_.length];
			for(int i = 0; i < bases_.length; i++){
				this.bases[i] = getNormalBase(bases_[i]);
			}
			this.score = d;
		}
		public String toString(){
			StringBuffer sb = new StringBuffer();
			for(int i = 0; i < bases.length; i++){
				sb.append(Std.monomer.displayBase(Std.monomer.getMonomers()[bases[i]]));
			}
			sb.append(" "+score);
			return sb.toString();
		}
		public final int[] bases;
		public final double score; 
	}
	public class InteriorLoop {
		public InteriorLoop(int[] basesTop_, int[] basesBottom_, double d){
			this.basesTop = new int[basesTop_.length];
			for(int i = 0; i < basesTop_.length; i++){
				this.basesTop[i] = getNormalBase(basesTop_[i]);
			}
			this.basesBottom = new int[basesBottom_.length];
			for(int i = 0; i < basesBottom_.length; i++){
				this.basesBottom[i] = getNormalBase(basesBottom_[i]);
			}
			this.score = d;
		}
		public final int[] basesTop;
		public final int[] basesBottom;
		public final double score; 
		public String toString(){
			StringBuffer sb = new StringBuffer();
			for(int i = 0; i < basesTop.length; i++){
				sb.append(Std.monomer.displayBase(Std.monomer.getMonomers()[basesTop[i]]));
			}
			sb.append(" ");
			for(int i = 0; i < basesBottom.length; i++){
				sb.append(Std.monomer.displayBase(Std.monomer.getMonomers()[basesBottom[i]]));
			}
			sb.append(" "+score);
			return sb.toString();
		}
	}
	
	public class DangleScore {
		private double score;
		private int Y;
		private int X1;
		private int X2;

		public DangleScore(int X1, int X2, int Y, double score) {
			this.X1 = getNormalBase(X1);
			this.X2 = getNormalBase(X2);
			this.Y = getNormalBase(Y);
			this.score = score;
		}
	}
	public class TerminalMismatchPairScore extends NearestNeighborScore{
		public TerminalMismatchPairScore(int W, int X, int Y, int Z,
				double score) {
			super(W, X, Y, Z, score);
		}
	}
	public class NearestNeighborScore{
		public NearestNeighborScore(int W, int X, int Y, int Z, double score){
			this.W = getNormalBase(W);
			this.X = getNormalBase(X);
			this.Y = getNormalBase(Y);
			this.Z = getNormalBase(Z);
			this.score = score;
		}
		public final int W,X,Y,Z;
		public final double score;
	}
	
	/**
	 * Taken from nupack!!!!
	 */
	public static final double WaterDensity(double T) {
		/* 
		     Calculates the number of moles of water per liter at T degrees

		     Density of water calculated using data from:
		     Tanaka M., Girard, G., Davis, R., Peuto A.,
		     Bignell, N.   Recommended table for the denisty
		     of water..., Metrologia, 2001, 38, 301-309
		 */
		
		double FreezingPointOfWater = 273.15;  
			
		T -= FreezingPointOfWater;

		double a1 = -3.983035;
		double a2 = 301.797;
		double a3 = 522528.9;
		double a4 = 69.34881;
		double a5 = 999.974950;


		return a5 * (1 - (T+a1)*(T+a1)*(T+a2)/a3/(T+a4)) / 18.0152;
	}
	
	public double getR_kcalmol(){
		return 1.3806504 * 6.02214179 / 4.184 / 1000;
	}
	public double getDeltaGAssoc(int numStrands, double T) {
		double H37Bimolecular = .2;
		double G37Bimolecular = 1.96;
		double S37Bimolecular = (G37Bimolecular - H37Bimolecular)/310.15;
		//in kcal / mol. Using the most recent experiments of this constant.
		double GBimolecular = H37Bimolecular + T * S37Bimolecular;
		return (numStrands-1)*(GBimolecular - getR_kcalmol() * T * Math.log(WaterDensity(T)));
	}
	
	//Old
	private void parseMFoldParamsFile(String prefix, String file, boolean dangleMode, String postum) throws FileNotFoundException {
		Scanner in = new Scanner(new File(file));
		int[] ordering = new int[]{A,C,G,T};
		int A = 0, B = 0;
		String prEnd = "\"3pr\"";
		while(in.hasNextLine()){
			String line = in.nextLine();
			if (line.contains("3' <-- 5'")){
				int X = 0, Y = 0;
				while(in.hasNextLine()){
					String line3 = in.nextLine();
					String[] line2 = line3.trim().split("\\s+");
					if (line2.length!=16){
						//throw new RuntimeException("Non-mult-16 line: "+line3);
						break;
					}
					String sep = ",";
					B = 0;
					for(Y = 0; Y < line2.length; Y++){
						if (!line2[Y].equals(".")){
							System.out.print(prefix);
							System.out.print(Std.monomer.displayBase(ordering[A])+sep+
									Std.monomer.displayBase(ordering[B])+sep);
							if (!dangleMode){
								System.out.print(
										Std.monomer.displayBase(ordering[X])+sep);
							} else {
								System.out.print(prEnd+sep);
							}
							System.out.print(
									Std.monomer.displayBase(ordering[Y%4])+sep);
							System.out.print(
									line2[Y].equals(".")?"0":line2[Y]);
							System.out.println(postum);
						}
						if (Y%4==3){
							B++;
						}
					}
					X++;
				}
				A++;
				if (dangleMode){
					if (A >= ordering.length){
						A = 0;
						prEnd = "\"5pr\"";
					}
				}
			}
		}
	}
}
