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
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import circdesigna.ZipExtractor;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.config.CircDesigNASystemElement;



/**
 * The actual Experimental Parameters Database.
 */
public class ExperimentalDuplexParamsImpl extends CircDesigNASystemElement implements NAExperimentDatabase{
	public ExperimentalDuplexParamsImpl(CircDesigNAConfig config){
		super(config);
		System.out.print("Unpacking Thermo Parameters ... ");
		
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
	private int getNormalBase(int nonnormalBase){
		return Std.monomer.getNormalBaseFromZero(nonnormalBase);
	}
	private static final int D2DECI(double value){
		return (int)(value * 100);
	}
	private int[][][][] getNNdeltaG_deci;
	public double getNNdeltaG(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaG_deci[W][X][Y][Z]/100.0;
	}
	public int getNNdeltaG_deci(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaG_deci[W][X][Y][Z];
	}
	//X,Y,D,primeEnd(0,1)
	private int[][][][] getDangle_deci;
	public double getDanglePenalty(int X, int Y, int D, boolean PrimeEnd3) {
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		D = getNormalBase(D);
		return getDangle_deci[X][Y][D][PrimeEnd3?1:0]/100.0;
	}
	public int getDanglePenalty_deci(int X, int Y, int D, boolean PrimeEnd3) {
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		D = getNormalBase(D);
		return getDangle_deci[X][Y][D][PrimeEnd3?1:0];
	}
	private int[][][][] getNNdeltaGterm_deci;
	public double getNNdeltaGterm(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaGterm_deci[W][X][Y][Z]/100.0;
	}
	public int getNNdeltaGterm_deci(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaGterm_deci[W][X][Y][Z];
	}
	
	
	public static void main(String[] args) throws Throwable{
		CircDesigNAConfig config = new CircDesigNAConfig();
		//config.setMode(CircDesigNAConfig.RNA_MODE);
		final ExperimentalDuplexParamsImpl x = new ExperimentalDuplexParamsImpl(config);

		//parseMFoldParamsFile("dangles.add(new DangleScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\dangle.dgd",true,"));");
		//parseMFoldParamsFile("nns.add(new NearestNeighborScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\stack.dgd",false,"));");
		//parseMFoldParamsFile("tns.add(new TerminalMismatchPairScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\tstackh.dgd",false,"));");

		//x.parseMFoldParamsFile("dangles.add(new DangleScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\dangle.dg",true,"));");
		//x.parseMFoldParamsFile("nns.add(new NearestNeighborScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\stack.dg",false,"));");
		//x.parseMFoldParamsFile("tns.add(new TerminalMismatchPairScore(","C:\\Users\\Benjamin\\PROGRAMMING\\Libraries\\Compiled-Proprietary\\unafold\\unafold-3.8\\data\\rules\\tstackh.dg",false,"));");
		
		//System.out.println(x.getDeltaGAssoc(2, 310.15));
		System.out.println(x.getNNdeltaG(A,T,G,T));
		System.out.println(x.getNNdeltaGterm(A,T,G,G));
		
		//OK!
	}

	/**
	 * Called after makeTable*
	 */
	public void reformScoreLists(ArrayList<NearestNeighborScore> nns,
			ArrayList<TerminalMismatchPairScore> tns,
			ArrayList<DangleScore> dangles) {
		//Ok! parse the arraylists to actual tables.
		getNNdeltaG_deci = new int[4][4][4][4];
		getNNdeltaGterm_deci = new int[4][4][4][4];
		getDangle_deci = new int[4][4][4][2];
		for(NearestNeighborScore nnsi : nns){
			getNNdeltaG_deci[nnsi.W]
			            [nnsi.X]
			             [nnsi.Y]
			              [nnsi.Z] = D2DECI(nnsi.score);
		}
		for(TerminalMismatchPairScore tnsi : tns){
			getNNdeltaGterm_deci[tnsi.W]
			                [tnsi.X]
			                 [tnsi.Y]
			                  [tnsi.Z] = D2DECI(tnsi.score);
		}
		for(DangleScore dang : dangles){
			getDangle_deci[dang.W]
			          [dang.X]
			           [dang.Z]
			            [dang.is3PrimeEnd?1:0] = D2DECI(dang.score);
		}
		
		//Print them out:
		//System.out.println("NNDelta: "+Arrays.deepToString(getNNdeltaG));
		//System.out.println("Terminating: "+Arrays.deepToString(getNNdeltaGterm));
		//System.out.println("Dangles: "+Arrays.deepToString(getDangle));
	}
	public class DangleScore extends NearestNeighborScore{
		private boolean is3PrimeEnd;
		public DangleScore(int W, int X, String prEnd, int Z, double score) {
			super(W, X, A, Z, score);
			is3PrimeEnd = prEnd.contains("3");
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
	
	public double getDeltaGAssoc(int numStrands, double T) {
		double H37Bimolecular = .2;
		double G37Bimolecular = 1.96;
		double S37Bimolecular = (G37Bimolecular - H37Bimolecular)/310.15;
		//in kcal / mol. Using the most recent experiments of this constant.
		double kB = 1.3806504 * 6.02214179 / 4.184 / 1000;
		double GBimolecular = H37Bimolecular + T * S37Bimolecular;
		return (numStrands-1)*(GBimolecular - kB * T * Math.log(WaterDensity(T)));
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
