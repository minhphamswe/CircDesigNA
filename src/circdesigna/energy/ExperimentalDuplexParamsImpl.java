package circdesigna.energy;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;


import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.Config.CircDesigNASystemElement;

/**
 * The actual Experimental Parameters Database.
 */
public class ExperimentalDuplexParamsImpl extends CircDesigNASystemElement implements NAExperimentDatabase{
	public ExperimentalDuplexParamsImpl(CircDesigNAConfig config){
		super(config);
		
		ZipInputStream paramZip = null;
		try {
			System.out.print("Unpacking Thermo Parameters ... ");
			paramZip = new ZipInputStream(ExperimentalDuplexParamsImpl.class.getResourceAsStream("/parameters.zip"));
			System.out.println("Done (1)");
		} catch (Throwable e){
			//Try loading it as a file.
			try {
				paramZip = new ZipInputStream(new FileInputStream("parameters.zip"));
				System.out.println("Done (2)");
				//System.out.println("Loaded parameters file from disk.");
			} catch (Throwable f){
				throw new RuntimeException("Could not load the parameters.zip file. Please include this file in the working directory!");
			}
		}
		ZipEntry nextEntry;
		String dG = null, dH = null;
		try {
			while((nextEntry= paramZip.getNextEntry())!=null){
				if (nextEntry.getName().startsWith(config.getParameterName())){
					ByteArrayOutputStream baos = new ByteArrayOutputStream();
					byte[] buf = new byte[1024];
					int read = -1;
					while((read=paramZip.read(buf))>0){
						baos.write(buf, 0, read);
					}
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
	private double[][][][] getNNdeltaG;
	public double getNNdeltaG(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaG[W][X][Y][Z];
	}
	//X,Y,D,primeEnd(0,1)
	private double[][][][] getDangle;
	public double getDanglePenalty(int X, int Y, int D, boolean PrimeEnd3) {
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		D = getNormalBase(D);
		return getDangle[X][Y][D][PrimeEnd3?1:0];
	}
	private double[][][][] getNNdeltaGterm;
	public double getNNdeltaGterm(int W, int X, int Y, int Z) {
		W = getNormalBase(W);
		X = getNormalBase(X);
		Y = getNormalBase(Y);
		Z = getNormalBase(Z);
		return getNNdeltaGterm[W][X][Y][Z];
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
		getNNdeltaG = new double[4][4][4][4];
		getNNdeltaGterm = new double[4][4][4][4];
		getDangle = new double[4][4][4][2];
		for(NearestNeighborScore nnsi : nns){
			getNNdeltaG[nnsi.W]
			            [nnsi.X]
			             [nnsi.Y]
			              [nnsi.Z] = nnsi.score;
		}
		for(TerminalMismatchPairScore tnsi : tns){
			getNNdeltaGterm[tnsi.W]
			                [tnsi.X]
			                 [tnsi.Y]
			                  [tnsi.Z] = tnsi.score;
		}
		for(DangleScore dang : dangles){
			getDangle[dang.W]
			          [dang.X]
			           [dang.Z]
			            [dang.is3PrimeEnd?1:0] = dang.score;
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
