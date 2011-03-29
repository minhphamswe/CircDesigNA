package DnaDesign.impl;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;
import static DnaDesign.AbstractPolymer.RnaDefinition.U;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import DnaDesign.ExperimentDatabase;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.Config.CircDesigNASystemElement;

/**
 * The actual Experimental Parameters Database.
 */
public class ExperimentalDuplexParamsImpl extends CircDesigNASystemElement implements ExperimentDatabase{
	public ExperimentalDuplexParamsImpl(CircDesigNAConfig config){
		super(config);
		
		ZipInputStream paramZip = new ZipInputStream(ExperimentalDuplexParamsImpl.class.getResourceAsStream("/parameters.zip")); 
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
		if (true){
			StandardizedThermoFileLoader.makeTable(this,dG,dH);
		} else {
			if (Std.isDNAMode()){
				makeTable_DNA();
			} else {
				makeTable_RNA();
			}
		}
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

	private void makeTable_RNA(){
		//****References:***//
		//Zucker's UNAFOLD 3.0 parameters.
		//The were taken from the following papers
		//Matthews 1999
		//Sort-of Zucker, 1999 (JMB) <--- this is RNA scores, but may have been extrapolated.
	
		ArrayList<NearestNeighborScore> nns = new ArrayList();
		{
			//AUTOWRITTEN
			nns.add(new NearestNeighborScore(A,U,A,U,-0.90));
			nns.add(new NearestNeighborScore(A,U,C,G,-2.10));
			nns.add(new NearestNeighborScore(A,U,G,C,-1.70));
			nns.add(new NearestNeighborScore(A,U,G,U,-0.50));
			nns.add(new NearestNeighborScore(A,U,U,A,-0.90));
			nns.add(new NearestNeighborScore(A,U,U,G,-1.00));
			nns.add(new NearestNeighborScore(C,G,A,U,-1.80));
			nns.add(new NearestNeighborScore(C,G,C,G,-2.90));
			nns.add(new NearestNeighborScore(C,G,G,C,-2.00));
			nns.add(new NearestNeighborScore(C,G,G,U,-1.20));
			nns.add(new NearestNeighborScore(C,G,U,A,-1.70));
			nns.add(new NearestNeighborScore(C,G,U,G,-1.90));
			nns.add(new NearestNeighborScore(G,C,A,U,-2.30));
			nns.add(new NearestNeighborScore(G,U,A,U,-1.10));
			nns.add(new NearestNeighborScore(G,C,C,G,-3.40));
			nns.add(new NearestNeighborScore(G,U,C,G,-2.10));
			nns.add(new NearestNeighborScore(G,C,G,C,-2.90));
			nns.add(new NearestNeighborScore(G,C,G,U,-1.40));
			nns.add(new NearestNeighborScore(G,U,G,C,-1.90));
			nns.add(new NearestNeighborScore(G,U,G,U,-0.40));
			nns.add(new NearestNeighborScore(G,C,U,A,-2.10));
			nns.add(new NearestNeighborScore(G,C,U,G,-2.10));
			nns.add(new NearestNeighborScore(G,U,U,A,-1.00));
			nns.add(new NearestNeighborScore(G,U,U,G,+1.50));
			nns.add(new NearestNeighborScore(U,A,A,U,-1.10));
			nns.add(new NearestNeighborScore(U,G,A,U,-0.80));
			nns.add(new NearestNeighborScore(U,A,C,G,-2.30));
			nns.add(new NearestNeighborScore(U,G,C,G,-1.40));
			nns.add(new NearestNeighborScore(U,A,G,C,-1.80));
			nns.add(new NearestNeighborScore(U,A,G,U,-0.80));
			nns.add(new NearestNeighborScore(U,G,G,C,-1.20));
			nns.add(new NearestNeighborScore(U,G,G,U,-0.20));
			nns.add(new NearestNeighborScore(U,A,U,A,-0.90));
			nns.add(new NearestNeighborScore(U,A,U,G,-1.10));
			nns.add(new NearestNeighborScore(U,G,U,A,-0.50));
			nns.add(new NearestNeighborScore(U,G,U,G,-0.40));
		}
		ArrayList<TerminalMismatchPairScore> tns = new ArrayList();
		{
			tns.add(new TerminalMismatchPairScore(A,U,A,A,-0.80));
			tns.add(new TerminalMismatchPairScore(A,U,A,C,-1.00));
			tns.add(new TerminalMismatchPairScore(A,U,A,G,-1.70));
			tns.add(new TerminalMismatchPairScore(A,U,A,U,-1.00));
			tns.add(new TerminalMismatchPairScore(A,U,C,A,-0.70));
			tns.add(new TerminalMismatchPairScore(A,U,C,C,-0.70));
			tns.add(new TerminalMismatchPairScore(A,U,C,G,-0.70));
			tns.add(new TerminalMismatchPairScore(A,U,C,U,-0.70));
			tns.add(new TerminalMismatchPairScore(A,U,G,A,-1.50));
			tns.add(new TerminalMismatchPairScore(A,U,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(A,U,G,G,-1.00));
			tns.add(new TerminalMismatchPairScore(A,U,G,U,-1.00));
			tns.add(new TerminalMismatchPairScore(A,U,U,A,-0.80));
			tns.add(new TerminalMismatchPairScore(A,U,U,C,-0.80));
			tns.add(new TerminalMismatchPairScore(A,U,U,G,-0.80));
			tns.add(new TerminalMismatchPairScore(A,U,U,U,-0.80));
			tns.add(new TerminalMismatchPairScore(C,G,A,A,-1.40));
			tns.add(new TerminalMismatchPairScore(C,G,A,C,-2.00));
			tns.add(new TerminalMismatchPairScore(C,G,A,G,-2.10));
			tns.add(new TerminalMismatchPairScore(C,G,A,U,-1.90));
			tns.add(new TerminalMismatchPairScore(C,G,C,A,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,C,C,-1.10));
			tns.add(new TerminalMismatchPairScore(C,G,C,G,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,C,U,-0.80));
			tns.add(new TerminalMismatchPairScore(C,G,G,A,-2.10));
			tns.add(new TerminalMismatchPairScore(C,G,G,C,-1.90));
			tns.add(new TerminalMismatchPairScore(C,G,G,G,-1.40));
			tns.add(new TerminalMismatchPairScore(C,G,G,U,-1.90));
			tns.add(new TerminalMismatchPairScore(C,G,U,A,-1.40));
			tns.add(new TerminalMismatchPairScore(C,G,U,C,-1.50));
			tns.add(new TerminalMismatchPairScore(C,G,U,G,-1.40));
			tns.add(new TerminalMismatchPairScore(C,G,U,U,-1.20));
			tns.add(new TerminalMismatchPairScore(G,C,A,A,-1.10));
			tns.add(new TerminalMismatchPairScore(G,C,A,C,-1.30));
			tns.add(new TerminalMismatchPairScore(G,C,A,G,-2.00));
			tns.add(new TerminalMismatchPairScore(G,C,A,U,-1.30));
			tns.add(new TerminalMismatchPairScore(G,U,A,A,-0.80));
			tns.add(new TerminalMismatchPairScore(G,U,A,C,-1.00));
			tns.add(new TerminalMismatchPairScore(G,U,A,G,-1.70));
			tns.add(new TerminalMismatchPairScore(G,U,A,U,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,C,A,-1.10));
			tns.add(new TerminalMismatchPairScore(G,C,C,C,-0.60));
			tns.add(new TerminalMismatchPairScore(G,C,C,G,-0.60));
			tns.add(new TerminalMismatchPairScore(G,C,C,U,-0.50));
			tns.add(new TerminalMismatchPairScore(G,U,C,A,-0.70));
			tns.add(new TerminalMismatchPairScore(G,U,C,C,-0.70));
			tns.add(new TerminalMismatchPairScore(G,U,C,G,-0.70));
			tns.add(new TerminalMismatchPairScore(G,U,C,U,-0.70));
			tns.add(new TerminalMismatchPairScore(G,C,G,A,-2.30));
			tns.add(new TerminalMismatchPairScore(G,C,G,C,-1.50));
			tns.add(new TerminalMismatchPairScore(G,C,G,G,-1.40));
			tns.add(new TerminalMismatchPairScore(G,C,G,U,-1.50));
			tns.add(new TerminalMismatchPairScore(G,U,G,A,-1.50));
			tns.add(new TerminalMismatchPairScore(G,U,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(G,U,G,G,-1.00));
			tns.add(new TerminalMismatchPairScore(G,U,G,U,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,U,A,-0.80));
			tns.add(new TerminalMismatchPairScore(G,C,U,C,-0.80));
			tns.add(new TerminalMismatchPairScore(G,C,U,G,-0.80));
			tns.add(new TerminalMismatchPairScore(G,C,U,U,-0.70));
			tns.add(new TerminalMismatchPairScore(G,U,U,A,-0.80));
			tns.add(new TerminalMismatchPairScore(G,U,U,C,-0.80));
			tns.add(new TerminalMismatchPairScore(G,U,U,G,-0.80));
			tns.add(new TerminalMismatchPairScore(G,U,U,U,-0.80));
			tns.add(new TerminalMismatchPairScore(U,A,A,A,-1.00));
			tns.add(new TerminalMismatchPairScore(U,A,A,C,-0.80));
			tns.add(new TerminalMismatchPairScore(U,A,A,G,-1.80));
			tns.add(new TerminalMismatchPairScore(U,A,A,U,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,A,A,-1.20));
			tns.add(new TerminalMismatchPairScore(U,G,A,C,-1.40));
			tns.add(new TerminalMismatchPairScore(U,G,A,G,-2.00));
			tns.add(new TerminalMismatchPairScore(U,G,A,U,-1.40));
			tns.add(new TerminalMismatchPairScore(U,A,C,A,-0.70));
			tns.add(new TerminalMismatchPairScore(U,A,C,C,-0.60));
			tns.add(new TerminalMismatchPairScore(U,A,C,G,-0.30));
			tns.add(new TerminalMismatchPairScore(U,A,C,U,-0.50));
			tns.add(new TerminalMismatchPairScore(U,G,C,A,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,C,C,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,C,G,-0.70));
			tns.add(new TerminalMismatchPairScore(U,G,C,U,-0.70));
			tns.add(new TerminalMismatchPairScore(U,A,G,A,-1.80));
			tns.add(new TerminalMismatchPairScore(U,A,G,C,-0.90));
			tns.add(new TerminalMismatchPairScore(U,A,G,G,-1.20));
			tns.add(new TerminalMismatchPairScore(U,A,G,U,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,G,A,-2.00));
			tns.add(new TerminalMismatchPairScore(U,G,G,C,-1.40));
			tns.add(new TerminalMismatchPairScore(U,G,G,G,-1.30));
			tns.add(new TerminalMismatchPairScore(U,G,G,U,-1.40));
			tns.add(new TerminalMismatchPairScore(U,A,U,A,-0.30));
			tns.add(new TerminalMismatchPairScore(U,A,U,C,-0.60));
			tns.add(new TerminalMismatchPairScore(U,A,U,G,-0.30));
			tns.add(new TerminalMismatchPairScore(U,A,U,U,-0.50));
			tns.add(new TerminalMismatchPairScore(U,G,U,A,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,U,C,-1.10));
			tns.add(new TerminalMismatchPairScore(U,G,U,G,-0.90));
			tns.add(new TerminalMismatchPairScore(U,G,U,U,-0.90));
		}
		ArrayList<DangleScore> dangles = new ArrayList();
		{
			dangles.add(new DangleScore(A,U,"3pr",A,-0.80));
			dangles.add(new DangleScore(A,U,"3pr",C,-0.50));
			dangles.add(new DangleScore(A,U,"3pr",G,-0.80));
			dangles.add(new DangleScore(A,U,"3pr",U,-0.60));
			dangles.add(new DangleScore(C,G,"3pr",A,-1.70));
			dangles.add(new DangleScore(C,G,"3pr",C,-0.80));
			dangles.add(new DangleScore(C,G,"3pr",G,-1.70));
			dangles.add(new DangleScore(C,G,"3pr",U,-1.20));
			dangles.add(new DangleScore(G,C,"3pr",A,-1.10));
			dangles.add(new DangleScore(G,C,"3pr",C,-0.40));
			dangles.add(new DangleScore(G,C,"3pr",G,-1.30));
			dangles.add(new DangleScore(G,C,"3pr",U,-0.60));
			dangles.add(new DangleScore(G,U,"3pr",A,-0.80));
			dangles.add(new DangleScore(G,U,"3pr",C,-0.50));
			dangles.add(new DangleScore(G,U,"3pr",G,-0.80));
			dangles.add(new DangleScore(G,U,"3pr",U,-0.60));
			dangles.add(new DangleScore(U,A,"3pr",A,-0.70));
			dangles.add(new DangleScore(U,A,"3pr",C,-0.10));
			dangles.add(new DangleScore(U,A,"3pr",G,-0.70));
			dangles.add(new DangleScore(U,A,"3pr",U,-0.10));
			dangles.add(new DangleScore(U,G,"3pr",A,-1.20));
			dangles.add(new DangleScore(U,G,"3pr",C,-0.50));
			dangles.add(new DangleScore(U,G,"3pr",G,-1.20));
			dangles.add(new DangleScore(U,G,"3pr",U,-0.70));
			dangles.add(new DangleScore(A,U,"5pr",A,-0.30));
			dangles.add(new DangleScore(A,U,"5pr",C,-0.10));
			dangles.add(new DangleScore(A,U,"5pr",G,-0.20));
			dangles.add(new DangleScore(A,U,"5pr",U,-0.20));
			dangles.add(new DangleScore(C,G,"5pr",A,-0.20));
			dangles.add(new DangleScore(C,G,"5pr",C,-0.30));
			dangles.add(new DangleScore(C,G,"5pr",G,-0.00));
			dangles.add(new DangleScore(C,G,"5pr",U,-0.00));
			dangles.add(new DangleScore(G,C,"5pr",A,-0.50));
			dangles.add(new DangleScore(G,C,"5pr",C,-0.30));
			dangles.add(new DangleScore(G,C,"5pr",G,-0.20));
			dangles.add(new DangleScore(G,C,"5pr",U,-0.10));
			dangles.add(new DangleScore(G,U,"5pr",A,-0.20));
			dangles.add(new DangleScore(G,U,"5pr",C,-0.20));
			dangles.add(new DangleScore(G,U,"5pr",G,-0.20));
			dangles.add(new DangleScore(G,U,"5pr",U,-0.20));
			dangles.add(new DangleScore(U,A,"5pr",A,-0.30));
			dangles.add(new DangleScore(U,A,"5pr",C,-0.30));
			dangles.add(new DangleScore(U,A,"5pr",G,-0.40));
			dangles.add(new DangleScore(U,A,"5pr",U,-0.20));
			dangles.add(new DangleScore(U,G,"5pr",A,-0.20));
			dangles.add(new DangleScore(U,G,"5pr",C,-0.20));
			dangles.add(new DangleScore(U,G,"5pr",G,-0.20));
			dangles.add(new DangleScore(U,G,"5pr",U,-0.20));		
		}
		
		reformScoreLists(nns,tns,dangles);
	}
	private void makeTable_DNA(){
		//****References:***//
		//Zucker's UNAFOLD parameters.
		//The were taken from the following papers
		//Santa Lucia 1998
		//Santa Lucia 1999
		//Sort-of Zucker, 1999 (JMB) <--- this is RNA scores, but may have been extrapolated.
	
		ArrayList<NearestNeighborScore> nns = new ArrayList();
		{
			//AUTOWRITTEN CODE.
			nns.add(new NearestNeighborScore(A,T,A,T,-1.00));
			nns.add(new NearestNeighborScore(A,T,C,G,-1.44));
			nns.add(new NearestNeighborScore(A,T,G,C,-1.28));
			nns.add(new NearestNeighborScore(A,T,G,T,+0.71));
			nns.add(new NearestNeighborScore(A,T,T,A,-0.88));
			nns.add(new NearestNeighborScore(A,T,T,G,+0.07));
			nns.add(new NearestNeighborScore(C,G,A,T,-1.45));
			nns.add(new NearestNeighborScore(C,G,C,G,-1.84));
			nns.add(new NearestNeighborScore(C,G,G,C,-2.17));
			nns.add(new NearestNeighborScore(C,G,G,T,-0.47));
			nns.add(new NearestNeighborScore(C,G,T,A,-1.28));
			nns.add(new NearestNeighborScore(C,G,T,G,-0.32));
			nns.add(new NearestNeighborScore(G,C,A,T,-1.30));
			nns.add(new NearestNeighborScore(G,T,A,T,+0.34));
			nns.add(new NearestNeighborScore(G,C,C,G,-2.24));
			nns.add(new NearestNeighborScore(G,T,C,G,-0.59));
			nns.add(new NearestNeighborScore(G,C,G,C,-1.84));
			nns.add(new NearestNeighborScore(G,C,G,T,+0.08));
			nns.add(new NearestNeighborScore(G,T,G,C,-0.32));
			nns.add(new NearestNeighborScore(G,T,G,T,+0.74));
			nns.add(new NearestNeighborScore(G,C,T,A,-1.44));
			nns.add(new NearestNeighborScore(G,C,T,G,-0.59));
			nns.add(new NearestNeighborScore(G,T,T,A,+0.07));
			nns.add(new NearestNeighborScore(G,T,T,G,+1.15));
			nns.add(new NearestNeighborScore(T,A,A,T,-0.58));
			nns.add(new NearestNeighborScore(T,G,A,T,+0.43));
			nns.add(new NearestNeighborScore(T,A,C,G,-1.30));
			nns.add(new NearestNeighborScore(T,G,C,G,+0.08));
			nns.add(new NearestNeighborScore(T,A,G,C,-1.45));
			nns.add(new NearestNeighborScore(T,A,G,T,+0.43));
			nns.add(new NearestNeighborScore(T,G,G,C,-0.47));
			nns.add(new NearestNeighborScore(T,G,G,T,+0.52));
			nns.add(new NearestNeighborScore(T,A,T,A,-1.00));
			nns.add(new NearestNeighborScore(T,A,T,G,+0.34));
			nns.add(new NearestNeighborScore(T,G,T,A,+0.71));
			nns.add(new NearestNeighborScore(T,G,T,G,+0.74));
		}
		ArrayList<TerminalMismatchPairScore> tns = new ArrayList();
		{
			//AUTOWRITTEN CODE
			tns.add(new TerminalMismatchPairScore(A,T,A,A,-0.70));
			tns.add(new TerminalMismatchPairScore(A,T,A,C,-0.30));
			tns.add(new TerminalMismatchPairScore(A,T,A,G,-0.50));
			tns.add(new TerminalMismatchPairScore(A,T,A,T,-0.90));
			tns.add(new TerminalMismatchPairScore(A,T,C,A,-0.60));
			tns.add(new TerminalMismatchPairScore(A,T,C,C,-0.20));
			tns.add(new TerminalMismatchPairScore(A,T,C,G,-1.00));
			tns.add(new TerminalMismatchPairScore(A,T,C,T,-0.30));
			tns.add(new TerminalMismatchPairScore(A,T,G,A,-0.60));
			tns.add(new TerminalMismatchPairScore(A,T,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(A,T,G,G,-0.40));
			tns.add(new TerminalMismatchPairScore(A,T,G,T,-0.50));
			tns.add(new TerminalMismatchPairScore(A,T,T,A,-0.78));
			tns.add(new TerminalMismatchPairScore(A,T,T,C,-0.30));
			tns.add(new TerminalMismatchPairScore(A,T,T,G,-0.50));
			tns.add(new TerminalMismatchPairScore(A,T,T,T,-0.40));
			tns.add(new TerminalMismatchPairScore(C,G,A,A,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,A,C,-0.80));
			tns.add(new TerminalMismatchPairScore(C,G,A,G,-0.90));
			tns.add(new TerminalMismatchPairScore(C,G,A,T,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,C,A,-0.80));
			tns.add(new TerminalMismatchPairScore(C,G,C,C,-0.50));
			tns.add(new TerminalMismatchPairScore(C,G,C,G,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,C,T,-0.70));
			tns.add(new TerminalMismatchPairScore(C,G,G,A,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,G,G,-0.90));
			tns.add(new TerminalMismatchPairScore(C,G,G,T,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,T,A,-1.00));
			tns.add(new TerminalMismatchPairScore(C,G,T,C,-0.60));
			tns.add(new TerminalMismatchPairScore(C,G,T,G,-0.90));
			tns.add(new TerminalMismatchPairScore(C,G,T,T,-0.90));
			tns.add(new TerminalMismatchPairScore(G,C,A,A,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,A,C,-0.70));
			tns.add(new TerminalMismatchPairScore(G,C,A,G,-0.80));
			tns.add(new TerminalMismatchPairScore(G,C,A,T,-1.00));
			tns.add(new TerminalMismatchPairScore(G,T,A,A,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,A,C,-0.20));
			tns.add(new TerminalMismatchPairScore(G,T,A,G,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,A,T,-0.60));
			tns.add(new TerminalMismatchPairScore(G,C,C,A,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,C,C,-0.60));
			tns.add(new TerminalMismatchPairScore(G,C,C,G,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,C,T,-0.70));
			tns.add(new TerminalMismatchPairScore(G,T,C,A,-0.20));
			tns.add(new TerminalMismatchPairScore(G,T,C,C,-0.20));
			tns.add(new TerminalMismatchPairScore(G,T,C,G,-0.90));
			tns.add(new TerminalMismatchPairScore(G,T,C,T,-0.20));
			tns.add(new TerminalMismatchPairScore(G,C,G,A,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,G,G,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,G,T,-0.80));
			tns.add(new TerminalMismatchPairScore(G,T,G,A,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,G,C,-0.90));
			tns.add(new TerminalMismatchPairScore(G,T,G,G,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,G,T,-0.50));
			tns.add(new TerminalMismatchPairScore(G,C,T,A,-1.00));
			tns.add(new TerminalMismatchPairScore(G,C,T,C,-0.60));
			tns.add(new TerminalMismatchPairScore(G,C,T,G,-0.90));
			tns.add(new TerminalMismatchPairScore(G,C,T,T,-0.90));
			tns.add(new TerminalMismatchPairScore(G,T,T,A,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,T,C,-0.20));
			tns.add(new TerminalMismatchPairScore(G,T,T,G,-0.50));
			tns.add(new TerminalMismatchPairScore(G,T,T,T,-0.20));
			tns.add(new TerminalMismatchPairScore(T,A,A,A,-0.60));
			tns.add(new TerminalMismatchPairScore(T,A,A,C,-0.40));
			tns.add(new TerminalMismatchPairScore(T,A,A,G,-0.50));
			tns.add(new TerminalMismatchPairScore(T,A,A,T,-0.60));
			tns.add(new TerminalMismatchPairScore(T,G,A,A,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,A,C,-0.20));
			tns.add(new TerminalMismatchPairScore(T,G,A,G,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,A,T,-0.50));
			tns.add(new TerminalMismatchPairScore(T,A,C,A,-0.50));
			tns.add(new TerminalMismatchPairScore(T,A,C,C,-0.20));
			tns.add(new TerminalMismatchPairScore(T,A,C,G,-1.00));
			tns.add(new TerminalMismatchPairScore(T,A,C,T,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,C,A,-0.20));
			tns.add(new TerminalMismatchPairScore(T,G,C,C,-0.20));
			tns.add(new TerminalMismatchPairScore(T,G,C,G,-0.80));
			tns.add(new TerminalMismatchPairScore(T,G,C,T,-0.20));
			tns.add(new TerminalMismatchPairScore(T,A,G,A,-0.60));
			tns.add(new TerminalMismatchPairScore(T,A,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(T,A,G,G,-0.40));
			tns.add(new TerminalMismatchPairScore(T,A,G,T,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,G,A,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,G,C,-1.00));
			tns.add(new TerminalMismatchPairScore(T,G,G,G,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,G,T,-0.50));
			tns.add(new TerminalMismatchPairScore(T,A,T,A,-0.80));
			tns.add(new TerminalMismatchPairScore(T,A,T,C,-0.30));
			tns.add(new TerminalMismatchPairScore(T,A,T,G,-0.60));
			tns.add(new TerminalMismatchPairScore(T,A,T,T,-0.30));
			tns.add(new TerminalMismatchPairScore(T,G,T,A,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,T,C,-0.20));
			tns.add(new TerminalMismatchPairScore(T,G,T,G,-0.50));
			tns.add(new TerminalMismatchPairScore(T,G,T,T,-0.20));
		}
		ArrayList<DangleScore> dangles = new ArrayList();
		{
			dangles.add(new DangleScore(A,T,"3pr",A,-0.12));
			dangles.add(new DangleScore(A,T,"3pr",C,0.28));
			dangles.add(new DangleScore(A,T,"3pr",G,-0.01));
			dangles.add(new DangleScore(A,T,"3pr",T,0.13));
			dangles.add(new DangleScore(C,G,"3pr",A,-0.82));
			dangles.add(new DangleScore(C,G,"3pr",C,-0.31));
			dangles.add(new DangleScore(C,G,"3pr",G,-0.01));
			dangles.add(new DangleScore(C,G,"3pr",T,-0.52));
			dangles.add(new DangleScore(G,C,"3pr",A,-0.92));
			dangles.add(new DangleScore(G,C,"3pr",C,-0.23));
			dangles.add(new DangleScore(G,C,"3pr",G,-0.44));
			dangles.add(new DangleScore(G,C,"3pr",T,-0.35));
			dangles.add(new DangleScore(G,T,"3pr",A,-0.20));
			dangles.add(new DangleScore(G,T,"3pr",C,-0.20));
			dangles.add(new DangleScore(G,T,"3pr",G,-0.20));
			dangles.add(new DangleScore(G,T,"3pr",T,-0.20));
			dangles.add(new DangleScore(T,A,"3pr",A,-0.48));
			dangles.add(new DangleScore(T,A,"3pr",C,-0.19));
			dangles.add(new DangleScore(T,A,"3pr",G,-0.50));
			dangles.add(new DangleScore(T,A,"3pr",T,-0.29));
			dangles.add(new DangleScore(T,G,"3pr",A,-0.20));
			dangles.add(new DangleScore(T,G,"3pr",C,-0.20));
			dangles.add(new DangleScore(T,G,"3pr",G,-0.20));
			dangles.add(new DangleScore(T,G,"3pr",T,-0.20));
			dangles.add(new DangleScore(A,T,"5pr",A,-0.50));
			dangles.add(new DangleScore(A,T,"5pr",C,-0.02));
			dangles.add(new DangleScore(A,T,"5pr",G,0.48));
			dangles.add(new DangleScore(A,T,"5pr",T,-0.10));
			dangles.add(new DangleScore(C,G,"5pr",A,-0.58));
			dangles.add(new DangleScore(C,G,"5pr",C,-0.34));
			dangles.add(new DangleScore(C,G,"5pr",G,-0.56));
			dangles.add(new DangleScore(C,G,"5pr",T,-0.61));
			dangles.add(new DangleScore(G,C,"5pr",A,-0.96));
			dangles.add(new DangleScore(G,C,"5pr",C,-0.52));
			dangles.add(new DangleScore(G,C,"5pr",G,-0.72));
			dangles.add(new DangleScore(G,C,"5pr",T,-0.58));
			dangles.add(new DangleScore(G,T,"5pr",A,-0.10));
			dangles.add(new DangleScore(G,T,"5pr",C,-0.10));
			dangles.add(new DangleScore(G,T,"5pr",G,-0.10));
			dangles.add(new DangleScore(G,T,"5pr",T,-0.10));
			dangles.add(new DangleScore(T,A,"5pr",A,-0.51));
			dangles.add(new DangleScore(T,A,"5pr",C,-0.42));
			dangles.add(new DangleScore(T,A,"5pr",G,-0.62));
			dangles.add(new DangleScore(T,A,"5pr",T,-0.71));
			dangles.add(new DangleScore(T,G,"5pr",A,-0.10));
			dangles.add(new DangleScore(T,G,"5pr",C,-0.10));
			dangles.add(new DangleScore(T,G,"5pr",G,-0.10));
			dangles.add(new DangleScore(T,G,"5pr",T,-0.10));
		}

		reformScoreLists(nns,tns,dangles);
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
}
