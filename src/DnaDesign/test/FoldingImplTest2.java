package DnaDesign.test;

import static DnaDesign.AbstractPolymer.DnaDefinition.A;
import static DnaDesign.AbstractPolymer.DnaDefinition.C;
import static DnaDesign.AbstractPolymer.DnaDefinition.G;
import static DnaDesign.AbstractPolymer.DnaDefinition.T;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane.RestoreAction;

import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureData;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.FoldingImpl;
public class FoldingImplTest2 {
	/**
	 * Tests base pairing observables calculations
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException{
		CircDesigNAConfig config = new CircDesigNAConfig();
		System.out.println("L1 L2 S1 S2 CA NU MF");
		for(int i = 0; i < 4000; i++){
			FoldingImpl fl = new FoldingImpl(config);
			DomainStructureData dsd = new DomainStructureData(config);
			String seq = null;// "AAAAATTTTTTGGGGGGCCCCCCCC";
			int seqLength;
			if (seq==null){
				seqLength = (int)(Math.random()*40+1);
			} else {
				seqLength = seq.length();
			}
			int max_did = 1;
			int[][] domain = new int[max_did+1][seqLength];
			int[][] domainMark= new int[max_did+1][seqLength];
			for(int did = 0; did <= max_did; did++){
				if (seq!=null){
					for(int k = 0; k < seqLength; k++){
						domain[did][k] = config.monomer.decodeConstraintChar(seq.charAt(k));
					}
				} else {
					for(int k = 0; k < seqLength; k++){
						domain[did][k] = randomChoice(A,C,G,T);
					}
				}
			}

			/*
			long now;
			double[][] basePairs;
			if (true){
				int N = seqLength + seqLength;
				basePairs = new double[N][N+1];
				now = System.nanoTime();
				fl.pairPrHybrid(basePairs, seqS, seqS, domain);
			} else {
				int N = seqLength;
				basePairs = new double[N][N+1];
				now = System.nanoTime();
				fl.pairPrSS(basePairs, seqS, domain);
			}
			long time = System.nanoTime()-now;
			System.out.println(seqLength+" "+(time/1e9));
			*/
			DomainSequence ds1 = new DomainSequence();
			ds1.setDomains(0, null);
			DomainSequence ds2 = new DomainSequence();
			ds2.setDomains(1, null);
			System.out.printf("%d %d ",seqLength,seqLength);
			double resultSelf = fl.mfeHybridDeltaG_viaMatrix(ds1, ds2, domain, domainMark);
			double resultNupack;
			{
				String prefix = "nupack0";
				StringBuffer seqs = new StringBuffer();
				seqs.append(max_did+1+"\n");
				for(int did = 0; did <= max_did; did++){
					for(int j = 0; j < seqLength; j++){
						seqs.append(dsd.Std.monomer.displayBase(domain[did][j]));
					}
					seqs.append("\n");
				}
				StringBuffer concs = new StringBuffer();
				File nupackDir = new File("nupackTest/");
				nupackDir.mkdir();
				RunNupackTool.runNupack(seqs.toString(), concs.toString(), 2, prefix, false, nupackDir);
				
				Scanner cxin = new Scanner(new File("nupackTest/"+prefix+".cx"));
				resultNupack = 0;
				while(cxin.hasNextLine()){
					String string = cxin.nextLine();
					String[] line = string.split("\\s+");
					if (line[0].startsWith("%")){
						continue;
					}
					if (line[1].equals("1") && line[2].equals("1")){
						resultNupack = new Double(line[3]);
					}
				}
			}
			for(int did = 0; did <= max_did; did++){
				for(int j = 0; j < seqLength; j++){
					System.out.print(dsd.Std.monomer.displayBase(domain[did][j]));
				}
				System.out.print(" ");
			}
			System.out.printf("%f %f",resultSelf, resultNupack);
			System.out.println();
		}
	}

	private static int randomChoice(int a, int b, int c, int d) {
		switch((int)(Math.random()*4)){
		case 0:return a;
		case 1:return b;
		case 2:return c;
		case 3:return d;
		}
		throw new RuntimeException("Assertion error: randomChoice");
	}
}
