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
package circdesigna.batch;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.TreeMap;


import circdesigna.AbstractDomainDesignTarget;
import circdesigna.CircDesigNAOptions;
import circdesigna.DesignIntermediateReporter;
import circdesigna.DomainDefinitions;
import circdesigna.CircDesigNA;
import circdesigna.CircDesigNA_SharedUtils;
import circdesigna.DomainSequence;
import circdesigna.SequenceDesigner;
import circdesigna.CircDesigNA.ScorePenalty;
import circdesigna.config.CircDesigNAConfig;
import circdesigna.energy.CircDesigNAMCSFolder;
import circdesigna.energy.NAFolding;
import circdesigna.impl.CircDesigNAImpl;
import circdesigna.impl.SequencePenaltiesImpl;
import circdesigna.impl.CircDesigNAImpl.HairpinOpening;
import circdesigna.impl.CircDesigNAImpl.MFEHybridScore;
import circdesigna.impl.CircDesigNAImpl.SelfFold;
import circdesigna.impl.CircDesigNAImpl.SelfSimilarityScore;
import circdesigna.impl.CircDesigNAImpl.VariousSequencePenalties;
import circdesigna.plugins.RunNupackTool;



public class DesignMultipleTimes {
	/**
	 * Input is a path to a file of the form
	 * NumTimesToRun
	 * TargetDirectory
	 * Prefix
	 * (Domain definitions)
	 * END
	 * (Molecule definitions)
	 * END
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException{
		Scanner in = new Scanner(System.in);
		int numTimesToRun = new Integer(readToken(in));
		int popSize = new Integer(readToken(in));
		int selfSimilarOpt = new Integer(readToken(in));
		double numSecondsEach = new Double(readToken(in));
		int maximumComplexSize = new Integer(readToken(in));
		double tournamentInterval = new Double(readToken(in));
		boolean runGA = new Integer(readToken(in))==1;
		String targetDir = in.nextLine();
		String prefix = in.nextLine();
		String Domains = readToEnd(in);
		String Molecules = readToEnd(in);
		in.close();
		
		CircDesigNAConfig config = new CircDesigNAConfig();

		
		DomainDefinitions dsd = new DomainDefinitions(config);
		DomainDefinitions.readDomainDefs(Domains, dsd);
		if (args.length==0){
			new File(targetDir).mkdir();
			System.out.println("Entering design mode");
			for(int i = 1; i <= numTimesToRun; i++){
				SequenceDesigner<CircDesigNAOptions> defaultDesigner = CircDesigNA.getDefaultDesigner(Molecules, Domains, config);
				CircDesigNAOptions options = defaultDesigner.getOptions();
				options.rule_ccend_option.setState(false);
				options.population_size.setState(popSize);
				options.selfSimilarityPenalty.setState(selfSimilarOpt);
				options.resourcePerMember.setState(tournamentInterval);
				options.standardUseGA.setState(runGA);

				File out2 = new File(targetDir+File.separator+prefix+i+".des");
				System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(out2))));
				designUsingCircDesigNA(defaultDesigner,numSecondsEach);
				System.err.println("Designed "+i);
				System.err.flush();
				System.out.close();
			}
			File out2 = new File(targetDir+File.separator+prefix+"rnd"+1+".des");
			System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(out2))));
			outputRandomDesigns(dsd);
			System.out.close();
			
			out2 = new File(targetDir+File.separator+prefix+"ra"+1+".des");
			System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(out2))));
			outputRationalDesign(dsd);
			System.out.close();
		} else {
			System.out.println("Entering Evaluate mode");
			if (args[0].equals("-evaluate")){
				boolean goNupack = true;
				if (args.length >= 2){
					if (args[1].equals("-nopack")){
						goNupack = false;
					}
				}
				for(int i = 1; i <= numTimesToRun; i++){
					for(String rprefix : new String[]{prefix,prefix+"rnd",prefix+"ra"}){
						File out2 = new File(targetDir+File.separator+rprefix+i+".des");
						if (!out2.exists()){
							System.err.println("Could not find file "+rprefix+i+".des");
							continue;
						}
						File out3 = new File(targetDir+File.separator+rprefix+i+".eval");
						System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(out3))));
						RunEvaluation(out2,Molecules,dsd,maximumComplexSize, goNupack);
						System.err.println("Evaluated "+i);
						System.out.close();
					}
				}
			} else if (args[0].equals("-evaluate_self")){
				CircDesigNAOptions options = CircDesigNAOptions.getDefaultOptions();
				options.selfSimilarityPenalty.setState(selfSimilarOpt);
				
				for(int i = 1; i <= numTimesToRun; i++){
					for(String rprefix : new String[]{prefix,prefix+"rnd",prefix+"ra"}){
						File out2 = new File(targetDir+File.separator+rprefix+i+".des");
						if (!out2.exists()){
							System.err.println("Could not find file "+rprefix+i+".des");
							continue;
						}
						File out3 = new File(targetDir+File.separator+rprefix+i+".seval");
						System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream(out3))));
						RunSelfEvaluation(out2,Molecules,dsd,options);
						System.out.close();
					}
				}
			}
		}
	}

	private static void RunSelfEvaluation(File out2, String molecules, DomainDefinitions dsd, CircDesigNAOptions options) throws FileNotFoundException {
		CircDesigNAConfig config = new CircDesigNAConfig();
		AbstractDomainDesignTarget target = new AbstractDomainDesignTarget(dsd, config);
		{
			for(String k : molecules.split("\n")){ //newlines have been sanitized
				target.addTargetStructure(k);
			}
		}
		
		CircDesigNAImpl ddi = new CircDesigNAImpl(new CircDesigNAMCSFolder(config),new SequencePenaltiesImpl(config),config);
		DesignIntermediateReporter DIR = new DesignIntermediateReporter();
		
		//generate domain from dsd.
		int[] domainLengths = dsd.domainLengths;
		int[][] domain = new int[domainLengths.length][];
		int[][] nulldomainMark = new int[domainLengths.length][];
		for(int i = 0; i < domain.length; i++){
			domain[i] = new int[domainLengths[i]];
			nulldomainMark[i] = new int[domainLengths[i]];
		}
		
		List<ScorePenalty> listPenalties = ddi.listPenalties(target, DIR, domain, options, dsd);
		
		//evaluate scores for all designed.
		Scanner in = new Scanner(out2);
		String[] Headers = new String[]{"ITR","Hybrid","SelfSimilar","HairpinOpen","SelfFold","Other","TOTAL"};
		System.out.printf("%-20s%-20s%-20s%-20s%-20s%-20s%-20s",(Object[])Headers);
		System.out.println();
		while(in.hasNextLine()){
			String line = in.nextLine();
			System.out.flush();
			if (line.startsWith("Iteration")){
				System.out.printf("%-20s",line.split("\\s+")[1]);

				GetSelfEvaluation(in,dsd,listPenalties,domain,nulldomainMark,molecules);
			}
		}
		in.close();
	}

	private static final int GS_MFEHyb = 0, GS_SS = GS_MFEHyb+1, GS_HP = GS_SS+1, GS_SF = GS_HP+1, GS_OTHER=GS_SF+1, GS_NUM=GS_OTHER+1;
	
	private static void GetSelfEvaluation(Scanner in, DomainDefinitions dsd,
			List<ScorePenalty> listPenalties, int[][] domain,
			int[][] nulldomainMark, String molecules) {
		
		parseDnaDesignOutput(dsd, in, molecules, domain);
		
		double[] score = new double[GS_NUM];
		double total = 0;
		for(ScorePenalty k : listPenalties){
			int num = -1;
			if (k instanceof MFEHybridScore){
				num = GS_MFEHyb;
			}
			if (k instanceof HairpinOpening){
				num = GS_HP;
			}
			if (k instanceof SelfSimilarityScore){
				num = GS_SS;
			}
			if (k instanceof SelfFold){
				num = GS_SF;
			}
			if (k instanceof VariousSequencePenalties){
				num = GS_OTHER;
			}
			if (num==-1){
				System.err.println(k.getClass());
			}
			double sc = k.getScore(domain, nulldomainMark);
			score[num] += sc;
			total += sc;
		}
		for(int k = 0; k < score.length; k++){
			System.out.printf("%-20.3f",score[k]);
		}
		System.out.printf("%-20.3f",total);
		System.out.println();
		
	}

	private static void outputRandomDesigns(DomainDefinitions dsd) {
		for(int itr = 1; itr <= 100; itr++){
			System.out.println("Iteration "+itr+" Score 0");
			for(int i = 0; i < dsd.domainLengths.length; i++){
				int len = dsd.domainLengths[i];
				System.out.println(dsd.getDomainName(i)+">");
				for(int j = 0; j < len; j++){
					System.out.print(dsd.Std.monomer.displayBase(1+(int)(Math.random()*4)));
				}
				System.out.println();
			}
			System.out.println();
		}
	}
	
	private static final String[] shortReedWords2005 = new String[]{
	     "CACATCACCAATATA", "TCCTCCAATTAATTA", "TCAATCTTTTTCAAT", "TATCTTTCCAATCTA", "TCATTCTTCTCTTAT", "ACTACAATCTCAATA", "TCCTAACCAAAAATT", "TTTACCTTTTCAAAT", "CAACTACATTTTCTA", "TCATCTTATCTCTCT", "TCATTAAATCCATCT", "CATAATACCTTCCTA", "ACTACTTACATTTCA", "TCTTCAATCTACTTA", "TTTACCTCTCTAATC", "ATCTCTTTTCTTTTC", "CTCTTTCACAAAAAT", "AACAAAACAACAAAT", "ACAAAATCTCTTACA", "ACCATTTTTTCACTA", "CTTAATTCTCTCACT", "CTCTCCTCAAATATT", "ATCATCCACATATTT", "TTCCACACTAAAATT", "ACCATCATTTATCTT", "TCTCCTTATTCATTT", "TACCATACCAATTTT", "AAACCACCATAATTA", "TCCTAATCCTCTTAA", "TCCTACTCATAACTA", "CTCAACTTTCAAATT", "AACTTTACTATCCAT", "ATTCACCAAACTTTA", "TATACCATCTTTCAA", "TCCACTAATAAAACT", "CAATAATCCAACATT", "CTCTAATTTCTTCCT", "CACCTACATCAAATA", "AACAAACCATAACTA", "TCTACCTATTCACTA",  
	};

	private static void outputRationalDesign(DomainDefinitions dsd) {
		for(int i = 0; i < dsd.domainLengths.length; i++){
			int len = dsd.domainLengths[i];
			if (len > 15){
				System.err.println("Aborted rational design");
				return;
			}
		}
		System.out.println("Iteration 1 Score 0");
		for(int i = 0; i < dsd.domainLengths.length; i++){
			int len = dsd.domainLengths[i];
			System.out.println(dsd.getDomainName(i)+">");
			//Output a domain of that length
			System.out.println(shortReedWords2005[i].substring(0,len));
		}
		System.out.println();
	}
	
	private static void designUsingCircDesigNA(
			SequenceDesigner<CircDesigNAOptions> defaultDesigner, double numSecondsEach) {
		long diff = 0;
		defaultDesigner.resume();
		long now = System.nanoTime();
		while(!defaultDesigner.isFinished()){
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			if ((diff=System.nanoTime()-now) > numSecondsEach*1e9){
				break;
			}
		}
		defaultDesigner.abort();
		//System.err.print("HE");
		while(!defaultDesigner.isFinished()){
			try {
				Thread.sleep(10);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		//System.err.println("Y!");
		System.out.println("END. TIME:");
		System.out.println(String.format("%.3f",diff/1e9));
		
	}

	public static void RunEvaluation(File out3, String molecules, DomainDefinitions dsd, int maximumComplexSize, boolean goNupack) throws IOException {
		Scanner in = new Scanner(out3);
		System.out.printf("%-20s%-20s","ITR","SCORE");
		if (goNupack){
			System.out.printf("%-20s","NUPACK");
		}
		System.out.println();
		ArrayList lastMoleculeParse = new ArrayList();
		while(in.hasNextLine()){
			String line = in.nextLine();
			System.out.flush();
			if (line.startsWith("Iteration")){
				int i = new Integer(line.split("\\s+")[1]);
				double score = new Double(line.split("\\s+")[3]);
				System.out.printf("%-20d%-20.3f",i,score);
				if (goNupack){
					double nupack = getDefectScore(in,dsd,molecules,maximumComplexSize, lastMoleculeParse);
					System.out.printf("%-20e",nupack);
				}
				System.out.println();
			}	
		}
		in.close();
	}
	

	private static double getDefectScore(Scanner in, DomainDefinitions dsd, String molecules, int maximumComplexSize, ArrayList<TestMolecule> lastMoleculeParse) throws IOException {
		int[][] domains = new int[dsd.domainLengths.length][];
		for(int i = 0; i < domains.length; i++){
			domains[i] = new int[dsd.domainLengths[i]];
		}
		
		ArrayList<TestMolecule> moleculeParse = parseDnaDesignOutput(dsd,in,molecules,domains);
		
		boolean equal = true;
		for(int i = 0; i < lastMoleculeParse.size(); i++){
			if (!lastMoleculeParse.get(i).equals(moleculeParse.get(i))){
				equal = false;
				break;
			}
		}
		if(lastMoleculeParse.size()>0){
			if (equal){
				return -1;
			}
		}
		lastMoleculeParse.clear();
		lastMoleculeParse.addAll(moleculeParse);
		
		
		//System.out.println(Arrays.deepToString(domains));

		String prefix = "nupack"+System.nanoTime();
		//Nupack score

		StringBuffer molecules1 = new StringBuffer();
		molecules1.append(moleculeParse.size()+"\n");
		for(TestMolecule k : moleculeParse){
			molecules1.append(k.seq+"\n");
		}
		StringBuffer concs = new StringBuffer();
		for(int i = 0; i < moleculeParse.size(); i++){
			concs.append("5e-6"+"\n");
		}
		
		File nupackDir = new File("nupackTest/");
		
		RunNupackTool.runNupack(molecules1.toString(), concs.toString(), maximumComplexSize, prefix, RunNupackTool.OPCODE_COMPLEXES_AND_CONC, nupackDir);
		
		//Ok, nupack ran. Analyze!
		return calculateDefectScoreFromNupackOutput(moleculeParse, dsd,"nupackTest/"+prefix);
	}
	private static ArrayList<TestMolecule> parseDnaDesignOutput(DomainDefinitions dsd, Scanner in, String molecules, int[][] domains) {
			String[] molecule = new String[1];
			StringBuffer seq = new StringBuffer();
			while(in.hasNextLine()){
				String line2 = in.nextLine();
				if (line2.trim().equals("")){
					break;
				}
				if (line2.contains(">")){
					addMolecule(molecule,dsd,seq,domains);
					seq = new StringBuffer();
					molecule[0] = line2.substring(0,line2.length()-1);
				} else {
					//Sequence containing line
					seq.append(line2);
				}
			}
			addMolecule(molecule,dsd,seq,domains);

			//Ok, put molecules together
			ArrayList<DomainSequence> alreadyPrintedSequences = new ArrayList();
			ArrayList<TestMolecule> moleculeParse = new ArrayList();
			for(String q : molecules.split("\n")){
				String[] a = q.split("\\s+",2);
				int ct = 0;
				splitLoop: for(String subStrand : a[1].split("}")){
					DomainSequence ds = new DomainSequence();
					ds.setDomains(subStrand,dsd,null);
					for(DomainSequence g : alreadyPrintedSequences){
						if (g.equals(ds)){
							continue splitLoop;
						}
					}
					alreadyPrintedSequences.add(ds);

					TestMolecule mol = new TestMolecule();
					mol.ID = a[0]+ct++;
					mol.molecule = subStrand;
					seq = new StringBuffer();
					int len = ds.length(domains);
					for(int i = 0; i < len; i++){
						seq.append(dsd.Std.monomer.displayBase(ds.base(i, domains, dsd.Std.monomer)));
					}
					mol.seq = seq.toString();
					moleculeParse.add(mol);
				}
			}
			return moleculeParse;
	}
	

	private static double calculateDefectScoreFromNupackOutput(
			ArrayList<TestMolecule> moleculeParse, DomainDefinitions dsd,
			String filePrefix) throws FileNotFoundException {
		LinkedList<TestComplex> complex = new LinkedList();
		//Get concentrations
		Scanner inConc = new Scanner(new File(filePrefix+".eq"));
		TreeMap<Integer, Double> concStructures = new TreeMap();
		while(inConc.hasNextLine()){
			String line2 = inConc.nextLine();
			if (line2.startsWith("%") || line2.trim().equals("")){
				continue;
			}
			String[] line = line2.split("\\s+");
			concStructures.put(new Integer(line[0]),new Double(line[line.length-1]));
		}
		//System.err.println(concStructures);
		if (false){
			Scanner inkey = new Scanner(new File(filePrefix+".ocx-key"));

			while(inkey.hasNextLine()){
				String line = inkey.nextLine();
				if (line.startsWith("%")){
					continue;
				}
				String[] line2 = line.split("\\s+");
				StringBuffer ds = new StringBuffer();
				//ds.append("[");
				
				for(int i = 2; i< line2.length; i++){
					ds.append(moleculeParse.get(new Integer(line2[i])-1).molecule);
					ds.append("|");
				}
				 
				/*
				for(int i = 0; i < moleculeParse.size(); i++){
					ds.append(moleculeParse.get(i).molecule);
					ds.append("|");
				}
				*/
				//ds.append("}");
				complex.add(new TestComplex(ds.toString(), dsd, line2[0]));
			}
			inkey.close();
		}
		Scanner inPpairs = new Scanner(new File(filePrefix+".cx-epairs"));
		//Generate pseudo domains from dsd
		int[][] domain = new int[dsd.domainLengths.length][];
		for(int i = 0; i < domain.length; i++){
			domain[i] = new int[dsd.domainLengths[i]];
		}
		double TOTALSCORE = 0;
		//Make a single giant molecule for all strands
		StringBuffer dsstr = new StringBuffer();
		for(int i = 0; i < moleculeParse.size(); i++){
			dsstr.append(moleculeParse.get(i).molecule);
			dsstr.append("|");
		}
		TestComplex toAnalyze = new TestComplex(dsstr.toString(),dsd,"1");
		int concID = 0;
		while(inPpairs.hasNextLine()){
			String line = inPpairs.nextLine();
			if (line.startsWith("%")){
				continue;
			}
			if (line.trim().equals("")){
				continue;
			}
			concID++;
			//int num = new Integer(line);
			//double fe = new Double(inPpairs.nextLine());
			if (false){
				toAnalyze = complex.removeFirst();
			}
			String struct = inPpairs.nextLine();
			double structureCont = 0;
			for(; inPpairs.hasNextLine(); ){
				String nextLine = inPpairs.nextLine();
				if (nextLine.startsWith("%")){
					break;
				}
				String[] line2 = nextLine.split("\\s+");
				int i = new Integer(line2[0])-1;
				int j = new Integer(line2[1])-1;
				double p = new Double(line2[2]);
				DomainSequence ds = toAnalyze.ds;
				int len = ds.length(domain);
				if (j < len){
					int domain1 = ds.domainAt(i, domain);
					int domain2 = ds.domainAt(j, domain);
						
					int offset1 = ds.offsetInto(i, domain, false); //Don't uncomplement - give offset 
					int offset2 = ds.offsetInto(j, domain, false);
					
					if (!CircDesigNA_SharedUtils.isAlignedAndShouldPair(ds, i, ds, j, domain)){
						//System.err.print("Shouldn't pair: ");
						//TOTALSCORE += 1;
						structureCont += p;
					} else {
						//System.err.print("Did pair: ");
					}
					//System.err.println(dsd.getDomainName(domain1)+" "+offset1+" "+dsd.getDomainName(domain2)+" "+offset2);
				}
			}
			//System.err.println(structureCont+" "+concStructures.get(concID));
			double conc = 0;
			if (concStructures.containsKey(concID)){
				conc = concStructures.get(concID);
			}
			TOTALSCORE += structureCont * conc;
			
		}
		inPpairs.close();
		return TOTALSCORE;
	}

	private static void addMolecule(String[] molecule, DomainDefinitions dsd, StringBuffer seq, int[][] domains) {
		if (molecule[0]==null){
			return;
		}
		int domainNum = dsd.lookupDomainName(molecule[0]);
		for(int k = 0; k < seq.length(); k++){
			domains[domainNum][k] = dsd.Std.monomer.decodeBaseChar(seq.charAt(k));
		}
		/*
		TestMolecule newM = new TestMolecule();
		newM.ID = ""+domains.size();
		newM.molecule = molecule[0];
		newM.seq = seq.toString();
		//domains.add(newM);
		 */
		molecule[0] = null;
	}
	private static class TestComplex{
		public TestComplex(String domains, DomainDefinitions dsd, String name){
			ds = new DomainSequence();
			ds.setDomains(domains, dsd, null);
			this.name = name;
		}
		public String name;
		private DomainSequence ds;
	}
	private static class TestMolecule{
		public String seq;
		public String molecule;
		public String ID;
		public boolean equals(Object o){
			TestMolecule b = (TestMolecule)o;
			return b.seq.equals(seq) && molecule.equals(b.molecule);
		}
	}

	public static String readToken(Scanner in) {
		String line = in.nextLine();
		return line.split("\\s+")[0];
	}

	public static String readToEnd(Scanner in) {
		StringBuffer toRet = new StringBuffer();
		while(in.hasNextLine()){
			String line = in.nextLine();
			if (line.equals("END")){
				break;
			}
			toRet.append(line+"\n");
		}
		return toRet.toString();
	}
}
