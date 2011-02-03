package DnaDesign.test;

import java.util.ArrayList;
import java.util.TreeMap;

import DnaDesign.DomainSequence;
import DnaDesign.Config.CircDesigNAConfig;
import DnaDesign.impl.CodonCode;

public class DomainDesignerTest1 {
	/*
	public static void main34(String[] args) throws IOException{
		DomainSequence ds = new DomainSequence();
		DomainSequence ds2 = new DomainSequence();
		double numSamples = 1e6;
		//String seq = "AAAAAAAAATTTTTTTTTTTT";
		File out = new File("C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\SSSscores\\xmer.txt");
		//File out = new File("C:\\Users\\Benjamin\\CLASSWORK\\002. UT UNDERGRADUATE GENERAL\\EllingtonLab\\AutoAmplifierDesign\\SSSscores\\xmer-pair.txt");
		PrintWriter o = new PrintWriter(new FileWriter(out));
		DomainDesigner_5_RandomDesigner2 ed = new DomainDesigner_5_RandomDesigner2();
		for(int i = 0; i < numSamples; i++){
			int merlen = (int)(Math.random()*1000) + 8;
			int[][] domain = new int[2][merlen];
			Arrays.fill(domain[0],1);
			int[][] domainMarkings = new int[2][merlen];
			deepFill(domainMarkings, DNAMARKER_DONTMUTATE);
			ds.setDomains(0);
			ds2.setDomains(1);
			for(int count = 0; count < 1; count++){
				for(int k = 0; k < domain[0].length; k++){
					domain[0][k] = (int)((Math.random()*4) + 1); 
				};
				for(int k = 0; k < domain[1].length; k++){
					domain[1][k] = (int)((Math.random()*4) + 1); 
				};
				//ed.displayDomains(domain);
				long now = System.nanoTime();
				double results_Matrix = ed.foldSingleStranded_viaMatrix(ds,domain,domainMarkings);
				//double results_Matrix = ed.pairscore_viaMatrix(ds, ds2, domain, domainMarkings);
				double results_Matrix_T = (System.nanoTime()-now)/1e9;
				now = System.nanoTime();
				double results_Unafold = ed.foldSingleStranded_viaUnafold(ds,domain,domainMarkings);
				//double results_Unafold = ed.pairscore_viaUnafold(ds, ds2, domain, domainMarkings);
				double results_Unafold_T = (System.nanoTime()-now)/1e9;
				
				o.println(ed.displaySequence(ds, domain)+" "+ed.displaySequence(ds2, domain)+" "+results_Matrix+" "+results_Unafold+" "+results_Matrix_T+" "+results_Unafold_T);
				o.flush();
			}
			System.out.println("merlen "+merlen);
		}
		o.close();
	}
	public static void main4(String[] args) throws Throwable{
		DomainSequence ds = new DomainSequence();
		String seq = "AGAACCGGCAACCAGACAGACAGAGTAACAGCCTAGAAATCAGCGCGCTCGGATCAGTGAATTCTCTTGAGCGTTGCAATACGACGTTGGCCAAGACACCCATTGTCTGTCCCTAGCCACCTATAGAGACTAAGAGGGAGAGATAATCTAACTGAGATACTACCCCAAAGGTCCAAAGTGAGACTTTTAAAAACGTGATGTGCCCGAATTTCTCCATGGGATTCGGCTTCAAAAAATGATACGTAGGCACGAATTGAGGAGCGCTAAGACCATTCTCTGAGACACGCTTCCGACACCAGACGCCCTACGGTGCGGTTACCCGGCGACCCTGTCATCATTGAGCCTACCCGGCCGTGACACATGCATCCCGCTATGGGCACGGATAGCTGACTATCCATTGTACTGCGTAAACGCAGAGACCGGAGAGTCCCGAGTTCCGGCTAGGCTGAAGGAGTCGGGTCAGAAACTCCGGATTTCCATGTTAAGCGAAACGCGTGCATGATGGCCCCTCAACCAACACCAACCTTGGGAATCTACAGTGACGTTTGCATTCTATGCTTACAATCTAAGGGTCTGATCGAAGTGCCGCCACCGTGTCCAAAGTATAGGTCTAGCAAGAACCGATCGCTAACCAGCCATTCAGGCCAGAAGACACCCGGTAGTGCTGCCCAGTCTACCACGGTGCTGGTCAGCGAGTTAAAATTGGCCTAAGAGAGGGCGACGGGGGAGCGATATGCGAGCATAACGTGAACGAAATCGTCGCGAGCGACAGACGCGCAGACCACGGTTGCGCCATGGGG";
		//String seq = "AAAAAAAAATTTTTTTTTTTT";
		int[][] domain = new int[1][seq.length()];
		int[][] domainMarkings = new int[1][seq.length()];
		int count = 0;
		DomainDesigner_5_RandomDesigner2 ed = new DomainDesigner_5_RandomDesigner2();
		for(char d : seq.toCharArray()){
			domain[0][count++]= ed.getLockedBase(d);
		}
		ds.setDomains(0);
		System.out.println(ed.foldSingleStranded(ds,domain,domainMarkings));
		System.out.println(Arrays.deepToString(domainMarkings));
	}
	*/
	public static void main42(String[] args) throws Throwable{
		TreeMap<Integer, String> lock = new TreeMap<Integer, String>();
		//lock.put(6,"GTTC");
		
		TreeMap<Integer, String> initial = new TreeMap<Integer, String>();

		String t7RnaPolymerase = "ATGGGGAGCTCGCATCACCATCACCATCACGGATCCAACACGATTAACATCGCTAAGAACGACTTCTCTGACATCGAACTGGCTGCTATCCCGTTCAACACTCTGGCTGACCATTACGGTGAGCGTTTAGCTCGCGAACAGTTGGCCCTTGAGCATGAGTCTTACGAGATGGGTGAAGCACGCTTCCGCAAGATGTTTGAGCGTCAACTTAAAGCTGGTGAGGTTGCGGATAACGCTGCCGCCAAGCCTCTCATCACTACCCTACTCCCTAAGATGATTGCACGCATCAACGACTGGTTTGAGGAAGTGAAAGCTAAGCGCGGCAAGCGCCCGACAGCCTTCCAGTTCCTGCAAGAAATCAAGCCGGAAGCCGTAGCGTACATCACCATTAAGACCACTCTGGCTTGCCTAACCAGTGCTGACAATACAACCGTTCAGGCTGTAGCAAGCGCAATCGGTCGGGCCATTGAGGACGAGGCTCGCTTCGGTCGTATCCGTGACCTTGAAGCTAAGCACTTCAAGAAAAACGTTGAGGAACAACTCAACAAGCGCGTAGGGCACGTCTACAAGAAAGCATTTATGCAAGTTGTCGAGGCTGACATGCTCTCTAAGGGTCTACTCGGTGGCGAGGCGTGGTCTTCGTGGCATAAGGAAGACTCTATTCATGTAGGAGTACGCTGCATCGAGATGCTCATTGAGTCAACCGGAATGGTTAGCTTACACCGCCAAAATGCTGGCGTAGTAGGTCAAGACTCTGAGACTATCGAACTCGCACCTGAATACGCTGAGGCTATCGCAACCCGTGCAGGTGCGCTGGCTGGCATCTCTCCGATGTTCCAACCTTGCGTAGTTCCTCCTAAGCCGTGGACTGGCATTACTGGTGGTGGCTATTGGGCTAACGGTCGTCGTCCTCTGGCGCTGGTGCGTACTCACAGTAAGAAAGCACTGATGCGCTACGAAGACGTTTACATGCCTGAGGTGTACAAAGCGATTAACATTGCGCAAAACACCGCATGGAAAATCAACAAGAAAGTCCTAGCGGTCGCCAACGTAATCACCAAGTGGAAGCATTGTCCGGTCGAGGACATCCCTGCGATTGAGCGTGAAGAACTCCCGATGAAACCGGAAGACATCGACATGAATCCTGAGGCTCTCACCGCGTGGAAACGTGCTGCCGCTGCTGTGTACCGCAAGGACAAGGCTCGCAAGTCTCGCCGTATCAGCCTTGAGTTCATGCTTGAGCAAGCCAATAAGTTTGCTAACCATAAGGCCATCTGGTTCCCTTACAACATGGACTGGCGCGGTCGTGTTTACGCTGTGTCAATGTTCAACCCGCAAGGTAACGATATGACCAAAGGACTGCTTACGCTGGCGAAAGGTAAACCAATCGGTAAGGAAGGTTACTACTGGCTGAAAATCCACGGTGCAAACTGTGCGGGTGTCGATAAGGTTCCGTTCCCTGAGCGCATCAAGTTCATTGAGGAAAACCACGAGAACATCATGGCTTGCGCTAAGTCTCCACTGGAGAACACTTGGTGGGCTGAGCAAGATTCTCCGTTCTGCTTCCTTGCGTTCTGCTTTGAGTACGCTGGGGTACAGCACCACGGCCTGAGCTATAACTGCTCCCTTCCGCTGGCGTTTGACGGGTCTTGCTCTGGCATCCAGCACTTCTCCGCGATGCTCCGAGATGAGGTAGGTGGTCGCGCGGTTAACTTGCTTCCTAGTGAAACCGTTCAGGACATCTACGGGATTGTTGCTAAGAAAGTCAACGAGATTCTACAAGCAGACGCAATCAATGGGACCGATAACGAAGTAGTTACCGTGACCGATGAGAACACTGGTGAAATCTCTGAGAAAGTCAAGCTGGGCACTAAGGCACTGGCTGGTCAATGGCTGGCTTACGGTGTTACTCGCAGTGTGACTAAGCGTTCAGTCATGACGCTGGCTTACGGGTCCAAAGAGTTCGGCTTCCGTCAACAAGTGCTGGAAGATACCATTCAGCCAGCTATTGATTCCGGCAAGGGTCTGATGTTCACTCAGCCGAATCAGGCTGCTGGATACATGGCTAAGCTGATTTGGGAATCTGTGAGCGTGACGGTGGTAGCTGCGGTTGAAGCAATGAACTGGCTTAAGTCTGCTGCTAAGCTGCTGGCTGCTGAGGTCAAAGATAAGAAGACTGGAGAGATTCTTCGCAAGCGTTGCGCTGTGCATTGGGTAACTCCTGATGGTTTCCCTGTGTGGCAGGAATACAAGAAGCCTATTCAGACGCGCTTGAACCTGATGTTCCTCGGTCAGTTCCGCTTACAGCCTACCATTAACACCAACAAAGATAGCGAGATTGATGCACACAAACAGGAGTCTGGTATCGCTCCTAACTTTGTACACAGCCAAGACGGTAGCCACCTTCGTAAGACTGTAGTGTGGGCACACGAGAAGTACGGAATCGAATCTTTTGCACTGATTCACGACTCCTTCGGTACCATTCCGGCTGACGCTGCGAACCTGTTCAAAGCAGTGCGCGAAACTATGGTTGACACATATGAGTCTTGTGATGTACTGGCTGATTTCTACGACCAGTTCGCTGACCAGTTGCACGAGTCTCAATTGGACAAAATGCCAGCACTTCCGGCTAAAGGTAACTTGAACCTCCGTGACATCTTAGAGTCGGACTTCGCGTTCGCGTAA"; 
		if (args.length > 0 && args[0].length()>0){
			t7RnaPolymerase = args[0].toUpperCase().replaceAll("\\s+","");
		}
		initial.put(0,t7RnaPolymerase);
		
		//Test
		CodonCode cc = null;
		cc = new CodonCode(new CircDesigNAConfig());
		
		/*
		int[] t7Dna = new int[t7RnaPolymerase.length()];
		for(int i = 0; i < t7Dna.length; i++){
			t7Dna[i] = getRegularBase(t7RnaPolymerase.charAt(i));
		}
		for(int i = 0; i < t7Dna.length; i += 3){
			cc.mutateToOther(t7Dna, i);
			System.out.print(cc.decodeAmino(t7Dna, i));
		}
		System.out.println();
		System.exit(0);
		*/
		
		
		ArrayList<DomainSequence> junctions = new ArrayList();
		//utilJunctionSplitter(junctions,"[1*|4*|7|5*|6*}");
		//utilJunctionSplitter(junctions,"[3*|2*|7*|1*}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[1|7|2|3|7*|1*|4*|7|3*|2*|7*}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[3|7*|4|1|7|3*|2*|7*|1*|4*|7|5*|6*}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[4|7*|5|6|7|4*|1*|7*|6*|5*|7}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions,"[6|7|1|4|7*|6*|5*|7|4*|1*|7*|2*|3*}");

		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions, "[1<8.>|2<8(>|3<8(>|4<11.>|3*<8)>|2*<8)>|5<8>}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions, "[3<8.>|4*<11(>|3*<8.>|2*<8.>|4<11)>}");
		//DomainDesigner_SharedUtils.utilJunctionSplitter(junctions, "[3*<8.>|2*<8.>|1*<8.>}");
		//new DomainDesigner_2_RandomFreshDomains().main(5,new int[]{8,8,8,11,8,8,8},Integer.MAX_VALUE,lock,junctions);

		ArrayList<String> strands = new ArrayList();
		
		strands.add("[1<"+t7RnaPolymerase.length()+".>");
		
		/*
		strands.add("[1<8.>|2<8(>|3<8(>|4<11.>|3*<8)>|2*<8)>|5<8>}");
		strands.add("[3<8.>|4*<11(>|3*<8.>|2*<8.>|4<11)>}");
		strands.add("[3*<8.>|2*<8.>|1*<8.>}");
		strands.add("[3*<8(>|2*<8(>|1*<8(>}[1<8)>|2<8)>|3<8)>|4<11.>|3*<8.>|2*<8.>|5<8.>}");
		 */
		
		
		/*
		strands.add("[1<8.>|2<8(>|3<8(>|4<11.>|3*<8)>|2*<8)>|5<8>}");
		strands.add("[3<8.>|4*<11(>|3*<8.>|2*<8.>|4<11)>}");
		strands.add("[3*<8.>|2*<8.>|1*<8.>}");
		strands.add("[1<8>|2<8(>|3<8(>|4<11(>|3*<(>|2*|5<8>}[3<)>|4*<)>|3*<)>|2*<)>|4}");
		strands.add("[3*<8(>|2*<8(>|1*<8(>}[1<8)>|2<8)>|3<8)>|4<11.>|3*<8.>|2*<8.>|5<8.>}");
		
		strands.add("[6<8.>|7<8(>|8<8(>|9<11.>|8*<8)>|7*<8)>|10<8>}");
		strands.add("[8<8.>|9*<11(>|8*<8.>|7*<8.>|9<11)>}");
		strands.add("[8*<8.>|7*<8.>|6*<8.>}");
		strands.add("[6<8>|7<8(>|8<8(>|9<11(>|8*<(>|7*|10<8>}[8<)>|9*<)>|8*<)>|7*<)>|9}");
		strands.add("[8*<8(>|7*<8(>|6*<8(>}[6<8)>|7<8)>|8<8)>|9<11.>|8*<8.>|7*<8.>|10<8.>}");
		*/
		
		//strands.add("[1<80.>|2*<80.>}");
		
		//strands.add("[1<8>|7<4(>|2<8(>|3<8(>|7*<(>|1*|4*<8>|7<)>|3*<)>|2*<)>|7*<)>");
		//strands.add("[3<8>|7*<4(>|4<8(>|1<8(>|7<(>|3*|2*<8>|7*<)>|1*<)>|4*<)>|7<)>|5*<8>|6*<8>");
		//strands.add("[4<8>|7*<4(>|5<8(>|6<8(>|7<(>|4*|1*<8>|7*<)>|6*<)>|5*<)>|7<)>");
		//strands.add("[6<8>|7<4(>|1<8(>|4<8(>|7*<(>|6*|5*<8>|7<)>|4*<)>|1*<)>|7*<)>|2*<8>|3*<8>");
		
		//TODO
		/*
		DDSeqDesigner dsd = makeDesigner(strands, lock, initial, cc);

		dsd.resume();
		Scanner in = new Scanner(System.in);
		in.nextLine();
		in.close();
		Thread.sleep(10);
		System.out.println(dsd.getResult());

		System.exit(0);
		*/
	}
}
