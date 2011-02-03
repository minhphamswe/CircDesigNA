package DnaDesign.test;

import java.util.ArrayList;
import java.util.Scanner;

import DnaDesign.DomainDesigner_SharedUtils;
import DnaDesign.DomainPolymerGraph;
import DnaDesign.DomainSequence;
import DnaDesign.DomainStructureBNFTree;
import DnaDesign.DomainStructureData;
import DnaDesign.Config.CircDesigNAConfig;

/**
 * Checks that the methods in DomainDesigner_SharedUtils are equally good at handling
 * the BNFTree data structure and the PolymerGraph datastructure.
 */
public class MoleculeRepresentationTest {
	public static void main(String[] args){		

		CircDesigNAConfig config = new CircDesigNAConfig();
		Scanner in = new Scanner(System.in);
		System.out.println("Enter in the domain defs, END when finished");
		StringBuffer domainDefs = new StringBuffer();
		while(in.hasNextLine()){
			String line = in.nextLine();
			if (line.equals("END")){
				break;
			}
			domainDefs.append(line);
			domainDefs.append("\n");
		}
		String domainDefsBlock = domainDefs.toString();
		DomainStructureData dsd = new DomainStructureData(config);
		DomainStructureData.readDomainDefs(domainDefsBlock, dsd);

		System.out.println("Enter in a molecule");
		String mol;
		while(true){
			mol = in.nextLine().trim();
			if (mol.length()!=0){
				break;
			}
		}
		String[] mol2 = mol.split("\\s+");
		mol = mol2[mol2.length-1];
		
		testSharedUtils(dsd,mol);
	}

	private static void testSharedUtils(DomainStructureData dsd, String mol) {
		DomainStructureBNFTree dsg1 = new DomainStructureBNFTree(dsd);
		DomainStructureBNFTree.readStructure("A", mol, dsg1);
		DomainPolymerGraph dsg2 = new DomainPolymerGraph(dsd);
		DomainPolymerGraph.readStructure("A", mol, dsg2);
		int k = 1;
		{
			System.out.println("TEST"+k+++" : SINGLE STRANDED");
			ArrayList<DomainSequence> test1 = new ArrayList<DomainSequence>(){
				public String toString(){
					return "Test BNFTree";
				}
			};
			ArrayList<DomainSequence> test2 = new ArrayList<DomainSequence>(){
				public String toString(){
					return "Test PolymerGraph";
				}
			};;
			DomainDesigner_SharedUtils.utilSingleStrandedFinder(dsg1, test1);
			DomainDesigner_SharedUtils.utilSingleStrandedFinder(dsg2, test2);
			for(ArrayList<DomainSequence> g : new ArrayList[]{test1,test2}){
				System.out.println("TEST:"+g);
				for(DomainSequence q : g){
					System.out.println(q.toString(dsd));
				}
			}
		}
		{
			System.out.println("TEST"+k+++" : HAIRPIN INTERNALS");
			ArrayList<DomainSequence> test1 = new ArrayList<DomainSequence>(){
				public String toString(){
					return "Test BNFTree";
				}
			};
			ArrayList<DomainSequence> test2 = new ArrayList<DomainSequence>(){
				public String toString(){
					return "Test PolymerGraph";
				}
			};;
			DomainDesigner_SharedUtils.utilHairpinInternalsFinder(dsg1, test1);
			DomainDesigner_SharedUtils.utilHairpinInternalsFinder(dsg2, test2);
			for(ArrayList<DomainSequence> g : new ArrayList[]{test1,test2}){
				System.out.println("TEST:"+g);
				for(DomainSequence q : g){
					System.out.println(q.toString(dsd));
				}
			}
		}
		{
			System.out.println("TEST"+k+++" : HAIRPIN CLOSING");
			ArrayList<DomainSequence[]> test1 = new ArrayList<DomainSequence[]>(){
				public String toString(){
					return "Test BNFTree";
				}
			};
			ArrayList<DomainSequence[]> test2 = new ArrayList<DomainSequence[]>(){
				public String toString(){
					return "Test PolymerGraph";
				}
			};;
			DomainDesigner_SharedUtils.utilHairpinClosingFinder(dsg1, test1);
			//DomainDesigner_SharedUtils.utilHairpinClosingFinder(dsg2, test2);
			for(ArrayList<DomainSequence[]> g : new ArrayList[]{test1,test2}){
				System.out.println("TEST:"+g);
				for(DomainSequence[] q : g){
					for(DomainSequence d : q){
						System.out.print(d.toString(dsd));
					}
					System.out.println();
				}
			}
		}
	}
}
