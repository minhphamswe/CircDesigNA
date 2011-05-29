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
package circdesigna.test;

import java.util.ArrayList;
import java.util.Scanner;

import circdesigna.DomainDefinitions;
import circdesigna.CircDesigNA_SharedUtils;
import circdesigna.DomainPolymerGraph;
import circdesigna.DomainSequence;
import circdesigna.DomainStructureBNFTree;
import circdesigna.config.CircDesigNAConfig;



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
		DomainDefinitions dsd = new DomainDefinitions(config);
		DomainDefinitions.readDomainDefs(domainDefsBlock, dsd);

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

	private static void testSharedUtils(DomainDefinitions dsd, String mol) {
		DomainStructureBNFTree dsg1 = new DomainStructureBNFTree(dsd);
		DomainStructureBNFTree.readStructure("A "+mol, dsg1);
		DomainPolymerGraph dsg2 = new DomainPolymerGraph(dsd);
		DomainPolymerGraph.readStructure("A "+ mol, dsg2);
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
			CircDesigNA_SharedUtils.utilSingleStrandedFinder(dsg1, test1);
			CircDesigNA_SharedUtils.utilSingleStrandedFinder(dsg2, test2);
			for(ArrayList<DomainSequence> g : new ArrayList[]{test1,test2}){
				System.out.println("TEST:"+g);
				for(DomainSequence q : g){
					System.out.println(q.toString(dsd));
				}
			}
		}
		/*
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
		*/
		{
			/*
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
			DomainDesigner_SharedUtils.utilHairpinInternalsFinder(dsg1, test1);
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
			*/
		}
	}
}
