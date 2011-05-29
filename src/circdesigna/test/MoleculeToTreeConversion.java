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
import java.util.Collection;
import java.util.Scanner;

import circdesigna.AbstractComplex;
import circdesigna.DomainDefinitions;
import circdesigna.DomainPolymerGraph;
import circdesigna.DomainStructureBNFTree;
import circdesigna.config.CircDesigNAConfig;



public class MoleculeToTreeConversion {
	public static void main(String[] args){
		Collection<AbstractComplex> dsd = getInputTree(2, GRAPH);
		for(AbstractComplex ds : dsd){
			System.out.println(ds.getStructureString()+"=");
			if (ds instanceof DomainPolymerGraph){
				DomainPolymerGraph ds2 = (DomainPolymerGraph)ds;
				Collection<Integer> getStrandRotations = ds2.getStrandRotations();
				for(int k : getStrandRotations){
					DomainPolymerGraph rotation = ds2.getRotation(ds2.getMoleculeName()+"x"+k,k);
					System.out.println(rotation.getStructureString());
				}
			}
			System.out.println(ds);
			System.out.println("////////");
		}
	}
	private static final int BNF = 0, GRAPH = 1, BOTH = 2;
	public static Collection<AbstractComplex> getInputTree(int numTrees, int structureForm){
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
		CircDesigNAConfig config = new CircDesigNAConfig();
		String domainDefsBlock = domainDefs.toString();
		ArrayList<AbstractComplex> list = new ArrayList();
		DomainDefinitions dsd = new DomainDefinitions(config);
		DomainDefinitions.readDomainDefs(domainDefsBlock, dsd);
		for(int i = 0 ;i < numTrees; i++){
			System.out.println("Enter the molecule to be converted to a tree");
			String mol;
			while(true){
				mol = in.nextLine().trim();
				if (mol.length()!=0){
					break;
				}
			}
			String[] mol2 = mol.split("\\s+");
			mol = mol2[mol2.length-1];
			AbstractComplex dsg;
			switch(structureForm){
			case BNF:
				dsg = new DomainStructureBNFTree(dsd);
				DomainStructureBNFTree.readStructure("A "+mol,(DomainStructureBNFTree) dsg);
				list.add(dsg);
				break;
			case GRAPH:
				dsg = new DomainPolymerGraph(dsd);
				DomainPolymerGraph.readStructure("A "+mol,(DomainPolymerGraph)dsg);
				list.add(dsg);
				break;
			case BOTH:
				dsg = new DomainStructureBNFTree(dsd);
				DomainStructureBNFTree.readStructure("A "+mol,(DomainStructureBNFTree) dsg);
				list.add(dsg);
				dsg = new DomainPolymerGraph(dsd);
				DomainPolymerGraph.readStructure("A "+mol,(DomainPolymerGraph)dsg);
				list.add(dsg);
			}
		}
		return list;
	}
}
