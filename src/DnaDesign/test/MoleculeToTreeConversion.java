package DnaDesign.test;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Scanner;

import DnaDesign.DomainStructureData;

public class MoleculeToTreeConversion {
	public static void main(String[] args){
		Collection<DomainStructureData> dsd = getInputTree(2);
		for(DomainStructureData ds : dsd){
			System.out.println(DomainStructureData.getStructureString(ds)+"=");
			System.out.println(ds);
			System.out.println("////////");
		}
	}
	public static Collection<DomainStructureData> getInputTree(int numTrees){
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
		ArrayList<DomainStructureData> list = new ArrayList();
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
			DomainStructureData dsd = new DomainStructureData();
			dsd.readDomainDefs(domainDefsBlock, dsd);
			DomainStructureData.readStructure("A",mol,dsd);
			list.add(dsd);
		}
		return list;
	}
}
