package DnaDesign;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import DnaDesign.DomainStructureData.DomainStructure;
import DnaDesign.DomainStructureData.HairpinStem;
import DnaDesign.DomainStructureData.SingleStranded;

public class DomainDesigner_SharedUtils {
	public static void utilJunctionSplitter(List<DomainSequence> out, String toJunctionize, DomainStructureData dsd){
		int[] bases = DomainDesigner_SharedUtils.utilReadSequence(toJunctionize,dsd);
		for(int b = 0; b < bases.length; b++){
			int nexB = b+1;
			if (nexB < bases.length){
				DomainSequence m = new DomainSequence();
				m.setDomains(bases[b],bases[nexB]);
				out.add(m);
				//System.out.println(m);
			}
		}
	}

	public static int[] utilReadSequence(String toJunctionize, DomainStructureData dsd) {
		toJunctionize = toJunctionize.replaceAll("[\\[}]","|").replaceAll("\\s+","");
		toJunctionize = toJunctionize.replaceAll("[|]+","|");
		if (toJunctionize.endsWith("|")){
			toJunctionize = toJunctionize.substring(0,toJunctionize.length()-1);
		}
		if (toJunctionize.startsWith("|")){
			toJunctionize = toJunctionize.substring(1);
		}
		String[] commands = toJunctionize.split("[|]");
		int[] bases = new int[commands.length];
		for(int i = 0; i < bases.length; i++){
			commands[i] = commands[i].replaceAll("[^\\w|*]","");
			if (commands[i].charAt(commands[i].length()-1)=='*'){
				bases[i] = dsd.nameMap.get(commands[i].substring(0,commands[i].length()-1));
				bases[i] |= DomainDesigner_ByRandomPartialMutations.DNA_COMPLEMENT_FLAG;
			} else {
				bases[i] = dsd.nameMap.get(commands[i]);
			}
		}
		return bases;
	}

	public static void utilRemoveDuplicateSequences(List<DomainSequence> seqToSynthesize) {
		ListIterator<DomainSequence> itr = seqToSynthesize.listIterator();
		big: while(itr.hasNext()){
			DomainSequence seq = itr.next();
			for(DomainSequence q : seqToSynthesize){
				if (q!=seq && q.equals(seq)){
					itr.remove();
					continue big;
				}
			}
		}
	}

	public static void utilVanillaTargetFinder(DomainStructureData dsd,
			ArrayList<DomainSequence> MustBeVanilla) {
		ArrayList<Integer> freeList = new ArrayList<Integer>();
		
		for(DomainStructure ds : dsd.structures){
			utilVanillaTargetFinder_R(dsd,ds,MustBeVanilla,freeList);
		}
		if (freeList.size()>0){
			DomainSequence freeListD = new DomainSequence();
			freeListD.setDomains(freeList);
			MustBeVanilla.add(freeListD);
		}
	}

	private static void utilVanillaTargetFinder_R(DomainStructureData dsd, DomainStructure ds,
			ArrayList<DomainSequence> MustBeVanilla, ArrayList<Integer> freeList) {

		//TODO: free strands need to have "hairpins" in the middle.
		//For now, we just add them all.
		if (ds instanceof SingleStranded){
			for(int k : ds.sequencePartsInvolved){
				freeList.add((dsd.domains[k]));
			}
		}
		if (ds instanceof HairpinStem){
			ArrayList<Integer> subFree = new ArrayList();
			for(DomainStructure g : ds.subStructure){
				utilVanillaTargetFinder_R(dsd,g,MustBeVanilla,subFree);
			}	
			if (subFree.size()>0){
				DomainSequence freeListD = new DomainSequence();
				freeListD.setDomains(subFree);
				MustBeVanilla.add(freeListD);
			}
		
		}
		
	}

}
