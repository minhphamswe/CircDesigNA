package DnaDesign;

import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import DnaDesign.DomainStructureData.DomainStructure;
import DnaDesign.DomainStructureData.HairpinStem;
import DnaDesign.DomainStructureData.SingleStranded;
import DnaDesign.DomainStructureData.ThreePFivePOpenJunc;

public class DomainDesigner_SharedUtils {
	public static void utilJunctionSplitter(List<DomainSequence> out, String toJunctionize, DomainStructureData dsd){
		int[] bases = DomainDesigner_SharedUtils.utilReadSequence(toJunctionize,dsd);
		for(int b = 0; b < bases.length; b++){
			int nexB = b+1;
			if (nexB < bases.length){
				DomainSequence m = new DomainSequence();
				m.setDomains(bases[b],bases[nexB],dsd);
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
				bases[i] |= DNA_COMPLEMENT_FLAG;
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
					q.appendMoleculeNames(seq);
					itr.remove();
					continue big;
				}
			}
		}
	}


	public static void utilHairpinClosingFinder(DomainStructureData dsd,
			ArrayList<DomainSequence[]> hairpins) {

		DomainStructure[] dsH = new DomainStructure[3];
		for(int i = 0; i < dsd.structures.length; i++){
			dsH[2] = dsd.structures[i];
			for(int k = 0; k < 2; k++){
				dsH[k] = dsH[k+1];
			}
			if (i+1 < dsd.structures.length){
				dsH[2] = dsd.structures[i+1];
			}
			utilHairpinClosingFinder_R(dsd,dsH,hairpins);
			//From now on, put in last.
		}
		//Get the last in there.
		if (dsd.structures.length>1){
			utilHairpinClosingFinder_R(dsd,dsH,hairpins);
		}
	}
	private static void utilHairpinClosingFinder_R(DomainStructureData dsd,
			DomainStructure ds[], ArrayList<DomainSequence[]> hairpins) {
		if (ds[1].subStructure.isEmpty()){
			return;
		}
		if (ds[1] instanceof HairpinStem){
			//null instanceof is always false
			if (ds[0] instanceof SingleStranded && ds[2] instanceof SingleStranded){
				//Surrounded by singles - add score.
				addHairpinLoopOpening(ds[2],ds[0],hairpins,dsd);
			}
			boolean isBroken = false;
			DomainStructure[] dsH = new DomainStructure[3];
			for(int i = 0; i < ds[1].subStructure.size(); i++){
				DomainStructure g = dsH[2] = ds[1].subStructure.get(i);
				for(int k = 0; k < 2; k++){
					dsH[k] = dsH[k+1];
				}
				if (i+1 < ds[1].subStructure.size()){
					dsH[2] = ds[1].subStructure.get(i+1);
				}
				if (g instanceof ThreePFivePOpenJunc){
					isBroken = true;
				}
				if (g instanceof HairpinStem){
					utilHairpinClosingFinder_R(dsd,dsH,hairpins);
				}
			}
			if (!isBroken){
				DomainStructure left = ds[1].subStructure.get(0);
				DomainStructure right = ds[1].subStructure.get(ds[1].subStructure.size()-1);
				if (left instanceof SingleStranded && right instanceof SingleStranded){
					addHairpinLoopOpening(left,right,hairpins,dsd);
				}
			}
		}
	}

	private static void addHairpinLoopOpening(DomainStructure left,
			DomainStructure right, ArrayList<DomainSequence[]> hairpins, DomainStructureData dsd) {
		DomainSequence leftS = new DomainSequence();
		DomainSequence rightS = new DomainSequence();
		leftS.setDomains(left,dsd);
		rightS.setDomains(right,dsd);
		hairpins.add(new DomainSequence[]{leftS,rightS});
	}

	public static void utilVanillaTargetFinder(DomainStructureData dsd,
			ArrayList<DomainSequence> MustBeVanilla) {
		ArrayList<Integer> freeList = new ArrayList<Integer>();
		
		for(DomainStructure ds : dsd.structures){
			utilVanillaTargetFinder_R(dsd,ds,MustBeVanilla,freeList);
		}
		if (freeList.size()>0){
			DomainSequence freeListD = new DomainSequence();
			freeListD.setDomains(freeList,dsd);
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
				freeListD.setDomains(subFree,dsd);
				MustBeVanilla.add(freeListD);
			}
		}
		
	}

	/**
	 * Returns both a new domainDefs and a Molecules text.
	 */
	public static String[] duplicateSystem(String domainDef, String molecul) {
		String[] toRet = new String[2];
		String LN = System.getProperty("line.separator");
		//First the domainDefs
		int whichV = 0;
		StringBuffer out;
		out = new StringBuffer();
		{
			Scanner domainDefIn = new Scanner(domainDef);
			while(domainDefIn.hasNextLine()){
				String line = domainDefIn.nextLine();
				if (line.trim().length()!=0){
					out.append(line);
					out.append(LN);
				}
			}
			whichV = 1;
			for(boolean actuallyAppend : new boolean[]{false,true}){
				Scanner in = new Scanner(domainDef);
				while(in.hasNextLine()){
					String line = in.nextLine();
					if (line.trim().length()==0){
						continue;
					}
					String[] linefake = line.split("\\s+");
					//Add "vX" to first object.
					int v = linefake[0].lastIndexOf('v');
					int vInd = linefake[0].length();
					int newV = -1;
					if (v!=-1){
						vInd = v;
						try {
							newV = new Integer(linefake[0].substring(v+1));
						} catch (NumberFormatException e){
							//e.printStackTrace();
						}
					} else {
						newV = 1;
					}
					whichV = Math.max(whichV,newV);
					if (newV < whichV){
						continue;
					}
					if (actuallyAppend){
						out.append(linefake[0].substring(0,vInd)+"v"+(whichV+1));
						out.append(line.substring(linefake[0].length()));
						out.append(LN);
					}
				}
				in.close();
			}
		}
		toRet[0] = out.toString();
		out = new StringBuffer();
		{
			Scanner moleculIN = new Scanner(molecul);
			while(moleculIN.hasNextLine()){
				String line = moleculIN.nextLine();
				out.append(line);
				out.append(LN);
			}
			moleculIN.close();
			Scanner in = new Scanner(molecul);
			while(in.hasNextLine()){
				String line = in.nextLine();
				String[] linefake = line.split("\\s+");
				if (linefake[0].endsWith("v"+whichV) || whichV == 1){
					//Translate this guy!
					int subLen = ("v"+whichV).length();
					if (whichV==1){
						subLen = 0;
					}
					out.append(linefake[0].substring(0,linefake[0].length()-subLen)+"v"+(whichV+1)+" ");
					Pattern pat = Pattern.compile("\\w+|\\W+");
					Matcher m = pat.matcher(linefake[1]);
					ArrayList<String> line2S = new ArrayList<String>();
					while(m.find()){
						line2S.add(m.group(0));
					}
					String[] line2 = new String[line2S.size()]; line2S.toArray(line2);
					for(int k = 0; k < line2.length; k++){
						if (line2[k].matches("\\W+")){
							out.append(line2[k]);
						} else {
							if (whichV!=1){
								int v = line2[k].lastIndexOf('v');
								if (v==-1){
									throw new RuntimeException("Invalid duplicated molecule: "+line);
								} else {
									try {
										if (new Integer(line2[k].substring(v+1))>=0){
											if (!line2[k].endsWith("v"+whichV)){
												throw new RuntimeException("Invalid duplicated molecule: "+line);
											}
											line2[k] = line2[k].substring(0,v);
										}
									} catch (NumberFormatException e){
										throw new RuntimeException("Invalid duplicated molecule: "+line);
									}
								}
							}
							out.append(line2[k]+"v"+(whichV+1));
						}
					}
					out.append(LN);
				}
			}
			in.close();
			
		}
		toRet[1] = out.toString();
		return toRet;
	}


}
