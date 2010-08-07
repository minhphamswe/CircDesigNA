package DNASim$DB;
import static DNASim$DB.DomainDesigner_ByRandomPartialMutations.DNA_COMPLEMENT_FLAG;
import static DNASim$DB.DomainDesigner_ByRandomPartialMutations.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DomainStructureData {
	public int[] domainLengths;
	public DomainStructure[] structures;
	public int[] domains;
	public static final int DEFAULT_TO_LENGTH = 8;
	
	private static final Pattern regexp = Pattern.compile("(\\d+\\*?)(<(.*?)>)?|}");
	public static void readStructure(String dnaString, DomainStructureData out){
		Matcher m = regexp.matcher(dnaString);
		int whichDomain = 0, seqId = 0;
		LinkedList<Integer> parens = new LinkedList();
		TreeMap<Integer, Integer> lengthMap = new TreeMap();
		TreeMap<Integer, Integer> lockMap = new TreeMap();
		TreeMap<Integer,DomainStructure> out2 = new TreeMap();
		TreeMap<Integer, Integer> whichDomainsInOrder = new TreeMap();
		int highestDomainUsed = 0;
		while(m.find()){
			try {
				if (m.group(0).equals("}")){
					//3- end. 
					whichDomain --; //No domain!
					out2.put(seqId,new ThreePFivePOpenJunc());
				} else {
					String domainName = m.group(1);
					//Decode which domain
					int domainNameL = domainName.length();
					if (domainName.endsWith("*")){
						domainNameL --;
					}
					int numberDomain = new Integer(domainName.substring(0,domainNameL));
					int numberDomain2 = numberDomain;
					if (numberDomain2 < 0){
						throw new RuntimeException("Invalid domain ID: "+numberDomain2);
					}
					if (domainName.endsWith("*")){
						numberDomain2 |= DNA_COMPLEMENT_FLAG;
					}
					highestDomainUsed = Math.max(highestDomainUsed,numberDomain);
					whichDomainsInOrder.put(whichDomain,numberDomain2);
					
					String match = m.group(3);
					if (match==null){
						match = "";
					}
					String lenStr = match.replaceAll("[^0-9]|\\s+", "");
					int length = lenStr.length()==0?-1:Integer.decode(lenStr);
					String struct = match.replaceAll("[0-9]|\\s+","");
					if (struct.length()==0 || struct.contains(".")){
						out2.put(seqId,new SingleStranded(whichDomain));
					} else if (struct.contains("(")){
						out2.put(seqId,new SingleStranded(whichDomain));
						parens.add(seqId);
					} else if (struct.contains(")")){
						if (parens.isEmpty()){
							throw new RuntimeException("Empty Stack: "+m.group(0));
						}
						int mP = parens.removeLast();
						DomainStructure remove = out2.remove(mP);
						if (!(remove instanceof SingleStranded)){
							throw new RuntimeException("?huh?");
						}
						//Replace with a hairpinStem.
						HairpinStem create = new HairpinStem(remove.sequencePartsInvolved[0],whichDomain);
						SortedMap<Integer, DomainStructure> subMap = out2.subMap(mP, seqId);
						ArrayList<Integer> holder = new ArrayList();
						for(DomainStructure q : subMap.values()){
							create.addSubStructure(q);
						}
						for(int q : subMap.keySet()){
							holder.add(q);
						}
						for(int p : holder){
							out2.remove(p);
						}
						out2.put(seqId, create);
					}
					if (length >= 0){
						if (length==0){
							throw new RuntimeException("Domain of length 0: "+numberDomain);
						}
						if (lengthMap.containsKey(numberDomain)){
							int was  = lengthMap.get(numberDomain);
							if (was!=length){
								throw new RuntimeException("Differing lengths for same domain: "+numberDomain);
							}
						}
						lengthMap.put(numberDomain,length);
					}
				}
			} finally {
				whichDomain++;
				seqId++;
			}
		}
		if (highestDomainUsed==0){
			throw new RuntimeException("Empty strand; no domains");
		}
		//Debug
		int numDomains = whichDomain--;
		out.domainLengths = new int[highestDomainUsed+1];
		Arrays.fill(out.domainLengths,-1);
		out.domains = new int[numDomains];
		for(Entry<Integer, Integer> q : lengthMap.entrySet()){
			out.domainLengths[q.getKey()] = q.getValue();
		}
		for(Entry<Integer, Integer> q : whichDomainsInOrder.entrySet()){
			out.domains[q.getKey()] = q.getValue();
		}
		out.structures = new DomainStructure[out2.size()];
		int i = 0;
		for(DomainStructure struct : out2.values()){
			out.structures[i++]=struct;
			if (struct instanceof HairpinStem){
				((HairpinStem)struct).handleSubConformation(out.domainLengths,out.domains);
			}
		}
	}
	
	public static class DomainStructure {
		public DomainStructure (int ... sequencePartsInvolved){
			this.sequencePartsInvolved = sequencePartsInvolved;
		}
		//Numbering is 0 - length of input structure sequence.
		public int[] sequencePartsInvolved;
		public ArrayList<DomainStructure> subStructure = new ArrayList();
		public float random0 = (float)Math.random();
		public void addSubStructure(DomainStructure q) {
			subStructure.add(q);
		}
	}
	public static class HairpinStem extends DomainStructure {
		public HairpinStem(int ... whichDomain) {
			super(whichDomain);
		}
		
		public int leftRightBreak = -1;
		public int innerCurveCircumference = 0;
		
		public void handleSubConformation(int[] domainLengths, int[] domains) {
			loop:for(DomainStructure q : subStructure){
				if (q instanceof HairpinStem){
					((HairpinStem)q).handleSubConformation(domainLengths,domains);
				}
			}
			boolean isConnected = true;
			if (subStructure.size()>1){
				int index = 0;
				loop:for(DomainStructure q : subStructure){
					if (q instanceof ThreePFivePOpenJunc){
						isConnected = false;
						leftRightBreak = index;
						break loop;
					}
					index++;
				};
			}
			if (isConnected){
				if (subStructure.size()==1 && (subStructure.get(0) instanceof HairpinStem)){
					innerCurveCircumference = 0; //Just continue the stem
				} else {
					innerCurveCircumference = 3; //The "opening" of the hairpin takes up some of the ring
					loop:for(DomainStructure q : subStructure){
						if (q instanceof HairpinStem){
							innerCurveCircumference += 4; //3 for hairpin space
						} else if (q instanceof SingleStranded){
							for(int p : q.sequencePartsInvolved){
								innerCurveCircumference += domainLengths[domains[p] & DNA_SEQ_FLAGSINVERSE];
							}

						}
					}
				}
			}
		}
		
		public String toString(){
			StringBuffer sb = new StringBuffer();
			String line = System.getProperty("line.separator");
			sb.append(super.toString());
			sb.append(" ");
			if (innerCurveCircumference > 0)
				sb.append(innerCurveCircumference);
			for(int i = 0; i < subStructure.size(); i++){
				String[] sub = subStructure.get(i).toString().split(line);
				for(String d : sub){
					if (leftRightBreak==-1){
						sb.append(line+">"+d);
					} else {
						if (i > leftRightBreak){
							sb.append(line+"R"+d);
						} else {
							sb.append(line+"L"+d);
						}
					}
				}
			}
			return sb.toString();
		}
	}
	public static class SingleStranded extends DomainStructure {
		public SingleStranded(int ... whichDomain) {
			super(whichDomain);
		}
		public void addSubStructure(DomainStructure q) {
			throw new RuntimeException("Cannot add to Single Stranded Structures");
		}
	}
	public static class ThreePFivePOpenJunc extends DomainStructure{
		public ThreePFivePOpenJunc() {
			super();
		}
		//When 2 dna strands are combined, don't connect them (i.e., don't make it look like ligation)

		public void addSubStructure(DomainStructure q) {
			throw new RuntimeException("Cannot add to 3'-5' junc");
		}
	}	
	public int getDomainLength(int i) {
		return domainLengths[domains[i]&DNA_SEQ_FLAGSINVERSE];
	}
}
