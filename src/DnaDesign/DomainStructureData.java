package DnaDesign;
import static DnaDesign.DomainDesigner_ByRandomPartialMutations.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainDesigner_ByRandomPartialMutations.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;
import java.util.Scanner;
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
	/**
	 * Invertible.
	 */
	public Map<String, Integer> nameMap = new TreeMap();
	public String getDomainName(int domain){
		String postpend = ((domain & DNA_COMPLEMENT_FLAG)!=0?"*":"");
		for(Entry<String, Integer> q : nameMap.entrySet()){
			if (q.getValue()==(domain & DNA_SEQ_FLAGSINVERSE)){
				return q.getKey()+postpend;
			}
		}
		return ""+(domain & DNA_SEQ_FLAGSINVERSE)+postpend;
	}
	private Map<Integer, String> domainConstraints = new TreeMap();
	/**
	 * Returns the constraint for the NONCOMPLEMENTED version of this domain.
	 * You must handle complementation yourself in handling the constraints!
	 */
	public String getConstraint(int domain){
		domain &= DNA_SEQ_FLAGSINVERSE;
		if (domainConstraints.containsKey(domain)){
			return domainConstraints.get(domain);
		}
		StringBuffer wT = new StringBuffer();
		for(int k = 0; k < domainLengths[domain]; k++){
			wT.append("-");
		}
		return wT.toString();
	}
	
	private static final Pattern regexp = Pattern.compile("(\\w+\\*?)(<?(.*?)>?)?|}");
	
	/** 
	 * Pass in whole block, please
	 **/
	public static void readDomainDefs(String domainDefsBlock, DomainStructureData out){
		out.nameMap.clear();
		out.domainConstraints.clear();
		ArrayList<Integer> domainLengths = new ArrayList();
		Scanner in = new Scanner(domainDefsBlock);
		int k = -1;
		while(in.hasNextLine()){
			String[] line = in.nextLine().split("\\s+");
			if (line.length<2){
				continue;
			}
			k++;
			domainLengths.add(-1);//initialize length
			out.nameMap.put(line[0],k);
			int seqIndex = 1;
			if (line[1].matches("\\d+")){
				//We have length map (optional)
				seqIndex = 2;
				domainLengths.set(k,new Integer(line[1]));
			}
			//Sequence constraints...
			if (line.length>seqIndex){
				//Regions of characters enclosed in square bracket will be lowercased, meaning "lock".
				StringBuffer parseBracketLowerCase = new StringBuffer();
				boolean inBracket = false;
				for(int e = 0; e < line[seqIndex].length(); e++){
					char kc = line[seqIndex].charAt(e);
					if (kc=='['){
						inBracket = true;
						continue;
					}
					if (kc==']'){
						if (!inBracket){
							//Oops, stack underflow.
							throw new RuntimeException("Stack Underflow of constraint bracket: "+line[seqIndex]);
						}
						inBracket = false;
						continue;
					}
					if (inBracket){
						kc = Character.toLowerCase(kc);
					}
					parseBracketLowerCase.append(kc);
				}
				line[seqIndex] = parseBracketLowerCase.toString();
				//Ok! load the constraint. Default to "unconstrained".
				if (line[seqIndex].equalsIgnoreCase("TBD")){
					//Lol, silly format.
				} else {
					out.domainConstraints.put(k,line[seqIndex]);
					domainLengths.set(k,line[seqIndex].length());
				}
			}
		}
		out.domainLengths = new int[domainLengths.size()];
		for(k = 0; k < domainLengths.size(); k++){
			out.domainLengths[k] = domainLengths.get(k);
		}
	}
	public static void readStructure(String dnaString, DomainStructureData out){
		out.domains = null;
		out.structures = null;
		Matcher m = regexp.matcher(dnaString);
		int whichDomain = 0, seqId = 0;
		LinkedList<Integer> parens = new LinkedList();
		TreeMap<Integer, Integer> lengthMap = new TreeMap();
		TreeMap<Integer, Integer> lockMap = new TreeMap();
		TreeMap<Integer,DomainStructure> out2 = new TreeMap();
		TreeMap<Integer, Integer> whichDomainsInOrder = new TreeMap();
		int highestDomainUsed = -1;
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
					int numberDomain = out.lookupDomainName(domainName.substring(0,domainNameL));
					int numberDomain2 = numberDomain;
					if (numberDomain2 < 0){
						//TODO: Domain targetted exceptions?
						throw new RuntimeException("Invalid domain: "+domainName);
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
		if (highestDomainUsed<0){
			throw new RuntimeException("Empty strand; no domains");
		}
		//Debug
		int numDomains = whichDomain--;
		out.domainLengths = expandToLength(out.domainLengths, highestDomainUsed+1,-1);
		for(Entry<Integer, Integer> q : lengthMap.entrySet()){
			out.domainLengths[q.getKey()] = q.getValue();
		}
		//This is dependent on the input sequence. Not on Domain Definitions.
		out.domains = new int[numDomains];
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
	private static int[] expandToLength(int[] old, int newSize, int fillValue){
		int[] newOld = new int[newSize];
		Arrays.fill(newOld,fillValue);
		if (old!=null){
			System.arraycopy(old,0,newOld,0,Math.min(old.length,newOld.length));
		}
		return newOld;
	}
	
	private int lookupDomainName(String substring) {
		if (!nameMap.containsKey(substring)){
			return -1;
		}
		return nameMap.get(substring);
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
		public int countDomainsRecursively() {
			int ret = 0;
			if (this instanceof SingleStranded){
				ret += ((SingleStranded)this).sequencePartsInvolved.length;
			}
			for(DomainStructure q : subStructure){
				ret += q.countDomainsRecursively();
			}
			return ret;
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
