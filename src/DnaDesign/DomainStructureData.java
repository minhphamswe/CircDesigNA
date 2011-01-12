package DnaDesign;
import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The structural tree derived from a molecule description string (see the BNF grammar)
 */
public class DomainStructureData {
	public int[] domainLengths;
	public DomainStructure[] structures;
	public int[] domains;
	public String moleculeName;
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
	private Map<Integer, String> compositionConstraints = new TreeMap();
	private static final int FLAG_CONSERVEAMINOS = 2;
	private Map<Integer, Integer> otherRuleFlags = new TreeMap();
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

	/** 
	 * Pass in whole block, please
	 **/
	public static void readDomainDefs(String domainDefsBlock, DomainStructureData out){
		out.structures = null;
		out.nameMap.clear();
		out.domainConstraints.clear();
		out.compositionConstraints.clear();
		out.otherRuleFlags.clear();
		ArrayList<Integer> domainLengths = new ArrayList<Integer>();
		Scanner in = new Scanner(domainDefsBlock);
		int k = -1;
		while(in.hasNextLine()){
			k = parseDomainDefLine(in.nextLine(),out,domainLengths,k);
		}
		out.domainLengths = new int[domainLengths.size()];
		for(k = 0; k < domainLengths.size(); k++){
			out.domainLengths[k] = domainLengths.get(k);
			//Test the seq constraints
			//Don't use Default, this prevents us from crashing on conflicting definition.
			DesignSequenceConstraints dsc = new DesignSequenceConstraints();
			out.loadConstraints(k, dsc, true);
		}
	}

	private static int parseDomainDefLine(String nextLine, DomainStructureData out, List<Integer> domainLengths, int k) {
		String[] line = nextLine.trim().split("\\s+");
		if (line.length<=1){ //0 is sort of impossible though.
			return k;
		}
		k++;
		domainLengths.add(-1);//initialize length
		String domainID = line[0];
		if (domainID.length()==0){
			throw new RuntimeException("'Empty' Domain ID: on line "+k);
		}
		if (Character.isWhitespace(nextLine.charAt(0))){
			throw new RuntimeException("'Empty' Domain ID: on line "+k+". (Leading whitespace?)");
		}
		out.nameMap.put(domainID,k);
		int seqIndex = 1;
		int seqLen = -1;
		if (line[1].matches("\\d+")){
			//We have length map (optional)
			seqIndex = 2;
			domainLengths.set(k,seqLen = new Integer(line[1]));
		}
		//Sequence constraints...
		if (line.length>seqIndex){
			//Regions of characters enclosed in square bracket will be lowercased, meaning "lock".
			if (line[seqIndex].charAt(0)!='-'){
				StringBuffer constraintParsed = new StringBuffer();
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
							throw new RuntimeException("Stack Underflow of nonconstraint bracket: "+line[seqIndex]);
						}
						inBracket = false;
						continue;
					}
					if (Character.toLowerCase(kc)=='u'){
						kc = 'T'; // rna / dna treated equally.
					}
					//The case of characters is significant to deeper stages of the pipeline;
					//Capital means mutable.
					kc = inBracket ? Character.toUpperCase(kc) : Character.toLowerCase(kc);
					constraintParsed.append(kc);
				}
				line[seqIndex] = constraintParsed.toString();
				//Ok! load the constraint. Default to "unconstrained".
				if (line[seqIndex].equalsIgnoreCase("TBD")){
					//Lol, silly format.
				} else {
					out.domainConstraints.put(k,line[seqIndex]);
					domainLengths.set(k,seqLen = line[seqIndex].length());
				}
				seqIndex++;
			}
			//Do we have flags?
			int flagSum = 0;
			Pattern decodeArg = Pattern.compile("\\-(\\w+)(\\((.*)\\))?");
			for(int flag = seqIndex; flag < line.length; flag++){
				Matcher m = decodeArg.matcher(line[flag]);
				if(m.find()){
					String paramName = m.group(1);
					try {
						String args = m.group(3);
						if (paramName.equalsIgnoreCase("p")){
							if (seqLen%3!=0){
								throw new RuntimeException("Domain "+domainID+" not a valid protein sequence - length not a multiple of 3");
							}
							flagSum |= FLAG_CONSERVEAMINOS;
						} 
						if (paramName.equalsIgnoreCase("seq")){
							String old2 = out.compositionConstraints.get(k);
							if (old2==null){
								old2 = "";
							}
							while(old2.endsWith(",")){
								old2 = old2.substring(0,old2.length()-1);
							}
							if (old2.contains(",")){
								old2+=",";
							}
							out.compositionConstraints.put(k, old2+args);
						}
					} catch (Throwable e){
						throw new RuntimeException("Invalid args to '-"+paramName+"': "+e.getMessage());
					}
				}
			}
			out.otherRuleFlags.put(k,flagSum);
		}
		if (seqLen==-1){
			throw new RuntimeException("Assertion error - no length for domain '"+domainID+"'");
		}
		if (seqLen < 2){
			throw new RuntimeException("1-base domains not allowed for now.");
		}
		return k;
	}

	
	public static void readStructure(String moleculeName, String dnaString, DomainStructureData out){
		//Tested and verified. Group 2 : Domain, with comp flag, Group 1 + Group 3: Structural flag
		final Pattern regexp = Pattern.compile("(.*?)(\\w+\\*?)(.*?)($|[\\|\\}\\[]+)");
		
		out.moleculeName = moleculeName;
		out.domains = null;
		out.structures = null;
		Matcher m = regexp.matcher(dnaString);
		int whichDomain = 0, seqId = 0;
		LinkedList<Integer> parens = new LinkedList();
		TreeMap<Integer, Integer> lockMap = new TreeMap();
		TreeMap<Integer,DomainStructure> out2 = new TreeMap();
		TreeMap<Integer, Integer> whichDomainsInOrder = new TreeMap();
		int highestDomainUsed = -1;
		while(m.find()){
			try {
				String domainName = m.group(2);
				//Decode which domain, extract the "star" character (reverse complement flag)
				int domainNameL = domainName.length();
				if (domainName.endsWith("*")){
					domainNameL --;
				}
				int numberDomain = out.lookupDomainName(domainName.substring(0,domainNameL));
				int numberDomain2 = numberDomain;
				if (numberDomain2 < 0){
					throw new RuntimeException("Invalid domain: "+domainName);
				}
				if (domainName.endsWith("*")){
					numberDomain2 |= DNA_COMPLEMENT_FLAG;
				}
				highestDomainUsed = Math.max(highestDomainUsed,numberDomain);
				whichDomainsInOrder.put(whichDomain,numberDomain2);

				String match = m.group(1)+m.group(3);
				if (match==null){
					match = "";
				}
				String struct = match;
				if (struct.contains("(")){
					out2.put(seqId,new SingleStranded(whichDomain));
					parens.add(seqId);
				} else if (struct.contains(")")){
					if (parens.isEmpty()){
						throw new RuntimeException("Empty Stack: "+m.group(0));
					}
					int mP = parens.removeLast();
					DomainStructure remove = out2.remove(mP);
					if (!(remove instanceof SingleStranded)){
						throw new RuntimeException("Assertion error.");
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
				} else {
					//if (struct.length()==0 || struct.contains(".")){
					out2.put(seqId,new SingleStranded(whichDomain));
				}
				final String delim = m.group(4);
				if (delim!=null){
					//Checking valid molecule "enders"
					if (delim.contains("[")){
						for(int k = 0; k < delim.length(); k++){
							//Only valid use of '[' is after '}', except for the first character (taken care of above)
							if (delim.charAt(k)=='[' && (k==0 || delim.charAt(k-1)!='}')){
								throw new RuntimeException("Invalid structure: Use } before [.");
							}
						}
					}
					if (delim.matches("(.*)\\[(.*)\\}(.*)")){
						throw new RuntimeException("Empty strand; no domains");
					}
					if (delim.contains("}")){
						//3- end. 
						out2.put(++seqId,new ThreePFivePOpenJunc());
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
		//This is dependent on the input sequence. Not on Domain Definitions.
		out.domains = new int[numDomains];
		for(Entry<Integer, Integer> q : whichDomainsInOrder.entrySet()){
			out.domains[q.getKey()] = q.getValue();
		}
		out.structures = new DomainStructure[out2.size()];
		int i = 0;
		for(DomainStructure struct : out2.values()){
			//Top level substructures - recursively handle subconformation
			out.structures[i]=struct;
			if (i>0 && out.structures[i-1] instanceof ThreePFivePOpenJunc){
				throw new RuntimeException("Invalid use of 5' end. Molecule is not connected.");
			}
			i++;
			struct.handleSubConformation(out.domainLengths,out.domains);
		}		
		return;
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
		/**
		 * +1 for each domain, +1 for each 3' end.
		 */
		public int countLabeledElements() {
			int ret = 0;
			ret += this.sequencePartsInvolved.length;
			for(DomainStructure q : subStructure){
				ret += q.countLabeledElements();
			}
			if (this instanceof ThreePFivePOpenJunc){
				ret++;
			}
			return ret;
		}
		public void handleSubConformation(int[] domainLengths, int[] domains) {
			
		}
		
		/**
		 * Returns the "width" of this domain structure as seen by its parent.
		 * For instance, a hairpin loop wants to know how wide a substructure is so it can calculate
		 * the total circumference of the loop.
		 */
		public static int getOuterLevelSpace(DomainStructure q,int[] domainLengths, int[] domains) {
			int circ = 0;
			if (q instanceof HairpinStem){
				circ+= 1; 
			} else if (q instanceof SingleStranded){
				SingleStranded s = (SingleStranded)q;
				for(int p : q.sequencePartsInvolved){
					circ += domainLengths[domains[p] & DNA_SEQ_FLAGSINVERSE];
				}
			} else if (q instanceof ThreePFivePOpenJunc){
				circ += 2;
			}
			return circ;
		}
	}
	public static class HairpinStem extends DomainStructure {
		public HairpinStem(int ... whichDomain) {
			super(whichDomain);
		}
		
		public int leftRightBreak = -1;
		public int innerCurveCircumference = 0;
		
		public void handleSubConformation(int[] domainLengths, int[] domains) {
			if (subStructure.size()==0){
				throw new RuntimeException("Hairpin with no loop: Invalid");
			}
			loop:for(DomainStructure q : subStructure){
				if (q instanceof HairpinStem){
					((HairpinStem)q).handleSubConformation(domainLengths,domains);
				}
			}
			if (subStructure.size()>1){
				int index = 0;
				loop:for(DomainStructure q : subStructure){
					if (q instanceof ThreePFivePOpenJunc){
						leftRightBreak = index;
						break loop;
					}
					index++;
				};
			}
			if (subStructure.size()==1 && (subStructure.get(0) instanceof HairpinStem)){
				innerCurveCircumference = 0; //Just continue the stem
			} else {
				innerCurveCircumference = 0; //The "opening" of the hairpin takes up some of the ring
				loop:for(DomainStructure q : subStructure){
					innerCurveCircumference += getOuterLevelSpace(q,domainLengths,domains);
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
			if (whichDomain.length>1){
				System.out.println("MultisingleStrand");
			}
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
	public String toString(){
		String lineS = System.getProperty("line.separator");
		StringBuffer sb = new StringBuffer();
		for(int k = 0; k < structures.length; k++){
			sb.append(structures[k].toString());
			if (k+1<structures.length){
				sb.append(lineS);
			}
		}
		return sb.toString();
	}
	
	private static class getStructureString_helperClass{
		private StringBuffer sb;
		private boolean lastWasThreePrime = true, lastWasItem = false;
		private DomainStructureData dsd;
		public getStructureString_helperClass(DomainStructureData dsd){
			this.dsd = dsd;
			sb = new StringBuffer();
			openStructure();
		}
		private void openStructure() {
			if (lastWasThreePrime){
				sb.append("[");
			}
			lastWasThreePrime = false;
			if (lastWasItem){
				sb.append("|");
			}
			lastWasItem = false;
		}
		public String toString(){
			return sb.toString();
		}
		public void previsit(DomainStructure q) {
			if (q instanceof SingleStranded){
				openStructure();
				sb.append(dsd.getDomainName(dsd.domains[q.sequencePartsInvolved[0]]));
				lastWasItem = true;
			}
			if (q instanceof HairpinStem){
				openStructure();
				sb.append("(");
				sb.append(dsd.getDomainName(dsd.domains[q.sequencePartsInvolved[0]]));
				lastWasItem = true;
			}
			if (q instanceof ThreePFivePOpenJunc){
				lastWasItem = false;
				openStructure();
				sb.append("}");
				lastWasThreePrime = true;
			}
		}
		public void postvisit(DomainStructure q) {
			//Only has effect on hairpin loops, which have a pre- and post- part.
			//That is, all other actions occur on the previsit.
			if (q instanceof HairpinStem){
				openStructure();
				sb.append(dsd.getDomainName(dsd.domains[q.sequencePartsInvolved[1]]));
				sb.append(")");
				lastWasItem = true;
			}
		}
	}
	public static String getStructureString(DomainStructureData in){
		getStructureString_helperClass out = new getStructureString_helperClass(in);
		for(DomainStructure q : in.structures){
			getStructureString_r(q,out);
		}
		return out.toString();
	}
	private static void getStructureString_r(DomainStructure q, getStructureString_helperClass out) {
		out.previsit(q);
		for(DomainStructure d : q.subStructure){
			getStructureString_r(d,out);
		}
		out.postvisit(q);
	}

	public int getDomainLength(int i) {
		return domainLengths[domains[i]&DNA_SEQ_FLAGSINVERSE];
	}

	boolean maintainAminos(int k) {
		Integer integer = otherRuleFlags.get(k);
		if (integer==null){
			return false;
		}
		return (integer & FLAG_CONSERVEAMINOS)!=0;
	}

	public void loadConstraints(int i, DesignSequenceConstraints dsc, boolean crashOnOverwrite) {
		String args = compositionConstraints.get(i);
		if (args==null){
			return;
		}
		String[] array = args.split(",");
		if (array.length%3!=0){
			throw new RuntimeException("Each contraint has 3 parts: base, min, and max");
		}
		
		for(int k = 0; k < array.length; k+=3){
			ArrayList<Integer> bases = new ArrayList();
			for(char q : array[k].toCharArray()){
				if (Character.isLetter(q)){
					q = Character.toUpperCase(q);
					int base = DnaDefinition.decodeBaseChar(q);
					if (base == 0 || base >= DNAFLAG_ADD){
						throw new RuntimeException("Invalid base: "+array[k]);
					}
					bases.add(base);
				}
			}
			//pure base.
			int num1 = parsePercent(array[k+1],domainLengths[i],false);
			int num2 = parsePercent(array[k+2],domainLengths[i],true);
			if (num1 < -1 || num2 < -1){
				throw new RuntimeException("Bound values must be >= -1. -1 means no bound.");
			}
			if (num2 !=-1 && num1 != -1 && num2 < num1){
				throw new RuntimeException("Invalid bound: max < min");
			}
			int[] bases2 = new int[bases.size()];
			int count = 0;
			for(int j : bases){
				bases2[count++]=j;
			}
			boolean hadMin = dsc.setMinConstraint(num1, bases2);
			boolean hadMax = dsc.setMaxConstraint(num2, bases2);
			if (crashOnOverwrite && (hadMax||hadMin)){
				throw new RuntimeException("Duplicate bounds for "+array[k]);
			}
		}
	}
	private int parsePercent(String percent, int of, boolean roundUp){
		boolean isP = percent.endsWith("%");
		boolean isF = percent.contains(".");
		if (isP || isF){
			double value;
			if (isP){
				value = new Double(percent.substring(0,percent.length()-1))/100.;	
			} else {
				//isF only
				value = new Double(percent);
			}
			value *= of;
			if (roundUp){
				return (int) Math.ceil(value);
			} else {
				return (int) Math.floor(value);
			}
		}
		return new Integer(percent);
	}
}
