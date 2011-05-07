package DnaDesign;

import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A more physical representation of a domain structure than the polymer graph, which already
 * splits the molecule input into a tree structure based on looping.
 * 
 * This parser does NOT generalize to pseudoknotted structures. It is also not obvious how to
 * perform "cyclic rotations" on this tree-type structure.
 * 
 * It is highly convenient for graphics routines, and it supports certain error-checking features
 * which DomainPolymerGraph may not. Specifically, this class ensures that any hairpin you are designing
 * is non-empty.
 * 
 * @author Benjamin
 */
public class DomainStructureBNFTree implements AbstractComplex{
	private DomainDefinitions domainDefs;
	public String moleculeName;
	
	//Indexed by position in the string.
	public int[] domains;
	//Tree view of the molecule (used only for displaying the molecule)
	public DomainStructure[] structures;
	public int outerCurveCircum = -1;
	
	public DomainStructureBNFTree(DomainDefinitions domainDefs){
		this.domainDefs = domainDefs;
	}

	private static class ParserState {
		private String structure;
		private DomainStructureBNFTree out;
		public int i, whichDomain, seqId;

		private LinkedList<Integer> parens;
		TreeMap<Integer,DomainStructure> out2 = new TreeMap();
		TreeMap<Integer, Integer> whichDomainsInOrder = new TreeMap();

		public ParserState(String structure, DomainStructureBNFTree out) {
			this.structure = structure;
			this.out = out;
			i = 0;

			parens = new LinkedList();
			out2 = new TreeMap();
			whichDomainsInOrder = new TreeMap();
			whichDomain = 0;
			seqId = 0;
		}
		private void err(String fmt, Object ... args) {
			throw new RuntimeException(String.format(fmt,args));
		}
		private String local(){
			int num = 8;
			if (structure.length()<=num){
				return "("+structure+")";
			}
			int begin = Math.max(0,i-num/2);
			int end = Math.min(i+num/2,structure.length());
			return "..."+structure.substring(begin,end)+"...";
		}
		private char charAt(int i) {
			return structure.charAt(i);
		}
		private boolean isElement(int i2, String string) {
			if (i2 >= structure.length()){
				return false;
			}
			return string.contains(charAt(i2)+"");
		}
		//begin published section
		void parse(){
			if (!isElement(i,"[")){
				parseV();
				if (i < structure.length()){
					err("Strands must begin with a '[' character.", local());
				}
			} else {
				parseU();
			}
		}
		void parseU() {
			if (!isElement(i,"["))
				err("Strands must begin with a 3' end ('[' character).", local());
			i++;
			parseV();
			if (i < structure.length()){
				if (!isElement(i,"}"))
					err("Strands must end with a 5' end ('}' character).", local());
				i++;
				if (out2.isEmpty() || out2.lastEntry().getValue() instanceof ThreePFivePOpenJunc){
					//Empty strand.
				} else {
					out2.put(seqId++,new ThreePFivePOpenJunc());
				}
			}
			//If there are more non-whitespace characters, parseU again.
			if (i < structure.length())
				parseU();
		}
		void parseV(){
			parseD();
			if (isElement(i,"|")){
				i++;
				parseV();
			}
		}
		private void parseD() {
			char form = '.';
			String domainName = null;
			if (isElement(i,"().")){
				//Case 1: SB.
				form = charAt(i);
				i++;
				domainName = parseB();
			} else {
				//Case 2: BS.
				domainName = parseB();
				if (isElement(i,"().")){
					form = charAt(i);
					i++;
				}
			}
			addDomain(form,domainName);
		}
		private void addDomain(char form, String domainName) {
			int domainNameL = domainName.length();
			if (domainName.endsWith("*")){
				domainNameL --;
			}
			int numberDomain = out.domainDefs.lookupDomainName(domainName.substring(0,domainNameL));
			int numberDomain2 = numberDomain;
			if (numberDomain2 < 0){
				throw new RuntimeException("Invalid domain: "+domainName);
			}

			if (out.domainDefs.domainLengths[numberDomain2] == 0){
				return; //Ignore domains of length 0 when parsing the structures.
			}
			
			if (domainName.endsWith("*")){
				numberDomain2 |= DNA_COMPLEMENT_FLAG;
			}
			//highestDomainUsed = Math.max(highestDomainUsed,numberDomain);
			whichDomainsInOrder.put(whichDomain,numberDomain2);
				
			switch(form){
			case '(':
				out2.put(seqId,new SingleStranded(whichDomain));
				parens.add(seqId);
				break;
			case ')':
				if (parens.isEmpty()){
					err("Empty Stack: %s",local());
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
				break;
			case '.':
				//if (struct.length()==0 || struct.contains(".")){
				out2.put(seqId,new SingleStranded(whichDomain));
				break;
			default:
				err("Invalid structural flag: %s",form+"");
			}
			
			seqId++;
			whichDomain++;
		}
		private String parseB() {
			StringBuffer toRet = new StringBuffer();
			while(true){
				if (isElement(i,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_")){
					toRet.append(charAt(i));
					i++;
				} else { 
					break;
				}
			}
			if (toRet.length()==0){
				err("Empty domain %s",local());
			}
			if (isElement(i,"*")){
				toRet.append(charAt(i));
				i++;
			}
			return toRet.toString();
		}
		public void output() {
			//Debug
			int numDomains = whichDomain--;
			//This is dependent on the input sequence. Not on Domain Definitions.
			out.domains = new int[numDomains];
			for(Entry<Integer, Integer> q : whichDomainsInOrder.entrySet()){
				out.domains[q.getKey()] = q.getValue();
			}
			out.structures = new DomainStructure[out2.size()];
			out.outerCurveCircum = 0;
			int i = 0;
			int numHairpinsOnOuterLoop = 0;
			for(DomainStructure struct : out2.values()){
				//Top level substructures - recursively handle subconformation
				out.structures[i]=struct;
				if (i>0 && out.structures[i-1] instanceof ThreePFivePOpenJunc){
					throw new RuntimeException("Invalid use of 5' end. Molecule is not connected.");
				}
				i++;
				struct.handleSubConformation(out.domainDefs.domainLengths,out.domains);
				
				if (struct instanceof HairpinStem){
					numHairpinsOnOuterLoop++;
				}
				out.outerCurveCircum += struct.getOuterLevelSpace(struct,out.domainDefs.domainLengths,out.domains);
			}	
			if (numHairpinsOnOuterLoop < 2){
				//Linear structure
				out.outerCurveCircum = -1;
			}
		}
	}
	public static void readStructure(String moleculeName, String structure, DomainStructureBNFTree out){
		out.moleculeName = moleculeName;
		out.domains = null;
		out.structures = null;
		
		//First disregard all whitespace.
		//This does not really require regular expressions.
		structure = structure.replaceAll("\\s","");
		//readStructure_(moleculeName,structure,out);

		ParserState p = new ParserState(structure,out);
		p.parse();
		p.output();
	}

	public static void readStructure_(String moleculeName, String dnaString, DomainStructureBNFTree out){
		//Tested and verified. Group 2 : Domain, with comp flag, Group 1 + Group 3: Structural flag
		final Pattern regexp = Pattern.compile("(.*?)(\\w+\\*?)(.*?)($|[\\|\\}\\[]+)");
		
		out.moleculeName = moleculeName;
		out.domains = null;
		out.structures = null;
		Matcher m = regexp.matcher(dnaString);
		int whichDomain = 0, seqId = 0;
		LinkedList<Integer> parens = new LinkedList();
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
				int numberDomain = out.domainDefs.lookupDomainName(domainName.substring(0,domainNameL));
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
		out.outerCurveCircum = 0;
		int i = 0;
		int numHairpinsOnOuterLoop = 0;
		for(DomainStructure struct : out2.values()){
			//Top level substructures - recursively handle subconformation
			out.structures[i]=struct;
			if (i>0 && out.structures[i-1] instanceof ThreePFivePOpenJunc){
				throw new RuntimeException("Invalid use of 5' end. Molecule is not connected.");
			}
			i++;
			struct.handleSubConformation(out.domainDefs.domainLengths,out.domains);
			
			if (struct instanceof HairpinStem){
				numHairpinsOnOuterLoop++;
			}
			out.outerCurveCircum += struct.getOuterLevelSpace(struct,out.domainDefs.domainLengths,out.domains);
		}	
		if (numHairpinsOnOuterLoop < 2){
			//Linear structure
			out.outerCurveCircum = -1;
		}
		return;
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
			//Check valid hairpin
			if (domainLengths[domains[sequencePartsInvolved[0]] & DNA_SEQ_FLAGSINVERSE]!=domainLengths[domains[sequencePartsInvolved[1]] & DNA_SEQ_FLAGSINVERSE]){
				throw new RuntimeException("Invalid duplex between domains of different lengths.");
			}
			loop:for(DomainStructure q : subStructure){
				if (q instanceof HairpinStem){
					((HairpinStem)q).handleSubConformation(domainLengths,domains);
				}
			}
			boolean isEmptyLoop = true;
			if (subStructure.size()>0){
				int index = 0;
				loop:for(DomainStructure q : subStructure){
					if (q instanceof SingleStranded && getOuterLevelSpace(q, domainLengths, domains) > 0){
						isEmptyLoop = false;
					}
					if (!(q instanceof SingleStranded)){
						isEmptyLoop = false;
					}
					if (q instanceof ThreePFivePOpenJunc){
						leftRightBreak = index;
						break loop;
					}
					index++;
				};
			}
			if (isEmptyLoop){
				throw new RuntimeException("Hairpin with no loop: Invalid");
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
		private DomainStructureBNFTree dsg;
		private DomainDefinitions dsd;
		public getStructureString_helperClass(DomainStructureBNFTree dsg){
			this.dsg = dsg;
			this.dsd = dsg.domainDefs;
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
				sb.append(dsd.getDomainName(dsg.domains[q.sequencePartsInvolved[0]]));
				lastWasItem = true;
			}
			if (q instanceof HairpinStem){
				openStructure();
				sb.append("(");
				sb.append(dsd.getDomainName(dsg.domains[q.sequencePartsInvolved[0]]));
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
				sb.append(dsd.getDomainName(dsg.domains[q.sequencePartsInvolved[1]]));
				sb.append(")");
				lastWasItem = true;
			}
		}
	}
	public String getStructureString(){
		getStructureString_helperClass out = new getStructureString_helperClass(this);
		for(DomainStructure q : structures){
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

	public String getMoleculeName() {
		return moleculeName;
	}

	public DomainDefinitions getDomainDefs() {
		return domainDefs;
	}
}
