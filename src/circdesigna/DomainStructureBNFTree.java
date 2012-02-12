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
package circdesigna;

import static circdesigna.DomainSequence.NA_COMPLEMENT_FLAG;
import static circdesigna.DomainSequence.NA_COMPLEMENT_FLAGINV;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.Map.Entry;

import circdesigna.parser.CDNA2PublicParser;
import circdesigna.parser.CDNA2Token;
import circdesigna.parser.CDNA2PublicParser.ParseResult;
import circdesignagui.math.CircumscribedPolygonTool;
import circdesignagui.math.CircumscribedPolygonTool.CircumscribedPolygon;

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
public class DomainStructureBNFTree extends AbstractComplex{
	private DomainDefinitions domainDefs;
	public String moleculeName;
	
	//Indexed by position in the string.
	public int[] domains;
	//Tree view of the molecule (used only for displaying the molecule)
	public DomainStructure[] structures;
	public CircumscribedPolygon outerCurve;
	
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
		//begin published section
		void parse(){
			ParseResult t = CDNA2PublicParser.parse(structure, out.domainDefs);
			out.moleculeName = t.moleculeName;
			
			for(Object token : t.parse){
				if (token instanceof CDNA2Token.Domain){
					addDomain((CDNA2Token.Domain)token);
				} else if (token instanceof CDNA2Token.FivePrimeEnd){
					//Ignore.
				} else if (token instanceof CDNA2Token.ThreePrimeEnd){
					out2.put(seqId++,new ThreePFivePOpenJunc());
					//XXX
					whichDomainsInOrder.put(whichDomain,-1);
					whichDomain++;
				}
			}
		}
		private void addDomain(CDNA2Token.Domain d) {
			String domainName = d.name;
			int domainNameL = domainName.length();
			if (domainName.endsWith("*")){
				domainNameL --;
			}
			int numberDomain = out.domainDefs.lookupDomainName(domainName.substring(0,domainNameL));
			int numberDomain2 = numberDomain;
			if (numberDomain2 < 0){
				throw new RuntimeException("Invalid domain: "+domainName);
			}
			
			/*
			if (out.domainDefs.domainLengths[numberDomain2] == 0){
				return; //Ignore domains of length 0 when parsing the structures.
			}
			*/
			
			if (domainName.endsWith("*")){
				numberDomain2 |= NA_COMPLEMENT_FLAG;
			}
			//highestDomainUsed = Math.max(highestDomainUsed,numberDomain);
			whichDomainsInOrder.put(whichDomain,numberDomain2);
				
			if (d.open){
				out2.put(seqId,new SingleStranded(whichDomain));
				parens.add(seqId);
			} else if (d.close){
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
			
			seqId++;
			whichDomain++;
		}
		public void output() {
			//Debug
			int numDomains = whichDomain;
			//This is dependent on the input sequence. Not on Domain Definitions.
			out.domains = new int[numDomains];
			for(Entry<Integer, Integer> q : whichDomainsInOrder.entrySet()){
				out.domains[q.getKey()] = q.getValue();
			}
			out.structures = new DomainStructure[out2.size()];
			//XXX
			//S.add(ThreePFivePOpenJunc.size);
			int i = 0;
			for(DomainStructure struct : out2.values()){
				//Top level substructures - recursively handle subconformation
				out.structures[i]=struct;
				if (i>0 && out.structures[i-1] instanceof ThreePFivePOpenJunc){
					throw new RuntimeException("Invalid use of 5' end. Molecule is not connected.");
				}
				i++;
				struct.handleSubConformation(out.domainDefs.domainLengths,out.domains);
			}
			
			out.buildOuterCurve();
		}
	}
	public static void readStructure(String moleculeDefinition, DomainStructureBNFTree out){
		out.domains = null;
		out.structures = null;
		
		//readStructure_(moleculeName,structure,out);

		ParserState p = new ParserState(moleculeDefinition,out);
		p.parse();
		p.output();
	}
	/*
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
	*/

	public void buildOuterCurve() {
		ArrayList<Float> S = new ArrayList();
		S.add(1.2f);
		int numHairpinsOnOuterLoop= 0;
		for(DomainStructure struct : structures){
			DomainStructure.getOuterLevelSpace(S,struct,domainDefs.domainLengths,domains);
			if (struct instanceof HairpinStem){
				numHairpinsOnOuterLoop++;
			}
		}
		
		outerCurve = new CircumscribedPolygon();
		outerCurve.S = new float[S.size()];
		for(int k = 0; k < S.size(); k++){
			outerCurve.S[k] = S.get(k);
		}
		CircumscribedPolygonTool.solvePolygonProblem(outerCurve);
		
		if (numHairpinsOnOuterLoop < 2){
			//Linear structure
			outerCurve = null;
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
		 * @param holder 
		 */
		public static boolean getOuterLevelSpace(List<Float> holder, DomainStructure q,int[] domainLengths, int[] domains) {
			if (q instanceof HairpinStem){
				if (holder!=null) holder.add(((HairpinStem)q).openingSize);
				return true;
			} else if (q instanceof SingleStranded){
				boolean nonEmpty = false;
				for(int p : q.sequencePartsInvolved){
					int len = domainLengths[domains[p] & NA_COMPLEMENT_FLAGINV];
					for(int k = 0; k < len; k++){
						if (holder!=null) holder.add(1f);
						nonEmpty= true;
					}
				}
				return nonEmpty;
			} else if (q instanceof ThreePFivePOpenJunc){
				if (holder!=null) holder.add(((ThreePFivePOpenJunc)q).size);
				return true;
			}
			throw new RuntimeException("Unrecognized structure: "+q);
		}
	}
	public static class HairpinStem extends DomainStructure {
		public float openingSize = 2;
		public float hydrogenBondStrength = 1; //Maximum alpha at 1
		
		public HairpinStem(int ... whichDomain) {
			super(whichDomain);
		}
		
		public int leftRightBreak;
		public CircumscribedPolygon innerCurve;
		
		public void handleSubConformation(int[] domainLengths, int[] domains) {
			//Check valid hairpin
			if (domainLengths[domains[sequencePartsInvolved[0]] & NA_COMPLEMENT_FLAGINV]!=domainLengths[domains[sequencePartsInvolved[1]] & NA_COMPLEMENT_FLAGINV]){
				throw new RuntimeException("Invalid duplex between domains of different lengths.");
			}
			loop:for(DomainStructure q : subStructure){
				if (q instanceof HairpinStem){
					((HairpinStem)q).handleSubConformation(domainLengths,domains);
				}
			}
			boolean isEmptyLoop = true;
			leftRightBreak = -1;
			for(int i = 0; i < subStructure.size(); i++){
				DomainStructure q = subStructure.get(i);
				if (getOuterLevelSpace(null, q, domainLengths, domains)){
					isEmptyLoop = false;
				}
				if (q instanceof ThreePFivePOpenJunc){
					if (leftRightBreak==-1){
						leftRightBreak = i;
					} else {
						//leftRightBreak = -2; //More than one break. Only nonnegative values are handled.
						throw new RuntimeException("Invalid use of 5' end. Molecule is not connected.");
					}
				}
			}
			if (isEmptyLoop){
				throw new RuntimeException("Hairpin with no loop: Invalid");
			}
			if (subStructure.size()==1 && (subStructure.get(0) instanceof HairpinStem)){
				innerCurve = null; //Just continue stem.
			} else {
				//Nonempty contents.
				innerCurve = new CircumscribedPolygon();
				ArrayList<Float> S = new ArrayList();
				//Add the opening.
				S.add(openingSize);
				for(DomainStructure q : subStructure){
					getOuterLevelSpace(S,q,domainLengths,domains);
				};
				innerCurve.S = new float[S.size()];
				for(int k = 0; k < S.size(); k++){
					innerCurve.S[k] = S.get(k);
				}
				CircumscribedPolygonTool.solvePolygonProblem(innerCurve);
			}
		}
		

		public String toString(){
			StringBuffer sb = new StringBuffer();
			String line = System.getProperty("line.separator");
			sb.append(super.toString());
			sb.append(" ");
			if (innerCurve != null)
				sb.append(Arrays.toString(innerCurve.S));
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
		public float size = 1.2f;

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

	public Collection<DomainStructure> listStructures(){
		ArrayList<DomainStructure> toRet = new ArrayList();
		for(DomainStructure p : structures){
			listStructures_r(p,toRet);
		}
		return toRet;
	}
	private void listStructures_r(DomainStructure parent, Collection<DomainStructure> list){
		list.add(parent);
		for(DomainStructure p : parent.subStructure){
			listStructures_r(p, list);
		}
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
