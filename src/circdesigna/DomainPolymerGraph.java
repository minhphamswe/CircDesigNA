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

import static circdesigna.GeneralizedInteractiveRegion.NA_COMPLEMENT_FLAG;

import java.awt.Color;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import circdesigna.parser.CDNA2PublicParser;
import circdesigna.parser.CDNA2PublicParser.ParseResult;
import circdesigna.parser.CDNA2Token;

/**
 * A data structure for representing DNA complexes. This one is based on the circular polymer graph.
 * @author Benjamin
 */
public class DomainPolymerGraph extends AbstractComplex{
	public DomainPolymerGraph(DomainDefinitions domainDefs){
		this.domainDefs = domainDefs;
		data = new BasePolymerGraph();
	}	
	private DomainPolymerGraph(DomainDefinitions domainDefs, AbstractPolymerGraph aa){
		this.domainDefs = domainDefs;
		data = aa;
	}	
	public DomainPolymerGraph getRotation(String moleculeName, int rotation){
		if (rotation == 0){
			return this; //Preserves things like being a CanonicalRotation.
		}
		return new DomainPolymerGraph(domainDefs, new RotatedPolymerGraph(data,rotation));
	}
	private DomainDefinitions domainDefs;
	private AbstractPolymerGraph data;
	public String moleculeName;
	public int length(){
		return data.length();
	}
	public int getDomain(int position){
		return data.getDomain(position);
	}
	public int getDomainPair(int position){
		return data.getDomainPair(position);
	}
	public CircDesigNAStyle getStyle(int position){
		return data.getStyle(position);
	}
	public boolean setDomainPair(int i, int j){
		return data.setDomainPair(i,j);
	}
	/**
	 * Assigns a number to each base, numbering it inside a region in the planar nonpseudoknotted graph.
	 */
	public int[] getDomainLevels() {
		int[] toRet = new int[length()];
		int[] levels = new int[length()+1];
		int levels_ptr = 0;
		for(int i =0; i < toRet.length; i++){
			if (getDomainPair(i) >= 0){
				if (getDomainPair(i) > i){
					levels_ptr++;
					levels[levels_ptr] = i+1; //Ensure that they are distinct.
				}	
			}
			if (getDomain(i)<0){
				toRet[i] = -1;
			} else {
				toRet[i] = levels[levels_ptr];
			}
			if (getDomainPair(i) >= 0){
				if (getDomainPair(i) < i){
					levels_ptr--;
				}	
			}
		}
		return toRet;
	}
	public boolean equals(Object other){
		if (other instanceof DomainPolymerGraph){
			DomainPolymerGraph apg = (DomainPolymerGraph)other;
			if (apg.length() != length()){
				return false;
			}
			for(int i = 0; i < length(); i++){
				if (getDomain(i)!=apg.getDomain(i))
					return false;
				if (getDomainPair(i)!=apg.getDomainPair(i))
					return false;
			}
			return true;
		}
		return false;
	}
	public boolean equalsRotation(DomainPolymerGraph other){
		for(int i : getStrandRotations()){
			if (other.equals(getRotation("A", i))){
				return true;
			}
		}
		return false;
	}
	
	private abstract static class AbstractPolymerGraph {
		public abstract BasePolymerGraph getCyclicIndependentForm();
		public abstract int length();
		/**
		 * -1 means 3' end
		 */
		public abstract int getDomain(int position);
		/**
		 * -1 if not paired
		 */
		public abstract int getDomainPair(int position);
		
		public abstract CircDesigNAStyle getStyle(int position);
		
		public abstract boolean setDomainPair(int i, int j);
	}
	public class BasePolymerGraph extends AbstractPolymerGraph{
		//Indexed by position in the string.
		//Domain -1 is a 3' end.
		private int[] domains;
		private CircDesigNAStyle[] styles;
		//-1 means "not paired".
		private int[] domain_pairs;
		public BasePolymerGraph getCyclicIndependentForm() {
			return this;
		}
		public int getDomain(int position) {
			return domains[position];
		}
		public int getDomainPair(int position) {
			return domain_pairs[position];
		}
		public CircDesigNAStyle getStyle(int position) {
			return styles[position];
		}
		public int length() {
			return domains.length;
		}
		public boolean setDomainPair(int i, int j) {
			if (j == i){
				throw new RuntimeException("Cannot pair a domain with itself. "+i);
			}
			int oldTarget = domain_pairs[i];
			if (oldTarget>=0){
				domain_pairs[oldTarget] = -1;
				domain_pairs[i] = -1;
			}
			if (j>=0){
				int oldTargetOfTarget = domain_pairs[j];
				if (oldTargetOfTarget>=0){
					domain_pairs[oldTargetOfTarget] = -1;
					domain_pairs[j] = -1;
				}
				if (getLevel_(j)!=getLevel_(i)){
					if (oldTargetOfTarget>=0){
						domain_pairs[oldTargetOfTarget] = j;
						domain_pairs[j] = oldTargetOfTarget;
					}
					if (oldTarget>=0){
						domain_pairs[oldTarget] = i;
						domain_pairs[i] = oldTarget;
					}
					//throw new RuntimeException("Pairing "+i+" to "+j+" would cause a pseudoknot.");
					return false;
				}
				domain_pairs[j] = i;
			}
			domain_pairs[i] = j;
			return true;
		}
		private int getLevel_(int i) {
			int pointer = 0;
			int[] levels = new int[i+1];
			for(int k = 0; k < i; k++){
				if (domain_pairs[k] >= 0){
					if (domain_pairs[k] < k){
						pointer--;
					} else {
						pointer++;
						levels[pointer] = k+1;
					}	
				}
			}
			return levels[pointer];
		}
	}
	public class RotatedPolymerGraph extends BasePolymerGraph{
		private BasePolymerGraph base;
		private int rotation;
		public RotatedPolymerGraph(AbstractPolymerGraph base, int rotation){
			this.base = base.getCyclicIndependentForm();
			if (base instanceof RotatedPolymerGraph){
				rotation += ((RotatedPolymerGraph)base).rotation;
				rotation = wrap(rotation, base.length());
			}
			this.rotation = rotation;
		}
		public BasePolymerGraph getCyclicIndependentForm() {
			return base;
		}
		public int getDomain(int position) {
			return base.domains[wrap(position+rotation,base.domains.length)];
		}
		public CircDesigNAStyle getStyle(int position){
			return base.styles[wrap(position+rotation,base.domains.length)];
		}
		private int wrap(int i, int length) {
			while (i < 0){
				i += length;
			}
			i %= length;
			return i;
		}
		public int getDomainPair(int position) {
			int originalPair = base.domain_pairs[wrap(position+rotation,base.domains.length)];
			if (originalPair<0){
				return originalPair;
			}
			int newPair = wrap(originalPair-rotation,base.domains.length);
			return newPair;
		}
		public boolean setDomainPair(int i, int j){
			if (j!=-1){
				j = wrap(j+rotation,base.domains.length);
			}
			if (i < 0){
				return false;
			}
			i = wrap(i+rotation,base.domains.length);
			return super.setDomainPair(i,j);
		}
		public int length() {
			return base.domains.length;
		}
	}

	public static void readStructure(String moleculeDefinition, DomainPolymerGraph out2){
		
		BasePolymerGraph out = out2.data.getCyclicIndependentForm();
		out.domains = null;
		out.domain_pairs = null;

		//Tested and verified. Group 2 : Domain, with comp flag, Group 1 + Group 3: Structural flag
		ParseResult t = CDNA2PublicParser.parse(moleculeDefinition, out2.domainDefs);
		out2.moleculeName = t.moleculeName;
		
		ArrayList<Integer> domains_tmp = new ArrayList();
		ArrayList<CircDesigNAStyle> styles_tmp = new ArrayList();
		ArrayList<Integer> domain_pairs_tmp = new ArrayList();
		LinkedList<Integer> stack = new LinkedList();

		CircDesigNAStyle style = CircDesigNAStyle.getDefaultStyle();
		for(Object token : t.parse){
			if (token instanceof CDNA2Token.Option){
				{
					CircDesigNAStyle newStyle = new CircDesigNAStyle();
					newStyle.color = style.color;
					//Swap pointers
					style = newStyle;
				}
				
				CDNA2Token.Option d = (CDNA2Token.Option)token;
				try {
				    Field field = Color.class.getField(d.optionWord);
				    style.color = (Color)field.get(null);
				} catch (Exception e) {
				    // Not defined
				}
				try {
				    Field field = CircDesigNAStyle.Colors.class.getField(d.optionWord);
				    style.color = (Color)field.get(null);
				} catch (Exception e) {
				    // Not defined
				}
			}
			
			if (token instanceof CDNA2Token.Domain){
				CDNA2Token.Domain d = (CDNA2Token.Domain)token;
				String domainName = d.name;
				//Decode which domain, extract the "star" character (reverse complement flag)
				int domainNameL = domainName.length();
				if (domainName.endsWith("*")){
					domainNameL --;
				}
				int numberDomain = out2.domainDefs.lookupDomainName(domainName.substring(0,domainNameL));
				int numberDomain2 = numberDomain;
				if (numberDomain2 < 0){
					throw new RuntimeException("Invalid domain: "+domainName);
				}
				
				if (out2.domainDefs.domainLengths[numberDomain2] == 0){
					//Ignore domains of length 0 when parsing the structures.
				} else {
					if (domainName.endsWith("*")){
						numberDomain2 |= NA_COMPLEMENT_FLAG;
					}

					int thisIndex = domains_tmp.size();
					domains_tmp.add(numberDomain2);
					styles_tmp.add(style);
					domain_pairs_tmp.add(-1); //fill with -1s, overwrite if adding pair.
					if (d.open){
						stack.push(thisIndex);
					} else if (d.close){
						int openParens = stack.pop();
						domain_pairs_tmp.set(openParens,thisIndex);
						domain_pairs_tmp.set(thisIndex,openParens);
					} else {
						//Do nothing.
					}
				}
			} else if (token instanceof CDNA2Token.FivePrimeEnd){
				//Ignore.
			} else if (token instanceof CDNA2Token.ThreePrimeEnd){
				if (domains_tmp.isEmpty() || domains_tmp.get(domains_tmp.size()-1)==-1){
					//Empty strand.
				} else {
					//3- end.
					styles_tmp.add(style);
					domains_tmp.add(-1);
					domain_pairs_tmp.add(-1); //fill with -1s, overwrite if adding pair.
				}
			}
		}

		//This is dependent on the input sequence. Not on Domain Definitions.
		out.domains = new int[domains_tmp.size()];
		out.domain_pairs = new int[domains_tmp.size()];
		out.styles = new CircDesigNAStyle[domains_tmp.size()];
		for(int k = 0; k < out.domains.length; k++){
			out.domains[k] = domains_tmp.get(k);
			out.domain_pairs[k] = domain_pairs_tmp.get(k);
			out.styles[k] = styles_tmp.get(k);
		}
		
		//Erase rotation afterwards.
		out2.data = out;
		return;
	}		

	public String toString(){
		return "Polymer graph of "+moleculeName;
	}

	public String getStructureString(){
		return getStructureString(0,length());
	}
	public String getStructureString(int start, int end){
		return getStructureString(data, start, end);
	}
	public String getStructureString(AbstractPolymerGraph apg, int start, int end){
		StringBuffer sb = new StringBuffer();
		int depth = 0;
		boolean needsOpenBracket = true, needsCloseItem = false;
		for(int k = start; k < end; k++){
			if (needsOpenBracket){
				sb.append("[");
				needsOpenBracket = false;
			}
			int domain = apg.getDomain(k);
			if (domain>=0){
				//getDomainName does the works - adds an asterix, looks up the name, etc.
				sb.append(domainDefs.getDomainName(domain));
				int pair = apg.getDomainPair(k);
				if (pair >= 0){
					if (pair < k){
						//the close parens
						sb.append(")");
					} else {
						sb.append("(");	
					}
				} else {
					sb.append(" ");
				}
				needsCloseItem = true;
			} else {
				sb.append("}");
				needsOpenBracket = true;
				needsCloseItem = false;
			}
		}
		if (needsCloseItem){
			sb.append("}");
		}
		return sb.toString();
	}

	/**
	 * Returns the positions of all the 3' ends.
	 */
	public List<Integer> getStrandRotations() {
		ArrayList<Integer> toRet = new ArrayList<Integer>();
		//Length -1, because the last one doesn't count.
		for(int k = 0; k < length()-1; k++){
			if (getDomain(k)==-1){
				toRet.add(k+1);
			}
		}
		return toRet;
	}

	/**
	 * Returns an array of tuples (a,b) where the characters a,b are a strand in the molecule.
	 * Essentially partitions the molecule into strands.
	 */
	public int[][] getStrands() {
		List<Integer> rots = getStrandRotations();
		rots.add(0,0);
		int[][] strands = new int[rots.size()][];
		rots.add(length());
		for(int i = 0; i < strands.length; i++){
			strands[i] = new int[]{rots.get(i),rots.get(i+1)};
		}
		return strands;
	}
	
	public String getMoleculeName() {
		return moleculeName;
	}

	public DomainDefinitions getDomainDefs() {
		return domainDefs;
	}
		
	public CanonicalDomainPolymerGraph getCanonicalForm(){
		if (this instanceof CanonicalDomainPolymerGraph){
			return (CanonicalDomainPolymerGraph)this;
		}
		AbstractPolymerGraph selected = data;
		String best = getStructureString(selected, 0, length());
		List<Integer> rotations = getStrandRotations();
		for(int k : rotations){
			AbstractPolymerGraph test = new RotatedPolymerGraph(data,k);
			String got = getStructureString(test, 0, length());
			if (got.compareTo(best) < 0){
				best = got;
				selected = test;
			}
		}
		CanonicalDomainPolymerGraph toRet = new CanonicalDomainPolymerGraph(this, selected);
		{
			DomainPolymerGraph test = new DomainPolymerGraph(toRet.getDomainDefs());
			DomainPolymerGraph.readStructure("A " +toRet.getStructureString(), test);
		}
		return toRet;
	}
	
	public static class CanonicalDomainPolymerGraph extends DomainPolymerGraph {
		public CanonicalDomainPolymerGraph(DomainPolymerGraph g, AbstractPolymerGraph canonical) {
			super(g.domainDefs,canonical);
		}
	}

}
