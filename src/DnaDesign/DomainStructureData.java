package DnaDesign;
import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DomainSequence.DNA_COMPLEMENT_FLAG;
import static DnaDesign.DomainSequence.DNA_SEQ_FLAGSINVERSE;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * The information parsed from the DomainDefs block is encoded in an object of this class.
 */
public class DomainStructureData {
	/**
	 * Invertible.
	 */
	public Map<String, Integer> nameMap = new TreeMap();
	//keys are domains converted to numbers.
	public String getDomainName(int domain){
		String postpend = ((domain & DNA_COMPLEMENT_FLAG)!=0?"*":"");
		for(Entry<String, Integer> q : nameMap.entrySet()){
			if (q.getValue()==(domain & DNA_SEQ_FLAGSINVERSE)){
				return q.getKey()+postpend;
			}
		}
		return ""+(domain & DNA_SEQ_FLAGSINVERSE)+postpend;
	}
	public int lookupDomainName(String substring) {
		if (!nameMap.containsKey(substring)){
			return -1;
		}
		return nameMap.get(substring);
	}
	public int[] domainLengths;
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

	
	private static int parsePercent(String percent, int of, boolean roundUp){
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
