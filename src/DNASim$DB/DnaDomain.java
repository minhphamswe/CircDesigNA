package DNASim$DB;

import TaiGameCore.GameDataBase;

public class DnaDomain extends GameDataBase{
	public DnaDomain(String hash) {
		super(hash);
		// TODO Auto-generated constructor stub
	}
	public static int readNA(char k, boolean flip){
		k = Character.toLowerCase(k);
		if (!flip){
			switch(k){
			case 'g':
				return G;
			case 'c':
				return C;
			case 'a':
				return A;
			case 't':
				return T;
			case 'u':
				return T;
			}
		} else { 
			switch(k){
			case 'g':
				return C;
			case 'c':
				return G;
			case 'a':
				return T;
			case 't':
				return A;
			case 'u':
				return A;
			}
		}
		throw new RuntimeException("Invalid na character: "+k);
	}
	public static int G = 0;
	public static int C = 1;
	public static int A = 2;
	public static int T = 3;
	/**
	 * Nucleoside list
	 */
	public int[] data;
	/**
	 * Metadata on nucleosides.
	 */
	public int[] dataMarking;

	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////
	///////////////////
	/////////////////

	public void autoWrittenDeSerializeCode(){
		data = ((IntArrayEntry)readField("data", new IntArrayEntry(new int[]{}))).getIntArray();
		dataMarking = ((IntArrayEntry)readField("dataMarking", new IntArrayEntry(new int[]{}))).getIntArray();
	}
	public void autoWrittenSerializeCode(){
		writeField("data", new IntArrayEntry(data));
		writeField("dataMarking", new IntArrayEntry(dataMarking));
	}
}
