package DNASim$DB;

import TaiGameCore.GameDataBase;

public class DnaStrand extends GameDataBase{
	public static void main(String[] args){
		//Probability that strand contains a 6 letter run.
		for(int k = 6; k < 300; k++){
			int num = 0;
			int p;
			for(p = 0; p < 100000; p++){
				int strCt = 0;
				for(int q = 0; q < k; q++){
					int got = (int)(Math.random()*4);
					if (got > 0){
						strCt = 0;
					} else {
						strCt++;
					}
					if(strCt>4){
						num++;
						break;
					}
				}
			}
			System.out.println(k+" "+num/(double)p);
		}
		/*
		for(int k = 6; k < 30; k++){
			System.out.print(k+" ");
			float ex = 0;
			float var = 0;
			for(int y = 6; y < k; y++){
				float instaEx = (float)(Math.pow(2/3.,y)*(k-y)); 
				ex += instaEx;
				var += u - instaEx*instaEx;
			}
			ex /= 3;
			var /= 9;
			System.out.println(ex+" "+var);
		}
		*/
	}
	public DnaStrand(String hash) {
		super(hash);
	}

	//In sequential order.
	public int[] domains;

	public int numCopies;
	
	@DefaultValue(value="3")
	public int numDimensions;
	
	public float[] positions;
	public float[] velocities;
	
	public int[] bindersStrand;
	public int[] bindersInd;
	
	public String name;
	
	
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
		domains = ((IntArrayEntry)readField("domains", new IntArrayEntry(new int[]{}))).getIntArray();
		numCopies = ((IntEntry)readField("numCopies", new IntEntry())).getInt();
		numDimensions = ((IntEntry)readField("numDimensions", new IntEntry(3))).getInt();
		positions = ((FloatArrayEntry)readField("positions", new FloatArrayEntry(new float[]{}))).getFloatArray();
		velocities = ((FloatArrayEntry)readField("velocities", new FloatArrayEntry(new float[]{}))).getFloatArray();
		bindersStrand = ((IntArrayEntry)readField("bindersStrand", new IntArrayEntry(new int[]{}))).getIntArray();
		bindersInd = ((IntArrayEntry)readField("bindersInd", new IntArrayEntry(new int[]{}))).getIntArray();
		name = ((StringEntry)readField("name", new StringEntry(""))).getString();
	}
	public void autoWrittenSerializeCode(){
		writeField("domains", new IntArrayEntry(domains));
		writeField("numCopies", new IntEntry(numCopies));
		writeField("numDimensions", new IntEntry(numDimensions));
		writeField("positions", new FloatArrayEntry(positions));
		writeField("velocities", new FloatArrayEntry(velocities));
		writeField("bindersStrand", new IntArrayEntry(bindersStrand));
		writeField("bindersInd", new IntArrayEntry(bindersInd));
		writeField("name", new StringEntry(name));
	}

}
