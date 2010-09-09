package DnaDesign;


/**
 * I needed a place to hold the various properties of DNA.
 */
public class DnaDefinition {
	//Explicit definition of DNA: Mapping from a subset of Z to all of DNA 
	public static final int A = 1, T = 2, G = 3, C = 4, D = 5, H = 6, P = 7, Z = 8;
	//A number guaranteed to be bigger than the largest base.
	public static final int DNAFLAG_ADD = 10;
	//Explicit definition of DNAxDNA -> Real . All other values are 0.
	public static final double GCstr = 3;
	public static final double ATstr = 2;
	public static final double DHstr = GCstr;
	public static final double PZstr = GCstr;
	public static final double GTstr = 0.1;
	//Binding score function.
	public static final double bindScore(int a, int b){
		a = noFlags(a);
		b = noFlags(b);
		if ((a==A && b==T) || (a==T && b==A)){
			return ATstr;
		} else
			if ((a==G && b==C) || (a==C && b==G)){
				return GCstr;
			} else
				if ((a==D && b==H) || (a==H && b==D)){
					return DHstr;
				} else 
					if ((a==P && b==Z) || (a==Z && b==P)){
						return DHstr;
					} else
						if ((a==G && b==T) || (a==T && b==G)){
							return GTstr;
						} else {
							return 0;
						}
	}
	public static final String DisplayBase(int base) {
		base = noFlags(base);
		switch(base){
		case G:
			return "G";
		case A:
			return "A";
		case C:
			return "C";
		case T:
			return "T";
		case D:
			return "D";
		case H:
			return "H";
		case P:
			return "P";
		case Z:
			return "Z";
		default:
			throw new IllegalArgumentException("Unrecognized Base: "+base);
		}
	}
	/**
	 * Returns a randomly chosen base, not equal to i, which is with probability 'f' an exotic base.
	 * @param maxPZ 
	 * @param mut_new 
	 */
	public static int pickNewBase(int[] mut_new, int i, int maxISO, int maxPZ) {
		if (i != noFlags(i)){
			throw new RuntimeException("Cannot mutate bases with flags: "+i);
		}
		if (maxISO==0 && maxPZ==0){
			
		} else {
			//Predefined: try exotic 1/10 of the time.
			if (Math.random()<.1f){
				//Select from EXOTIC
				int actualIso = 0, actualPZ = 0;
				boolean canISO = false, canPZ = false;
				if (actualIso < maxISO){
					actualIso = countDH(mut_new);
					if (actualIso < maxISO){
						canISO = true;
					}
				}
				if (actualPZ < maxPZ){
					actualPZ = countPZ(mut_new);
					if (actualPZ < maxPZ){
						canPZ = true;
					}
				}
				if (canISO || canPZ){
					double doesIso = Math.random()*(canISO?1:0);
					double doesPZ = Math.random()*(canPZ?1:0);
					//best
					int comp = (int)(Math.random()*2);
					if (doesPZ >= doesIso){
						if (P+comp==i){
							comp = (comp+1)%2;
						}
						return P+comp;
					} else {
						if (D+comp==i){
							comp = (comp+1)%2;
						}
						return D+comp;
					}
				}
			}
		}
		//Change to normal 
		return (int)((i-A)+(Math.random()*3+1))%4+A;
	}
	public static int countDH(int[] d){
		int sum = 0;
		for(int q : d){
			if (q==D || q==H){
				sum++;
			}
		}
		return sum;
	}
	public static int countPZ(int[] d){
		int sum = 0;
		for(int q : d){
			if (q==P || q==Z){
				sum++;
			}
		}
		return sum;
	}
	public static int complement(int i) {
		i = noFlags(i);
		if ((i&1)!=0){
			return i+1;
		}
		return i-1;
	}
	public static int noFlags(int i) {
		return i%DNAFLAG_ADD;
	}
}
