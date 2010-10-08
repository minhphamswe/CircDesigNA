package DnaDesign;

import static DnaDesign.DesignSequenceConstraints.LOCK_FLAG;
import static DnaDesign.DnaDefinition.A;
import static DnaDesign.DnaDefinition.C;
import static DnaDesign.DnaDefinition.D;
import static DnaDesign.DnaDefinition.DNAFLAG_ADD;
import static DnaDesign.DnaDefinition.G;
import static DnaDesign.DnaDefinition.H;
import static DnaDesign.DnaDefinition.P;
import static DnaDesign.DnaDefinition.T;


/**
 * I needed a place to hold the various properties of DNA.
 */
public class DnaDefinition {
	//Explicit definition of DNA: Mapping from a subset of Z to all of DNA 
	public static final int A = 1, T = 2, G = 3, C = 4, D = 5, H = 6, P = 7, Z = 8;
	//A number guaranteed to be bigger than the largest base.
	public static final int DNAFLAG_ADD = Z+2;
	//Blech. to get rid of.
	public static final double GCstr = 2;
	public static final double ATstr = 1;
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
	 * See getNormalBase, but the number returned will be in [0,3].
	 */
	public static int getNormalBaseFromZero(int nonnormalBase) {
		return getNormalBase(nonnormalBase)-1;
	}
	/**
	 * Returns a base from {A,C,T,G} which is most like 'x'.
	 */
	public static int getNormalBase(int x){
		switch(x){
		case G:
			return G;
		case A:
			return A;
		case C:
			return C;
		case T:
			return T;
		case D:
			return C;
		case H:
			return G;
		case P:
			return C;
		case Z:
			return G;
		default:
			throw new IllegalArgumentException("Unrecognized Base: "+x);
		}
	} 

	public static int decodeBaseChar(char charAt){
		charAt = Character.toUpperCase(charAt);
		switch(charAt){
		case 'G':
			return G;
		case 'A':
			return A;
		case 'C':
			return C;
		case 'T':
			return T;
		case 'D':
			return D;
		case 'H':
			return H;
		case 'P':
			return P;
		case 'Z':
			return Z;
		default:
			throw new RuntimeException("Invalid char: "+charAt);
		}
	}

	/**
	 * Returns the unflagged complement of base i
	 */
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
