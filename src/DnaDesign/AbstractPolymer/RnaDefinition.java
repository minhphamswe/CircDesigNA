package DnaDesign.AbstractPolymer;


/**
 * This class defines the syntax / encoding of RNA
 */
public class RnaDefinition extends DnaDefinition{
	//Explicit definition of DNA: Mapping from a subset of Z to all of DNA
	public static final int A = DnaDefinition.A, U = DnaDefinition.T, 
		G = DnaDefinition.G, C = DnaDefinition.C, 
		//D = DnaDefinition.D, H = DnaDefinition.H, 
		P = DnaDefinition.P, Z = DnaDefinition.Z;

	public int decodeBaseChar(char charAt) {
		charAt = Character.toUpperCase(charAt);
		if (charAt=='U' || charAt=='T'){
			return U;
		}
		return super.decodeBaseChar(charAt);
	}
	public String displayBase(int base) {
		if (base==U){
			return "U";
		} else {
			return super.displayBase(base);
		}
	}
}
