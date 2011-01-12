package DnaDesign.Exception;

/**
 * A runtime exception due to invalidity of a certain input molecule
 */
public class InvalidDNAMoleculeException extends RuntimeException{
	public final int myStrand;
	public InvalidDNAMoleculeException(String string, int kq) {
		super(string);
		myStrand = kq;
	}
}
