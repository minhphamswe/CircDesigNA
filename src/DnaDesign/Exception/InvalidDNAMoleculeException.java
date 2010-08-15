package DnaDesign.Exception;

public class InvalidDNAMoleculeException extends RuntimeException{
	/**
	 * A runtime exception due to a certain input molecule
	 */
	public final int myStrand;
	public InvalidDNAMoleculeException(String string, int kq) {
		super(string);
		myStrand = kq;
	}
}
