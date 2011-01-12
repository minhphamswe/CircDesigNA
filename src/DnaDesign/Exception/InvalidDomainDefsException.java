package DnaDesign.Exception;

/**
 * Exception due to the syntax of the DomainDefs block.
 */
public class InvalidDomainDefsException extends RuntimeException{
	public InvalidDomainDefsException(String string) {
		super(string);
	}
}
