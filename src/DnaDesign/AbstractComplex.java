package DnaDesign;

/**
 * An abstraction of a molecular complex. It combines a structural representation with a reference
 * to a Domain Definition set necessary for comprehending the structure. 
 */
public interface AbstractComplex {
	public String getMoleculeName();
	public DomainDefinitions getDomainDefs();
	/**
	 * Should return the same string that was parsed to create this abstract complex.
	 */
	public String getStructureString();
}
