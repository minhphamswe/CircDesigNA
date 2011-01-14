package DnaDesign;

public interface AbstractComplex {
	public String getMoleculeName();
	public DomainStructureData getDomainDefs();
	/**
	 * Should return the same string that was parsed to create this abstract complex.
	 */
	public String getStructureString();
}
