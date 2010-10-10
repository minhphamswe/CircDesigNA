package DnaDesign;

public interface DesignerCode {
	public boolean mutateToOther(int[][] domain, int whichDomain, int i);
	public boolean isValid(int[][] domain, int whichDomain);
}
