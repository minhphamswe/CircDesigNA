package DnaDesign.AbstractDesigner;

public abstract class SingleMemberDesigner  <T extends PopulationDesignMember<T>> {
	public abstract double getOverallScore(T q);
	public abstract boolean mutateAndTestAndBackup(T q);
}
