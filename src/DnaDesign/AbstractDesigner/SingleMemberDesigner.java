package DnaDesign.AbstractDesigner;

public abstract class SingleMemberDesigner  <T extends PopulationDesignMember<T>> {
	public abstract double getOverallScore(T q);
	public abstract boolean mutateAndTestAndBackup(T q);
	public abstract boolean mutateAndTest(T original, T into);
	public abstract boolean fourPtCrossoverAndTest(T a, T b, T into);
}
