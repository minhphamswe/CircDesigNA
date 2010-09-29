package DnaDesign.AbstractDesigner;

public abstract class PopulationDesignMember<T extends PopulationDesignMember> implements Comparable<T> {
	//Population design members are sorted by number in the population.
	protected int myID = 0;
	public final int compareTo(T o) {
		return myID - o.myID;
	}
	protected final PopulationDesignMember<T> designerCopyConstructor(int myID){
		PopulationDesignMember<T> toRet = designerCopyConstructor();
		toRet.myID = myID;
		toRet.seedFromOther(this);
		return toRet;
	}
	/**
	 * Creates a new instance of this design member. A deep copy is not required, as
	 * this call will always be followed up with a call to "seed".
	 */
	public abstract PopulationDesignMember designerCopyConstructor();
	public abstract void seedFromOther(PopulationDesignMember pdm);
}
