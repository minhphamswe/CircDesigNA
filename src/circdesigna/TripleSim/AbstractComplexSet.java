package circdesigna.TripleSim;

import java.util.ArrayList;

import circdesigna.AbstractComplex;

public class AbstractComplexSet <E extends AbstractComplex> extends ArrayList<E>{
	public boolean add(E c){
		for(E q : this){
			if (q.equals(c)){
				q.annotate(c.getAnnotationTree());
				return true;
			}
		}
		super.add(c);
		return true;
	}
}
