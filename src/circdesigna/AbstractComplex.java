/*
  Part of the CircDesigNA Project - http://cssb.utexas.edu/circdesigna
  
  Copyright (c) 2010-11 Ben Braun
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation, version 2.1.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General
  Public License along with this library; if not, write to the
  Free Software Foundation, Inc., 59 Temple Place, Suite 330,
  Boston, MA  02111-1307  USA
*/
package circdesigna;

import java.util.ArrayList;
import java.util.Collection;

/**
 * An abstraction of a molecular complex. It combines a structural representation with a reference
 * to a Domain Definition set necessary for comprehending the structure. 
 */
public abstract class AbstractComplex {
	public abstract String getMoleculeName();
	public abstract DomainDefinitions getDomainDefs();
	/**
	 * Should return a string equivalent to one that was parsed to create this abstract complex.
	 */
	public abstract String getStructureString();
	
	public void annotate(String message){
		annotation.add(message);	
	}
	public void annotate(Annotation k){
		annotation.add(k);
	}
	public Annotation getAnnotationTree(){
		return annotation;
	}
	private Annotation annotation = new Annotation();

	public Collection<Annotation> getAnnotationLevel(String match) {
		Collection<Annotation> collect = new ArrayList();
		getAnnotationLevel(collect, annotation, match);
		return collect;
	}
	private void getAnnotationLevel(Collection<Annotation> collect, Annotation a, String match) {
		for(String q : a.annotations){
			if (q.matches(match)){
				collect.add(a);
				return; //!!!!!! We do not go past the first level where a match is found.
			}
		}
		for(Annotation q : a.others){
			getAnnotationLevel(collect, q, match);
		}
	}
	public static class Annotation{
		public ArrayList<Annotation> others = new ArrayList();
		public ArrayList<String> annotations = new ArrayList();
		
		public void clear(){
			others.clear();
			annotations.clear();
		}
		
		public Annotation add(Annotation other){
			if (!annotations.isEmpty()){
				Annotation mine = new Annotation();
				mine.add(this);
				mine.add(other);
				return mine;
			} else {
				others.add(other);
				return this;
			}
		}
		public Annotation add(String annotation){
			if (others.size() > 1){
				Annotation mine = new Annotation();
				mine.add(this);
				mine.add(annotation);
				return mine;
			} else {
				annotations.add(annotation);
				return this;
			}
		}
		private void print(int depth){
			for(Annotation o : others){
				o.print(depth+1);
			}
			for(String k : annotations){
				for(int u = 0; u < depth; u++){
					System.out.print(">");
				}
				System.out.println(k);
			}
		}
		public void print() {
			print(0);
		}
	}
}
