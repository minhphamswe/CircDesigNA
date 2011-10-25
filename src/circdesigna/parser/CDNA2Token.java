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
package circdesigna.parser;

import static circdesigna.CircDesigNA_SharedUtils.isComplements;

import java.util.ArrayList;

import circdesigna.DomainDefinitions;

public class CDNA2Token {
	public static class Domain {
		public boolean open, close, ss;
		public String name;
		public Domain pair;
		public Domain(String name){
			this.name = name;
		}
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(name);
			if (open) sb.append("(");
			if (close) sb.append(")");

			return sb.toString();
		}
		public void setOpen(ArrayList parenStack){
			if (open){
				throw new RuntimeException(String.format("Double open parenthesis on %s",this.toString()));
			}
			open = true;
			if (ss || close){
				throw new RuntimeException(String.format("Conflicting structure on %s",this.toString()));
			}
			parenStack.add(this);
		}
		public void setClosed(ArrayList parenStack) {
			if (parenStack.isEmpty()){
				throw new RuntimeException(String.format("Stack underflow on %s",this.toString()));
			}
			pair = (Domain) parenStack.remove(parenStack.size()-1);
			if (close){
				throw new RuntimeException(String.format("Double close parenthesis on %s",this.toString()));
			}
			close = true;
			if (ss || open){
				throw new RuntimeException(String.format("Conflicting structure on %s",this.toString()));
			}
		}
		public void setSingleStranded() {
			if (close || open){
				throw new RuntimeException(String.format("Conflicting structure on %s",this.toString()));
			}
			if (ss){
				throw new RuntimeException(String.format("Double open parenthesis on %s",this.toString()));
			}
			ss = true;
		}
		public void setComplement(){
			if (name.endsWith("*")){
				throw new RuntimeException(String.format("Double complement (*) on %s",this.toString()));
			}
			name+="*";
		}
		public void validate(DomainDefinitions domains){
			if (close){
				boolean validPair = false;
				if (pair.name.endsWith("*")){
					if (pair.name.equals(name+"*")){
						validPair = true;
					}
				} else {
					if ((pair.name+"*").equals(name)){
						validPair = true;
					}
				}
				if (!validPair){
					//Are they reverse complementary?
					boolean areRevCompFixed = false;
					try {
						int i = domains.lookupDomainName(name);
						int j = domains.lookupDomainName(pair.name);
						areRevCompFixed = isComplements(i,j,domains);
					} catch (Throwable e){
						areRevCompFixed = false;
					}
					if (!areRevCompFixed){
						throw new RuntimeException(String.format("Not complementary: %s and %s",name,pair.name));
					}
				}
			}
		}
	}
	public static class FivePrimeEnd{
		public FivePrimeEnd(ArrayList braceStack) {
			if (!braceStack.isEmpty()){
				throw new RuntimeException("Nesting [ brackets.");
			}
			braceStack.add(this);
		}

		public String toString(){
			return "[";
		}
	}
	public static class ThreePrimeEnd{
		public ThreePrimeEnd(ArrayList braceStack) {
			if (braceStack.isEmpty()){
				throw new RuntimeException("Stack underflow on }. Remember, correct molecule format is <name> <molecule>");
			}
			braceStack.clear();
		}

		public String toString(){
			return "}";
		}
	}
}
