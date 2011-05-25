package DnaDesign.parser;

import java.util.ArrayList;

public class CDNA2Token {
	public static class Domain {
		public boolean open, close, ss;
		public String name;
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
			Domain pair = (Domain) parenStack.remove(parenStack.size()-1);
			if (close){
				throw new RuntimeException(String.format("Double close parenthesis on %s",this.toString()));
			}
			close = true;
			if (ss || open){
				throw new RuntimeException(String.format("Conflicting structure on %s",this.toString()));
			}
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
				throw new RuntimeException(String.format("Not complementary: %s and %s",name,pair.name));
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
				throw new RuntimeException("Stack underflow on }.");
			}
			braceStack.clear();
		}

		public String toString(){
			return "}";
		}
	}
}
