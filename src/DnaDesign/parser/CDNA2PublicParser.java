package DnaDesign.parser;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;

import beaver.Parser.Exception;

public class CDNA2PublicParser {
	public static class Result {
		public String moleculeName;
		public ArrayList parse;
		public String toString(){
			StringBuffer sb = new StringBuffer();
			sb.append(moleculeName);
			sb.append(" ");
			for(Object k : parse){
				sb.append(k.toString()+" ");
			}
			return sb.toString();
		}
	}
	public static Result parse(String info){
		CDNA2Parser cp = new CDNA2Parser();
		Result res = new Result();
		try {
			ArrayList parse = (ArrayList) cp.parse(new CDNA2Scanner(new StringReader(info)));
			res.moleculeName = (String) parse.remove(0);
			res.parse = new ArrayList();
			res.parse.addAll(parse);
			
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		} catch (Exception e) {
			if (e.getCause().equals("")){
				throw new RuntimeException(e.getMessage());
			}
			throw new RuntimeException(e.getCause());
		}
		return res;
	}
}
