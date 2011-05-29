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

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;

import beaver.Parser.Exception;

public class CDNA2PublicParser {
	public static class ParseResult {
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
	public static ParseResult parse(String info){
		CDNA2Parser cp = new CDNA2Parser();
		ParseResult res = new ParseResult();
		try {
			ArrayList parse = (ArrayList) cp.parse(new CDNA2Scanner(new StringReader(info)));
			res.moleculeName = (String) parse.remove(0);
			for(Object o : parse){
				if (o instanceof CDNA2Token.Domain){
					((CDNA2Token.Domain)o).validate();
				}
			}
			
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
