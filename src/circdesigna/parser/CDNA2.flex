package circdesigna.parser;

import beaver.Symbol;
import beaver.Scanner;

import circdesigna.parser.CDNA2Parser.Terminals;

%%

%class CDNA2Scanner
%extends Scanner
%function nextToken
%type Symbol
%yylexthrow Scanner.Exception
%eofval{
	return newToken(Terminals.EOF, "end-of-file");
%eofval}
%unicode
%line
%column
%{
	private Symbol newToken(short id)
	{
		return new Symbol(id, yyline + 1, yycolumn + 1, yylength());
	}

	private Symbol newToken(short id, Object value)
	{
		return new Symbol(id, yyline + 1, yycolumn + 1, yylength(), value);
	}
%}
LineTerminator = \r|\n|\r\n
WhiteSpace     = {LineTerminator} | [ \t\f]

DomainName = [:jletterdigit:] [:jletterdigit:]*

%%

{WhiteSpace}+   { /* ignore */ }

<YYINITIAL> {
	{DomainName}    { return newToken(Terminals.DOMAINNAME, yytext()); }

	"("         { return newToken(Terminals.LPAREN, yytext()); }
	")"         { return newToken(Terminals.RPAREN, yytext()); }
	"*"         { return newToken(Terminals.MULT,   yytext()); }
	"["			{ return newToken(Terminals.LSQBRACE,   yytext()); }
	"}"			{ return newToken(Terminals.RCBRACE,   yytext()); }
	"."			{ return newToken(Terminals.DOT,   yytext()); }
}

.|\n            { throw new Scanner.Exception("unexpected character '" + yytext() + "'"); }
