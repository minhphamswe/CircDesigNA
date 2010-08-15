package DnaDesign;

public class CodonCode {
	/**
	 * Returns which amino acid is being encoded at [i,i+2]
	 */
	public char decodeAmino(int[] out, int i){
		int result = pow8(out,i);
		return forwardTable[result];
	}
	/**
	 * Mutates the codons [i,i+2] and [j,j+2] to different codons, both of the original amino acid, but
	 * with complementarity reduced. 
	 */
	public void mutateToFixComplement(int[] is, int i, int j) {
		final int pow8i = pow8(is,i);
		final int pow8j = pow8(is,j);
		int allowableComplement = 0;
		if (getComplementNum(pow8i,pow8j) <= allowableComplement){
			return; //Nothing to change!
		}
		char aminoI = forwardTable[pow8i];
		char aminoJ = forwardTable[pow8j];
		int[] possibleI = reverseTable[aminoI];
		int[] possibleJ = reverseTable[aminoJ];
		if (possibleI.length*possibleJ.length==1){
			//throw new RuntimeException("Amino "+amino+" has just one codon!");
			return; //Unfortunately, this is all we can do...
		}
		boolean fixed = false;
		int bestScore = 50;
		int bestI = -1, bestJ = -1;
		int rOffI = (int) (Math.random()*possibleI.length);
		int rOffj = (int) (Math.random()*possibleJ.length);
		biggloop: for(int tryPossibleIr = 0; tryPossibleIr < possibleI.length; tryPossibleIr++){
			for(int tryPossibleJr = 0; tryPossibleJr < possibleJ.length; tryPossibleJr++){
				int tryPossibleI = (tryPossibleIr + rOffI)%possibleI.length;
				int tryPossibleJ = (tryPossibleIr + rOffj)%possibleJ.length;
				/*
				if (possibleI[tryPossibleI] == pow8i && possibleJ[tryPossibleJ] == pow8j){
					if (Math.random()<.5){
						tryPossibleI = (tryPossibleI+1)%possibleI.length;
					} else {
						tryPossibleJ = (tryPossibleJ+1)%possibleJ.length;
					}
				}
				*/
				int newI = possibleI[tryPossibleI];
				int newJ = possibleJ[tryPossibleJ];
				int score = getComplementNum(newI, newJ);
				if (score <= allowableComplement){
					bestI = newI;
					bestJ = newJ;
					fixed = true;
					break biggloop;
				}
				if (score <= bestScore){
					bestScore = score;
					bestI = newI;
					bestJ = newJ;
				}
			}
		}
		if (!fixed){
			//System.out.println(bestScore);
			//System.err.println("Couldn't fix "+aminoI+" "+aminoJ);
		}

		is[i] = (bestI>>6)&7;
		is[i+1] = (bestI>>3)&7;
		is[i+2] = (bestI)&7;

		is[j] = (bestJ>>6)&7;
		is[j+1] = (bestJ>>3)&7;
		is[j+2] = (bestJ)&7;
	}
	public static void main(String[] args){
		CodonCode c = new CodonCode();
		System.out.println(getComplementNum(c.reverseTable['A'][0],c.reverseTable['M'][0]));
	}
	private static final int getComplementNum(int pow81, int pow82){
		int comp = 0;
		for(int offset = -2; offset <= 2; offset++){
			int subComp = 0;
			for(int y = 0; y < 3; y++){
				int ny = 2-y+offset;
				if (ny < 0 || ny >= 3){
					continue;
				}
				if (((pow81>>(y*3))&7) == 5 - ((pow82>>(ny*3))&7)){
					subComp++;
				}
			}
			comp = Math.max(comp,subComp);
		}
		return comp;
	}
	private static final int pow8(int[] domain, int codonStart){
		return domain[codonStart]*64+domain[codonStart+1]*8+domain[codonStart+2];
	}
	/**
	 * Mutates the codon at [i,i+2] to a different codon of the same amino acid.
	 */
	public void mutateToOther(int[] out, int i){
		//Code duplicated from decodeAmino
		int pow8 = pow8(out,i);
		char amino = forwardTable[pow8];
		
		//k.
		int[] possible = reverseTable[amino];
		if (possible.length==1){
			//throw new RuntimeException("Amino "+amino+" has just one codon!");
			return; //Unfortunately, this is all we can do...
		}
		int which = (int)(Math.random()*possible.length);
		int newPow8 = possible[which];
		//Guarantee we picked something different.
		if (pow8==newPow8){
			which = (which+1)%possible.length;
			newPow8 = possible[which];
		}
		
		//Ok, dedicate newPow8
		out[i] = (newPow8>>6)&7;
		out[i+1] = (newPow8>>3)&7;
		out[i+2] = (newPow8)&7;
	}
	private char[] forwardTable = new char[512];
	private int[][] reverseTable = new int['Z'][0];
	private void writeCodon(int power8Rep, char base){
		int[] newReverseTable = new int[reverseTable[base].length+1];
		System.arraycopy(reverseTable[base],0,newReverseTable,0,newReverseTable.length-1);
		newReverseTable[newReverseTable.length-1] = power8Rep;
		reverseTable[base] = newReverseTable;
		//Forward
		forwardTable[power8Rep] = base;
	}
	{
		writeCodon(98,'A');
		writeCodon(100,'A');
		writeCodon(97,'A');
		writeCodon(99,'A');
		writeCodon(203,'C');
		writeCodon(204,'C');
		writeCodon(84,'D');
		writeCodon(83,'D');
		writeCodon(82,'E');
		writeCodon(81,'E');
		writeCodon(220,'F');
		writeCodon(219,'F');
		writeCodon(74,'G');
		writeCodon(76,'G');
		writeCodon(73,'G');
		writeCodon(75,'G');
		writeCodon(276,'H');
		writeCodon(275,'H');
		writeCodon(210,'*');
		writeCodon(209,'*');
		writeCodon(202,'*');
		writeCodon(154,'I');
		writeCodon(156,'I');
		writeCodon(155,'I');
		writeCodon(146,'K');
		writeCodon(145,'K');
		writeCodon(282,'L');
		writeCodon(284,'L');
		writeCodon(281,'L');
		writeCodon(283,'L');
		writeCodon(218,'L');
		writeCodon(217,'L');
		writeCodon(153,'M');
		writeCodon(148,'N');
		writeCodon(147,'N');
		writeCodon(290,'P');
		writeCodon(292,'P');
		writeCodon(289,'P');
		writeCodon(291,'P');
		writeCodon(274,'Q');
		writeCodon(273,'Q');
		writeCodon(138,'R');
		writeCodon(266,'R');
		writeCodon(137,'R');
		writeCodon(268,'R');
		writeCodon(265,'R');
		writeCodon(267,'R');
		writeCodon(140,'S');
		writeCodon(226,'S');
		writeCodon(228,'S');
		writeCodon(139,'S');
		writeCodon(225,'S');
		writeCodon(227,'S');
		writeCodon(162,'T');
		writeCodon(164,'T');
		writeCodon(161,'T');
		writeCodon(163,'T');
		writeCodon(90,'V');
		writeCodon(92,'V');
		writeCodon(89,'V');
		writeCodon(91,'V');
		writeCodon(201,'W');
		writeCodon(212,'Y');
		writeCodon(211,'Y');
	}
}
