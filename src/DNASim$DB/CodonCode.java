package DNASim$DB;

public class CodonCode {
	/**
	 * Returns which amino acid is being encoded at [i,i+2]
	 */
	public char decodeAmino(int[] out, int i){
		int result = out[i]*64+out[i+1]*8+out[i+2];
		return forwardTable[result];
	}
	/**
	 * Mutates the codon at [i,i+2] to a different codon of the same amino acid.
	 */
	public void mutateToOther(int[] out, int i){
		//Code duplicated from decodeAmino
		int pow8 = out[i]*64+out[i+1]*8+out[i+2];
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
