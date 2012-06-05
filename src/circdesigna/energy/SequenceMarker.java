package circdesigna.energy;

import circdesigna.DomainSequence;

public class SequenceMarker {
	private int[][] seq_origin;
	private int[][] domain_markings;
	private int N;
	public SequenceMarker(int N, int[][] seq_origin, int[][] domain_markings){
		this.N = N;
		this.seq_origin = seq_origin;
		this.domain_markings = domain_markings;
	}
	public void mark(int i){
		if (i < 0 || i >= N){
			throw new ArrayIndexOutOfBoundsException();
		}
		int whichDomain = seq_origin[i][0] & DomainSequence.NA_COMPLEMENT_FLAGINV;
		int whichBase = seq_origin[i][1]; //Index is into uncomplemented form.
		domain_markings[whichDomain][whichBase] ++;
	}
}
