package interactome.input;

import interactome.data.Refseq;

public class RefseqInput {
	public Refseq refseq;
	
	public long[] overwrap_counts;
	public long[] starting_counts;
	
	public double coverage(int depth) {
		int covered = 0;
		for (int i=0; i<this.refseq.length; i++) {
			if (this.overwrap_counts[i] >= depth) covered++;
		}
		
		return (double)covered / this.refseq.length;
	}
	
	public RefseqInput(Refseq refseq) {
		this.refseq = refseq;
		if (!refseq.is_invalid) {
			this.overwrap_counts = new long[refseq.length];
			this.starting_counts = new long[refseq.length];
		}
	}
}
