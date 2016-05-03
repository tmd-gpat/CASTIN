package interactome.input;

import interactome.data.Refseq;

public class RefseqInput {
	public Refseq refseq;
	
	public long rawCount = 0;
	public long[] overlap_counts;
	public long[] starting_counts;

	public double true_expression;
	
	public double coverage(int depth) {
		int covered = 0;
		for (int i=0; i<this.refseq.length; i++) {
			if (this.overlap_counts[i] >= depth) covered++;
		}
		
		return (double)covered / this.refseq.length;
	}
	
	public RefseqInput(Refseq refseq) {
		this.refseq = refseq;
		if (!refseq.is_invalid) {
			this.overlap_counts = new long[refseq.length];
			this.starting_counts = new long[refseq.length];
		}
	}
}
