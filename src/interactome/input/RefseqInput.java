package interactome.input;

import interactome.data.Refseq;

public class RefseqInput {
	public Refseq refseq;
	
	public long rawCount = 0;
	public long[] overlap_counts;
	public long[] starting_counts;

	public double true_expression;
	
	public double v_init; // initial value of correction
	public double v_end; // ending value of correction
	public float[] residue_init;
	public float[] residue_end;
	
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

			this.residue_init = new float[refseq.length];
			this.residue_end = new float[refseq.length];
		}
	}
}
