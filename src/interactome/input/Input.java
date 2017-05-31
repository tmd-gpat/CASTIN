package interactome.input;

import interactome.data.BioDB;
import interactome.data.Gene;
import interactome.data.Refseq;

import java.util.HashMap;

public abstract class Input {
	public HashMap<String, RefseqInput> refseq_inputs;
	public HashMap<String, GeneInput> gene_inputs;
	
	// statistics of match-length
	public long[] match_length_distribution = new long[1001];
	
	public Input() {
		BioDB biodb = BioDB.getInstance();
		
		refseq_inputs = new HashMap<String, RefseqInput>();
		for (String refseq_id : biodb.all_refseq_ids) {
			Refseq refseq = biodb.refseq_db.get(refseq_id);
			if (!refseq.is_invalid)
				refseq_inputs.put(refseq_id, new RefseqInput(refseq));
		}
		
		gene_inputs = new HashMap<String, GeneInput>();
		for (String entrez_id : biodb.all_entrez_ids) {
			Gene gene = biodb.gene_db.get(entrez_id);
			gene_inputs.put(entrez_id, new GeneInput(gene));
		}
		
		for (int i=0; i<match_length_distribution.length; i++) match_length_distribution[i] = 0;
	}
	
	// abstract method for loading sequence files.
	public abstract boolean loadFile();	
}
