package interactome.analysis;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Collections;

import interactome.Logger;
import interactome.data.BioDB;
import interactome.data.Refseq;
import interactome.input.Input;

public class DynamicParameters {
	static DynamicParameters _instance;
	
	// constants
	final double minimum_coverage = 0.8;
	final int maximum_parameters_genes = 200;
	
	// instance fields
	Input input;
	Refseq[] parameterRefseqs;
	
	public static DynamicParameters getInstance() {
		return _instance;
	}
	
	public static DynamicParameters createInstance(Input input) {
		_instance = new DynamicParameters(input);
		return _instance;
	}
	
	private DynamicParameters(Input input) {
		this.input = input;
	}
	
	public void selectParameterRefseqs() {
		BioDB biodb = BioDB.getInstance();
		
		Logger.logf("\nselecting parameters genes.");
		
		Set<Refseq> candidate_refseqs = new HashSet<Refseq>();
		
		// accept refseq without a splice variant
		for (Map.Entry<String, Refseq> entry : biodb.refseq_db.entrySet()) {
			Refseq refseq = entry.getValue();
			if (refseq.gene != null && refseq.gene.variants.size() == 1) {
				candidate_refseqs.add(refseq);
			}
		}
		
		// accept refseq with adequate coverage
		for (Iterator<Refseq> i = candidate_refseqs.iterator(); i.hasNext(); ) {
			Refseq refseq = i.next();
			if (this.input.refseq_inputs.get(refseq.refseq_id).coverage(1) < minimum_coverage) {
				i.remove();
			}
		}
		
		// sort by length
		ArrayList<Refseq> candidate_array = new ArrayList<Refseq>(candidate_refseqs);
		Collections.sort(candidate_array, new Comparator<Refseq>() {
			public int compare(Refseq a, Refseq b) {
				return b.length - a.length;
			}
		});
		
		this.parameterRefseqs = new Refseq[Math.min(candidate_array.size(), maximum_parameters_genes)];
		for (int i=0; i<this.parameterRefseqs.length; i++) {
			this.parameterRefseqs[i] = candidate_array.get(i);
		}
		
		return;
	}
	
}
