package interactome.analysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Collections;

import interactome.Logger;
import interactome.Option;
import interactome.Pair;
import interactome.data.BioDB;
import interactome.data.Refseq;
import interactome.input.Input;

public class DynamicParameters {
	static DynamicParameters _instance;
	
	// constants
	final double minimum_coverage = 0.8;
	final int minimum_coverage_depth = 1;
	final int maximum_parameters_genes = 200;
	final int minimum_length = 1000;
	
	final int mappability_min_match_length = 50;
	final int mappability_max_match_length = 250;
	
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
		
		// accept refseq with adequate coverage and length
		for (Iterator<Refseq> i = candidate_refseqs.iterator(); i.hasNext(); ) {
			Refseq refseq = i.next();
			if (this.input.refseq_inputs.get(refseq.refseq_id).coverage(minimum_coverage_depth) < minimum_coverage ||
				refseq.length < minimum_length) {
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
		Logger.logf("parameter estimation genes:");
		for (int i=0; i<this.parameterRefseqs.length && i<maximum_parameters_genes; i++) {
			Logger.logf(candidate_array.get(i).refseq_id);
			this.parameterRefseqs[i] = candidate_array.get(i);
		}
		
		return;
	}
	
	public boolean calculateWeightedMappability() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
		
		Logger.logf("\ncalculating weighted mappabilities.");
		
		// load self-matching table
		try {
			FileReader fr = new FileReader(
				option.directional_mode == 0 ? option.settings.get("self_matching") : option.settings.get("self_matching_direct")
						);
			BufferedReader br = new BufferedReader(fr);
			
			String l;
			Refseq refseq = null;
			while ((l = br.readLine()) != null) {
				if (l.charAt(0) == '>') {
					refseq = biodb.refseq_db.get(l.substring(1));
				} else {
					String[] row = l.split("\t");
					// self-matching entry is given by 0-order
					if (refseq != null) {
						refseq.unmappable_regions.add(new Pair<Integer, Integer>(Integer.valueOf(row[1]), Integer.valueOf(row[2])));
					}
				}
			}
			
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		// calculate match-length-ratio
		double[] ratios = new double[mappability_max_match_length+1];
		long dist_sum = 0;
//		Logger.logf("match length distribution:");
		for (int i=mappability_min_match_length; i<=mappability_max_match_length; i++) {
			dist_sum += input.match_length_distribution[i];
//			Logger.logf("%d: %d", i, input.match_length_distribution[i]);
		}
		for (int i=mappability_min_match_length; i<=mappability_max_match_length; i++) {
			ratios[i] = (double)input.match_length_distribution[i] / dist_sum;
//			Logger.logf("%f", ratios[i]);
		}
		
		// construct weight-buffer
		double[][][] weight_buffer = new double[mappability_max_match_length+1][mappability_max_match_length+1][mappability_max_match_length];
		for (int mappable_min=mappability_min_match_length; mappable_min<=mappability_max_match_length; mappable_min++) {
			for (int mappable_max=mappable_min; mappable_max<=mappability_max_match_length; mappable_max++) {
				for (int match_length=mappable_min; match_length<=mappable_max; match_length++) {
					for (int pos=0; pos<match_length; pos++) {
						weight_buffer[mappable_min][mappable_max][pos] += ratios[match_length];
					}
				}
			}
		}
		
		int n = 0;
		for (String refseq_id : biodb.all_refseq_ids) {
			if ((n++) % 10000 == 0) {
				Logger.logf("%d", n-1);
			}
			
			Refseq refseq = biodb.refseq_db.get(refseq_id);
			refseq.mappability = new float[refseq.length];
			for (int i=0; i<refseq.length; i++) refseq.mappability[i] = 0;
			
			// for each start position of a matching
			for (int p=0; p<refseq.length; p++) {
				// check the minimum- and maximum- unique matching length
				int min_mappable = mappability_min_match_length;
				for (Pair<Integer, Integer> region : refseq.unmappable_regions) {
					if (region.first <= p) {
						min_mappable = Math.max(region.second - p + 2, min_mappable);
					}
				}
				int max_mappable = Math.min(mappability_max_match_length, refseq.length-p);
				if (min_mappable > max_mappable) continue;
				
				// e.g., for a read-starting position p,
				// min_mappable = 50, max_mappable = 51,
				// ratios = 90% @50bp, 5% @51bp, 5% for others (global distribution)
				// mappability[p -- p+49] get (0.9 + 0.05) from p
				// mappability[p+50] get (0 + 0.05) from p
				
				// O(N) using weight-buffer
				for (int l=0; l<max_mappable; l++) {
					refseq.mappability[p+l] += weight_buffer[min_mappable][max_mappable][l];
				}
			}
		}
		
		Logger.logf("\ndone.");
		
		return true;
	}

	public void checkLongUnmappableRefseqs() {
		BioDB biodb = BioDB.getInstance();
		
		// check long (>=10000bp), unmappable (<50%) refseqs.
		for (String refseq_id : biodb.all_refseq_ids) {
			Refseq refseq = biodb.refseq_db.get(refseq_id);
			if (refseq.length > 10000) {
				int mappable_position_count = 0;
				for (int i=0; i<refseq.length; i++)
					if (refseq.mappability[i] > 0) mappable_position_count++;
				if (mappable_position_count < refseq.length/2)
					refseq.is_long_and_unmappable = true;
			}
		}
	}
}
