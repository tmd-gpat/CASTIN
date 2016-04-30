package interactome.input;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import interactome.Logger;
import interactome.Option;
import interactome.data.BioDB;
import interactome.data.Gene;
import interactome.data.Refseq;

public class PairedEndInput extends Input {
	BioDB biodb;
	Option option;
	
	@Override
	public boolean loadFile() {
		this.option = Option.getInstance();
		this.biodb = BioDB.getInstance();
		
		long i;
		long cancer_read_count = 0;
		long stromal_read_count = 0;
		
		Logger.logf("\nstart loading RNA-seq file (paired-ended)");
		try {
			FileReader fr_1 = new FileReader(option.input_prefix_paired + "_1_" + option.input_paired_length + ".sam");
			FileReader fr_2 = new FileReader(option.input_prefix_paired + "_2_" + option.input_paired_length + ".sam");
			BufferedReader br_1 = new BufferedReader(fr_1);
			BufferedReader br_2 = new BufferedReader(fr_2);
			
			// next row data
			String[] next_1, next_2;
			
			// skip header lines
			while (true) {
				next_1 = br_1.readLine().split("\t");
				if (next_1[0].charAt(0) != '@') break;
			}
			while (true) {
				next_2 = br_2.readLine().split("\t");
				if (next_2[0].charAt(0) != '@') break;
			}
			
			ArrayList<String[]> data_1 = new ArrayList<String[]>();
			ArrayList<String[]> data_2 = new ArrayList<String[]>();
			
			// multiple-multiple hits are stored
			ArrayList<DetectedData> reservation = new ArrayList<PairedEndInput.DetectedData>();
			long[] length_stat = new long[20000];
			for (int j=0; j<length_stat.length; j++) length_stat[j] = 0;
			
			for (i=0; ; i++) {
				if (i%1000000 == 0) {
					Logger.logf("processed %d reads.", i);
				}
				
				data_1.clear();
				data_2.clear();
				boolean end = false;
				
				String read_name = next_1[0].split("#")[0];
				Set<Refseq> found_refseqs = new HashSet<Refseq>();
				
				// accept sam entries appropriate for the current directional mode
				// 2nd column of sam file means...
				//     0:  forward mapping
				//     16: reverse mapping
				//     4:  not mapped
				
				// load valid rows from [input 1]
				while (true) {
					int map_flag = Integer.valueOf(next_1[1]);
					if (option.directional_mode == 0 ||
						(option.directional_mode == 1 && (map_flag & 16) == 0) ||
						(option.directional_mode == 2 && (map_flag & 16) != 0)) {
						Refseq refseq = biodb.refseq_db.get(next_1[2]);
						if (refseq != null && !refseq.is_invalid &&
							next_1[5].equals(option.input_paired_length + "M")) {
							found_refseqs.add(refseq);
							data_1.add(next_1);
						}
					}
					
					String line = br_1.readLine();
					if (line == null) {
						end = true;
						break;
					}
					next_1 = line.split("\t");
					if (!next_1[0].split("#")[0].equals(read_name)) {
						break;
					}
				}
				// load valid rows from [input 2]
				while (true) {
					int map_flag = Integer.valueOf(next_2[1]);
					if (option.directional_mode == 0 ||
						(option.directional_mode == 1 && (map_flag & 16) != 0) ||
						(option.directional_mode == 2 && (map_flag & 16) == 0)) {
						Refseq refseq = biodb.refseq_db.get(next_2[2]);
						if (refseq != null && !refseq.is_invalid &&
							next_2[5].equals(option.input_paired_length + "M")) {
							found_refseqs.add(refseq);
							data_2.add(next_2);
						}
					}
					
					String line = br_2.readLine();
					if (line == null) {
						end = true;
						break;
					}
					next_2 = line.split("\t");
					if (!next_2[0].split("#")[0].equals(read_name)) {
						break;
					}
				}
				
				if (data_1.size() == 0 || data_2.size() == 0) continue;
				
				// check unique-gene-hit condition
				Gene gene_1 = biodb.refseq_db.get(data_1.get(0)[2]).gene;
				Gene gene_2 = biodb.refseq_db.get(data_2.get(0)[2]).gene;
				boolean multiple_gene_hit = false;
				for (int j=1; j<data_1.size(); j++) {
					Gene gene_j = biodb.refseq_db.get(data_1.get(j)[2]).gene;
					if (gene_1 != gene_j) {
						multiple_gene_hit = true;
						break;
					}
				}
				for (int j=1; j<data_2.size(); j++) {
					Gene gene_j = biodb.refseq_db.get(data_2.get(j)[2]).gene;
					if (gene_2 != gene_j) {
						multiple_gene_hit = true;
						break;
					}
				}
				if (multiple_gene_hit) continue;
				
				// check the gene correspondence among the pair
				if (gene_1 != gene_2) continue;
				
				// find mapping-pairs
				HashMap<Refseq, DetectedData> detections = new HashMap<Refseq, PairedEndInput.DetectedData>();
				for (Refseq found_refseq : found_refseqs) {
					detections.put(found_refseq, new DetectedData(found_refseq));
				}
				for (String[] data : data_1) {
					Refseq refseq = biodb.refseq_db.get(data[2]);
					DetectedData detection = detections.get(refseq);
					detection.pos_a.add(Integer.valueOf(data[3]));
					detection.reverse_a.add((Integer.valueOf(data[1]) & 16) != 0);
				}
				for (String[] data : data_2) {
					Refseq refseq = biodb.refseq_db.get(data[2]);
					DetectedData detection = detections.get(refseq);
					detection.pos_b.add(Integer.valueOf(data[3]));
					detection.reverse_b.add((Integer.valueOf(data[1]) & 16) != 0);
				}
				
				// count valid mapping-pair
				boolean valid_read = false;
				for (Map.Entry<Refseq, DetectedData> entry : detections.entrySet()) {
					DetectedData data = entry.getValue();
					if (data.pos_a.size() == 1 && data.pos_b.size() == 1 &&
						data.reverse_a.get(0) != data.reverse_b.get(0)) {
						
						int pos_a = data.pos_a.get(0);
						int pos_b = data.pos_b.get(0);
						int wrapped_length = Math.abs(pos_a - pos_b) + option.input_paired_length;
						this.incrementPair(entry.getKey(), pos_a, pos_b);
						if (wrapped_length < length_stat.length) {
							length_stat[wrapped_length]++;
						}
						valid_read = true;
					} else if (data.pos_a.size() > 0 && data.pos_b.size() > 0) {
						reservation.add(data);
					}
				}
				
				if (valid_read) {
					if (gene_1.tax_id.equals(this.option.settings.get("cancer_taxonomy"))) {
						cancer_read_count++;
					} else {
						stromal_read_count++;
					}
				}
				
				if (end) break;
			}
			
			// calculate wrapped length average & sd
			double average = 0;
			long count = 0;
			for (int j=0; j<length_stat.length; j++) {
				average += j * length_stat[j];
				count += length_stat[j];
			}
			average /= count;
			
			double sd = 0;
			for (int j=0; j<length_stat.length; j++) {
				sd += (average - j) * (average - j) * length_stat[j];
			}
			sd /= count;
			sd = Math.sqrt(sd);
			
			Logger.logf("average wrapped length = %f (sd = %f)", average, sd);
			
			// process reserved multiple-multiple detections
			// positions whose wrap-length is most close to the average are accepted
			for (DetectedData detect : reservation) {
				int best_wrap_length = Integer.MAX_VALUE;
				int final_pos_a = 0;
				int final_pos_b = 0;
				for (int x=0; x<detect.pos_a.size(); x++) {
					for (int y=0; y<detect.pos_b.size(); y++) {
						int pos_a = detect.pos_a.get(x);
						int pos_b = detect.pos_b.get(y);
						boolean rev_a = detect.reverse_a.get(x);
						boolean rev_b = detect.reverse_b.get(y);
						if (rev_a == rev_b) continue;
						
						int wrap_length = Math.abs(pos_a - pos_b) + option.input_paired_length;
						if (Math.abs(wrap_length - average) < Math.abs(best_wrap_length - average)
							&& Math.abs(wrap_length - average) < sd*2) {
							best_wrap_length = wrap_length;
							final_pos_a = pos_a;
							final_pos_b = pos_b;
						}
					}
				}
				
				if (best_wrap_length < Integer.MAX_VALUE) {
					this.incrementPair(detect.refseq, final_pos_a, final_pos_b);
					if (detect.refseq.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
						cancer_read_count++;
					} else {
						stromal_read_count++;
					}
				}
			}
			
			br_1.close();
			br_2.close();
			fr_1.close();
			fr_2.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		Logger.logf("%d reads are loaded.", i);
		Logger.logf("%d reads are mapped to unique genes.", cancer_read_count + stromal_read_count);
		Logger.logf("(cancer: %d, stroma: %d)", cancer_read_count, stromal_read_count);
		
		return true;
	}

	void incrementPair(Refseq refseq, int pos_a, int pos_b) {
		RefseqInput entry = refseq_inputs.get(refseq.refseq_id);
		
		entry.rawCount++;
		entry.starting_counts[pos_a]++;
		entry.starting_counts[pos_b]++;
		
		for (int i=pos_a; i<pos_a+option.input_paired_length && i<refseq.length; i++) {
			entry.overwrap_counts[i]++;
		}
		for (int i=pos_b; i<pos_b+option.input_paired_length && i<refseq.length; i++) {
			entry.overwrap_counts[i]++;
		}
	}

	class DetectedData {
		public Refseq refseq;
		public ArrayList<Integer> pos_a;
		public ArrayList<Integer> pos_b;
		public ArrayList<Boolean> reverse_a;
		public ArrayList<Boolean> reverse_b;
		public ArrayList<Integer> len_a;
		public ArrayList<Integer> len_b;
		
		public DetectedData(Refseq refseq) {
			this.refseq = refseq;
			this.pos_a = new ArrayList<Integer>();
			this.pos_b = new ArrayList<Integer>();
			this.reverse_a = new ArrayList<Boolean>();
			this.reverse_b = new ArrayList<Boolean>();
			this.len_a = new ArrayList<Integer>();
			this.len_b = new ArrayList<Integer>();
		}
	}
}
