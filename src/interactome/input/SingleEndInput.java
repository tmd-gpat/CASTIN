package interactome.input;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Random;

import interactome.Logger;
import interactome.Option;
import interactome.data.BioDB;
import interactome.data.Gene;
import interactome.data.Refseq;

public class SingleEndInput extends Input {
	BioDB biodb;
	final int MIN_MATCH_LENGTH = 50;
	Random r;
	
	public SingleEndInput() {
		super();
		this.biodb = BioDB.getInstance();
		this.r = new Random();
	}

	@Override
	public boolean loadFile() {
		Option option = Option.getInstance();
		
		Logger.logf("\nstart loading RNA-seq file (single-ended)");
		try {
			FileReader fr = new FileReader(option.input_file_single);
			BufferedReader br = new BufferedReader(fr);
			
			long row_count = 0; // number of sam entry
			long map_count = 0; // number of mapped sam entry
			long accepted_read_count = 0;
			long cancer_read_count = 0;
			long stromal_read_count = 0;
			
			String current_read = "";
			ArrayList<String> mapped_refseq_ids = new ArrayList<String>();
			ArrayList<Integer> mapped_positions = new ArrayList<Integer>();
			ArrayList<Integer> mapped_lengthes = new ArrayList<Integer>(); 
			
			String line = null;
			while ((line = br.readLine()) != null) {
				if (line.charAt(0) == '@') continue;
				if (row_count++ % 1_000_000 == 0) Logger.logf("processed %d rows", row_count);
				
				String[] row = line.split("\t");
				String read_id = row[0];
				String refseq_id = row[2];
				
				// no hit
				if (refseq_id.contains("*")) continue;
				map_count++;
				
				// calculate match-length
				int match_length = 0;
				String[] length_chunks = row[5].split("[A-Z]");
				
				int j = 0;
				for (int i=0; i<row[5].length(); i++) {
					char c = row[5].charAt(i);
					if (c == 'M' || c == 'D' || c == 'N') {
						match_length += Integer.valueOf(length_chunks[j]);
					}
					if (Character.isLetter(c)) j++;
				}
				if (match_length > this.match_length_distribution.length) {
					Logger.errorf("too long match length! (%d bp, readID = %s)", match_length, read_id);
					continue;
				}
				this.match_length_distribution[match_length]++;
				
				int match_position = Integer.valueOf(row[3]);
				
				// ignore too short match length
				if (match_length < MIN_MATCH_LENGTH) continue;
				
				Refseq refseq = biodb.refseq_db.get(refseq_id);
				if (refseq == null) continue;
				
				// when current read was changed
				if (!current_read.equals(read_id)) {
					// if the read was mapped to an unique gene
					if (isSingleGeneHit(mapped_refseq_ids)) {
						// count up
						countRead(
							mapped_refseq_ids.toArray(new String[]{}),
							mapped_positions.toArray(new Integer[]{}),
							mapped_lengthes.toArray(new Integer[]{}));
						accepted_read_count++;
						
						if (biodb.refseq_db.get(mapped_refseq_ids.get(0)).tax_id.equals(option.settings.get("cancer_taxonomy"))) {
							cancer_read_count++;
						} else {
							stromal_read_count++;
						}
					}
					// clear buffers
					mapped_refseq_ids.clear();
					mapped_positions.clear();
					mapped_lengthes.clear();
				}
				
				mapped_refseq_ids.add(refseq_id);
				mapped_positions.add(match_position);
				mapped_lengthes.add(match_length);
			}
			// process final read
			if (isSingleGeneHit(mapped_refseq_ids)) {
				// count up
				countRead(
					mapped_refseq_ids.toArray(new String[]{}),
					mapped_positions.toArray(new Integer[]{}),
					mapped_lengthes.toArray(new Integer[]{}));
				accepted_read_count++;
				
				if (biodb.refseq_db.get(mapped_refseq_ids.get(0)).tax_id.equals(option.settings.get("cancer_taxonomy"))) {
					cancer_read_count++;
				} else {
					stromal_read_count++;
				}
			}
			
			Logger.logf("%d SAM rows is loaded.", row_count);
			Logger.logf("%d map is contained in the SAM file.", map_count);
			Logger.logf("%d reads is mapped to a unique gene.", accepted_read_count);
			Logger.logf("(cancer: %d, stroma: %d)", cancer_read_count, stromal_read_count);
			
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		return true;
	}

	boolean isSingleGeneHit(ArrayList<String> refseq_ids) {
		if (refseq_ids.size() == 0) return false;

		Refseq refseq = null;
		refseq = biodb.refseq_db.get(refseq_ids.get(0));
		if (refseq.is_invalid) return false;
		
		Gene gene = refseq.gene;
		for (int i=1; i<refseq_ids.size(); i++) {
			refseq = biodb.refseq_db.get(refseq_ids.get(i));
			if (refseq.is_invalid) return false;
			if (refseq.gene != gene) return false;
		}
		
		return true;
	}

	void countRead(String[] refseq_ids, Integer[] positions, Integer[] lengthes) {
		boolean[] used = new boolean[refseq_ids.length];
		for (int i=0; i<refseq_ids.length; i++) used[i] = false;
		
		ArrayList<Integer>indexes = new ArrayList<Integer>();
		for (int i=0; i<refseq_ids.length; i++) {
			if (used[i]) continue;
			indexes.clear();
			
			for (int j=i; j<refseq_ids.length; j++) {
				if (refseq_ids[i].equals(refseq_ids[j])) {
					indexes.add(j);
					used[j] = true;
				}
			}
			
			int indexForCount = r.nextInt(indexes.size());
			this.incrementSingle(refseq_ids[indexForCount], positions[indexForCount], lengthes[indexForCount]);
		}
	}

	void incrementSingle(String refseq_id, int position, int length) {
		RefseqInput entry = refseq_inputs.get(refseq_id);
		
		entry.rawCount++;
		entry.starting_counts[position]++;
		
		// 端を飛び出るタイプのマッチは切り捨てる
		for (int i=0; i<length && (position+i)<entry.overwrap_counts.length; i++) {
			entry.overwrap_counts[position+i]++;
		}
	}
}
