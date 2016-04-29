package interactome.data;

import interactome.Pair;

import java.util.ArrayList;

public class Refseq {
	public String refseq_id; // e.g., NM_******
	
	public Gene gene;
	public int length;
	public String chr; // chromosome
	public int tax_id; // e.g., 9606
	
	public boolean is_invalid; // lacking refLink entry

	public short[] gc_percent;
	public double[] mappability;

	public ArrayList<Pair> unmappable_regions;
	
	public Refseq(String refseq_id) {
		this.refseq_id = refseq_id;
		this.is_invalid = false;
	}
	
	// number of mappable positions where n-bp read starting from there is uniquely mappable
	public int mappable_position_count_for_specific_read_length(int n) {
		int count = 0;
		for (int i=0; i<=length-n; i++) {
			boolean is_mappable = true;
			for (Pair pair : unmappable_regions) {
				if ((Integer)pair.first <= i && (i+n) < (Integer)pair.second) {
					is_mappable = false;
					break;
				}
			}
			if (is_mappable) count++;
		}
		return count;
	}
}
