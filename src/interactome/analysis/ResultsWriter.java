package interactome.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Map;

import interactome.Logger;
import interactome.Option;
import interactome.input.Input;
import interactome.input.RefseqInput;

public class ResultsWriter {
	public static ResultsWriter _instance;
	public static ResultsWriter getInstance() {
		return _instance;
	}
	public static ResultsWriter createInstance(Input input, DynamicParameters dynamicParameters) {
		_instance = new ResultsWriter(input, dynamicParameters);
		return _instance;
	}
	
	Input input;
	DynamicParameters dynamicParameters;
	
	public ResultsWriter(Input input, DynamicParameters dynamicParameters) {
		this.input = input;
		this.dynamicParameters = dynamicParameters;
	}
	
	// Refseq_cancer / Refseq_stroma
	public void writeRefseqFiles() {
		Option option = Option.getInstance();
		
		Logger.logf("writing Refseq Files.");
		
		// cancer
		ArrayList<RefseqInput> inputs = new ArrayList<RefseqInput>();
		for (Map.Entry<String, RefseqInput> entry : this.input.refseq_inputs.entrySet()) {
			if (entry.getValue().refseq.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
				inputs.add(entry.getValue());
			}
		}
		this.writeRefseqFile("Refseq_cancer.txt", inputs.toArray(new RefseqInput[]{}));
		
		// stroma
		inputs.clear();
		for (Map.Entry<String, RefseqInput> entry : this.input.refseq_inputs.entrySet()) {
			if (entry.getValue().refseq.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
				inputs.add(entry.getValue());
			}
		}
		this.writeRefseqFile("Refseq_stroma.txt", inputs.toArray(new RefseqInput[]{}));
	}
	
	private void writeRefseqFile(String filename, RefseqInput[] rows) {
		Option option = Option.getInstance();
		
		try {
			FileWriter fw = new FileWriter(option.output_path + "/" + filename);
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"Refseq ID",
				"Gene Symbol",
				"raw count",
				"count (from Poly-A, GC corrected)",
				"length",
				"mappable reads (" + option.input_paired_length + "bp)",
				"no cover",
				"10 cover",
				"20 cover",
				"30 cover",
				"long unmappable refseq"
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}
			
			for (RefseqInput row : rows) {
				Object[] data = new Object[] {
					row.refseq.refseq_id,
					row.refseq.gene.symbol,
					row.rawCount,
					row.true_expression,
					row.refseq.length,
					row.refseq.mappable_position_count_for_specific_read_length(option.input_paired_length),
					1.0 - row.coverage(1),
					row.coverage(10),
					row.coverage(20),
					row.coverage(30),
					row.refseq.is_long_and_unmappable ? "Y" : "N",
				};
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					if (i < data.length-1) bw.write("\t");
					else bw.write("\n");
				}
			}
			
			bw.close();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	// KEGGHPRD_result.txt
	public void writeKEGGHPRDresult() {
		
	}
}
