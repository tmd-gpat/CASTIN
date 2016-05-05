package interactome.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;

import interactome.Logger;
import interactome.Option;
import interactome.analysis.Analysis.InteractionResult;
import interactome.data.BioDB;
import interactome.data.Interaction;
import interactome.input.GeneInput;
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
			if (entry.getValue().refseq.tax_id.equals(option.settings.get("stromal_taxonomy"))) {
				inputs.add(entry.getValue());
			}
		}
		this.writeRefseqFile("Refseq_stroma.txt", inputs.toArray(new RefseqInput[]{}));
		
		Logger.logf("wrote Refseq files.");
	}
	
	public void writeSymbolFiles() {
		Option option = Option.getInstance();
		
		// cancer
		ArrayList<GeneInput> inputs = new ArrayList<GeneInput>();
		for (Map.Entry<String, GeneInput> entry : this.input.gene_inputs.entrySet()) {
			if (entry.getValue().gene.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
				inputs.add(entry.getValue());
			}
		}
		this.writeSymbolFile("Symbol_cancer.txt", inputs.toArray(new GeneInput[]{}));
		
		// stroma
		inputs.clear();
		for (Map.Entry<String, GeneInput> entry : this.input.gene_inputs.entrySet()) {
			if (entry.getValue().gene.tax_id.equals(option.settings.get("stromal_taxonomy"))) {
				inputs.add(entry.getValue());
			}
		}
		this.writeSymbolFile("Symbol_stroma.txt", inputs.toArray(new GeneInput[]{}));
		
		Logger.logf("wrote Symbol files.");
	}
	
	public void writeKEGGHPRDFiles() {
		writeKEGGHPRDResult();
		
		writeKEGGHPRDResultCancerLigand();
		writeKEGGHPRDResultCancerReceptor();
		writeKEGGHPRDResultStromaLigand();
		writeKEGGHPRDResultStromaReceptor();
		
		Logger.logf("wrote KEGGHPRD files.");
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
	
	private void writeSymbolFile(String filename, GeneInput[] rows) {
		Option option = Option.getInstance();
		
		try {
			FileWriter fw = new FileWriter(option.output_path + "/" + filename);
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"Gene Symbol",
				"Refseq ID",
				"raw count",
				"count (from Poly-A, GC corrected)",
				"count (from Poly-A, GC corrected, Total=300k)",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}
			
			for (GeneInput row : rows) {
				Object[] data = new Object[] {
					row.gene.symbol,
					row.representativeRefseq.refseq_id,
					row.representativeRefseqInput.rawCount,
					row.representativeRefseqInput.true_expression,
					row.normalizedExpression,
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
	private void writeKEGGHPRDResult() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
	
		try {
			FileWriter fw = new FileWriter(option.output_path + "/KEGGHPRD_result.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"ligand",
				"receptor",
				"count(ligand)",
				"count(receptor)",
				"raw count(ligand)",
				"raw count(receptor)",
				"average",
				"ligand ratio",
				"count of the other ligands",
				"receptor ratio",
				"interaction type",
				"pathway",
				"link",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}
			
			for (Interaction interaction : biodb.interactions) {
				// cancer -> stroma
				Object[] data;
				if (interaction.valid_cancer_to_stroma) {
					data = new Object[] {
						interaction.ligand_symbol + "(cancer)",
						interaction.receptor_symbol + "(stroma)",
						interaction.ginput_ligand_cancer.normalizedExpression,
						interaction.ginput_receptor_stroma.normalizedExpression,
						interaction.ginput_ligand_cancer.representativeRefseqInput.rawCount,
						interaction.ginput_receptor_stroma.representativeRefseqInput.rawCount,
						interaction.average_cancer2stroma,
						interaction.ligand_ratio_cancer == -1 ? "NA" : interaction.ligand_ratio_cancer,
						interaction.ligand_posession_for_same_receptor,
						interaction.receptor_ratio_stroma == -1 ? "NA" : interaction.receptor_ratio_stroma,
						interaction.type,
						interaction.kegg,
						interaction.url,
					};
				} else {
					data = new Object[] {
						interaction.ligand_symbol + "(cancer)",
						interaction.receptor_symbol + "(stroma/NA)",
						interaction.ginput_ligand_cancer.normalizedExpression,
						"NA",
					};
				}
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					if (i < data.length-1) bw.write("\t");
					else bw.write("\n");
				}
				// stroma -> cancer
				if (interaction.valid_stroma_to_cancer) {
					data = new Object[] {
						interaction.ligand_symbol + "(stroma)",
						interaction.receptor_symbol + "(cancer)",
						interaction.ginput_ligand_stroma.normalizedExpression,
						interaction.ginput_receptor_cancer.normalizedExpression,
						interaction.ginput_ligand_stroma.representativeRefseqInput.rawCount,
						interaction.ginput_receptor_cancer.representativeRefseqInput.rawCount,
						interaction.average_stroma2cancer,
						interaction.ligand_ratio_stroma == -1 ? "NA" : interaction.ligand_ratio_stroma,
						interaction.ligand_posession_for_same_receptor,
						interaction.receptor_ratio_cancer == -1 ? "NA" : interaction.receptor_ratio_cancer,
						interaction.type,
						interaction.kegg,
						interaction.url,
					};
				} else {
					data = new Object[] {
						interaction.ligand_symbol + "(stroma/NA)",
						interaction.receptor_symbol + "(cancer)",
						"NA",
						interaction.ginput_receptor_cancer.normalizedExpression,
					};
				}
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
	
	private void writeKEGGHPRDResultCancerLigand() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
	
		try {
			FileWriter fw = new FileWriter(option.output_path + "/KEGGHPRD_result_cancer_ligand.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"ligand",
				"average count",
				"receptor ratio",
				"ligand ratio",
				"activate or inhibit",
				"pathway",
				"url",
				"receptor (normalized count)",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}

			HashSet<String> used = new HashSet<String>();
			for (Interaction interaction : biodb.interactions) {
				if (used.contains(interaction.ligand_symbol)) continue;
				used.add(interaction.ligand_symbol);

				Object[] data;
				data = new Object[] {
					interaction.ligand_symbol,
					interaction.cancer_ligand_average,
					interaction.cancer_ligand_receptor_ratio == -1 ? "NA" : interaction.cancer_ligand_receptor_ratio,
					interaction.ligand_ratio_cancer == -1 ? "NA" : interaction.ligand_ratio_cancer,
					interaction.type,
					interaction.kegg,
					interaction.url,
				};
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					bw.write("\t");
				}
				// ranking
				for (InteractionResult res : interaction.cancer_ligand_receptor_ranking) {
					bw.write(res.symbol + "(" + res.expression + ")\t");
				}
				bw.write("\n");
			}
			
			bw.close();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	private void writeKEGGHPRDResultCancerReceptor() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
	
		try {
			FileWriter fw = new FileWriter(option.output_path + "/KEGGHPRD_result_cancer_receptor.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"ligand",
				"average count",
				"ligand ratio",
				"receptor ratio",
				"activate or inhibit",
				"pathway",
				"url",
				"ligand (normalized count)",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}
			
			HashSet<String> used = new HashSet<String>();
			for (Interaction interaction : biodb.interactions) {
				if (used.contains(interaction.receptor_symbol)) continue;
				used.add(interaction.receptor_symbol);
				Object[] data;
				data = new Object[] {
					interaction.receptor_symbol,
					interaction.cancer_receptor_average,
					interaction.cancer_receptor_ligand_ratio == -1 ? "NA" : interaction.cancer_receptor_ligand_ratio,
					interaction.receptor_ratio_cancer == -1 ? "NA" : interaction.receptor_ratio_cancer,
					interaction.type,
					interaction.kegg,
					interaction.url,
				};
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					bw.write("\t");
				}
				// ranking
				for (InteractionResult res : interaction.cancer_receptor_ligand_ranking) {
					bw.write(res.symbol + "(" + res.expression + ")\t");
				}
				bw.write("\n");
			}
			
			bw.close();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	private void writeKEGGHPRDResultStromaLigand() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
	
		try {
			FileWriter fw = new FileWriter(option.output_path + "/KEGGHPRD_result_stroma_ligand.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"ligand",
				"average count",
				"receptor ratio",
				"ligand ratio",
				"activate or inhibit",
				"pathway",
				"url",
				"receptor (normalized count)",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}

			HashSet<String> used = new HashSet<String>();
			for (Interaction interaction : biodb.interactions) {
				if (!interaction.valid_stroma_to_cancer || used.contains(interaction.ligand_symbol)) continue;
				used.add(interaction.ligand_symbol);
				
				Object[] data;
				data = new Object[] {
					interaction.ligand_symbol,
					interaction.stroma_ligand_average,
					interaction.stroma_ligand_receptor_ratio == -1 ? "NA" : interaction.stroma_ligand_receptor_ratio,
					interaction.ligand_ratio_stroma == -1 ? "NA" : interaction.ligand_ratio_stroma,
					interaction.type,
					interaction.kegg,
					interaction.url,
				};
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					bw.write("\t");
				}
				// ranking
				for (InteractionResult res : interaction.stroma_ligand_receptor_ranking) {
					bw.write(res.symbol + "(" + res.expression + ")\t");
				}
				bw.write("\n");
			}
			
			bw.close();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
	
	private void writeKEGGHPRDResultStromaReceptor() {
		Option option = Option.getInstance();
		BioDB biodb = BioDB.getInstance();
	
		try {
			FileWriter fw = new FileWriter(option.output_path + "/KEGGHPRD_result_stroma_receptor.txt");
			BufferedWriter bw = new BufferedWriter(fw);
			
			// write header
			String[] column_names = new String[] {
				"ligand",
				"average count",
				"ligand ratio",
				"receptor ratio",
				"activate or inhibit",
				"pathway",
				"url",
				"ligand (normalized count)",
			};
			for (int i=0; i<column_names.length; i++) {
				bw.write(column_names[i]);
				if (i<column_names.length-1) bw.write("\t");
				else bw.write("\n");
			}

			HashSet<String> used = new HashSet<String>();
			for (Interaction interaction : biodb.interactions) {
				if (!interaction.valid_cancer_to_stroma || used.contains(interaction.receptor_symbol)) continue;
				used.add(interaction.receptor_symbol);
				
				Object[] data;
				data = new Object[] {
					interaction.receptor_symbol,
					interaction.stroma_receptor_average,
					interaction.stroma_receptor_ligand_ratio == -1 ? "NA" : interaction.stroma_receptor_ligand_ratio,
					interaction.receptor_ratio_stroma == -1 ? "NA" : interaction.receptor_ratio_stroma,
					interaction.type,
					interaction.kegg,
					interaction.url,
				};
				for (int i=0; i<data.length; i++) {
					bw.write(String.valueOf(data[i]));
					bw.write("\t");
				}
				// ranking
				for (InteractionResult res : interaction.cancer_receptor_ligand_ranking) {
					bw.write(res.symbol + "(" + res.expression + ")\t");
				}
				bw.write("\n");
			}
			
			bw.close();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}
	}
}
