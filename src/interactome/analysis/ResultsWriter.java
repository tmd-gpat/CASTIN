package interactome.analysis;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;

import org.tc33.jheatchart.HeatChart;

import interactome.Logger;
import interactome.Option;
import interactome.analysis.Analysis.InteractionResult;
import interactome.data.BioDB;
import interactome.data.Interaction;
import interactome.data.Refseq;
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
	
	public void drawHeatmaps() {
		// all
		drawPolyAHeatMap("from-poly-A", input.refseq_inputs.values(), false);
		drawPolyAHeatMap("from-poly-A", input.refseq_inputs.values(), true);
		drawGCHeatMap("GC", input.refseq_inputs.values(), false);
		drawGCHeatMap("GC", input.refseq_inputs.values(), true);
		
		// estimator
		ArrayList<RefseqInput> estim_rinputs = new ArrayList<RefseqInput>();
		for (int i=0; i<dynamicParameters.parameterRefseqs.length; i++) {
			estim_rinputs.add(input.refseq_inputs.get(dynamicParameters.parameterRefseqs[i].refseq_id));
		}
		drawPolyAHeatMap("from-poly-A_estim", estim_rinputs, false);
		drawPolyAHeatMap("from-poly-A_estim", estim_rinputs, true);
		drawGCHeatMap("GC_estim", estim_rinputs, false);
		drawGCHeatMap("GC_estim", estim_rinputs, true);
		
		Logger.logf("wrote heatmaps.");
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
	
	// KEGGHPRD_result(_cancer/stroma_ligand/receptor).txt
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

	/*
	 * draw heat maps where X-axis = poly-a dist, Y-axis = residue of regression equation for before/after correction:
	 * log(C_i_j) = log( mappability_ratio * v_i ) + (gc %) * alpha + (poly-A) * beta + [Residue]
	 * 	 C_i_j = count of reads overlapping position j of gene i
	 *   mappability_ratio = refer to the original paper
	 *   v_i = expression
	 * 
	 * before: v = v^1_i, alpha = beta = 0
	 * after : v = true_expression, alpha = alpha', beta = beta'
	 */
	private void drawPolyAHeatMap(String prefix, Collection<RefseqInput> rinputs, boolean after) {
		Option option = Option.getInstance();
		
		double residue_min = -15;
		double residue_max = 15;
		int polya_min = 0;
		int polya_max = 20000;
		
		int x_split = 200;
		int y_split = 200;
		
		double[][] data = new double[y_split+1][x_split+1];
		
		for (RefseqInput rinput : input.refseq_inputs.values()) {
			Refseq refseq = rinput.refseq;
			for (int i=0; i<rinput.refseq.length; i++) {
				if ((!after && Float.isNaN(rinput.residue_init[i])) || (after && Float.isNaN(rinput.residue_end[i])) ) continue;
				int x = (int)Math.round((double)(refseq.length - i) / (polya_max - polya_min) * x_split);
				if (x >= x_split) continue;
				
				float left = after ? rinput.residue_end[i] : rinput.residue_init[i];
				int y = (int)Math.round((residue_max - left) / (residue_max - residue_min) * y_split);
				if (y < 0) continue;
				if (y >= y_split) continue;
				
				data[y][x]++;
			}
		}
		
		try {
			HeatChart heatmap = new HeatChart(data);
			heatmap.setBackgroundColour(new Color(0, 0, 0));
			heatmap.setHighValueColour(new Color(0, 255, 0));
			heatmap.setLowValueColour(new Color(0, 0, 0));
			heatmap.setAxisColour(new Color(255, 255, 255));
			heatmap.setAxisLabelColour(new Color(255, 255, 255));
			heatmap.setAxisValuesColour(new Color(255, 255, 255));
			heatmap.setAxisLabelsFont(new Font("SansSerif", Font.PLAIN, 72));
			heatmap.setAxisValuesFont(new Font("SansSerif", Font.PLAIN, 64));
			heatmap.setColourScale(0.8);
			heatmap.setXAxisLabel("from poly-A");
			heatmap.setYAxisLabel("regression residue");
			heatmap.setChartMargin(32);
			
			Object[] x_values = new Object[x_split+1];
			Object[] y_values = new Object[y_split+1];
			for (int i=0; i<=x_split; i++) {
				x_values[i] = polya_min + (polya_max - polya_min) / x_split * i;
			}
			for (int i=0; i<=y_split; i++) {
				y_values[i] = residue_max - (residue_max - residue_min) / y_split * i;
			}
			heatmap.setXValues(x_values);
			heatmap.setYValues(y_values);
			heatmap.setXAxisValuesFrequency(20);
			heatmap.setYAxisValuesFrequency(20);
			heatmap.setShowXAxisValues(true);
			heatmap.setShowYAxisValues(true);
			
			heatmap.saveToFile(
				new File(option.output_path + "/" + prefix + "_" + (after ? "after" : "before") + "_correction.png"));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/*
	 * draw heat maps where X-axis = gc%, Y-axis = residue of regression equation for before/after correction:
	 * log(C_i_j) = log( mappability_ratio * v_i ) + (gc %) * alpha + (poly-A) * beta + [Residue]
	 * 	 C_i_j = count of reads overlapping position j of gene i
	 *   mappability_ratio = refer to the original paper
	 *   v_i = expression
	 * 
	 * before: v = v^1_i, alpha = beta = 0
	 * after : v = true_expression, alpha = alpha', beta = beta'
	 */
	private void drawGCHeatMap(String prefix, Collection<RefseqInput> rinputs, boolean after) {
		Option option = Option.getInstance();
		
		double residue_min = -15;
		double residue_max = 15;
		int gc_min = 0;
		int gc_max = 100;
		
		int x_split = 100;
		int y_split = 100;
		
		double[][] data = new double[y_split+1][x_split+1];
		
		for (RefseqInput rinput : input.refseq_inputs.values()) {
			Refseq refseq = rinput.refseq;
			for (int i=0; i<rinput.refseq.length; i++) {
				if ((!after && Float.isNaN(rinput.residue_init[i])) || (after && Float.isNaN(rinput.residue_end[i])) ) continue;
				int x = (int)Math.round((double)refseq.gc_percent[i] / (gc_max - gc_min) * x_split);
				if (x >= x_split) continue;
				
				float left = after ? rinput.residue_end[i] : rinput.residue_init[i];
				int y = (int)Math.round((residue_max - left) / (residue_max - residue_min) * y_split);
				if (y < 0) continue;
				if (y >= y_split) continue;
				
				data[y][x]++;
			}
		}
		
		try {
			HeatChart heatmap = new HeatChart(data);
			heatmap.setBackgroundColour(new Color(0, 0, 0));
			heatmap.setHighValueColour(new Color(0, 255, 0));
			heatmap.setLowValueColour(new Color(0, 0, 0));
			heatmap.setAxisColour(new Color(255, 255, 255));
			heatmap.setAxisLabelColour(new Color(255, 255, 255));
			heatmap.setAxisValuesColour(new Color(255, 255, 255));
			heatmap.setAxisLabelsFont(new Font("SansSerif", Font.PLAIN, 36));
			heatmap.setAxisValuesFont(new Font("SansSerif", Font.PLAIN, 32));
			heatmap.setColourScale(0.8);
			heatmap.setXAxisLabel("GC %");
			heatmap.setYAxisLabel("regression residue");
			heatmap.setChartMargin(32);
			
			Object[] x_values = new Object[x_split+1];
			Object[] y_values = new Object[y_split+1];
			for (int i=0; i<=x_split; i++) {
				x_values[i] = gc_min + (gc_max - gc_min) / x_split * i;
			}
			for (int i=0; i<=y_split; i++) {
				y_values[i] = residue_max - (residue_max - residue_min) / y_split * i;
			}
			heatmap.setXValues(x_values);
			heatmap.setYValues(y_values);
			heatmap.setXAxisValuesFrequency(10);
			heatmap.setYAxisValuesFrequency(10);
			heatmap.setShowXAxisValues(true);
			heatmap.setShowYAxisValues(true);
			
			heatmap.saveToFile(
				new File(option.output_path + "/" + prefix + "_" + (after ? "after" : "before") + "_correction.png"));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
