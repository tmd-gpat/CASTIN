package interactome.data;

import interactome.analysis.Analysis;

public class Interaction {
	public int id;
	public String type, kegg, url;
	public String ligand_symbol, receptor_symbol;
	public Gene[] ligand_cancer, receptor_cancer, ligand_stroma, receptor_stroma;
	public boolean valid_cancer_to_stroma = true;
	public boolean valid_stroma_to_cancer = true;
	
	// analyze results
	public double expression_ligand_cancer = 0;
	public double expression_ligand_stroma = 0;
	public double expression_receptor_cancer = 0;
	public double expression_receptor_stroma = 0;
	
	public long raw_ligand_cancer = 0;
	public long raw_ligand_stroma = 0;
	public long raw_receptor_cancer = 0;
	public long raw_receptor_stroma = 0;

	public double average_cancer2stroma;
	public double average_stroma2cancer;
	public double ligand_ratio_cancer;
	public double ligand_ratio_stroma;
	public double ligand_posession_for_same_receptor = 0;
	public double receptor_ratio_stroma;
	public double receptor_ratio_cancer;
	
	public double cancer_ligand_average;
	public double cancer_ligand_receptor_ratio;
	public Analysis.InteractionResult[] cancer_ligand_receptor_ranking;
	
	public double cancer_receptor_average;
	public double cancer_receptor_ligand_ratio;
	public Analysis.InteractionResult[] cancer_receptor_ligand_ranking;
	
	public double stroma_ligand_average;
	public double stroma_ligand_receptor_ratio;
	public Analysis.InteractionResult[] stroma_ligand_receptor_ranking;
	
	public double stroma_receptor_average;
	public double stroma_receptor_ligand_ratio;
	public Analysis.InteractionResult[] stroma_receptor_ligand_ranking;
	
	public Interaction(int id, String type, String kegg, String kegg_id,
					   String ligand_symbol, String receptor_symbol,
					   Gene[] ligand_cancer, Gene[] receptor_cancer,
					   Gene[] ligand_stroma, Gene[] receptor_stroma) {
		this.id = id;
		this.type = type;
		this.kegg = kegg;
		
		this.url = "http://www.genome.jp/kegg/pathway/hsa/hsa" + kegg_id + ".html";
		
		this.ligand_cancer = ligand_cancer;
		this.receptor_cancer = receptor_cancer;
		this.ligand_stroma = ligand_stroma;
		this.receptor_stroma = receptor_stroma;
		
		this.ligand_symbol = ligand_symbol;
		this.receptor_symbol = receptor_symbol;
		
		if (this.ligand_cancer.length == 0 || this.receptor_stroma.length == 0)
			this.valid_cancer_to_stroma = false;
		if (this.ligand_stroma.length == 0 || this.receptor_cancer.length == 0)
			this.valid_stroma_to_cancer = false;
	}
}
