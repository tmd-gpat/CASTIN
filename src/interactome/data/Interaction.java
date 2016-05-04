package interactome.data;

import interactome.input.GeneInput;

public class Interaction {
	public int id;
	public String type, kegg, url;
	public String ligand_symbol, receptor_symbol;
	public Gene[] ligand_cancer, receptor_cancer, ligand_stroma, receptor_stroma;
	public boolean valid_cancer_to_stroma = true;
	public boolean valid_stroma_to_cancer = true;
	
	// analyze results
	public GeneInput ginput_ligand_cancer = null, ginput_ligand_stroma = null;
	public GeneInput ginput_receptor_cancer = null, ginput_receptor_stroma = null;
	public double average_cancer2stroma;
	public double average_stroma2cancer;
	public double ligand_ratio_cancer;
	public double ligand_ratio_stroma;
	public double ligand_posession_for_same_receptor = 0;
	public double receptor_ratio_stroma;
	public double receptor_ratio_cancer;
	
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
		
		if (this.ligand_cancer.length == 0 || this.receptor_stroma.length == 0)
			this.valid_cancer_to_stroma = false;
		if (this.ligand_stroma.length == 0 || this.receptor_cancer.length == 0)
			this.valid_stroma_to_cancer = false;
	}
}
