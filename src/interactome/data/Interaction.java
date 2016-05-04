package interactome.data;

public class Interaction {
	int id;
	String type, kegg, url;
	Gene[] ligand_cancer, receptor_cancer, ligand_stroma, receptor_stroma;
	boolean valid_cancer_to_stroma = true;
	boolean valid_stroma_to_cancer = true;
	
	public Interaction(int id, String type, String kegg, String kegg_id,
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
