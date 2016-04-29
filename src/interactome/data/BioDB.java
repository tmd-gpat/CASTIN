package interactome.data;

import interactome.Logger;
import interactome.Option;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class BioDB {
	private static BioDB _instance;
	public static BioDB getInstance() { return _instance; }
	
	// databases
	public HashMap<String, Refseq> refseq_db;
	public HashMap<Integer, Gene> gene_db;
	public HashMap<Integer, Integer> homolog_cancer2stroma_db;
	public HashMap<Integer, Integer> homolog_stroma2cancer_db;
	
	public String[] cancer_refseq_ids;
	public String[] stromal_refseq_ids;
	
	public Integer[] cancer_entrez_ids;
	public Integer[] stromal_entrez_ids;
	
	private BioDB() {
		this.refseq_db = new HashMap<String, Refseq>();
		this.gene_db = new HashMap<Integer, Gene>();
		this.homolog_cancer2stroma_db = new HashMap<Integer, Integer>();
		this.homolog_stroma2cancer_db = new HashMap<Integer, Integer>();
	}
	
	public static BioDB createInstance() {
		BioDB instance = new BioDB();
		
		if (!instance.construct()) {
			return null;
		}
		
		BioDB._instance = instance;
		return BioDB._instance;
	}
	
	private boolean construct() {
		Option option = Option.getInstance();
		Logger.logf("\nstart loading BioDB.");
		
		// 1. load refseq list from refNames
		ArrayList<String> cancer_refseq_id_list = new ArrayList<String>();
		ArrayList<String> stromal_refseq_id_list = new ArrayList<String>();
		try {
			FileReader fr = null;
			BufferedReader br = null;
			String line = null;

			fr = new FileReader(option.settings.get("cancer_refNames"));
			br = new BufferedReader(fr);
			while ((line = br.readLine()) != null) {
				Refseq refseq = new Refseq(line.trim());
				refseq.tax_id = option.settings.get("cancer_taxonomy");
				this.refseq_db.put(refseq.refseq_id, refseq);
				cancer_refseq_id_list.add(refseq.refseq_id);
			}			
			br.close();
			fr.close();
			
			fr = new FileReader(option.settings.get("stromal_refNames"));
			br = new BufferedReader(fr);
			while ((line = br.readLine()) != null) {
				Refseq refseq = new Refseq(line.trim());
				refseq.tax_id = option.settings.get("stromal_taxonomy");
				this.refseq_db.put(refseq.refseq_id, refseq);
				stromal_refseq_id_list.add(refseq.refseq_id);
			}			
			br.close();
			fr.close();			
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		Logger.logf("%d cancer refseq & %d stromal refseq names are loaded from refNames.",
					cancer_refseq_id_list.size(), stromal_refseq_id_list.size());
		
		// 2. load refLink and relate refseq to genes
		ArrayList<Integer> cancer_entrez_id_list = new ArrayList<Integer>();
		ArrayList<Integer> stromal_entrez_id_list = new ArrayList<Integer>();
		
		try {
			FileReader fr = new FileReader(option.settings.get("refLink"));
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				String symbol = row[0];
				String refseq_id = row[2];
				int entrez_id = Integer.valueOf(row[6]);
				
				Refseq refseq = this.refseq_db.get(refseq_id);
				if (refseq != null) {
					Gene gene = this.gene_db.get(entrez_id);
					if (gene == null) {
						gene = new Gene(entrez_id, symbol);
						gene.tax_id = refseq.tax_id;
						this.gene_db.put(entrez_id, gene);

						if (refseq.tax_id == option.settings.get("cancer_taxonomy")) {
							cancer_entrez_id_list.add(entrez_id);
						} else {
							stromal_entrez_id_list.add(entrez_id);
						}
					}
					gene.variants.add(refseq);
					refseq.gene = gene;
				}
			}
			
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		Logger.logf("%d cancer genes & %d stromal genes are loaded from refLink.",
					cancer_entrez_id_list.size(), stromal_entrez_id_list.size());
		
		// 3. find refseqs whose gene information is missing in refLink.
		int missing_refseqs_cancer = 0;
		int missing_refseqs_stroma = 0;
		for (Map.Entry<String, Refseq> entry : this.refseq_db.entrySet()) {
			Refseq refseq = entry.getValue();
			if (refseq.gene == null) {
				refseq.is_invalid = true;
//				Logger.errorf("refseq %s is missing in refLink.", entry.getKey());
				if (refseq.tax_id == option.settings.get("cancer_taxonomy")) {
					missing_refseqs_cancer++;
				} else {
					missing_refseqs_stroma++;
				}
			}
		}
		Logger.logf("%d cancer refseqs & %d stromal refseqs was missing in refLink.",
					missing_refseqs_cancer, missing_refseqs_stroma);

		Logger.logf("loading BioDB done.");
		return true;
	}
}
