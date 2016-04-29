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
	public HashMap<String, Gene> gene_db;
	public HashMap<String, String[]> homolog_cancer2stroma_db;
	public HashMap<String, String[]> homolog_stroma2cancer_db;
	
	public String[] all_refseq_ids;
	public String[] cancer_refseq_ids;
	public String[] stromal_refseq_ids;
	
	public String[] all_entrez_ids;
	public String[] cancer_entrez_ids;
	public String[] stromal_entrez_ids;
	
	private BioDB() {
		this.refseq_db = new HashMap<String, Refseq>();
		this.gene_db = new HashMap<String, Gene>();
		this.homolog_cancer2stroma_db = new HashMap<String, String[]>();
		this.homolog_stroma2cancer_db = new HashMap<String, String[]>();
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
		Logger.logf("\nstart loading BioDB.");
		
		if (!this.loadRefNames()) return false;
		if (!this.loadRefLink()) return false;
		this.checkMissingRefseqs();
		if (!this.loadRefSeqLen()) return false;
		if (!this.loadHomologene()) return false;

		Logger.logf("loading BioDB done.");
		return true;
	}

	// 1. load refseq list from refNames
	private boolean loadRefNames() {
		Option option = Option.getInstance();
		
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
		
		this.cancer_refseq_ids = cancer_refseq_id_list.toArray(new String[]{});
		this.stromal_refseq_ids = stromal_refseq_id_list.toArray(new String[]{});
		this.all_refseq_ids = new String[this.cancer_refseq_ids.length + this.stromal_refseq_ids.length];
		System.arraycopy(this.cancer_refseq_ids, 0, this.all_refseq_ids, 0, this.cancer_refseq_ids.length);
		System.arraycopy(this.stromal_refseq_ids, 0, this.all_refseq_ids, this.cancer_refseq_ids.length, this.stromal_refseq_ids.length);
		
		Logger.logf("%d cancer refseq & %d stromal refseq names are loaded from refNames.",
					cancer_refseq_id_list.size(), stromal_refseq_id_list.size());
		
		return true;
	}
	
	// 2. load refLink and relate refseq to genes
	private boolean loadRefLink() {
		Option option = Option.getInstance();
		
		ArrayList<String> cancer_entrez_id_list = new ArrayList<String>();
		ArrayList<String> stromal_entrez_id_list = new ArrayList<String>();
		
		try {
			FileReader fr = new FileReader(option.settings.get("refLink"));
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				String symbol = row[0];
				String refseq_id = row[2];
				String entrez_id = row[6];
				
				Refseq refseq = this.refseq_db.get(refseq_id);
				if (refseq != null) {
					Gene gene = this.gene_db.get(entrez_id);
					if (gene == null) {
						gene = new Gene(entrez_id, symbol);
						gene.tax_id = refseq.tax_id;
						this.gene_db.put(entrez_id, gene);

						if (refseq.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
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
		this.cancer_entrez_ids = cancer_entrez_id_list.toArray(new String[]{});
		this.stromal_entrez_ids = stromal_entrez_id_list.toArray(new String[]{});
		this.all_entrez_ids = new String[this.cancer_entrez_ids.length + this.stromal_entrez_ids.length];
		System.arraycopy(this.cancer_entrez_ids, 0, this.all_entrez_ids, 0, this.cancer_entrez_ids.length);
		System.arraycopy(this.stromal_entrez_ids, 0, this.all_entrez_ids, this.cancer_entrez_ids.length, this.stromal_entrez_ids.length);

		
		Logger.logf("%d cancer genes & %d stromal genes are loaded from refLink.",
					cancer_entrez_id_list.size(), stromal_entrez_id_list.size());
		
		return true;
	}
	
	// 3. find refseqs whose gene information is missing in refLink.
	private void checkMissingRefseqs() {
		Option option = Option.getInstance();
		
		int missing_refseqs_cancer = 0;
		int missing_refseqs_stroma = 0;
		for (Map.Entry<String, Refseq> entry : this.refseq_db.entrySet()) {
			Refseq refseq = entry.getValue();
			if (refseq.gene == null) {
				refseq.is_invalid = true;
//				Logger.errorf("refseq %s is missing in refLink.", entry.getKey());
				if (refseq.tax_id.equals(option.settings.get("cancer_taxonomy"))) {
					missing_refseqs_cancer++;
				} else {
					missing_refseqs_stroma++;
				}
			}
		}
		Logger.errorf("%d cancer refseqs & %d stromal refseqs was missing in refLink.",
					  missing_refseqs_cancer, missing_refseqs_stroma);
	}
	
	// 4. load refSeqLen
	private boolean loadRefSeqLen() {
		Option option = Option.getInstance();
		try {
			FileReader fr = new FileReader(option.settings.get("cancer_refSeqLen"));
			BufferedReader br = new BufferedReader(fr);

			String line = null;
			while ((line = br.readLine()) != null) {
				String[] data = line.split("\t");
				
				String refseq_id = data[0];
				String length = data[1];
				
				Refseq refseq = this.refseq_db.get(refseq_id);
				if (refseq != null) {
					refseq.length = Integer.valueOf(length);
				}
			}
			
			br.close();
			fr.close();
			
			fr = new FileReader(option.settings.get("stromal_refSeqLen"));
			br = new BufferedReader(fr);

			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				
				String refseq_id = row[0];
				String length = row[1];
				
				Refseq refseq = this.refseq_db.get(refseq_id);
				if (refseq != null) {
					refseq.length = Integer.valueOf(length);
				}
			}
			
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		Logger.logf("refSeqLen is loaded.");
		
		return true;
	}
	
	// 5. homologene
	private boolean loadHomologene() {
		Option option = Option.getInstance();
		int used_entry_count = 0;
		
		try {
			FileReader fr = new FileReader(option.settings.get("homologene"));
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			String current_homologene_id = null;
			ArrayList<String> h_cancer_entrez_ids = new ArrayList<String>();
			ArrayList<String> h_stromal_entrez_ids = new ArrayList<String>();
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				String homologene_id = row[0];
				String tax_id = row[1];
				String entrez_id = row[2];
				
				if (!homologene_id.equals(current_homologene_id)) {
					// create multiple-genes to multiple-genes map
					if (current_homologene_id != null && h_cancer_entrez_ids.size() > 0 && h_stromal_entrez_ids.size() > 0) {
						String[] ha_cancer_entrez_ids = h_cancer_entrez_ids.toArray(new String[]{});
						String[] ha_stromal_entrez_ids = h_stromal_entrez_ids.toArray(new String[]{});
						
						for (String cancer_entrez_id : h_cancer_entrez_ids) {
							homolog_cancer2stroma_db.put(cancer_entrez_id, ha_stromal_entrez_ids);
						}
						for (String stromal_entrez_id : h_stromal_entrez_ids) {
							homolog_stroma2cancer_db.put(stromal_entrez_id, ha_cancer_entrez_ids);
						}
						used_entry_count++;
					}
					
					h_cancer_entrez_ids.clear();
					h_stromal_entrez_ids.clear();
					current_homologene_id = homologene_id;
				}
				if (tax_id.equals(option.settings.get("cancer_taxonomy"))) {
					h_cancer_entrez_ids.add(entrez_id);
				}
				if (tax_id.equals(option.settings.get("stromal_taxonomy"))) {
					h_stromal_entrez_ids.add(entrez_id);
				}
			}
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		Logger.logf("%d homologene entries for cancer/stromal genes are loaded.", used_entry_count);
		
		return true;
	}
}
