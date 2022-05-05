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
	public HashMap<String, Gene> gene_cancer_symbol_db;
	public HashMap<String, Gene> gene_stromal_symbol_db;
	public HashMap<Gene, Gene[]> homolog_cancer2stroma_db;
	public HashMap<Gene, Gene[]> homolog_cancer2cancer_db;
	public HashMap<Gene, Gene[]> homolog_stroma2cancer_db;
	public HashMap<Gene, Gene[]> homolog_stroma2stroma_db;
	
	public String[] all_refseq_ids;
	public String[] cancer_refseq_ids;
	public String[] stromal_refseq_ids;
	
	public String[] all_entrez_ids;
	public String[] cancer_entrez_ids;
	public String[] stromal_entrez_ids;
	
	public Interaction[] interactions;
	
	private BioDB() {
		this.refseq_db = new HashMap<String, Refseq>();
		this.gene_db = new HashMap<String, Gene>();
		this.gene_cancer_symbol_db = new HashMap<String, Gene>();
		this.gene_stromal_symbol_db = new HashMap<String, Gene>();
		this.homolog_cancer2stroma_db = new HashMap<Gene, Gene[]>();
		this.homolog_stroma2cancer_db = new HashMap<Gene, Gene[]>();
		this.homolog_cancer2cancer_db = new HashMap<Gene, Gene[]>();
		this.homolog_stroma2stroma_db = new HashMap<Gene, Gene[]>();
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
		if (!this.loadInteractions()) return false;

		Logger.logf("loading BioDB done.");

		Option option = Option.getInstance();
		if (!this.calculateGCPercent(option.settings.get("cancer_refMrna"))) return false;
		if (!this.calculateGCPercent(option.settings.get("stromal_refMrna"))) return false;
		
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
							this.gene_cancer_symbol_db.put(symbol, gene);
							this.homolog_cancer2cancer_db.put(gene, new Gene[]{gene});
							this.homolog_cancer2stroma_db.put(gene, new Gene[]{});
						} else {
							stromal_entrez_id_list.add(entrez_id);
							this.gene_stromal_symbol_db.put(symbol, gene);
							this.homolog_stroma2stroma_db.put(gene, new Gene[]{gene});
							this.homolog_stroma2cancer_db.put(gene, new Gene[]{});
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
		Logger.errorf("%d cancer refseqs & %d stromal refseqs were missing in refLink.",
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
			ArrayList<Gene> h_cancer_genes = new ArrayList<Gene>();
			ArrayList<Gene> h_stromal_genes = new ArrayList<Gene>();
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				String homologene_id = row[0];
				String tax_id = row[1];
				String entrez_id = row[2];
				
				if (!homologene_id.equals(current_homologene_id)) {
					// create multiple-genes to multiple-genes map
					if (current_homologene_id != null &&
						h_cancer_genes.size() > 0 || h_stromal_genes.size() > 0) {
						Gene[] ha_cancer_genes = h_cancer_genes.toArray(new Gene[]{});
						Gene[] ha_stromal_genes = h_stromal_genes.toArray(new Gene[]{});
						
						for (Gene cancer_gene : h_cancer_genes) {
							homolog_cancer2stroma_db.put(cancer_gene, ha_stromal_genes);
							homolog_cancer2cancer_db.put(cancer_gene, ha_cancer_genes);
						}
						for (Gene stromal_gene : h_stromal_genes) {
							homolog_stroma2cancer_db.put(stromal_gene, ha_cancer_genes);
							homolog_stroma2stroma_db.put(stromal_gene, ha_stromal_genes);
						}
						used_entry_count++;
					}
					
					h_cancer_genes.clear();
					h_stromal_genes.clear();
					current_homologene_id = homologene_id;
				}
				if (tax_id.equals(option.settings.get("cancer_taxonomy"))) {
					Gene gene = gene_db.get(entrez_id);
					if (gene != null) h_cancer_genes.add(gene);
				}
				if (tax_id.equals(option.settings.get("stromal_taxonomy"))) {
					Gene gene = gene_db.get(entrez_id);
					if (gene != null) h_stromal_genes.add(gene);
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

	// 6. KEGG_HPRD
	private boolean loadInteractions() {
		Option option = Option.getInstance();
		ArrayList<Interaction> interactions = new ArrayList<Interaction>();
		
		try {
			FileReader fr = new FileReader(option.settings.get("curated_HPRD_KEGG"));
			BufferedReader br = new BufferedReader(fr);

			br.readLine(); // skip header line
			String line = null;
			while ((line = br.readLine()) != null) {
				String[] row = line.split("\t");
				if (row.length < 12) {
					Logger.errorf("invalid row of curated_HPRD_KEGG: %s", line);
					continue;
				}
				int interaction_id = Integer.valueOf(row[0]);
				String ligand = row[1];
				String receptor = row[2];
				String kegg = "NA", kegg_id = "NA";
				if (row[13].split("_").length >= 2) {
					kegg = row[13].split("_")[1];
					kegg_id = row[13].split("_")[0];
				}
				String type = row[14].trim();
				
				// find genes
				Gene ligand_gene = gene_cancer_symbol_db.get(ligand);
				Gene receptor_gene = gene_cancer_symbol_db.get(receptor);
				if (ligand_gene != null && receptor_gene != null) {
					Gene[] ligand_cancer = this.homolog_cancer2cancer_db.get(ligand_gene);
					Gene[] ligand_stroma = this.homolog_cancer2stroma_db.get(ligand_gene);
					Gene[] receptor_cancer = this.homolog_cancer2cancer_db.get(receptor_gene);
					Gene[] receptor_stroma = this.homolog_cancer2stroma_db.get(receptor_gene);
					
					if (ligand_stroma.length == 0)
						Logger.errorf("ligand %s does not have stromal homologue gene. (interaction %d)", ligand, interaction_id);
					if (receptor_stroma.length == 0)
						Logger.errorf("receptor %s does not have stromal homologue gene. (interaction %d)", receptor, interaction_id);
					
					interactions.add(
						new Interaction(
							interaction_id, type, kegg, kegg_id, ligand, receptor,
							ligand_cancer, receptor_cancer,ligand_stroma, receptor_stroma));
				} else {
					if (ligand_gene == null) Logger.errorf("ligand %s in interaction %d is missing in refLink.", ligand, interaction_id);
					if (receptor_gene == null) Logger.errorf("receptor %s in interaction %d is missing in refLink.", receptor, interaction_id);
				}
			}
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		
		this.interactions = interactions.toArray(new Interaction[]{});
		Logger.logf("%d interactions are loaded.", this.interactions.length);		

		return true;
	}
	
	private boolean calculateGCPercent(String refMrna_filename) {
		Logger.logf("\ncalculating GC percents for " + refMrna_filename);
		Option option = Option.getInstance();

		// DEBUG
		long[] gc_dist = new long[101];
		
		try {
			// load cancer refMrna
			FileReader fr = new FileReader(refMrna_filename);
			BufferedReader br = new BufferedReader(fr);
			
			String l;
			String sequence = "";
			String current_refseq = "";
			int processed_refseqs = 0;
			
			while (true) {
				l = br.readLine();
				if (l == null || l.charAt(0) == '>') {
					if (!current_refseq.equals("")) {
						int n = option.read_length; // side range
						
						int[] is_gc_array = new int[sequence.length() + 2*n - 2];
						for (int i=0; i<sequence.length(); i++) {
							if (sequence.charAt(i) == 'g' || sequence.charAt(i) == 'c') {
								is_gc_array[i + n - 1] = 1;
							}
						}

						int[] buffer = new int[n*2-1]; // ring buffer
						int buffer_start = 0;
						
						// initialize buffer
						int[] gc_counts = new int[sequence.length()];
						int sum = 0;
						for (int i=0; i<buffer.length; i++) {					
							buffer[i] = is_gc_array[i];
							sum += buffer[i];
						}
						for (int scan_pos = n-1; scan_pos<n+sequence.length()-1; scan_pos++) {
							gc_counts[scan_pos-(n-1)] = sum;
							if (scan_pos == n+sequence.length()-2) break;
							
							sum -= buffer[buffer_start];
							buffer[buffer_start] = is_gc_array[scan_pos + n];
							sum += buffer[buffer_start];
							buffer_start = (buffer_start+1) % buffer.length;
						}

						Refseq refseq = this.refseq_db.get(current_refseq);
						refseq.gc_percent = new short[refseq.length];
						for (int i=0; i<sequence.length(); i++) {
							int count = Math.min(i+n-1, sequence.length()-1) - Math.max(i-n+1, 0) + 1;
							double gc_ratio = (double) gc_counts[i] / count;
							refseq.gc_percent[i] = (short)(gc_ratio * 100);
							gc_dist[refseq.gc_percent[i]]++;
						}
						
						processed_refseqs++;
						if (processed_refseqs % 10000 == 0) {
							Logger.logf("%d", processed_refseqs);
						}
					}
					if (l != null) {
						current_refseq = l.split(" ")[0].substring(1);
						sequence = "";
					}
				} else {
					sequence = sequence + l.trim().toLowerCase();
				}
				if (l == null) break;
			}
			
			br.close();
			fr.close();
		} catch (Exception e) {
			e.printStackTrace();
			return false;
			
		}
		
		Logger.logf("done.");
		return true;
	}
}
