package interactome.data;

import interactome.Logger;
import interactome.Option;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

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
				refseq.tax_id = Integer.valueOf(option.settings.get("cancer_taxonomy"));
				this.refseq_db.put(refseq.refseq_id, refseq);
				cancer_refseq_id_list.add(refseq.refseq_id);
			}			
			br.close();
			fr.close();
			
			fr = new FileReader(option.settings.get("stromal_refNames"));
			br = new BufferedReader(fr);
			while ((line = br.readLine()) != null) {
				Refseq refseq = new Refseq(line.trim());
				refseq.tax_id = Integer.valueOf(option.settings.get("stromal_taxonomy"));
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

		Logger.logf("loading BioDB done.");
		return true;
	}
}
