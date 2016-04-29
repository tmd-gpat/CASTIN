package interactome.data;

import java.util.ArrayList;

public class Gene {
	public String entrez_id;
	public String symbol;
	
	public ArrayList<Refseq> variants;
	public String tax_id;
	
	public Gene(String entrez_id, String symbol) {
		this.entrez_id = entrez_id;
		this.symbol = symbol;
		this.variants = new ArrayList<Refseq>();
	}
}
