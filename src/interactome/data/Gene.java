package interactome.data;

import java.util.ArrayList;

public class Gene {
	public int entrez_id;
	public String symbol;
	
	public ArrayList<Refseq> variants;
	
	public Gene(int entrez_id, String symbol) {
		this.entrez_id = entrez_id;
		this.symbol = symbol;
		this.variants = new ArrayList<Refseq>();
	}
}
