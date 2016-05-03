package interactome.analysis;

import java.util.Arrays;
import java.util.Map;

import interactome.Logger;
import interactome.Option;
import interactome.data.BioDB;
import interactome.data.Gene;
import interactome.input.GeneInput;
import interactome.input.Input;
import interactome.input.PairedEndInput;
import interactome.input.RefseqInput;
import interactome.input.SingleEndInput;

public class Analysis {
	private static Analysis _instance;
	public static Analysis createInstance() {
		Analysis._instance = new Analysis();
		return Analysis._instance;
	}
	public static Analysis getInstance() {
		return Analysis._instance;
	}
	
	public Input input;
	public DynamicParameters dp;
	
	private Analysis() {
		this.input = null;
	}
	
	public boolean run() {
		// load input reads
		this.input = this.loadInput();
		if (this.input == null) {
			return false;
		}
		
		// calculate dynamic parameters
		this.dp = DynamicParameters.createInstance(this.input);
		dp.selectParameterRefseqs();
		if (!dp.calculateWeightedMappability()) {
			return false;
		}
		
		// execute bias correction
		BiasCorrector bc = BiasCorrector.createInstance(this.input, this.dp);
		if (!bc.correctBias()) {
			return false;
		}
		
		// normalize expression sums
		this.normalize();
		
		return true;
	}
	
	public boolean outputResults() {
		return true;
	}
	
	private Input loadInput() {
		Option option = Option.getInstance();
		
		Input input = null;
		if (option.input_type_paired) {
			input = new PairedEndInput();
		} else {
			input = new SingleEndInput();
		}
		
		if (!input.loadFile()) {
			return null;
		}
		
		return input;
	}
	
	public void normalize() {
		BioDB biodb = BioDB.getInstance();
		
		Logger.logf("normalizing.");
		
		for (Map.Entry<String, Gene> entry : biodb.gene_db.entrySet()) {
			String entrez_id = entry.getKey();
			Gene gene = entry.getValue();
			
			RefseqInput[] refinputs = new RefseqInput[gene.variants.size()];
			for (int i=0; i<gene.variants.size(); i++) {
				refinputs[i] = input.refseq_inputs.get(gene.variants.get(i).refseq_id);
			}
			
			// accept refseq with maximum true_expression
			int maximum = 0;
			for (int i=1; i<gene.variants.size(); i++) {
				if (refinputs[maximum].true_expression < refinputs[i].true_expression) {
					maximum = i;
				}
			}
			
			GeneInput geneinput = input.gene_inputs.get(entrez_id);
			geneinput.representativeRefseq = gene.variants.get(maximum);
			geneinput.representativeExpression = refinputs[maximum].true_expression;
		}
		
		// normalize of cancer gene expressions
		// make sum of the middle 90% genes to 300,000 reads
		double[] expressions = new double[biodb.cancer_entrez_ids.length];
		for (int i=0; i<biodb.cancer_entrez_ids.length; i++) {
			expressions[i] = input.gene_inputs.get(biodb.cancer_entrez_ids[i]).representativeExpression;
		}
		Arrays.sort(expressions);
		double sum = 0;
		for (int i=(int)(expressions.length*0.95); i>expressions.length*0.05; i--) {
			sum += expressions[i];
		}
		if (sum == 0) {
			for (int i=0; i<expressions.length; i++) {
				sum += expressions[i];
			}
		}
		for (String entrez_id : biodb.cancer_entrez_ids) {
			GeneInput ginput = input.gene_inputs.get(entrez_id);
			ginput.normalizedExpression = ginput.representativeExpression * 300000.0 / sum;
		}

		// normalize of stromal gene expressions
		expressions = new double[biodb.stromal_refseq_ids.length];
		for (int i=0; i<biodb.stromal_entrez_ids.length; i++) {
			expressions[i] = input.gene_inputs.get(biodb.stromal_entrez_ids[i]).representativeExpression;
		}
		Arrays.sort(expressions);
		sum = 0;
		for (int i=(int)(expressions.length*0.95); i>expressions.length*0.05; i--) {
			sum += expressions[i];
		}
		if (sum == 0) {
			for (int i=0; i<expressions.length; i++) {
				sum += expressions[i];
			}
		}
		for (String entrez_id : biodb.stromal_entrez_ids) {
			GeneInput ginput = input.gene_inputs.get(entrez_id);
			ginput.normalizedExpression = ginput.representativeExpression * 300000.0 / sum;
		}
		Logger.logf("normalizing done.");
	}
}
