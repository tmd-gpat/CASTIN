package interactome.analysis;

import interactome.Option;
import interactome.input.Input;
import interactome.input.PairedEndInput;
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
		DynamicParameters dp = DynamicParameters.createInstance(this.input);
		dp.selectParameterRefseqs();
		
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
}
