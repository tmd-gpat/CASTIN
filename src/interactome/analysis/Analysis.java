package interactome.analysis;

import interactome.Logger;
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
		this.input = this.loadInput();
		if (this.input == null) {
			return false;
		}
		return true;
	}
	
	public boolean outputResults() {
		return true;
	}
	
	private Input loadInput() {
		Option option = Option.getInstance();
		
		Input input = null;
		if (option.input_type_paired) {
			Logger.logf("\nloading input file (paired-ended).");
			input = new PairedEndInput();
		} else {
			Logger.logf("\nloading input file (single-ended).");
			input = new SingleEndInput();
		}
		
		if (!input.loadFile()) {
			return null;
		}
		
		return input;
	}
}
