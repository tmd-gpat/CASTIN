package interactome.analysis;

import interactome.input.Input;

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
		return null;
	}
}
