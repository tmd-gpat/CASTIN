package interactome.analysis;

import interactome.input.Input;

public class ResultsWriter {
	public static ResultsWriter _instance;
	public static ResultsWriter getInstance() {
		return _instance;
	}
	public static ResultsWriter createInstance(Input input, DynamicParameters dynamicParameters) {
		_instance = new ResultsWriter(input, dynamicParameters);
		return _instance;
	}
	
	Input input;
	DynamicParameters dynamicParameters;
	
	public ResultsWriter(Input input, DynamicParameters dynamicParameters) {
		this.input = input;
		this.dynamicParameters = dynamicParameters;
	}
	
	public void writeKEGGHPRDresult() {
		
	}
}
