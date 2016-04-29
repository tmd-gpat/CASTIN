/*
 * Interaction Analysis v3
 * written by Ryohei Suzuki (rsuz.gpat@mri.tmd.ac.jp)
 *
 */

package interactome;

public class Main {
	public static void main(String[] args) {
		// load options
		Option option = Option.createInstance(args);
		if (option == null) {
			return;
		}
		
		// initialize logger
		Logger.initialize();
		Logger.logf("logger initialized.");
		
		// load settings.properties
		if (!option.loadSettingFile()) {
			Logger.close();
			return;
		}
		
		// initialize BioDB
		
	}
}
