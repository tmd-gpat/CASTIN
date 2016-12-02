package interactome;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Properties;

import gnu.getopt.*;

public class Option {
	private static Option _instance;
	public static Option getInstance() { return Option._instance; }
	
	// command-line options
	public String output_path = null;
	
	public String input_file_single = null;
	public String input_prefix_paired = null;
	public boolean input_type_paired = false;
	public int input_paired_length = 0;

	public int read_length = 50; 		// only for paired-ended mode
	public int directional_mode = 0;	// 0: count all, 1: Fwd-Rev, 2: Rev-Fwd

	public boolean output_mappability = false;
	public boolean output_gcpercent = false;
	
	public boolean output_cancer_fastq = false;
	public boolean output_stromal_fastq = false;
	
	// analysis settings (e.g., parameter files)
	public String setting_file_name = "settings.properties";
	public HashMap<String, String> settings;
	
	private Option() {
		this.settings = new HashMap<String, String>();
	}
	
	/*
	 * parse command-line arguments.
	 */
	public static Option createInstance(String[] args) {
		Option instance = new Option();
		
		Getopt options = new Getopt("Interaction", args, "s:p:o:mgd:l:yz");
		int c;
		
		while ((c = options.getopt()) != -1) {
			switch (c) {
			case 's':	// single-end input prefix
				instance.input_file_single = options.getOptarg();
				instance.input_type_paired = false;
				Logger.logf("input mode: single-ended");
				Logger.logf("input file: %s", instance.input_file_single);
				break;
			case 'p':	// paired-end input prefix
				instance.input_prefix_paired = options.getOptarg();
				instance.input_type_paired = true;
				Logger.logf("input mode: paired-ended");
				Logger.logf("input prefix: %s", instance.input_prefix_paired);
				break;
			case 'o':	// output dir
				instance.output_path = options.getOptarg();
				Logger.logf("output path: %s", instance.output_path);
				break;
			case 'm':	// mappability dump
				instance.output_mappability = true;
				Logger.logf("mappability: output");
				break;
			case 'g':	// gc-percent dump
				instance.output_gcpercent = true;
				Logger.logf("gc%: output");
				break;
			case 'd':	// directional mode
				instance.directional_mode = Integer.valueOf(options.getOptarg());
				String strand_str = "no";
				if (instance.directional_mode == 1) strand_str = "(forward, reverse)";
				if (instance.directional_mode == 2) strand_str = "(reverse, forward)";
				Logger.logf("\tdirectional mode:\t%s\n", strand_str);
				break;
			case 'l':	// read length (paired-end input only) ... to be deprecated
				instance.input_paired_length = Integer.valueOf(options.getOptarg());
				Logger.logf("read length (paired-end): %d", instance.input_paired_length);
				break;
			case 'y':   // output sorted tumor fastq
				instance.output_cancer_fastq = true;
				Logger.logf("output tumor fastq");
				break;
			case 'z':	// output sorted stromal fastq
				instance.output_stromal_fastq = true;
				Logger.logf("output stromal fastq");
				break;
			}
		}
		
		// error if both or none of (-p, -s) was specified
		if (!((instance.input_file_single == null) ^ (instance.input_prefix_paired == null))) {
			Logger.errorf("you should specify one of -s or -p option.");
			return null;
		}
		
		// error if both -s and -l were specified
		if (instance.input_file_single != null && instance.input_paired_length != 0) {
			Logger.errorf("read length (-l) cannot be specified in single-end mode.");
			return null;
		}
		/*
		if (instance.input_prefix_paired != null && instance.input_paired_length == 0) {
			Logger.errorf("read length (-l) should be specified in paired-end mode.");
			return null;
		}*/
		
		// error if output path was not specified
		if (instance.output_path == null) {
			Logger.errorf("you should specify output path by -o option.");
			return null;
		}
		
		Option._instance = instance;
		return Option._instance;
	}
	
	/*
	 * load settings.parameters
	 */
	public boolean loadSettingFile() {
		Logger.logf("\nstart loading parameters.");
		
		Properties config = new Properties();
		try {
			InputStream istream = new FileInputStream(new File(this.setting_file_name));
			config.load(istream);
			
			String[] required_items = {
				"refLink",
				"cancer_taxonomy",
				"cancer_refNames",
				"cancer_refSeqLen",
				"cancer_refMrna",
				"stromal_taxonomy",
				"stromal_refNames",
				"stromal_refSeqLen",
				"stromal_refMrna",
				"homologene",
				"self_matching",
				"self_matching_direct",
				"curated_HPRD_KEGG",
			};
			HashMap<String, String> optional_items = new HashMap<String, String>();
			optional_items.put("parameter_gene_minimum_length", "8000");
			
			for (String item : required_items) {
				String value = config.getProperty(item);
				if (value == null) {
					Logger.errorf("%s is missing in property file", item);
					return false;
				}
				this.settings.put(item, value);
				Logger.logf("%s: %s", item, value);
			}
			
			for (String key : optional_items.keySet()) {
				String value = config.getProperty(key);
				if (value == null) { // default value
					this.settings.put(key, optional_items.get(key));
				} else { // user-specified value
					this.settings.put(key, value);
					Logger.logf("%s: %s", key, value);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}

		Logger.logf("loading parameters done.");
		return true;
	}
}

