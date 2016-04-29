package interactome;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Logger {
	private static FileWriter _fw;
	private static BufferedWriter _bw;
	private static boolean _active = false;
	
	public static void initialize() {
		String logFileName = Option.getInstance().output_path + "/"
				+ (new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss")).format(new Date()) + ".log";
		try {
			_fw = new FileWriter(logFileName);
			_bw = new BufferedWriter(_fw);
			_active = true;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void logf(String format, Object... args) {
		System.out.printf(format + "\n", args);

		if (!_active) return;
		try {
			_bw.write(String.format(format + "\n", args));
			_bw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void errorf(String format, Object... args) {
		System.err.printf(format + "\n", args);

		if (!_active) return;
		try {
			_bw.write(String.format(format + "\n", args));
			_bw.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
	
	public static void close() {
		try {
			_bw.close();
			_fw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
