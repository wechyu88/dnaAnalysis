package analysisLammpsDump;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

public class autoCorScalarData {
	ArrayList<Double> data;
	
	
	
	public void read_data(String path) throws FileNotFoundException {
		File file = new File(path);
		FileReader fr = new FileReader(file);
		
	}
	

}
