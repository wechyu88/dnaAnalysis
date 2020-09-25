package analysisLammpsDump;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class waterListResTime {
	ArrayList<waterResTime> ListWaterRes;
	HashMap<Integer,Double> PtList;
	long [] sumBList;
	int totalWater;
	public static void main(String [] argv) throws IOException {
		String pathInput = argv[0];
		String pathOut = argv[1];
		int timesteps = Integer.valueOf(argv[2]);
		waterListResTime wlrt = new waterListResTime();
		wlrt.run(pathInput,pathOut,timesteps);
		
	}
	public void run(String pathInput,String pathOut, int timesteps) throws NumberFormatException, IOException {
		//String pathInput = "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\h120a\\residenceTime\\T298\\waterInMajor.txt";
		File infile = new File(pathInput);
		FileReader fr = new FileReader(infile);
		BufferedReader br = new BufferedReader(fr);
		String line;
		int temptime=0;
		HashMap<Integer,waterResTime> waterIdToRes = new HashMap<Integer,waterResTime>();
		while((line=br.readLine())!=null) {
			temptime++;
			if(temptime>timesteps) {
				break;
			}
			String [] molids = line.split(" ");
			for(int i=2;i<molids.length;i++) {
				int id = Integer.valueOf(molids[i]);

				if(waterIdToRes.containsKey(id)) {
					waterResTime tempRes = waterIdToRes.get(id);
					tempRes.AList[temptime-1]=1;
				}else {
					waterResTime tempRes = new waterResTime(timesteps,id);
					tempRes.AList[temptime-1]=1;
					waterIdToRes.put(id,tempRes);
				}
			}
		}
		ArrayList<waterResTime> wrlist = new ArrayList<waterResTime>(); 
		for(int key: waterIdToRes.keySet()) {
			waterIdToRes.get(key).setAllList();
			wrlist.add(waterIdToRes.get(key));
		}
		waterListResTime wlrt = new waterListResTime(wrlist);
		//String pathOut = "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\h120a\\residenceTime\\T298\\resInMajor.txt";
		File outfile = new File(pathOut);
		PrintWriter pout = new PrintWriter(outfile);
		
		for(int i=0;i<timesteps-1;i++) {
			long tempnum = wlrt.sumBList[i]-wlrt.sumBList[i+1];
			pout.println(i+1+" "+wlrt.PtList.get(i+1)+" "+tempnum+" "+wlrt.totalWater);
		}		
		pout.close();
		
	}
	public waterListResTime() {
		
	}
	public waterListResTime(ArrayList<waterResTime> waterRes) {
		this.ListWaterRes = waterRes;
		PtList = new HashMap<Integer,Double>();
		this.sumBList = new long[waterRes.get(0).BList.length];
		this.totalWater = waterRes.size();
		for(int m=1;m<ListWaterRes.get(0).BList.length;m++) {
			double tempPt = 0.0;
			for(waterResTime tempWater:ListWaterRes) {
				tempPt+=tempWater.BList[m-1]*1.0/(ListWaterRes.get(0).BList.length-m);
				this.sumBList[m-1]+=tempWater.BList[m-1];
				//System.out.println(tempPt);
			}
			PtList.put(m,tempPt);
		}
	}
	
}
