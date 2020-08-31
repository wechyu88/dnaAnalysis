package analysisLammpsDump;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class glycerolGr extends waterOtherGr{
	HashMap<Integer,HashMap<Integer,atomGroup>> timegly = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	public static void main(String[] argv) throws IOException {	
/*		String pathInput = argv[0];
		int numberCores= Integer.valueOf(argv[1]);
		int numberAverageStep = Integer.valueOf(argv[2]);
		String pathOutput = System.getProperty("user.dir");
		String pathwaterCor = pathOutput+"/autoWaterCor.txt";
		String pathglyCor = pathOutput+"/autoglyCor.txt";*/

	    //int numberStepBin=2;
		// int numberStepBin = Integer.valueOf(argv[2]);
		int numberStepBin=1;
		int numberAverage = 2;
		double binsize = 0.4;
	    String pathInput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\dump.lammpstrj2";
	    String pathOutput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\";
		String pathglygr = pathOutput+"glygr.txt";
	    
		String path = pathInput;
		PrintWriter pwglygr = new PrintWriter(new File(pathglygr));

		readdumpCustom rdc = new readdumpCustom(path);
		//AutoCorrelation acl = new AutoCorrelation(rdc,numberStepBin);
		glycerolGr glyGr =new glycerolGr(rdc,numberStepBin,binsize);
		glyGr.CalNSnap(numberAverage);
		int counttime=0;
		double sysden=0;
		for(Integer t:glyGr.timeDen.keySet()) {
			counttime++;
			sysden+=glyGr.timeDen.get(t);
		}
		sysden/=counttime;
		double sumcount=0;
		for(Integer i =0;i<10*glyGr.disCount.keySet().size();i++) {
			if(glyGr.disCount.containsKey(i)) {
				double d = (i+0.5)*binsize;
				double count = glyGr.disCount.get(i);
				double vol=4/3*Math.PI*((d+0.5)*(d+0.5)*(d+0.5)-(d-0.5)*(d-0.5)*(d-0.5))*binsize*binsize*binsize;
				double dens = count/vol/sysden;
				sumcount+=count;
				pwglygr.printf("%10.4f \t %10.4f \t %10.8f \t %10.4f \t%10.4f\n",d,count,dens,vol,sumcount);				
				
			}			
		}
		pwglygr.close();
		
	}
	public glycerolGr(readdumpCustom rdc,int nstep,double binsize) throws IOException{
		dumpOneStep dos= rdc.readnextNstep(nstep);
		
		while(dos.timestep!=-1) {
			
			System.out.println("Time: "+dos.timestep+"\n");
			HashMap<Integer,water> waterList = new HashMap<Integer,water>();
			HashMap<Integer,atomGroup> glyList = new HashMap<Integer,atomGroup>();
			
			//ArrayList<atomGroup> waterList = new atomGroup();
			
			for(Atom a:dos.atomlist) {
				if(a.type!=70&&a.type!=71) {
					if(glyList.containsKey(a.moelculeid)) {
						glyList.get(a.moelculeid).addAtom(a);
					}else {
						atomGroup newgroup = new atomGroup();
						newgroup.addAtom(a);
						glyList.put(a.moelculeid, newgroup );
					}						
				}else {
					if(waterList.containsKey(a.moelculeid)) {
						waterList.get(a.moelculeid).addAtom(a);
					}else {
						water newgroup = new water();
						newgroup.addAtom(a);
						waterList.put(a.moelculeid, newgroup );
					}					
				}
				
				
			}
			for(Integer id:waterList.keySet()) {
				waterList.get(id).calcom(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			for(Integer id:glyList.keySet()) {
				glyList.get(id).calcom(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			timelist.add(dos.timestep);
			System.out.println("time is "+dos.timestep);
			timeWater.put(dos.timestep,waterList);
			timegly.put(dos.timestep,glyList);
			double den = glyList.size()/(dos.xhi-dos.xlo)/(dos.yhi-dos.ylo)/(dos.zhi-dos.zlo);
			System.out.println("den is  "+den+" count is "+glyList.size()+" v is "+(dos.xhi-dos.xlo)*(dos.yhi-dos.ylo)*(dos.zhi-dos.zlo)+" xhi "+dos.zhi+" xlo "+dos.zlo);
			timeDen.put(dos.timestep,den);
			grTime.add(CalOneSnap(binsize,dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,glyList));
			dos = rdc.readnextNstep(nstep);
	
		}
	}
	
	

	


}
