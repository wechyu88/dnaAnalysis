package analysisLammpsDump;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

public class AutoCorrelationCharmDNAWATER extends AutoCorrelation{
	HashMap<Integer,HashMap<Integer,atomGroup>> timeDNADipole;
	public static void main(String[] argv) throws IOException, InterruptedException {	
		
		
		  String pathInput = argv[0]; int numberCores= Integer.valueOf(argv[1]); int
		  numberAverageStep = Integer.valueOf(argv[2]); String pathOutput =
		  System.getProperty("user.dir"); String pathwaterCor =
		  pathOutput+"/autoWaterCor.txt"; String pathDNACor =
		  pathOutput+"/autoDNACor.txt";
		 

		/*
		 * 
		 * int numberCores=1; int numberAverageStep=1; String pathInput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\snapshots"; String pathOutput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\"; String pathwaterCor =
		 * pathOutput+"autoWaterCor.txt"; String pathDNACor =
		 * pathOutput+"autoglyCor.txt";
		 */
		 double averageRatio =  Double.valueOf(argv[3]);
		String path = pathInput;
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationCharmDNAWATER acl = new AutoCorrelationCharmDNAWATER(numberAverageStep,rdc);
		acl.init();
		HashMap<Integer,Double> waterDipole =  acl.calTimeCor(averageRatio,numberCores,acl.timelist, acl.timeWaterDipole);
		HashMap<Integer,Double> DNADipole = acl.calTimeCor2(numberCores,acl.timelist, acl.timeDNADipole);

		PrintWriter pwwaterCor = new PrintWriter(new File(pathwaterCor));
		PrintWriter pwglyCor = new PrintWriter(new File(pathDNACor));
		for(Integer interval:waterDipole.keySet()) {
			pwwaterCor.println(interval+" "+waterDipole.get(interval));
			pwglyCor.println(interval+" "+DNADipole.get(interval));
		}
		pwwaterCor.close();
		pwglyCor.close();
		
	}
	
	public AutoCorrelationCharmDNAWATER(int n, readdumpCustom rdc) throws IOException {
		super(n, rdc);
		this.timeDNADipole= new HashMap<Integer,HashMap<Integer,atomGroup>>();
		// TODO Auto-generated constructor stub
		//this.init();
	}
	
	@Override
	public void init() throws IOException {
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		while(dos.timestep!=-1) {
			
			System.out.println("Time: "+dos.timestep+"\n");
			HashMap<Integer,water> waterList = new HashMap<Integer,water>();
			HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
		//	HashMap<Integer,atomGroup> saltlist = new HashMap<Integer,atomGroup>();
			
			//ArrayList<atomGroup> waterList = new atomGroup();
			atomGroup dnagroup = new atomGroup();
			DNAlist.put(1, dnagroup );	
		
			for(Atom a:dos.atomlist) {
				if(a.type!=70&&a.type!=71) {
					DNAlist.get(1).addAtom(a);				
				}else if(a.type==70||a.type==71){
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
				waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			for(Integer id:DNAlist.keySet()) {
				DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			timelist.add(dos.timestep);
			System.out.println("time is "+dos.timestep);
			timeWaterDipole.put(dos.timestep,waterList);
			for(int a:DNAlist.keySet()) {
				System.out.println(" dna "+a+" "+DNAlist.get(a).dipolex+" "+DNAlist.get(a).dipoley+" "+DNAlist.get(a).dipolez);
			}
			this.timeDNADipole.put(dos.timestep,DNAlist);

			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}

		
	}
	
	
	

}
