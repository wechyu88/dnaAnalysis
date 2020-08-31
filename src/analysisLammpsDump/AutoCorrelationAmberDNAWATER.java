package analysisLammpsDump;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class AutoCorrelationAmberDNAWATER extends AutoCorrelation{
	HashMap<Integer,HashMap<Integer,atomGroup>> timeDNADipole;
	HashMap<Integer,HashMap<Integer,water>> timeWaterDipole;
	public HashMap<Integer, Integer> idMajor = new HashMap<Integer, Integer>();
	public HashMap<Integer, Integer> idMinor = new HashMap<Integer, Integer>();
	public static void main(String[] argv) throws IOException, InterruptedException {	
		
		
		
		
		
		
		
		  String pathInput = argv[0]; int numberCores= Integer.valueOf(argv[1]); int
		  numberAverageStep = Integer.valueOf(argv[2]); String pathOutput =
		  System.getProperty("user.dir"); double rad = Double.valueOf(argv[3]);
		  
		  double averageRatio = Double.valueOf(argv[4]); String outWater = argv[5];
		  String outDNA = argv[6];
		  
		  String pathwaterCor = pathOutput+"/"+outWater; String pathDNACor =
		  pathOutput+"/"+outDNA; String minorMajorPath = argv[7]; int majorOrMinor =
		  Integer.valueOf(argv[8]);
		 
		  AutoCorrelationAmberDNAWATER acDNAWATER = new AutoCorrelationAmberDNAWATER();

		
		/*
		 * String pathInput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\h120a\\dump.head"; String
		 * pathOutput = "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\h120a\\"; int
		 * numberCores=1; int numberAverageStep=1; String outWater = "waterCor.txt";
		 * String outDNA = "DNAcor.txt"; String pathwaterCor = pathOutput+"/"+outWater;
		 * String pathDNACor = pathOutput+"/"+outDNA; String minorMajorPath =
		 * pathOutput+"minorOrMajor.txt"; double rad = 3.6; double averageRatio = 0.8;
		 * int majorOrMinor = 0;
		 */
		  
		// runWaterDNAShell(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,rad,averageRatio);
		  //runWaterDNAShell(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,rad,averageRatio);
		  acDNAWATER.countWaterInShellMinorMajor(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,rad,averageRatio,minorMajorPath,majorOrMinor);
		 // CalWaterInShellTime(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,rad,averageRatio);
		  //runWaterCor(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,averageRatio);
		  
		 acDNAWATER.CalWaterInShellTimeGrooves(pathInput,numberCores,numberAverageStep,outWater,outDNA,pathOutput,rad,averageRatio,minorMajorPath,majorOrMinor);
		
		  /*
		 * String path = pathInput; readdumpCustom rdc = new readdumpCustom(path);
		 * AutoCorrelationAmberDNAWATER acl = new
		 * AutoCorrelationAmberDNAWATER(numberAverageStep,rdc); acl.init2(rad);
		 * HashMap<Integer,Double> waterDipole =
		 * acl.calTimeCor(averageRatio,numberCores,acl.timelist, acl.timeWaterDipole);
		 * HashMap<Integer,Double> DNADipole = acl.calTimeCor2(numberCores,acl.timelist,
		 * acl.timeDNADipole);
		 * 
		 * 
		 * PrintWriter pwwaterCor = new PrintWriter(new File(pathwaterCor)); PrintWriter
		 * pwglyCor = new PrintWriter(new File(pathDNACor)); for(Integer
		 * interval:waterDipole.keySet()) {
		 * pwwaterCor.println(interval+" "+waterDipole.get(interval));
		 * pwglyCor.println(interval+" "+DNADipole.get(interval)); } pwwaterCor.close();
		 * pwglyCor.close();
		 */
		
	}
	public static void countWaterInShell(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double rad,double averageRatio) throws IOException, InterruptedException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		acl.init2(rad);
		
	}
	public  void countWaterInShellMinorMajor(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double rad,double averageRatio,String minorMajorPath,int MajorOrMinor) throws IOException, InterruptedException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;
		
		

		
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		acl.init5(rad,MajorOrMinor,minorMajorPath);
		
	}
	
	public void readMinorMajorAtom(String path) throws IOException {
		File file = new File(path);
		FileReader fr = new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String major = br.readLine();
		String [] majorAtom = major.split(" ");
		for(int i=1;i<majorAtom.length;i++) {
			idMajor.put(Integer.valueOf(majorAtom[i]),1);
		}
		
		String minor = br.readLine();
		String [] minorAtom = minor.split(" ");
		for(int i=1;i<minorAtom.length;i++) {
			idMinor.put(Integer.valueOf(minorAtom[i]),1);
		}		
		
	}
	
	public void CalWaterInShellTimeGrooves(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double rad,double averageRatio,String minorMajorPath,int MajorOrMinor) throws IOException, InterruptedException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;

		readdumpCustom rdc = new readdumpCustom(path);
		
		
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		acl.init4(rad,MajorOrMinor,minorMajorPath);		
		

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
	public static void CalWaterInShellTime(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double rad,double averageRatio) throws IOException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		HashMap<Integer,ArrayList<Integer>> timedis = acl.init3(rad);		
		

		PrintWriter pwwaterCor = new PrintWriter(new File(pathwaterCor));
		PrintWriter pwglyCor = new PrintWriter(new File(pathDNACor));
		for(int key:timedis.keySet()) {
			pwwaterCor.print(key+" ");
			for(int id:timedis.get(key)) {
				pwwaterCor.print(id+" ");
			}
			pwwaterCor.print("\n");
		}
		
		pwwaterCor.close();
		pwglyCor.close();
	}
	public static void runWaterDNAShell(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double rad,double averageRatio) throws IOException, InterruptedException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		acl.init2(rad);
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
	
	public static void runWaterCor(String pathInput, int numberCores,int numberAverageStep,String outWater,String outDNA,String pathOutput,double averageRatio) throws IOException, InterruptedException {
		String path = pathInput;
		String pathwaterCor = pathOutput+"/"+outWater; 
		String pathDNACor = pathOutput+"/"+outDNA;
		readdumpCustom rdc = new readdumpCustom(path);
		AutoCorrelationAmberDNAWATER acl = new AutoCorrelationAmberDNAWATER(numberAverageStep,rdc);
		acl.initAllWater();
		HashMap<Integer,Double> waterDipole =  acl.calTimeCor(averageRatio,numberCores,acl.timelist, acl.timeWaterDipole);
		//HashMap<Integer,Double> DNADipole = acl.calTimeCor2(numberCores,acl.timelist, acl.timeDNADipole);


		PrintWriter pwwaterCor = new PrintWriter(new File(pathwaterCor));
		PrintWriter pwglyCor = new PrintWriter(new File(pathDNACor));
		for(Integer interval:waterDipole.keySet()) {
			pwwaterCor.println(interval+" "+waterDipole.get(interval));
		//	pwglyCor.println(interval+" "+DNADipole.get(interval));
		}
		pwwaterCor.close();
		pwglyCor.close();
		
	}
	
	public AutoCorrelationAmberDNAWATER(int n, readdumpCustom rdc) throws IOException {
		super(n, rdc);
		this.timeDNADipole= new HashMap<Integer,HashMap<Integer,atomGroup>>();
		// TODO Auto-generated constructor stub
		//this.init();
	}
	
	public AutoCorrelationAmberDNAWATER() {
		// TODO Auto-generated constructor stub
		super();
	}
	//@Override
	public void init(double rad) throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		int init = 0;
		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {	
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
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
				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
					waterID.add(id);
				}
				this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
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
				//timelist.add(dos.timestep);
		
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				System.out.println(waterID.size()+" "+newWaterList.size());
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,waterList);

				this.timeDNADipole.put(dos.timestep,DNAlist);
				
				
			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}

		
	}
	public void init2(double rad) throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		int init = 0;
		String outputdir =  System.getProperty("user.dir")+"/";

		String pathwaterCor = outputdir+"countWaterInShell.txt";

		PrintWriter countWaterInShell = new PrintWriter(new File(pathwaterCor));

		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {
					
	
				
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				 
				
				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());
				 
				 
				 
				 
				 
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}		
				System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				
				  for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
				  waterID.add(id); }
				 
				this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());

				
				System.out.println(waterID.size()+" "+newWaterList.size());
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}			
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,newWaterList);

				this.timeDNADipole.put(dos.timestep,DNAlist);
				
				
			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}
		countWaterInShell.close();
		
	}

	public void init5(double rad,int MajorOrMinor,String minorMajorPath) throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		readMinorMajorAtom(minorMajorPath);

		int init = 0;
		String outputdir =  System.getProperty("user.dir")+"/";

		String pathwaterCor = outputdir+"countWaterInShell.txt";

		PrintWriter countWaterInShell = new PrintWriter(new File(pathwaterCor));

		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {
					
	
				
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				//HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				HashMap<Integer,water> newWaterList = waterWithinShell2(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad,MajorOrMinor);

				 
				
				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());
				 
				 
				 
				 
				 
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}		
				System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				
				  for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
				  waterID.add(id); }
				 
				this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				//HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				HashMap<Integer,water> newWaterList = waterWithinShell2(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad,MajorOrMinor);

				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());

				
				System.out.println(waterID.size()+" "+newWaterList.size());
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}			
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,newWaterList);

				this.timeDNADipole.put(dos.timestep,DNAlist);
				
				
			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}
		countWaterInShell.close();
		
	}

	public HashMap<Integer,ArrayList<Integer>> init3(double rad) throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		int init = 0;
		String outputdir =  System.getProperty("user.dir")+"/";
		//HashMap<Integer,ArrayList<Integer>> waterIdToTimes = new HashMap<Integer,ArrayList<Integer>>(); 
		//String pathwaterCor = outputdir+"countWaterInShell.txt";
		HashMap<Integer,Integer> waterIdToTimes = new HashMap<Integer,Integer>(); 
		HashMap<Integer,ArrayList<Integer>> timedis= new HashMap<Integer,ArrayList<Integer>>() ;

		//PrintWriter countWaterInShell = new PrintWriter(new File(pathwaterCor));
		int deltaT=0;
		int firstStep=0;
		int secondStep=1;
		int length=0;
		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {
					
	
				
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				 
				
				//countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());
				 
				 
				 
				 
				 
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);

					waterIdToTimes.put(id, 1);
					
				}
				length++;
				
				ArrayList<Integer> waterIdThisTime = new ArrayList<Integer>();
				waterIdThisTime.add(newWaterList.size());
				for(int tempid:waterIdToTimes.keySet()) {
					waterIdThisTime.add(tempid);
				}
				timedis.put(length,waterIdThisTime);

				
				//timedis.put(length,);
				
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}		
				firstStep = dos.timestep;
				System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				
				  for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
				  waterID.add(id); }
				 
				this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {

				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				//countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());

				
				System.out.println(waterID.size()+" "+newWaterList.size());
				//int overlap=0;
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
					//if(waterIdToTimes.containsKey(id)) {
						//overlap++;
					//}
					/*
					 * if(waterIdToTimes.containsKey(id)) {
					 * waterIdToTimes.get(id).add(dos.timestep); }else { ArrayList<Integer>
					 * timesteps = new ArrayList<Integer>(); timesteps.add(dos.timestep);
					 * waterIdToTimes.put(id,timesteps); }
					 */
				}
				int overlap=0;
				HashMap<Integer, Integer> waterIdToTimesClone = (HashMap<Integer, Integer>) waterIdToTimes.clone(); 
				
				for(Integer id:waterIdToTimes.keySet()) {
					if(newWaterList.containsKey(id)) {
						overlap++;
					}else {
						waterIdToTimesClone.remove(id);
					}
					/*
					 * if(waterIdToTimes.containsKey(id)) {
					 * waterIdToTimes.get(id).add(dos.timestep); }else { ArrayList<Integer>
					 * timesteps = new ArrayList<Integer>(); timesteps.add(dos.timestep);
					 * waterIdToTimes.put(id,timesteps); }
					 */
				}
				waterIdToTimes = waterIdToTimesClone;
				
				length++;
				ArrayList<Integer> waterIdThisTime = new ArrayList<Integer>();
				waterIdThisTime.add(overlap);
				for(int tempid:waterIdToTimesClone.keySet()) {
					waterIdThisTime.add(tempid);
				}
				timedis.put(length,waterIdThisTime);
				if(overlap==0) {
					break;
				}
				if(init==1) {
					secondStep = dos.timestep;
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,newWaterList);

				this.timeDNADipole.put(dos.timestep,DNAlist);
				
				init++;

			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}
		//countWaterInShell.close();
		/*
		 * deltaT = secondStep - firstStep; for(int id:waterIdToTimes.keySet()) {
		 * calculateStayTimeDis(timedis,waterIdToTimes.get(id),deltaT); }
		 */
		return timedis;
		
	}
	
	public void init4(double rad,int MajorOrMinor,String minorMajorPath ) throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		int init = 0;
		readMinorMajorAtom(minorMajorPath);

		String outputdir =  System.getProperty("user.dir")+"/";

		String pathwaterCor = outputdir+"countWaterInShell.txt";

		PrintWriter countWaterInShell = new PrintWriter(new File(pathwaterCor));

		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {
					
	
				
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
				HashMap<Integer,water> newWaterList = waterWithinShell2(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad,MajorOrMinor);
				
				 
				
				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());
				 
				 
				 
				 
				 
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}		
				System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				
				  for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
				  waterID.add(id); }
				 
				this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				HashMap<Integer,water> newWaterList = waterWithinShell2(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad,MajorOrMinor);
				
				countWaterInShell.println("timestep "+dos.timestep+" water number in shell "+rad+" is " +newWaterList.keySet().size());

				
				System.out.println(waterID.size()+" "+newWaterList.size());
				
				for(Integer id:newWaterList.keySet()) {
					newWaterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}			
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,newWaterList);

				this.timeDNADipole.put(dos.timestep,DNAlist);
				
				
			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}
		countWaterInShell.close();
		
	}

	
	
	public static void calculateStayTimeDis(HashMap<Integer,Integer> timedis,ArrayList<Integer> times,int interval) {// HashMap<> key:continuous time length, value, count,
		int length = 1;
		for(int i=1;i<times.size();i++) {
			if(times.get(i)-times.get(0)==interval) {
				length++;
			}else {
				timedis.put(length, timedis.getOrDefault(length, 0)+1);
				length=1;
			}
		}
		
	}

	public void initAllWater() throws IOException {//rad is setting water molecules within radius rad from DNA
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
		this.timeWaterDipole = new 	HashMap<Integer,HashMap<Integer,water>> ();
		int init = 0;
		while(dos.timestep!=-1) {
			System.out.println(dos.xhi);
			if(init==0) {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if(a.type==17||a.type==18) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				
				
				
				////////////////////
	
				timelist.add(dos.timestep);
			//	HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				for(Integer id:waterList.keySet()) {
					waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}		
				//System.out.println(dos.timestep+" "+waterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				this.timeWaterDipole.put(dos.timestep,waterList);
				//this.timeDNADipole.put(dos.timestep,DNAlist);
				init++;
			}else {
				HashMap<Integer,water> waterList = new HashMap<Integer,water>();
				HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
				for(Atom a:dos.atomlist) {
					if(a.type!=17&&a.type!=18&&a.type!=16) {
						if(DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						}else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup );
						}						
					}else if((a.type==17||a.type==18)) {
						if(waterList.containsKey(a.moelculeid)) {
							waterList.get(a.moelculeid).addAtom(a);
							if(a.type==17) {
								waterList.get(a.moelculeid).setOxy(a);
	
							}
							
						}else {
							water newgroup = new water();
							newgroup.addAtom(a);
							if(a.type==17) {
								newgroup.setOxy(a);
							}						
	
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}

				//timelist.add(dos.timestep);
		
				//HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				
				for(Integer id:waterList.keySet()) {
					waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:DNAlist.keySet()) {
					DNAlist.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}			
				////////////////////
	
				timelist.add(dos.timestep);
				this.timeWaterDipole.put(dos.timestep,waterList);

			//	this.timeDNADipole.put(dos.timestep,DNAlist);
				
				
			}
			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}

		
	}

	public HashMap<Integer,water> waterWithinShell(double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> DNAMol, double rad) {
		HashMap<Integer,water> newWaterList = new HashMap<Integer,water> ();
		ArrayList<Double> dislist = new ArrayList<Double>();
		for(Integer gid:DNAMol.keySet()) {//DNAMol
			atomGroup twoDna = DNAMol.get(gid);

			for(Integer wid:waterMol.keySet()) {
				
				water oneWater=waterMol.get(wid);
				

				double dis = waterDNAgr.caldis(oneWater,twoDna,xl,xh,yl,yh,zl,zh);
				dislist.add(dis);
				if(dis<rad) {
					newWaterList.put(wid,oneWater);
				}
			}
		}
		return newWaterList;
		
	}
	public HashMap<Integer,water> waterWithinShell2(double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> DNAMol, double rad,int majorOrMinor) {
		HashMap<Integer,water> newWaterList = new HashMap<Integer,water> ();
		for(Integer gid:DNAMol.keySet()) {//DNAMol
			atomGroup twoDna = DNAMol.get(gid);

			for(Integer wid:waterMol.keySet()) {
				
				water oneWater=waterMol.get(wid);
				

				double dis = caldis(oneWater,twoDna,xl,xh,yl,yh,zl,zh,majorOrMinor);

				if(dis<rad&&dis>0) {
					newWaterList.put(wid,oneWater);
				}
			}
		}
		return newWaterList;
		
	}	
	
	public double caldis(water a,atomGroup b,double xl,double xh,double yl,double yh,double zl,double zh,int majorOrMinor) {
		double dis=0.0;
		Atom mindatom = new Atom ();
		double mindis = 0.5*Integer.MAX_VALUE;
		int count=0;
		int closestId=0;

		for(Atom batm:b.atomlist) {
			//if(batm.mass<13&&batm.type<14) {
			//	continue;
			//}
			if(idMajor.containsKey(batm.id)||idMinor.containsKey(batm.id)) {
				
			}else {
				continue;
			}
			count++;
			/*
			 * if(majorOrMinor==1) {// major groove if(idMajor.containsKey(batm.id)) {
			 * 
			 * }else { continue; } }else { if(idMinor.containsKey(batm.id)) {
			 * 
			 * }else { continue; }
			 * 
			 * }
			 */
			double dx = a.Oxygen.x-batm.x;
			double dy = a.Oxygen.y-batm.y;
			double dz = a.Oxygen.z-batm.z;
			dx = dx - (xh-xl)*Math.round(dx/(xh-xl));
			dy = dy - (yh-yl)*Math.round(dy/(yh-yl));
			dz = dz - (zh-zl)*Math.round(dz/(zh-zl));
			dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	
			if(dis<mindis&&dis>0) {
				mindis = dis;
				closestId = batm.id;
			}
			
			if(dis<0.4) {
				System.out.println("find wierd");
			}
		}
		if(majorOrMinor==1) {
			if(idMajor.containsKey(closestId)) {
				
			}else {
				mindis = -1;
			}
			
		}else {
			if(idMinor.containsKey(closestId)) {
				
			}else {
				mindis = -1;
			}
			
		}


		return mindis;
	}
	

}
