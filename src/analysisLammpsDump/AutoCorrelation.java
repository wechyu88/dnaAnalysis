package analysisLammpsDump;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;



public class AutoCorrelation {
	ArrayList<Integer> timelist = new ArrayList<Integer>();
	HashMap<Integer,Double[]> timeDipole = new HashMap<Integer,Double[]>();
	HashMap<Integer,HashMap<Integer,water>> timeWaterDipole = new HashMap<Integer,HashMap<Integer,water>>();
	HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	readdumpCustom rdc;
	int averageN;// combine several(averageN) snapshots and average the coordinates to be one snapshots
	int everyK;
	class MyThread implements Callable<Double>{
		 	ArrayList<Integer> timelist;
		 	HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole;
		 	int start;
		 	int end;
		 	int deltaT;
			public MyThread(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole,int start,int end,int deltaT){
	
			}
			@Override
			public Double call() throws IOException   {
				Double cor = calCor2(timelist, timeDipole,start,end, deltaT);
				return cor;
			}
			
		}	
	
	public static void main(String[] argv) throws IOException, InterruptedException {	
		
		String pathInput = argv[0];
		int numberCores= Integer.valueOf(argv[1]);
		int numberAverageStep = Integer.valueOf(argv[2]);
		String pathOutput = System.getProperty("user.dir");
		String pathwaterCor = pathOutput+"/autoWaterCor.txt";
		String pathglyCor = pathOutput+"/autoglyCor.txt";
		double averageRatio = Double.valueOf(argv[3]);
	    //int numberStepBin=2;
		// int numberStepBin = Integer.valueOf(argv[2]);
/*		int numberCores=1;
		int numberAverageStep=1;
	    String pathInput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\dump.lammpstrj2";
	    String pathOutput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\";
		String pathwaterCor = pathOutput+"autoWaterCor.txt";
		String pathglyCor = pathOutput+"autoglyCor.txt";*/
	    
		String path = pathInput;
		readdumpCustom rdc = new readdumpCustom(path);
		//AutoCorrelation acl = new AutoCorrelation(rdc,numberStepBin);
		AutoCorrelation acl = new AutoCorrelation(numberAverageStep,rdc);
		acl.init();
		HashMap<Integer,Double> waterDipole =  acl.calTimeCor(averageRatio,numberCores,acl.timelist, acl.timeWaterDipole);

		HashMap<Integer,Double> glyDipole = acl.calTimeCor2(numberCores,acl.timelist, acl.timeglyDipole);

		PrintWriter pwwaterCor = new PrintWriter(new File(pathwaterCor));
		PrintWriter pwglyCor = new PrintWriter(new File(pathglyCor));
		for(Integer interval:waterDipole.keySet()) {
			pwwaterCor.println(interval+" "+waterDipole.get(interval));
			pwglyCor.println(interval+" "+glyDipole.get(interval));
		}
		pwwaterCor.close();
		pwglyCor.close();
		
	}
	public AutoCorrelation() {
		
	}
	public AutoCorrelation(readdumpCustom rdc,int nstep) throws IOException {
       // BufferedWriter pw1 = new BufferedWriter(new FileWriter(path5+"/dataStockMayer.txt")); 
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
				waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			for(Integer id:glyList.keySet()) {
				glyList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			timelist.add(dos.timestep);
			System.out.println("time is "+dos.timestep);
			timeWaterDipole.put(dos.timestep,waterList);
			timeglyDipole.put(dos.timestep,glyList);

			dos = rdc.readnextNstep(nstep);
	
		}
	}
	public void setEveryK(int k) {
		this.everyK = k;
	}
	public AutoCorrelation(int n,readdumpCustom rdc) throws IOException {
	       // BufferedWriter pw1 = new BufferedWriter(new FileWriter(path5+"/dataStockMayer.txt")); 
			this.averageN=n;
			this.rdc=rdc;
			//init();
	}
	
	
	public void init() throws IOException {
		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
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
				waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			for(Integer id:glyList.keySet()) {
				glyList.get(id).calglydipole();
			}
			timelist.add(dos.timestep);
			System.out.println("time is "+dos.timestep);
			timeWaterDipole.put(dos.timestep,waterList);
			timeglyDipole.put(dos.timestep,glyList);

			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
	
		}
		
	}
	
	

	public AutoCorrelation(String path) throws IOException {
		File file = new File(path);
		FileReader fr = new FileReader(file);
		BufferedReader br = new BufferedReader(fr);
		String line;
		line = br.readLine();
		while(line != null ) {
			String [] fields = line.split(" ");
			Double [] dipole = new Double [3];
			dipole[0]=Double.valueOf(fields[1]);
			dipole[1]=Double.valueOf(fields[2]);
			dipole[2]=Double.valueOf(fields[3]);
			int timestep = Integer.valueOf(fields[0]);
			timeDipole.put(timestep, dipole);
			timelist.add(timestep);
			line=br.readLine();
		}
	}

	
	  public HashMap<Integer,Double> calTimeCor(double averageRatio,int nprocForLoop,ArrayList<Integer>  timeStep,HashMap<Integer,HashMap<Integer,water>> timeDipole) throws
	  InterruptedException { 
			int length = timeStep.size();
			int start = 0;
			int end = start+(int)(averageRatio*length)-1;
			HashMap<Integer,Double> timeCor= new HashMap<Integer,Double>();
			//int nprocForLoop=2;
			//ExecutorService execs2 = Executors.newFixedThreadPool(nprocForLoop);

			for(int interval=0;interval<length-end;interval++) {
				//execs2.submit(new MyThread(timeStep,timeDipole,start,end,interval));

				timeCor.put(interval, calCor(timeStep,timeDipole,start,end,interval));
				System.out.println(" interval "+interval +" finish ");
			}
			//execs2.shutdown();
			//while (!execs2.awaitTermination(24L, TimeUnit.HOURS)) {
			   // System.out.println("Not yet. Still waiting for termination");
			//}
			return timeCor;
		}
	  public HashMap<Integer,Double> calVelTimeCor(double averageRatio,int nprocForLoop,ArrayList<Integer>  timeStep,HashMap<Integer,HashMap<Integer,water>> timeDipole) throws
	  InterruptedException { 
			int length = timeStep.size();
			int start = 0;
			int end = start+(int)(averageRatio*length)-1;
			HashMap<Integer,Double> timeCor= new HashMap<Integer,Double>();
			//int nprocForLoop=2;
			//ExecutorService execs2 = Executors.newFixedThreadPool(nprocForLoop);

			for(int interval=0;interval<length-end;interval++) {
				//execs2.submit(new MyThread(timeStep,timeDipole,start,end,interval));

				timeCor.put(interval, calVelCor(timeStep,timeDipole,start,end,interval));
				System.out.println(" interval "+interval +" finish ");
			}
			//execs2.shutdown();
			//while (!execs2.awaitTermination(24L, TimeUnit.HOURS)) {
			   // System.out.println("Not yet. Still waiting for termination");
			//}
			return timeCor;
		}
	  public HashMap<Integer,Double> calVelTimeCorHatom(double averageRatio,int nprocForLoop,ArrayList<Integer>  timeStep,HashMap<Integer,HashMap<Integer,water>> timeDipole) throws
	  InterruptedException { 
			int length = timeStep.size();
			int start = 0;
			int end = start+(int)(averageRatio*length)-1;
			HashMap<Integer,Double> timeCor= new HashMap<Integer,Double>();
			//int nprocForLoop=2;
			//ExecutorService execs2 = Executors.newFixedThreadPool(nprocForLoop);

			for(int interval=0;interval<length-end;interval++) {
				//execs2.submit(new MyThread(timeStep,timeDipole,start,end,interval));

				timeCor.put(interval, calVelCorHatom(timeStep,timeDipole,start,end,interval));
				System.out.println(" interval "+interval +" finish ");
			}
			//execs2.shutdown();
			//while (!execs2.awaitTermination(24L, TimeUnit.HOURS)) {
			   // System.out.println("Not yet. Still waiting for termination");
			//}
			return timeCor;
		}
	  public HashMap<Integer,Double> calVelTimeCorDNA(int sideChainOrBack,double averageRatio,int nprocForLoop,ArrayList<Integer>  timeStep,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole) throws
	  InterruptedException { 
			int length = timeStep.size();
			int start = 0;
			int end = start+(int)(averageRatio*length)-1;
			HashMap<Integer,Double> timeCor= new HashMap<Integer,Double>();
			//int nprocForLoop=2;
			//ExecutorService execs2 = Executors.newFixedThreadPool(nprocForLoop);

			for(int interval=0;interval<length-end;interval++) {
				//execs2.submit(new MyThread(timeStep,timeDipole,start,end,interval));

				timeCor.put(interval, calVelCorDNA(sideChainOrBack,timeStep,timeDipole,start,end,interval));
				System.out.println(" interval "+interval +" finish ");
			}
			//execs2.shutdown();
			//while (!execs2.awaitTermination(24L, TimeUnit.HOURS)) {
			   // System.out.println("Not yet. Still waiting for termination");
			//}
			return timeCor;
		}
	  
	public HashMap<Integer,Double> calTimeCor2(int nprocForLoop,ArrayList<Integer> timeStep,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole) throws InterruptedException {
		int length = timeStep.size();
		int start = 0;
		int end = start+(int)(0.5*length)-1;
		HashMap<Integer,Double> timeCor= new HashMap<Integer,Double>();
		//int nprocForLoop=2;
		//ExecutorService execs2 = Executors.newFixedThreadPool(nprocForLoop);

		for(int interval=0;interval<length-end;interval++) {
			//execs2.submit(new MyThread(timeStep,timeDipole,start,end,interval));

			timeCor.put(interval, calCor2(timeStep,timeDipole,start,end,interval));
			System.out.println(" interval "+interval +" finish ");
		}
		//execs2.shutdown();
		//while (!execs2.awaitTermination(24L, TimeUnit.HOURS)) {
		   // System.out.println("Not yet. Still waiting for termination");
		//}
		return timeCor;
	} 
	public double calCor(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,water>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorDeltT(timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			System.out.println(" dt= "+deltaT+" overlap molecule "+corlist.size());

			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	
	public double calVelCor(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,water>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorVelDeltT(timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	public double calVelCorHatom(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,water>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorVelDeltTHatom(timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	public double calVelCorDNA(int sideChainOrBack,ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		
		
		
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorVelDeltTDNA(sideChainOrBack,timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	
	public double calCor2(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorDeltT2(timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	
	public HashMap<Integer,Double> calcorDeltT(HashMap<Integer,water> time1group,HashMap<Integer,water> time2group){
		HashMap<Integer,Double>  corlist = new HashMap<Integer,Double> ();
		for (Integer molid:time1group.keySet()) {
			if(time2group.containsKey(molid)) {
				
			}else {
				continue;
			}
			double tempCor=dipoleMultiply(time1group.get(molid),time2group.get(molid));
			
			corlist.put(molid, tempCor);
		}
		return corlist;
	}
	public HashMap<Integer,Double> calcorVelDeltT(HashMap<Integer,water> time1group,HashMap<Integer,water> time2group){
		HashMap<Integer,Double>  corlist = new HashMap<Integer,Double> ();
		for (Integer molid:time1group.keySet()) {
			double tempCor=VelMultiply(time1group.get(molid),time2group.get(molid));
			corlist.put(molid, tempCor);
		}
		return corlist;
	}
	public HashMap<Integer,Double> calcorVelDeltTHatom(HashMap<Integer,water> time1group,HashMap<Integer,water> time2group){
		HashMap<Integer,Double>  corlist = new HashMap<Integer,Double> ();
		for (Integer molid:time1group.keySet()) {
			double tempCor=VelMultiplyHatom(time1group.get(molid),time2group.get(molid));
			corlist.put(molid, tempCor);
		}
		return corlist;
	}
	
	public HashMap<Integer,Double> calcorVelDeltTDNA(int sideChainOrBack,HashMap<Integer,atomGroup> time1group,HashMap<Integer,atomGroup> time2group){
		HashMap<Integer,Double>  corlist = new HashMap<Integer,Double> ();

		for (Integer molid:time1group.keySet()) {
			//double tempCor=VelMultiplyDNA(time1group.get(molid),time2group.get(molid));
			
			if(time2group.containsKey(molid)) {
				
			}else {
				continue;
			}
			double tempCor=0.0;
			
			if(time1group.get(molid).atomlist.size()==3) {
				
				tempCor=VelMultiplyDNA(time1group.get(molid),time2group.get(molid));

			}else {
				if(sideChainOrBack==0) {// 0 is side chain
					tempCor=VelMultiplyDNASideChain(time1group.get(molid),time2group.get(molid));
					
				}else if(sideChainOrBack==1) {// 1 is backbone
					tempCor=VelMultiplyDNABackBone(time1group.get(molid),time2group.get(molid));
					
				}else if(sideChainOrBack==2) {// 2 is all atom on DNA
					tempCor=VelMultiplyDNA(time1group.get(molid),time2group.get(molid));
					
				}
			}

			
			corlist.put(molid, tempCor);
		}
		return corlist;
	}
	public HashMap<Integer,Double> calcorDeltT2(HashMap<Integer,atomGroup> time1group,HashMap<Integer,atomGroup> time2group){
		HashMap<Integer,Double>  corlist = new HashMap<Integer,Double> ();
		for (Integer molid:time1group.keySet()) {
			double tempCor=dipoleMultiply(time1group.get(molid),time2group.get(molid));
			corlist.put(molid, tempCor);
		}
		return corlist;
	}
	public double dipoleMultiply(atomGroup g1,atomGroup g2) {
		double result=0.0;
		result = g1.dipolex*g2.dipolex+g1.dipoley*g2.dipoley+g1.dipolez*g2.dipolez;
		return result;
	}
	public double VelMultiply(atomGroup g1,atomGroup g2) {
		double result=0.0;
		/*
		 * g1.setIdToAtom(); g2.setIdToAtom(); for(int id:g1.idToAtom.keySet()) {
		 * 
		 * Atom a1=g1.idToAtom.get(id); Atom a2=g2.idToAtom.get(id); if(a1.type==17) {
		 * 
		 * }else { continue; }
		 * 
		 * result +=
		 * (a1.vx*a2.vx+a1.vy*a2.vy+a1.vz*a2.vz)/(lengthOfVec(a1.vx,a1.vy,a1.vz)*
		 * lengthOfVec(a2.vx,a2.vy,a2.vz));
		 * 
		 * }
		 */
		result = (g1.coVx*g2.coVx+g1.coVy*g2.coVy+g1.coVz*g2.coVz)/(lengthOfVec(g1.coVx,g1.coVy,g1.coVz)*lengthOfVec(g2.coVx,g2.coVy,g2.coVz));
		return result;
	}
	
	public double VelMultiplyHatom(atomGroup g1,atomGroup g2) {
		double result=0.0;
		
		  g1.setIdToAtom(); 
		  g2.setIdToAtom(); 
		  for(int id:g1.idToAtom.keySet()) {
		  
			  Atom a1=g1.idToAtom.get(id);
			  Atom a2=g2.idToAtom.get(id);
			  if(a1.type==18) {
		  
			  }else { 
				  continue; 
			  
			  }
		  
			  result +=
			  (a1.vx*a2.vx+a1.vy*a2.vy+a1.vz*a2.vz)/(lengthOfVec(a1.vx,a1.vy,a1.vz)*
			  lengthOfVec(a2.vx,a2.vy,a2.vz));
		  
		  }
		 
		//result = (g1.coVx*g2.coVx+g1.coVy*g2.coVy+g1.coVz*g2.coVz)/(lengthOfVec(g1.coVx,g1.coVy,g1.coVz)*lengthOfVec(g2.coVx,g2.coVy,g2.coVz));
		return result/2;
	}
	public double VelMultiplyDNABackBone(atomGroup g1,atomGroup g2) {
		double result=0.0;
		
		  g1.setIdToAtom(); 
		  g2.setIdToAtom(); 
		  int atomCount=0;
		  for(int id:g1.idToAtom.keySet()) {
		  
			  Atom a1=g1.idToAtom.get(id);
			  Atom a2=g2.idToAtom.get(id);
			  if(a1.type!=6&&a1.type!=7&&a1.type!=8&&a1.type!=10&&a1.type!=11&&a1.type!=15&&a1.type!=16) {//a1 is not side chain
		  
			  }else { 
				  continue; 
			  
			  }
			  atomCount++;
		  
			  result +=
			  (a1.vx*a2.vx+a1.vy*a2.vy+a1.vz*a2.vz)/(lengthOfVec(a1.vx,a1.vy,a1.vz)*
			  lengthOfVec(a2.vx,a2.vy,a2.vz));
		  
		  }
		 
		//result = (g1.coVx*g2.coVx+g1.coVy*g2.coVy+g1.coVz*g2.coVz)/(lengthOfVec(g1.coVx,g1.coVy,g1.coVz)*lengthOfVec(g2.coVx,g2.coVy,g2.coVz));
		return result/atomCount;
	}
	public double VelMultiplyDNASideChain(atomGroup g1,atomGroup g2) {
		double result=0.0;
		
		  g1.setIdToAtom(); 
		  g2.setIdToAtom(); 
		  int atomCount=0;
		  for(int id:g1.idToAtom.keySet()) {
		  
			  Atom a1=g1.idToAtom.get(id);
			  Atom a2=g2.idToAtom.get(id);

			  
			  if(a1.type==6||a1.type==7||a1.type==8||a1.type==10||a1.type==11||a1.type==15) {//a1 is side chain
				  
				  
			  }else { continue; }
			 
			  atomCount++;
		  
			  result +=
			  (a1.vx*a2.vx+a1.vy*a2.vy+a1.vz*a2.vz)/(lengthOfVec(a1.vx,a1.vy,a1.vz)*
			  lengthOfVec(a2.vx,a2.vy,a2.vz));
		  
		  }
		 
		//result = (g1.coVx*g2.coVx+g1.coVy*g2.coVy+g1.coVz*g2.coVz)/(lengthOfVec(g1.coVx,g1.coVy,g1.coVz)*lengthOfVec(g2.coVx,g2.coVy,g2.coVz));
		return result/atomCount;
	}	
	public double VelMultiplyDNA(atomGroup g1,atomGroup g2) {
		double result=0.0;
		
		  g1.setIdToAtom(); 
		  g2.setIdToAtom(); 
		  int atomCount=0;
		  for(int id:g1.idToAtom.keySet()) {
		  
			  Atom a1=g1.idToAtom.get(id);
			  Atom a2=g2.idToAtom.get(id);
			  if(a1.type!=16) {
		  
			  }else { 
				  continue; 
			  
			  }
			  atomCount++;
		  
			  result +=
			  (a1.vx*a2.vx+a1.vy*a2.vy+a1.vz*a2.vz)/(lengthOfVec(a1.vx,a1.vy,a1.vz)*
			  lengthOfVec(a2.vx,a2.vy,a2.vz));
		  
		  }
		 
		//result = (g1.coVx*g2.coVx+g1.coVy*g2.coVy+g1.coVz*g2.coVz)/(lengthOfVec(g1.coVx,g1.coVy,g1.coVz)*lengthOfVec(g2.coVx,g2.coVy,g2.coVz));
		return result/atomCount;
	}
	public double lengthOfVec(double x,double y,double z) {
		return Math.sqrt(x*x+y*y+z*z);
	}
	
	public double calCor(HashMap<Integer,Double[]> timeDipole,ArrayList<Integer> timelist,int start,int end,int deltT) {
		//start is the start timestep of dipole to calculate correlation
		double sumMdipole=0.0;//sumation of dipole(t_1)*dipole(t_1+deltaT)
		int numdipoles=0;
		for(int i =start;i<=end;i++) {
			numdipoles++;
			Double [] dipoleS=timeDipole.get(timelist.get(start));
			Double [] dipoleE=timeDipole.get(timelist.get(start+deltT));
			double dipoleSM=Math.sqrt(Math.pow(dipoleS[0], 2)+Math.pow(dipoleS[1], 2)+Math.pow(dipoleS[2], 2));
			double dipoleEM=Math.sqrt(Math.pow(dipoleE[0], 2)+Math.pow(dipoleE[1], 2)+Math.pow(dipoleE[2], 2));

			sumMdipole+=(dipoleS[0]*dipoleE[0]+dipoleS[1]*dipoleE[1]+dipoleS[2]*dipoleE[2])/(dipoleSM*dipoleEM);
			
			
		}
		double averMdipole = sumMdipole/numdipoles;
		
		
		return averMdipole;
		
	}

	
}
