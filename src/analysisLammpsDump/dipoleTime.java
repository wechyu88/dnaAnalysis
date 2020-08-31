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



public class dipoleTime {
	ArrayList<Integer> timelist = new ArrayList<Integer>();
	HashMap<Integer,Double[]> timeDipole = new HashMap<Integer,Double[]>();
	HashMap<Integer,HashMap<Integer,atomGroup>> timeWaterDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
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
				Double cor = calCor(timelist, timeDipole,start,end, deltaT);
				return cor;
			}
			
		}	
	public static void main(String[] argv) throws IOException, InterruptedException {	
		
		String pathInput = argv[0];
		int numberCores= Integer.valueOf(argv[1]);// set to 1
		int numberAverageStep = Integer.valueOf(argv[2]);//set to 1
		int moleculenum = Integer.valueOf(argv[3]);
		String pathOutput = System.getProperty("user.dir");
		String pathwaterTheta = pathOutput+"/autoWaterAngle.txt";
		String pathglyTheta = pathOutput+"/autoglyAngle.txt";
		

	    //int numberStepBin=2;
		// int numberStepBin = Integer.valueOf(argv[2]);
/*		int numberCores=1;
		int numberAverageStep=1;
		int moleculenum = 2;

	    String pathInput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\dump.lammpstrj2";
	    String pathOutput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\";
		String pathwaterTheta = pathOutput+"autoWaterAngle.txt";
		String pathglyTheta = pathOutput+"autoglyAngle.txt";*/
    
		String path = pathInput;
		readdumpCustom rdc = new readdumpCustom(path);
		//AutoCorrelation acl = new AutoCorrelation(rdc,numberStepBin);
		dipoleTime dipoleTime = new dipoleTime(numberAverageStep,rdc);

		//HashMap<Integer,Double> waterDipole =  acl.calTimeCor(numberCores,acl.timelist, acl.timeWaterDipole);
		//HashMap<Integer,Double> glyDipole = 
		ArrayList<dipole> dipolelistWater =	dipoleTime.calglyDipoleTime(dipoleTime.timeWaterDipole,moleculenum);
		ArrayList<dipole> dipolelistGly =	dipoleTime.calglyDipoleTime(dipoleTime.timeglyDipole,moleculenum);
		
		PrintWriter pwglyTheta = new PrintWriter(new File(pathglyTheta));
		PrintWriter pwWaterTheta = new PrintWriter(new File(pathwaterTheta));

		
		

		dipole initdipoleWater = dipolelistWater.get(0);
		dipole initdipolegly = dipolelistGly.get(0);
		double lengthdipoleinitW = Math.sqrt(initdipoleWater.dipoleX*initdipoleWater.dipoleX+initdipoleWater.dipoleY*initdipoleWater.dipoleY+initdipoleWater.dipoleZ*initdipoleWater.dipoleZ);
		double lengthdipoleinitG = Math.sqrt(initdipolegly.dipoleX*initdipolegly.dipoleX+initdipolegly.dipoleY*initdipolegly.dipoleY+initdipolegly.dipoleZ*initdipolegly.dipoleZ);
		int timecount=0;
		for(dipole dipolew: dipolelistWater) {
			//pwwaterCor.println(interval+" "+waterDipole.get(interval));
			double lengthdipole = Math.sqrt(dipolew.dipoleX*dipolew.dipoleX+dipolew.dipoleY*dipolew.dipoleY+dipolew.dipoleZ*dipolew.dipoleZ);
			double cosine = (dipolew.dipoleX*initdipoleWater.dipoleX+dipolew.dipoleY*initdipoleWater.dipoleY+dipolew.dipoleZ*initdipoleWater.dipoleZ)/lengthdipoleinitW/lengthdipole;
			pwWaterTheta.println(timecount+" "+dipolew.theta+" "+dipolew.phi);
			timecount++;
			//pwglyCor.println(interval+" "+glyDipole.get(interval));
		}
		timecount=0;
		for(dipole dipoleg: dipolelistGly) {
			//pwwaterCor.println(interval+" "+waterDipole.get(interval));
			double lengthdipole = Math.sqrt(dipoleg.dipoleX*dipoleg.dipoleX+dipoleg.dipoleY*dipoleg.dipoleY+dipoleg.dipoleZ*dipoleg.dipoleZ);
			double cosine = (dipoleg.dipoleX*initdipolegly.dipoleX+dipoleg.dipoleY*initdipolegly.dipoleY+dipoleg.dipoleZ*initdipolegly.dipoleZ)/lengthdipoleinitG/lengthdipole;
			pwglyTheta.println(timecount+" "+dipoleg.theta+" "+dipoleg.phi);
			timecount++;
			//pwglyCor.println(interval+" "+glyDipole.get(interval));
		}
		pwWaterTheta.close();
		pwglyTheta.close();
		
	}

	public dipoleTime(int n,readdumpCustom rdc) throws IOException {
	       // BufferedWriter pw1 = new BufferedWriter(new FileWriter(path5+"/dataStockMayer.txt")); 
			dumpOneStep dos= readdumpCustom.readAverageNStep(rdc.br,n);//n=1
			while(dos.timestep!=-1) {
				
				System.out.println("Time: "+dos.timestep+"\n");
				HashMap<Integer,atomGroup> waterList = new HashMap<Integer,atomGroup>();
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
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							waterList.put(a.moelculeid, newgroup );
						}					
					}
					
					
				}
				for(Integer id:waterList.keySet()) {
					waterList.get(id).caldipoleV();;
				}
				for(Integer id:glyList.keySet()) {
					glyList.get(id).caldipoleV();
				}
				timelist.add(dos.timestep);
				System.out.println("time is "+dos.timestep);
				timeWaterDipole.put(dos.timestep,waterList);
				timeglyDipole.put(dos.timestep,glyList);

				dos = readdumpCustom.readAverageNStep(rdc.br,n);
		
			}
		}
	public dipoleTime(String path) throws IOException {
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
	public ArrayList<dipole> calglyAVDipoleTime(HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole) {
		//HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
		//return two things, 
		//return vector length of the dipole vs time
		//return vector(t).*vector(0) the projection on the init dipole vs time
		// just return dipole list vs time can get this two things we need.
		ArrayList<dipole> dipolelist = new ArrayList<dipole>();
		for(Integer timelist:timeglyDipole.keySet()) {
			HashMap<Integer, atomGroup> glylist =timeglyDipole.get(timelist);
			double dx=0;
			double dy=0;
			double dz=0;
			double avlength=0;
			for(Integer mid:glylist.keySet()) {
				double lengthdipole = Math.sqrt(Math.pow(glylist.get(mid).dipoleVx, 2)+Math.pow(glylist.get(mid).dipoleVy, 2)+Math.pow(glylist.get(mid).dipoleVz, 2));
				avlength+=lengthdipole;
				dx+=glylist.get(mid).dipoleVx;
				dy+=glylist.get(mid).dipoleVy;
				dz+=glylist.get(mid).dipoleVz;
			}
			avlength/=glylist.size();
			dx/=glylist.size();
			dy/=glylist.size();
			dz/=glylist.size();
			dipole avdipole = new dipole(dx,dy,dz);
			dipolelist.add(avdipole);
		}

		return dipolelist;
		
	}
	public ArrayList<dipole> calglyDipoleTime(HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole,int n) {//n is the n_th moleclule to extract 
		//HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
		//return two things, 
		//return vector length of the dipole vs time
		//return vector(t).*vector(0) the projection on the init dipole vs time
		// just return dipole list vs time can get this two things we need.
		ArrayList<dipole> dipolelist = new ArrayList<dipole>();
		int molid=0;
		int num=0;
		for(Integer timelist:timeglyDipole.keySet()) {
			HashMap<Integer, atomGroup> glylist =timeglyDipole.get(timelist);
			for(Integer mid:glylist.keySet()) {
				num++;
				if(num==n) {
					molid=mid;
					break;					
				}else {
					continue;
				}

			}
			break;
		}

		
		for(Integer timelist:timeglyDipole.keySet()) {
			
			HashMap<Integer, atomGroup> glylist =timeglyDipole.get(timelist);
			atomGroup glygroup=glylist.get(molid);
			dipole newdipole = new dipole(glygroup.dipoleVx,glygroup.dipoleVy,glygroup.dipoleVz);
			dipolelist.add(newdipole);
		}
		return dipolelist;
		
	}
	public HashMap<Integer,Double> calTimeCor(int nprocForLoop,ArrayList<Integer> timeStep,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole) throws InterruptedException {
		int length = timeStep.size();
		int start = 0;
		int end = start+(int)(0.5*length)-1;
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
	public double calCor(ArrayList<Integer> timelist,HashMap<Integer,HashMap<Integer,atomGroup>> timeDipole,int start,int end,int deltaT) {
		//double cor=0.0;
		double sumCor=0.0;
		int totalnum=0;
		for(int i=start;i<=end;i++) {
			HashMap<Integer,Double> corlist =  calcorDeltT(timeDipole.get(timelist.get(i)),timeDipole.get(timelist.get(i+deltaT)));
			for(Integer molid:corlist.keySet()) {
				sumCor+=corlist.get(molid);
				totalnum++;
			}
		}
		return sumCor/totalnum;
	}
	
	public HashMap<Integer,Double> calcorDeltT(HashMap<Integer,atomGroup> time1group,HashMap<Integer,atomGroup> time2group){
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
