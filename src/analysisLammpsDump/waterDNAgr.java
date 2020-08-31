package analysisLammpsDump;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

public class waterDNAgr extends waterOtherGr{
	HashMap<Integer,HashMap<Integer,atomGroup>> timeDNADipole;
	HashMap<Integer,water> waterList = new HashMap<Integer,water>();
	HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
	HashMap<Integer,HashMap<Integer,water>> timeWaterDipole = new HashMap<Integer,HashMap<Integer,water>>();

	HashMap<Integer,HashMap<Integer,atomGroup>> timeDNA = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	HashMap<Integer, Integer> idMajor = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> idMinor = new HashMap<Integer, Integer>();
	
	public static void main(String[] argv) throws IOException, InterruptedException {	
			
		  String pathInput = argv[0]; int numberCores= Integer.valueOf(argv[1]); int
		  numberAverageStep = Integer.valueOf(argv[2]); double binsize =
		  Double.valueOf(argv[3]); int numberAverage = Integer.valueOf(argv[4]); String
		  pathOutput = System.getProperty("user.dir"); String pathwaterGr =
		  pathOutput+"/autoWaterCor.txt"; String pathDNACor =
		  pathOutput+"/autoDNACor.txt"; int DNAATOMTYPE = Integer.valueOf(argv[5]);
		  String minorMajorPath = argv[6]; int majorOrMinor = Integer.valueOf(argv[7]);
		  int numSnaps = Integer.valueOf(argv[8]);
		
		
		
		/*
		 * int numberCores=1; int numberAverageStep=1; int numberAverage = 10; double
		 * binsize = 0.005;
		 * 
		 * String pathInput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\dump10\\T278\\dump.head";
		 * String pathOutput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\dump10\\T278\\"; String
		 * pathwaterGr = pathOutput+"autoWaterCor.txt";
		 * 
		 * int DNAATOMTYPE = 7; String minorMajorPath =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\amberDNAwater\\dump10\\T278\\minorOrMajor.txt";
		 * int majorOrMinor = 1;
		 */
		 
		//  String pathDNACor =  pathOutput+"autoglyCor.txt";
		 
		 
	    
		String path = pathInput;

		readdumpCustom rdc = new readdumpCustom(path);
		waterDNAgr wdgr = new waterDNAgr(rdc,numberAverage,binsize,DNAATOMTYPE,minorMajorPath,majorOrMinor);
		//wdgr.readMinorMajorAtom(minorMajorPath);

	//	waterDNAgr acl = new waterDNAgr(numberAverageStep,rdc);
	//	acl.init();
		PrintWriter pwWDgr = new PrintWriter(new File(pathwaterGr));
	//	PrintWriter pwglyCor = new PrintWriter(new File(pathDNACor));
		wdgr.CalNSnap(numSnaps);
		int counttime=0;
		double sysden=0;
		for(Integer t:wdgr.timeDen.keySet()) {
			counttime++;
			
			sysden+=wdgr.timeDen.get(t);
			if(counttime==10) {
				break;
			}
		}
		sysden/=counttime;
		double sumcount=0;
		for(Integer i =0;i<100*wdgr.DistD.size();i++) {
			if(wdgr.disCount.containsKey(i)) {
				double d = (i+0.5)*binsize;
				double count = wdgr.disCount.get(i);
				double vol=4/3*Math.PI*((i+0.5)*(i+0.5)*(i+0.5)-(i-0.5)*(i-0.5)*(i-0.5))*binsize*binsize*binsize;
				double dens = count/vol/sysden;
				sumcount+=count;
				pwWDgr.printf("%10.4f \t %10.4f \t %10.8f \t %10.4f \t%10.4f\n",d,count,dens,vol,sumcount);				
				
			}			
		}
		pwWDgr.close();
		



	}
	
	

	
//	@Override
//	public void init() throws IOException {
//		dumpOneStep dos= readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
//		while(dos.timestep!=-1) {
//			
//			System.out.println("Time: "+dos.timestep+"\n");
//		//	HashMap<Integer,atomGroup> waterList = new HashMap<Integer,atomGroup>();
//		//	HashMap<Integer,atomGroup> DNAlist = new HashMap<Integer,atomGroup>();
//		//	HashMap<Integer,atomGroup> saltlist = new HashMap<Integer,atomGroup>();
//			
//			//ArrayList<atomGroup> waterList = new atomGroup();
//			atomGroup dnagroup = new atomGroup();
//			DNAlist.put(1, dnagroup );	
//		
//			for(Atom a:dos.atomlist) {
//				if(a.type!=17&&a.type!=18&&a.type!=16) {
//					DNAlist.get(1).addAtom(a);				
//				}else if(a.type==17||a.type==18){
//					if(waterList.containsKey(a.moelculeid)) {
//						waterList.get(a.moelculeid).addAtom(a);
//						if(a.type==17) {
//							waterList.get(a.moelculeid).setOxy(a);
//						}
//					}else {
//						water newgroup = new water();
//						newgroup.addAtom(a);
//						if(a.type==17) {
//							newgroup.setOxy(a);
//						}
//						waterList.put(a.moelculeid, newgroup );
//					}					
//				}
//					
//				
//				
//				
//			}
//			for(Integer id:waterList.keySet()) {
//				waterList.get(id).caldipole();
//			}
//			for(Integer id:DNAlist.keySet()) {
//				DNAlist.get(id).caldipole();
//			}
//			timelist.add(dos.timestep);
//			System.out.println("time is "+dos.timestep);
//			timeWaterDipole.put(dos.timestep,waterList);
//			for(int a:DNAlist.keySet()) {
//				System.out.println(" dna "+a+" "+DNAlist.get(a).dipolex+" "+DNAlist.get(a).dipoley+" "+DNAlist.get(a).dipolez);
//			}
//			this.timeDNADipole.put(dos.timestep,DNAlist);
//
//			dos = readdumpCustom.readAverageNStep(this.rdc.br,this.averageN);
//	
//		}
//
//		
//	}
//	
	
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
	public waterDNAgr(readdumpCustom rdc,int nstep,double binsize,int DNATypeAtom,String minorMajorPath,int majorOrMinor) throws IOException{
		dumpOneStep dos= rdc.readnextNstep(nstep);
		readMinorMajorAtom(minorMajorPath);
		
		while(dos.timestep!=-1) {
			waterList.clear();
			DNAlist.clear();
			System.out.println("Time: "+dos.timestep+"\n");

			
			//ArrayList<atomGroup> waterList = new atomGroup();
			
			for(Atom a:dos.atomlist) {
				if(a.type==DNATypeAtom) {
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
				waterList.get(id).calcom(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			for(Integer id:DNAlist.keySet()) {
				DNAlist.get(id).calcom(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
			}
			timelist.add(dos.timestep);
			System.out.println("time is "+dos.timestep);
			timeWater.put(dos.timestep,waterList);
			timeDNA.put(dos.timestep,DNAlist);
			double den = waterList.size()/(dos.xhi-dos.xlo)/(dos.yhi-dos.ylo)/(dos.zhi-dos.zlo);
			System.out.println("den is  "+den+" count is "+waterList.size()+" v is "+(dos.xhi-dos.xlo)*(dos.yhi-dos.ylo)*(dos.zhi-dos.zlo)+" xhi "+dos.zhi+" xlo "+dos.zlo);
			timeDen.put(dos.timestep,den);
			gr tempgr = CalOneSnap4(binsize,dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,majorOrMinor);
			grTime.add(tempgr);
			dos = rdc.readnextNstep(nstep);
	
		}
	}
	public gr CalOneSnap2(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> DNAMol) {
		gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
		for(Integer gid:DNAMol.keySet()) {//DNAMol
		//	System.out.println("gly size is "+glyMol.size());
			atomGroup twoDna = DNAMol.get(gid);

			oneTimeGr.addAverNum();
			for(Integer wid:waterMol.keySet()) {
				
				water oneWater=waterMol.get(wid);
				

				double dis = caldis(oneWater,twoDna,xl,xh,yl,yh,zl,zh);
				if(dis!=0) {
					oneTimeGr.addDis(dis);									
				}
			}
		}
		oneTimeGr.CalculateDist();
		return oneTimeGr;
		
	}
	
	public gr CalOneSnap3(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> DNAMol) {
		gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
		int dnaatom = 0;
		for(Integer gid:DNAMol.keySet()) {//DNAMol
		//	System.out.println("gly size is "+glyMol.size());
			atomGroup twoDna = DNAMol.get(gid);
			dnaatom = 0;
			oneTimeGr.addAverNum();
			for(Atom a:twoDna.atomlist) {
				//if(a.mass<13||a.type==16||a.type==14) {
				//	continue;
				//}

				
				dnaatom++;
				for(Integer wid:waterMol.keySet()) {
					
					water oneWater=waterMol.get(wid);
					

					double dis = caldis(oneWater,a,xl,xh,yl,yh,zl,zh);
					if(dis!=0) {
						oneTimeGr.addDis(dis);									
					}
				}				
			}

		}
		oneTimeGr.Calgr(dnaatom);
		return oneTimeGr;
		
	}// without major and minor
	
	
	////
	//// CalOneSnap4 include major and minor distinction
	public gr CalOneSnap4(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> DNAMol,int majorOfMinor) {// major 1,minor 0
		HashMap<Integer, Integer> idMajorOrMinor = new HashMap<Integer,Integer> ();
		if(majorOfMinor==1) {
			idMajorOrMinor = idMajor;
		}else {
			idMajorOrMinor = idMinor;
			
		}
		gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
		int dnaatom = 0;
		for(Integer gid:DNAMol.keySet()) {//DNAMol
		//	System.out.println("gly size is "+glyMol.size());
			atomGroup twoDna = DNAMol.get(gid);
			dnaatom = 0;
			oneTimeGr.addAverNum();
			for(Atom a:twoDna.atomlist) {
				//if(a.mass<13||a.type==16||a.type==14) {
				//	continue;
				//}
				if(idMajorOrMinor.containsKey(a.id)) {
					
				}else {
					continue;
				}
				
				dnaatom++;
				for(Integer wid:waterMol.keySet()) {
					
					water oneWater=waterMol.get(wid);
					

					double dis = caldis(oneWater,a,xl,xh,yl,yh,zl,zh);
					if(dis!=0) {
						oneTimeGr.addDis(dis);									
					}
				}				
			}
		}
		oneTimeGr.Calgr(dnaatom);
		return oneTimeGr;
		
	}
	

	//@Override
	public static double caldis(water a,atomGroup b,double xl,double xh,double yl,double yh,double zl,double zh) {
		double dis=0.0;
		Atom mindatom = new Atom ();
		double mindis = 0.5*Integer.MAX_VALUE;
		int count=0;
		for(Atom batm:b.atomlist) {
			if(batm.mass<13&&batm.type<14) {
				continue;
			}
			count++;
			double dx = a.Oxygen.x-batm.x;
			double dy = a.Oxygen.y-batm.y;
			double dz = a.Oxygen.z-batm.z;
			dx = dx - (xh-xl)*Math.round(dx/(xh-xl));
			dy = dy - (yh-yl)*Math.round(dy/(yh-yl));
			dz = dz - (zh-zl)*Math.round(dz/(zh-zl));
			dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	
			if(dis<mindis&&dis>0) {
				mindis = dis;
			}
			if(dis<0.4) {
				System.out.println("find wierd");
			}
		}

		return mindis;
	}
	public double caldis(water a,Atom b,double xl,double xh,double yl,double yh,double zl,double zh) {
		double dis=0.0;
		Atom mindatom = new Atom ();
		double mindis = 1.0*Integer.MAX_VALUE;
//		for(Atom batm:b.atomlist) {

		double dx = a.Oxygen.x-b.x;
		double dy = a.Oxygen.y-b.y;
		double dz = a.Oxygen.z-b.z;
		dx = dx - (xh-xl)*Math.round(dx/(xh-xl));
		dy = dy - (yh-yl)*Math.round(dy/(yh-yl));
		dz = dz - (zh-zl)*Math.round(dz/(zh-zl));
		dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	


//		}

		return dis;
	}
	

}
