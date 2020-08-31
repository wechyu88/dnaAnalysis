package analysisLammpsDump;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class DisDistr {
	ArrayList<Integer> timelist = new ArrayList<Integer>();
	HashMap<Integer,HashMap<Integer,atomGroup>> timeWater = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	HashMap<Integer,HashMap<Integer,atomGroup>> timegly = new HashMap<Integer,HashMap<Integer,atomGroup>>();
	ArrayList<gr> grTime = new ArrayList<gr>();
	ArrayList<Double> DistD= new ArrayList<Double>();// dis for all the bins 
	HashMap<Integer,Double> disCount= new HashMap<Integer,Double>();
	public static void main(String[] argv) throws IOException, InterruptedException {			
		String pathInput = argv[0];
		int numberStepBin = Integer.valueOf(argv[1]);//dumpfile every numberStepbin timestep read once
		int numberAverage = Integer.valueOf(argv[2]);
		double binsize = Double.valueOf(argv[3]);

	    String pathOutput = System.getProperty("user.dir");
		String pathDisD = pathOutput+"/DisD.txt";

		//String path = pathInput;

/*		int numberStepBin=2;
		int numberAverage = 4;
		double binsize = 0.2;
	    String pathInput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\dump.lammpstrj2";
	    String pathOutput = "C:\\cygwin64\\home\\wechy\\vinh\\0.0415\\";
		String pathDisD = pathOutput+"DisD2.txt";*/
		PrintWriter pwDisD = new PrintWriter(new File(pathDisD));
		readdumpCustom rdc = new readdumpCustom(pathInput);
		DisDistr DisD = new DisDistr(rdc,numberStepBin,binsize);
		DisD.CalNSnap(numberAverage);
		for(Integer i =0;i<2*DisD.disCount.keySet().size();i++) {
			if(DisD.disCount.containsKey(i)) {
				double d = (i+0.5)*binsize;
				double count = DisD.disCount.get(i);
				double dens = count/(d*d*binsize);
				pwDisD.printf("%10.4f \t %10.4f \t %10.8f\n",d,count,dens);				
				
				
				
			}			
		}
/*		for(int i =0;i<DisD.DistD.size();i++) {
			if(DisD.disCount.containsKey(DisD.DistD.get(i))) {
				double d = DisD.DistD.get(i);
				double count = DisD.disCount.get(d);
				double dens = count/(4*Math.PI*d*d*(d-DisD.DistD.get(i-1)));
				pwDisD.println(DisD.DistD.get(i)+"\t"+DisD.disCount.get(DisD.DistD.get(i))+"\t"+dens);				
			}
		}*/

		pwDisD.close();
		
		
	}
	
	public DisDistr(readdumpCustom rdc,int nstep,double binsize) throws IOException {
	       // BufferedWriter pw1 = new BufferedWriter(new FileWriter(path5+"/dataStockMayer.txt")); 
			dumpOneStep dos= rdc.readnextNstep(nstep);
		
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
					waterList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				for(Integer id:glyList.keySet()) {
					glyList.get(id).caldipole(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi);
				}
				timelist.add(dos.timestep);
				System.out.println("time is "+dos.timestep);
				timeWater.put(dos.timestep,waterList);
				timegly.put(dos.timestep,glyList);
				grTime.add(CalOneSnap(binsize,dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,glyList));
				dos = rdc.readnextNstep(nstep);
		
			}
		}
		public void CalNSnap(int n) {
			//ArrayList<Double> DistD= new ArrayList<Double>();// dis for all the bins 
			HashMap<Integer,Double> disCountn= new HashMap<Integer,Double>();
			DistD=grTime.get(0).DistD;
			
			for(int i=0;i<n;i++) {
				System.out.println(" i is "+i);
				for(Integer dis:grTime.get(i).disCount.keySet()) {
					if(disCountn.containsKey(dis)) {
						//System.out.println(" i is " +i);
						disCountn.put(dis, disCountn.get(dis)+grTime.get(i).disCount.get(dis));
					}else {
						disCountn.put(dis, grTime.get(i).disCount.get(dis));
					}
				}

			}
			for(Integer dis :disCountn.keySet()) {
				disCountn.put(dis, disCountn.get(dis)/(n*1.0));
			}
			this.disCount=disCountn;

		}
		public gr CalOneSnap2(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,atomGroup> waterMol,HashMap<Integer,atomGroup> glyMol) {
			gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
			oneTimeGr.addAverNum();
			for(Integer wid:waterMol.keySet()) {
				atomGroup onewater=waterMol.get(wid);
				//oneTimeGr.addAverNum();
				double dis = 100000.00;
				for(Integer gid:glyMol.keySet()) {
					atomGroup onegly = glyMol.get(gid);
					for(Atom watom:onewater.atomlist) {
						if(watom.type==70) {//64,12 are carbon, 67 are oxygen of gly
							//oneTimeGr.addDis(oglydis(onegly,watom,xl,xh,yl,yh,zl,zh));
							double tempdis = oglydis(onegly,watom,xl,xh,yl,yh,zl,zh);
							if(tempdis<dis) {
								dis = tempdis;
							}
						}
					}
					
				}
				oneTimeGr.addDis(dis);
				
			}
			oneTimeGr.CalculateDist();
			return oneTimeGr;
			
		}
		public gr CalOneSnap(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,atomGroup> waterMol,HashMap<Integer,atomGroup> glyMol) {
			gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
			for(Integer gid:glyMol.keySet()) {
				atomGroup onegly=glyMol.get(gid);
				oneTimeGr.addAverNum();
				
				for(Integer wid:waterMol.keySet()) {
					atomGroup onewater = waterMol.get(wid);
					for(Atom watom:onewater.atomlist) {
						if(watom.type==70) {//64,12 are carbon, 67 are oxygen of gly
							oneTimeGr.addDis(oglydis(onegly,watom,xl,xh,yl,yh,zl,zh));
						}
					}
					
				}
			}
			oneTimeGr.CalculateDist();
			return oneTimeGr;
			
		}
		
		public double oglydis(atomGroup gly,Atom oxy,double xl,double xh,double yl,double yh,double zl,double zh) {
			double dis=xh-xl;
			for(Atom gl:gly.atomlist) {
				if(gl.type==64||gl.type==12||gl.type==67) {
					double tempdis=caldis(gl,oxy,xl,xh,yl,yh,zl,zh);
					if(dis>tempdis) {
						dis=tempdis;
					}
				}
			}
			
			return dis;
		}
		public double oglydis2(atomGroup gly,Atom oxy,double xl,double xh,double yl,double yh,double zl,double zh) {
			double dis=xh-xl;
			for(Atom gl:gly.atomlist) {
				//if(gl.type==64||gl.type==12||gl.type==67) {
					double tempdis=caldis(gl,oxy,xl,xh,yl,yh,zl,zh);
					if(dis>tempdis) {
						dis=tempdis;
					}
				//}
			}
			
			return dis;
		}
		public double caldis(Atom a,Atom b,double xl,double xh,double yl,double yh,double zl,double zh) {
			double dis=0.0;
			double dx = a.x-b.x;
			double dy = a.y-b.y;
			double dz = a.z-b.z;
			dx = dx - (xh-xl)*Math.round(dx/(xh-xl));
			dy = dy - (yh-yl)*Math.round(dy/(yh-yl));
			dz = dz - (zh-zl)*Math.round(dz/(zh-zl));
			dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	
			if(dis>0.9*(xh-xl)) {
				System.out.println("this dis not correct");
				System.out.println("dx is "+dx+"ax is "+a.x+" bx is "+b.x+" xh is "+xh+" xl is "+xl);
				System.out.println("dy is "+dy+"ay is "+a.y+" by is "+b.y+" yh is "+yh+" yl is "+yl);
				System.out.println("dz is "+dz+"az is "+a.z+" bz is "+b.z+" zh is "+zh+" zl is "+zl);

			}
			return dis;
		}
}
