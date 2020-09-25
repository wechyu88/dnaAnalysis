package analysisLammpsDump;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class readdumpCustom {
	BufferedReader br;
	public static void main(String[] argv) throws IOException {
		String pathInput = argv[0];
		
	    String pathOutput = System.getProperty("user.dir");
       // BufferedWriter pw1 = new BufferedWriter(new FileWriter(path5+"/dataStockMayer.txt")); 
		String path = pathInput;
		readdumpCustom rdc = new readdumpCustom(path);
		dumpOneStep dos= readnextStep(rdc.br);
		HashMap<Integer,HashMap<Integer,atomGroup>> timeWaterDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
		HashMap<Integer,HashMap<Integer,atomGroup>> timeglyDipole = new HashMap<Integer,HashMap<Integer,atomGroup>>();
		
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
			timeWaterDipole.put(dos.timestep,waterList);
			timeglyDipole.put(dos.timestep,glyList);

			dos = readnextStep(rdc.br);
	
		}

	}

	
	
	public readdumpCustom(String path) throws FileNotFoundException {
		File file = new File(path);
		FileReader fr = new FileReader(file);
		this.br = new BufferedReader(fr);
	}
	public dumpOneStep readnextNstep(int n) throws IOException {
		dumpOneStep oneStep = new dumpOneStep();
		String line="";
		for(int j =0;j<n-1;j++) {
			line=br.readLine();

			int timestep=-1;
			int numAtom=0;
			if(line!=null&&line.matches("ITEM: TIMESTEP")) {
				line=br.readLine();
			}else {
				oneStep.setTime(timestep);
			//	System.out.println("first line is not timestep");
				return oneStep;
			}
			line=br.readLine();
			if(line.matches("ITEM: NUMBER OF ATOMS")) {
				line=br.readLine();
				numAtom = Integer.valueOf(line);
			}else {
				System.out.println("second line is not number of atom");
			}	
			line= br.readLine();
			if(line.matches("ITEM: BOX.*")) {
				line=br.readLine();
				line=br.readLine();
				line=br.readLine();
				//numAtom = Integer.valueOf(line);
			}else {
				System.out.println("third line is not box");
			}
			line=br.readLine();
			if(line.matches("ITEM: ATOMS.*")) {
				for(int i =1;i<=numAtom;i++) {
					line=br.readLine();
				}
				//System.out.println(line);
			}
			//return oneStep;
		}
		oneStep = new dumpOneStep();
		int timestep=-1;
		int numAtom=0;
		line=br.readLine();
		if(line!=null&&line.matches("ITEM: TIMESTEP")) {
			line=br.readLine();
			timestep = Integer.valueOf(line);
			oneStep.setTime(timestep);
		}else {
			oneStep.setTime(timestep);
		//	System.out.println("first line is not timestep");
			return oneStep;
		}
		line=br.readLine();
		if(line.matches("ITEM: NUMBER OF ATOMS")) {
			line=br.readLine();
			numAtom = Integer.valueOf(line);
			oneStep.setnumAtom(numAtom);
		}else {
			System.out.println("second line is not number of atom");
		}	
		line= br.readLine();
		if(line.matches("ITEM: BOX.*")) {
			line=br.readLine();
			String [] x = line.split(" ");
			float xlow = Float.valueOf(x[0]);
			float xhigh=Float.valueOf(x[1]);
			line=br.readLine();
			String [] y = line.split(" ");
			float ylow = Float.valueOf(y[0]);
			float yhigh=Float.valueOf(y[1]);
			line=br.readLine();
			String [] z = line.split(" ");
			float zlow = Float.valueOf(z[0]);
			float zhigh=Float.valueOf(z[1]);
			oneStep.setboundary(xlow, xhigh, ylow, yhigh, zlow, zhigh);
			//numAtom = Integer.valueOf(line);
		}else {
			System.out.println("third line is not box");
		}
		ArrayList<Atom>atmlist = new ArrayList<Atom>();
		line=br.readLine();
		if(line.matches("ITEM: ATOMS.*")) {
			for(int i =1;i<=numAtom;i++) {
				line=br.readLine();
				String [] fields = line.split(" ");
				Atom atm = new Atom(Integer.valueOf(fields[0]));
				atm.setCharge(Double.valueOf(fields[4]));
				atm.setMass(Double.valueOf(fields[3]));
				atm.setMoleculeId(Integer.valueOf(fields[1]));
				atm.setType(Integer.valueOf(fields[2]));
				atm.setXYZ(Double.valueOf(fields[5]), Double.valueOf(fields[6]), Double.valueOf(fields[7]));				
				atmlist.add(atm);
			}
		}
		oneStep.setAtomlist(atmlist);
		return oneStep;
	}
	public void close() throws IOException {
		this.br.close();
		
	}
	public static dumpOneStep readnextStep(BufferedReader br) throws IOException {
		String line=br.readLine();
		dumpOneStep oneStep = new dumpOneStep();
		int timestep=-1;
		int numAtom=0;
		if(line!=null&&line.matches("ITEM: TIMESTEP")) {
			line=br.readLine();
			timestep = Integer.valueOf(line);
			oneStep.setTime(timestep);
		}else {
			oneStep.setTime(timestep);
		//	System.out.println("first line is not timestep");
			return oneStep;
		}
		line=br.readLine();
		if(line.matches("ITEM: NUMBER OF ATOMS")) {
			line=br.readLine();
			numAtom = Integer.valueOf(line);
			oneStep.setnumAtom(numAtom);
		}else {
			System.out.println("second line is not number of atom");
		}	
		line= br.readLine();
		if(line.matches("ITEM: BOX.*")) {
			line=br.readLine();
			String [] x = line.split(" ");
			float xlow = Float.valueOf(x[0]);
			float xhigh=Float.valueOf(x[1]);
			line=br.readLine();
			String [] y = line.split(" ");
			float ylow = Float.valueOf(y[0]);
			float yhigh=Float.valueOf(y[1]);
			line=br.readLine();
			String [] z = line.split(" ");
			float zlow = Float.valueOf(z[0]);
			float zhigh=Float.valueOf(z[1]);
			oneStep.setboundary(xlow, xhigh, ylow, yhigh, zlow, zhigh);
			//numAtom = Integer.valueOf(line);
		}else {
			System.out.println("third line is not box");
		}
		ArrayList<Atom>atmlist = new ArrayList<Atom>();
		line=br.readLine();
		System.out.println(oneStep.timestep+" time ");
		if(line.matches("ITEM: ATOMS.*")) {
			for(int i =1;i<=numAtom;i++) {
				line=br.readLine();
				//System.out.println(line);

				String [] fields = line.split(" ");
				Atom atm = new Atom(Integer.valueOf(fields[0]));
				atm.setCharge(Double.valueOf(fields[4]));
				atm.setMass(Double.valueOf(fields[3]));
				atm.setMoleculeId(Integer.valueOf(fields[1]));
				atm.setType(Integer.valueOf(fields[2]));
				atm.setXYZ(Double.valueOf(fields[5]), Double.valueOf(fields[6]), Double.valueOf(fields[7]));	
				if(fields.length<=8) {
					atm.setV(0.0, 0.0,0.0);
					
				}else {
					atm.setV(Double.valueOf(fields[8]), Double.valueOf(fields[9]), Double.valueOf(fields[10]));
				}
				atmlist.add(atm);
			}
		}
		oneStep.setAtomlist(atmlist);
		return oneStep;
	}
	public static dumpOneStep readAverageNStep(BufferedReader br,int n) throws IOException {
		ArrayList<dumpOneStep> oneStepList = new ArrayList<dumpOneStep>();
		for(int i=0;i<n;i++) {
			dumpOneStep tempOneStep =  readnextStep(br);
		//	System.out.println(tempOneStep.timestep+" readaverage");

			oneStepList.add(tempOneStep);
			if(tempOneStep.timestep==-1) {
				return tempOneStep;
			}
		}
		
		
		dumpOneStep oneStep = new dumpOneStep();
		oneStep.setTime(oneStepList.get(oneStepList.size()-1).timestep);
		ArrayList<Atom> averageList = new ArrayList<Atom>();
		oneStep.setboundary(oneStepList.get(0).xlo,oneStepList.get(0).xhi,oneStepList.get(0).ylo,oneStepList.get(0).yhi,oneStepList.get(0).zlo,oneStepList.get(0).zhi);
		for(Atom a:oneStepList.get(0).atomlist) {
			Atom  tempatom = new Atom(a.id);
			tempatom.setCharge(a.charge);
			tempatom.setMass(a.mass);
			tempatom.setMoleculeId(a.moelculeid);
			tempatom.setType(a.type);
			double x=0.0,y=0.0,z=0.0;
			double vx=0.0,vy=0.0,vz=0.0;
			for(dumpOneStep os : oneStepList) {
				Atom atomb = os.atomMap.get(a.id);
				x+=atomb.x;
				vx+=atomb.vx;
				y+=atomb.y;
				vy+=atomb.vy;
				z+=atomb.z;
				vz+=atomb.vz;
			}
			x/=oneStepList.size();
			y/=oneStepList.size();
			z/=oneStepList.size();
			vx/=oneStepList.size();
			vy/=oneStepList.size();
			vz/=oneStepList.size();
			tempatom.setXYZ(x, y, z);		
			tempatom.setV(a.vx, a.vy, a.vz);

			averageList.add(tempatom);
		}

		
		
		oneStep.setAtomlist(averageList);
		
		return oneStep;
	}
	
	public static dumpOneStep readAverageNEveryKStep(BufferedReader br,int n, int k) throws IOException {
		// average n every k step such as , 1 2 3 4 5 6, average 2 every 2 step means average 1 and 2, then average 5 and 6 then average 9 and 10
		ArrayList<dumpOneStep> oneStepList = new ArrayList<dumpOneStep>();
		for(int i=0;i<n;i++) {
			dumpOneStep tempOneStep =  readnextStep(br);
			System.out.println(tempOneStep.timestep+" readaverage");

			oneStepList.add(tempOneStep);
			if(tempOneStep.timestep==-1) {
				return tempOneStep;
			}
		}

		
		
		dumpOneStep oneStep = new dumpOneStep();
		oneStep.setTime(oneStepList.get(oneStepList.size()-1).timestep);
		ArrayList<Atom> averageList = new ArrayList<Atom>();
		oneStep.setboundary(oneStepList.get(0).xlo,oneStepList.get(0).xhi,oneStepList.get(0).ylo,oneStepList.get(0).yhi,oneStepList.get(0).zlo,oneStepList.get(0).zhi);
		for(Atom a:oneStepList.get(0).atomlist) {
			Atom  tempatom = new Atom(a.id);
			tempatom.setCharge(a.charge);
			tempatom.setMass(a.mass);
			tempatom.setMoleculeId(a.moelculeid);
			tempatom.setType(a.type);
			float x=(float) 0.0,y=(float) 0.0,z=(float) 0.0;
			float vx=(float) 0.0,vy=(float) 0.0,vz=(float) 0.0;
			for(dumpOneStep os : oneStepList) {
				Atom atomb = os.atomMap.get(a.id);
				x+=atomb.x;
				vx+=atomb.vx;
				y+=atomb.y;
				vy+=atomb.vy;
				z+=atomb.z;
				vz+=atomb.vz;
			}
			x/=oneStepList.size();
			y/=oneStepList.size();
			z/=oneStepList.size();
			vx/=oneStepList.size();
			vy/=oneStepList.size();
			vz/=oneStepList.size();
			tempatom.setXYZ(x, y, z);		
			tempatom.setV(a.vx, a.vy, a.vz);

			averageList.add(tempatom);
		}

		
		
		oneStep.setAtomlist(averageList);
		for(int i=0;i<k;i++) {
			dumpOneStep tempOneStep =  readnextStep(br);

			if(tempOneStep.timestep==-1) {
				return oneStep;
			}
		}
		return oneStep;
	}
	
}
