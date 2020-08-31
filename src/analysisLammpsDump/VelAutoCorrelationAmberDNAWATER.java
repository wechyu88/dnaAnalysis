package analysisLammpsDump;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class VelAutoCorrelationAmberDNAWATER extends AutoCorrelation {
	HashMap<Integer, HashMap<Integer, atomGroup>> timeDNADipole;
	HashMap<Integer, HashMap<Integer, water>> timeWaterDipole;
	HashMap<Integer, HashMap<Integer, water>> timeWaterVel;

	public static void main(String[] argv) throws IOException, InterruptedException {
		  String pathInput = argv[0]; int numberCores= Integer.valueOf(argv[1]); int
		  numberAverageStep = Integer.valueOf(argv[2]); 
		  int everyK = Integer.valueOf(argv[3]);
		  String pathOutput =
		  System.getProperty("user.dir"); double rad = Double.valueOf(argv[4]);		  
		  double averageRatio = Double.valueOf(argv[5]);
		 // String waterOut = argv[5];
		  String DNAOUT = argv[6];
		  int sideChainOrBack = Integer.valueOf(argv[7]);//0 side, 1 backbone, 2 all
		//  String pathwaterCor =
		//		  pathOutput+"/"+waterOut;
		  String pathDNACor =
				  pathOutput+"/"+DNAOUT;
		/*
		 * int numberCores=1; int numberAverageStep=1; String pathInput =
		 * "C:\\cygwin64\\home\\wechy\\DNA\\VelAutoCor\\dump.lammpstrj"; String
		 * pathOutput = "C:\\cygwin64\\home\\wechy\\DNA\\VelAutoCor\\"; String
		 * pathwaterCor = pathOutput+"velDNA.txt"; //String pathDNACor =
		 * //pathOutput+"velDNA.txt"; double rad = 3; double averageRatio =0.9;
		 * 
		 */
		 

		String path = pathInput;
		readdumpCustom rdc = new readdumpCustom(path);
		VelAutoCorrelationAmberDNAWATER acl = new VelAutoCorrelationAmberDNAWATER(numberAverageStep, rdc);
		acl.setEveryK(everyK);
		acl.init(rad);		
		HashMap<Integer, Double> waterVel = acl.calVelTimeCorDNA(sideChainOrBack,averageRatio, numberCores, acl.timelist,
				acl.timeDNADipole);

		PrintWriter pwwaterCor = new PrintWriter(new File(pathDNACor));
		for (Integer interval : waterVel.keySet()) {
			pwwaterCor.println(interval + " " + waterVel.get(interval));
		}
		pwwaterCor.close();
	}

	public VelAutoCorrelationAmberDNAWATER(int n, readdumpCustom rdc) throws IOException {
		super(n, rdc);
		this.timeDNADipole = new HashMap<Integer, HashMap<Integer, atomGroup>>();
		this.timeWaterVel = new HashMap<Integer, HashMap<Integer, water>>();

		// TODO Auto-generated constructor stub
		// this.init();
	}

	// @Override
	public void init(double rad) throws IOException {// rad is setting water molecules within radius rad from DNA
		dumpOneStep dos = readdumpCustom.readAverageNEveryKStep(this.rdc.br, this.averageN,this.everyK);
		ArrayList<Integer> waterID = new ArrayList<Integer>();
		this.timeWaterDipole = new HashMap<Integer, HashMap<Integer, water>>();
		int init = 0;
		
		while (dos.timestep != -1) {
			System.out.println(dos.xhi);
			if (init == 0) {

				HashMap<Integer, water> waterList = new HashMap<Integer, water>();

				HashMap<Integer, atomGroup> DNAlist = new HashMap<Integer, atomGroup>();
				for (Atom a : dos.atomlist) {
					if (a.type != 17 && a.type != 18 && a.type != 16) {
						if (DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						} else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup);
						}
					} else if (a.type == 17 || a.type == 18) {
						//continue;
						
						  if (waterList.containsKey(a.moelculeid)) {
						  waterList.get(a.moelculeid).addAtom(a); if (a.type == 17) {
						  waterList.get(a.moelculeid).setOxy(a);
						  
						  }
						  
						  } else { water newgroup = new water(); newgroup.addAtom(a); if (a.type == 17)
						  { newgroup.setOxy(a); }
						  
						  waterList.put(a.moelculeid, newgroup); }
						 
					}

				}
				

				 
				
				
			//	timelist.add(dos.timestep);
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);


				

				//System.out.println(dos.timestep+" "+newWaterList.size()+" "+DNAlist.size()+" "+rad+" "+dos.xlo+" "+dos.xhi+" "+dos.ylo+" "+dos.yhi+" "+dos.zlo+" "+dos.zhi);
				//this.timeWaterDipole.put(dos.timestep,newWaterList);
				
				
				//  for(int id:this.timeWaterDipole.get(dos.timestep).keySet()) {
				 // waterID.add(id); }

				
				for (Integer id : DNAlist.keySet()) {
					DNAlist.get(id).calCoV();
				}
				  for (Integer id : newWaterList.keySet()) { 
					  newWaterList.get(id).calCoV();
	
					  DNAlist.put(id, newWaterList.get(id));
				  }
				////////////////////

				timelist.add(dos.timestep);

				this.timeDNADipole.put(dos.timestep, DNAlist);
				init++;
			} else {
				HashMap<Integer, water> waterList = new HashMap<Integer, water>();

				HashMap<Integer, atomGroup> DNAlist = new HashMap<Integer, atomGroup>();
				for (Atom a : dos.atomlist) {
					if (a.type != 17 && a.type != 18 && a.type != 16) {
						if (DNAlist.containsKey(a.moelculeid)) {
							DNAlist.get(a.moelculeid).addAtom(a);
						} else {
							atomGroup newgroup = new atomGroup();
							newgroup.addAtom(a);
							DNAlist.put(a.moelculeid, newgroup);
						}
					} else if ((a.type == 17 || a.type == 18) ) {
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
						
						//continue;
					}
				}
				
				//timelist.add(dos.timestep);
				
				HashMap<Integer,water> newWaterList = waterWithinShell(dos.xlo,dos.xhi,dos.ylo,dos.yhi,dos.zlo,dos.zhi,waterList,DNAlist,rad);
				


				////////////////////
				
				for (Integer id : DNAlist.keySet()) {
					DNAlist.get(id).calCoV();
				}

				  for (Integer id : newWaterList.keySet()) { 
					  newWaterList.get(id).calCoV();
	
					  DNAlist.put(id, newWaterList.get(id));
				  }
				

				timelist.add(dos.timestep);
				this.timeDNADipole.put(dos.timestep, DNAlist);

			}
			dos = readdumpCustom.readAverageNEveryKStep(this.rdc.br, this.averageN,this.everyK);

		}

	}

	public HashMap<Integer, water> waterWithinShell(double xl, double xh, double yl, double yh, double zl, double zh,
			HashMap<Integer, water> waterMol, HashMap<Integer, atomGroup> DNAMol, double rad) {
		HashMap<Integer, water> newWaterList = new HashMap<Integer, water>();
		ArrayList<Double> dislist = new ArrayList<Double>();
		for (Integer gid : DNAMol.keySet()) {// DNAMol
			atomGroup twoDna = DNAMol.get(gid);

			for (Integer wid : waterMol.keySet()) {

				water oneWater = waterMol.get(wid);

				double dis = waterDNAgr.caldis(oneWater, twoDna, xl, xh, yl, yh, zl, zh);
				dislist.add(dis);
				if (dis < rad) {
					newWaterList.put(wid, oneWater);
				}
			}
		}
		return newWaterList;

	}

}
