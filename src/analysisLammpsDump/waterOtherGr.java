package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class waterOtherGr {
	ArrayList<Double> DistD= new ArrayList<Double>();// dis for all the bins 
	ArrayList<gr> grTime = new ArrayList<gr>();
	HashMap<Integer,Double> disCount= new HashMap<Integer,Double>();
	HashMap<Integer,HashMap<Integer,water>> timeWater = new HashMap<Integer,HashMap<Integer,water>>();
	HashMap<Integer,Double> timeDen = new HashMap<Integer,Double> (); //time - system density

	ArrayList<Integer> timelist = new ArrayList<Integer>();
	public waterOtherGr() {
		
	}
	public gr CalOneSnap(double binsize,double xl,double xh,double yl,double yh,double zl,double zh,HashMap<Integer,water> waterMol,HashMap<Integer,atomGroup> glyMol) {
		gr oneTimeGr = new gr(binsize,(int)(0.5*(xh-xl)/binsize));
		for(Integer gid:glyMol.keySet()) {
		//	System.out.println("gly size is "+glyMol.size());
			atomGroup onegly=glyMol.get(gid);
			oneTimeGr.addAverNum();
			for(Integer gid2:glyMol.keySet()) {
				atomGroup twogly = glyMol.get(gid2);
				double dis = caldis(onegly,twogly,xl,xh,yl,yh,zl,zh);
				if(dis!=0) {
					oneTimeGr.addDis(dis);									
				}
			}
		}
		oneTimeGr.CalculateDist();
		return oneTimeGr;
		
	}
	
	public double caldis(atomGroup a,atomGroup b,double xl,double xh,double yl,double yh,double zl,double zh) {
		double dis=0.0;
		double dx = a.comx-b.comx;
		double dy = a.comy-b.comy;
		double dz = a.comz-b.comz;
		dx = dx - (xh-xl)*Math.round(dx/(xh-xl));
		dy = dy - (yh-yl)*Math.round(dy/(yh-yl));
		dz = dz - (zh-zl)*Math.round(dx/(zh-zl));
		dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	
		return dis;
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

}
