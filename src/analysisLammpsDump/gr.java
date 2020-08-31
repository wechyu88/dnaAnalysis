package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class gr {
	double binsize;
	int maxbinnum;
	ArrayList<Double> dis = new ArrayList<Double>();
	int numAver=0;
	ArrayList<Double> DistD= new ArrayList<Double>();// dis for all the bins 
	HashMap<Integer,Double> disCount= new HashMap<Integer,Double>();
	public gr(double binsize,int maxbinnum) {
		this.binsize=binsize;
		this.maxbinnum=maxbinnum;
		for(int i=0;i<maxbinnum;i++) {
			DistD.add((i+0.5)*binsize);
		}
	}

	public void addDis(Double dis) {
		this.dis.add(dis);
	}
	public void addAverNum() {
		// such as count oxygen around gly, 
		//for each gly we can count once, add one gly, this number add one
		numAver++;
	}
	public void CalculateDist() {
		int count=0;
		for(Double onedis:dis) {
			count++;
			int binnum =(int) (onedis/binsize);
			//if(binnum<maxbinnum) {
				//double d = DistD.get(binnum);
				if(disCount.containsKey(binnum)) {
					disCount.put(binnum,disCount.get(binnum)+1.0);
				}else {
					disCount.put(binnum, 1.0);
				}				
			//}
		}
		System.out.println(" glysize here is "+numAver+" disnum is "+count);
		double numcount=0;
		for(Integer dis :disCount.keySet()) {
			numcount+=disCount.get(dis)/(numAver*1.0);
			disCount.put(dis, disCount.get(dis)/(numAver*1.0));
		}
		System.out.println(" glysize here is "+numAver+" disnum is "+numcount);

	}
	public void Calgr(int n) {
		int count=0;
		for(Double onedis:dis) {
			count++;
			int binnum =(int) (onedis/binsize);
			//if(binnum<maxbinnum) {
				//double d = DistD.get(binnum);
				if(disCount.containsKey(binnum)) {
					disCount.put(binnum,disCount.get(binnum)+1.0);
				}else {
					disCount.put(binnum, 1.0);
				}				
			//}
		}
		System.out.println(" glysize here is "+numAver+" disnum is "+count);
		double numcount=0;
		for(Integer dis :disCount.keySet()) {
			numcount+=disCount.get(dis)/(numAver*1.0*n);
			disCount.put(dis, disCount.get(dis)/(numAver*1.0*n));
		}
		System.out.println(" glysize here is "+numAver+" disnum is "+numcount);

	}
}
