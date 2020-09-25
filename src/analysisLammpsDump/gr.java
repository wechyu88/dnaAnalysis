package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class gr {
	double binsize;
	int maxbinnum;
	ArrayList<Double> dis = new ArrayList<Double>();
	int numAver=0;
	ArrayList<Double> DistD= new ArrayList<Double>();// dis for all the bins 
	ArrayList<prdfBin> DistBin = new ArrayList<prdfBin>();
	int cubicNum=0;// used for prdf, number of cubic
	HashMap<Integer,Integer> disCountBinNum = new HashMap<Integer,Integer>();//in each dis, how many cubics at this dis
	HashMap<Integer,Integer> disCountWaterNum = new HashMap<Integer,Integer>();//in each dis, how many water molecules at this dis
	double cubicV;
	double waterDen;
	HashMap<Integer,Double> disCount= new HashMap<Integer,Double>();
	public gr(double binsize,int maxbinnum) {
		this.binsize=binsize;
		this.maxbinnum=maxbinnum;
		for(int i=0;i<maxbinnum;i++) {
			DistD.add((i+0.5)*binsize);
		}
	}
	public gr(double binsize,int maxbinnum,double cubicV,double waterDen) {
		this.binsize=binsize;
		this.maxbinnum=maxbinnum;
		for(int i=0;i<maxbinnum;i++) {
			DistD.add((i+0.5)*binsize);
		}
		this.cubicV = cubicV;
		this.waterDen=waterDen;
	}
	public void addprdfBin(prdfBin prdfb) {
		DistBin.add(prdfb);
		
	}
	public void addDis(Double dis) {
		this.dis.add(dis);
	}
	public void addAverNum() {
		// such as count oxygen around gly, 
		//for each gly we can count once, add one gly, this number add one
		numAver++;
	}
	public void CalculatePrdf() {
		for(prdfBin prdfb : DistBin) {
			int binnum = (int) (prdfb.dis/binsize);
			if(disCountBinNum.containsKey(binnum)) {
				disCountBinNum.put(binnum,disCountBinNum.get(binnum)+1);
				disCountWaterNum.put(binnum,disCountWaterNum.get(binnum)+prdfb.inSideWater.size());
				if(prdfb.inSideWater.size()>1) {
					//System.out.println("find one cube more than one water");
				}

			}else {
				
				disCountBinNum.put(binnum,1);
				disCountWaterNum.put(binnum,prdfb.inSideWater.size());				
			}
		}
		for(int key:disCountBinNum.keySet()) {
			if(disCountWaterNum.get(key)==0) {
				continue;
			}
			if(disCountWaterNum.get(key)!=disCountBinNum.get(key)) {
				//System.out.println("find a layer bin not same with water");
			}
			double prob = disCountWaterNum.get(key)*1.0/(disCountBinNum.get(key)*cubicV)/waterDen;
			disCount.put(key,prob);
		}
		DistBin.clear();
		
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
