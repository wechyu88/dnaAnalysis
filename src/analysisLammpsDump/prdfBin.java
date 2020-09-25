package analysisLammpsDump;

import java.util.ArrayList;

public class prdfBin {
	double cx, cy,cz;//center coordinates x,y,z
	ArrayList<water> inSideWater;
	double dis;	
	int x,y,z;//indexes of this bin in the whole cubic box
	public prdfBin(int x,int y,int z, double cx, double cy,double cz) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		this.inSideWater=new ArrayList<water>();
	}
	public void addAtom(water a) {
		this.inSideWater.add(a);
	}
	public void setDis(double dis) {
		this.dis=dis;
	}
	

}
