package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class prdfDNAWater {
	prdfBin [][][] bins;
	double xlo,xhi;
	double ylo,yhi;
	double zlo,zhi;
	int xn,yn,zn;// number of bins on each direction
	double xbin,ybin,zbin;//binsize of x ,y,z
	public prdfDNAWater(double xlo,double xhi,double ylo,double yhi,double zlo,double zhi,int xn,int yn,int zn) {
		this.xlo=xlo;this.xhi=xhi;
		this.ylo=ylo;this.yhi=yhi;
		this.zlo=zlo;this.zhi=zhi;
		this.xn=xn;this.yn=yn;this.zn=zn;
		xbin=(xhi-xlo)/((double)1.0*xn);
		ybin=(yhi-ylo)/((double)1.0*yn);
		zbin=(zhi-zlo)/((double)1.0*zn);
		this.bins = new prdfBin[xn][yn][zn];
		for(int i=0;i<xn;i++) {
			for(int j=0;j<yn;j++) {
				for(int k=0;k<zn;k++) {
					bins[i][j][k]=generatePrdfBin(i,j,k);
					
				}
				
			}
			
		}
		
	}
	
	public prdfBin generatePrdfBin(int x,int y, int z) {//x from 0 to xn-1;
		double cx = xlo+(x+(double)0.5)*(xhi-xlo)/((double)1.0*xn);
		double cy = ylo+(y+(double)0.5)*(yhi-ylo)/((double)1.0*yn);
		double cz = zlo+(z+(double)0.5)*(zhi-zlo)/((double)1.0*zn);
		prdfBin oneBin = new prdfBin(x, y,z,cx,cy,cz);
		return oneBin;
		
	}
	public void setWater(HashMap<Integer,water> waterList) {
		for(int i:waterList.keySet()) {
			water w = waterList.get(i);

			int xIndex;
			int yIndex;
			int zIndex;
			
			
			double dis=0.0;
			double dx = w.Oxygen.x-xlo;
			double dy = w.Oxygen.y-ylo;
			double dz = w.Oxygen.z-zlo;
			dx = dx - (xhi-xlo)*Math.round(dx/(xhi-xlo));
			dy = dy - (yhi-ylo)*Math.round(dy/(yhi-ylo));
			dz = dz - (zhi-zlo)*Math.round(dz/(zhi-zlo));
			while(dx<0) dx+=(xhi-xlo);
			while(dy<0) dy+=(yhi-ylo);
			while(dz<0) dz+=(zhi-zlo);
			//dis = Math.sqrt(dx*dx+dy*dy+dz*dz);	
			xIndex = (int)(dx/xbin);
			yIndex = (int)(dy/ybin);
			zIndex = (int)(dz/zbin);
			//System.out.println(xIndex+" "+yIndex+" "+zIndex+"\n");
			if(xIndex>=xn||xIndex<0) {
				System.out.println("xIndex is not in bound "+xIndex+" xn "+xn+" dx "+dx+" bin "+xbin);

			}
			if(yIndex>=yn||yIndex<0) {
				System.out.println("yIndex is not in bound "+yIndex+" yn "+yn+" dy "+dy+" bin "+zbin);

			}
			if(zIndex>zn||zIndex<0) {
				System.out.println("zIndex is not in bound "+zIndex+" zn "+zn+" dz "+dz+" bin "+zbin);

			}
			
			bins[xIndex][yIndex][zIndex].addAtom(w);
		}
		
	}
		
	
}
