package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class dumpOneStep {
	int timestep;
	int numAtom;
	ArrayList<Atom> atomlist = new ArrayList<Atom>();
	float xlo;
	float xhi;
	float ylo;
	float yhi;
	float zlo;
	float zhi;
	HashMap<Integer,Atom> atomMap = new HashMap<Integer,Atom>();
	public dumpOneStep() {
		
	}
	public void setTime(int time) {
		this.timestep = time;
	}
	public void setnumAtom(int numA) {
		this.numAtom = numA;
	}
	public void setAtomlist(ArrayList<Atom> atomlist) {
		this.atomlist=atomlist;
		for(Atom a:atomlist) {
			atomMap.put(a.id, a);
		}
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((atomMap == null) ? 0 : atomMap.hashCode());
		result = prime * result + ((atomlist == null) ? 0 : atomlist.hashCode());
		result = prime * result + numAtom;
		result = prime * result + timestep;
		long temp;
		temp = Double.doubleToLongBits(xhi);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(xlo);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(yhi);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(ylo);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(zhi);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(zlo);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		dumpOneStep other = (dumpOneStep) obj;
		if (atomMap == null) {
			if (other.atomMap != null)
				return false;
		} else if (!atomMap.equals(other.atomMap))
			return false;
		if (atomlist == null) {
			if (other.atomlist != null)
				return false;
		} else if (!atomlist.equals(other.atomlist))
			return false;
		if (numAtom != other.numAtom)
			return false;
		if (timestep != other.timestep)
			return false;
		if (Double.doubleToLongBits(xhi) != Double.doubleToLongBits(other.xhi))
			return false;
		if (Double.doubleToLongBits(xlo) != Double.doubleToLongBits(other.xlo))
			return false;
		if (Double.doubleToLongBits(yhi) != Double.doubleToLongBits(other.yhi))
			return false;
		if (Double.doubleToLongBits(ylo) != Double.doubleToLongBits(other.ylo))
			return false;
		if (Double.doubleToLongBits(zhi) != Double.doubleToLongBits(other.zhi))
			return false;
		if (Double.doubleToLongBits(zlo) != Double.doubleToLongBits(other.zlo))
			return false;
		return true;
	}
	public void setboundary(float xl,float xh,float yl,float yh,float zl,float zh) {
		xlo=xl;
		xhi=xh;
		ylo=yl;
		yhi=yh;
		zlo=zl;
		zhi=zh;
	}
}
