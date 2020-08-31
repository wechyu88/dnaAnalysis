package analysisLammpsDump;
import java.util.ArrayList;

public class Atom {
	int id;
	double mass;
	ArrayList<Atom> connections;
	double x;
	double y;
	double z;
	double charge=0;
	double vx;
	double vy;
	double vz;
	boolean visited = false;
	int type;
	int moelculeid;
	public Atom(int id) {
		this.id=id;
	}
	public Atom() {
		
	}
	public void setMass(double m) {
		this.mass=m;
	}
	public void setCharge(double charge) {
		this.charge=charge;
	}
	public void setXYZ(double x,double y, double z) {
		this.x=x;
		this.y=y;
		this.z=z;
	};
	public void setV(double vx,double vy,double vz) {
		this.vx=vx;
		this.vy=vy;
		this.vz=vz;
	}
	public void setVisit(boolean visited) {
		this.visited=visited;
	}
	public void setType(int type) {
		this.type = type;
	}
	public void setMoleculeId(int mid) {
		this.moelculeid=mid;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(charge);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((connections == null) ? 0 : connections.hashCode());
		result = prime * result + id;
		temp = Double.doubleToLongBits(mass);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + moelculeid;
		result = prime * result + type;
		result = prime * result + (visited ? 1231 : 1237);
		temp = Double.doubleToLongBits(x);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(y);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(z);
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
		Atom other = (Atom) obj;
		if (Double.doubleToLongBits(charge) != Double.doubleToLongBits(other.charge))
			return false;
		if (connections == null) {
			if (other.connections != null)
				return false;
		} else if (!connections.equals(other.connections))
			return false;
		if (id != other.id)
			return false;
		if (Double.doubleToLongBits(mass) != Double.doubleToLongBits(other.mass))
			return false;
		if (moelculeid != other.moelculeid)
			return false;
		if (type != other.type)
			return false;
		if (visited != other.visited)
			return false;
		if (Double.doubleToLongBits(x) != Double.doubleToLongBits(other.x))
			return false;
		if (Double.doubleToLongBits(y) != Double.doubleToLongBits(other.y))
			return false;
		if (Double.doubleToLongBits(z) != Double.doubleToLongBits(other.z))
			return false;
		return true;
	}
	public void setConnections(ArrayList<Atom> atoms) {
		this.connections=atoms;
	}
	
}
