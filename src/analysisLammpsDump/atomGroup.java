package analysisLammpsDump;

import java.util.ArrayList;
import java.util.HashMap;

public class atomGroup {
	ArrayList<Atom> atomlist = new ArrayList<Atom>();
	HashMap<Integer,Atom> idToAtom = new HashMap<Integer,Atom>();
	float comx;//center of mass of x
	float comy;//center of mass of y
	float comz;
	float cocx;//center of charge of x
	float cocy;
	float cocz;
	float totalC;
	float totalM;
	float dipolex;//unit dipole
	float dipoley;
	float dipolez;
	float dipoleVx;//dipole without normalize
	float dipoleVy;
	float dipoleVz;
	float theta;//the angle between vector and z axis
	float phi;// the angle between vector projection on xy plane to x axis
	float coVx;
	float coVy;
	float coVz;
	public atomGroup() {
		
	}
	public void addAtom(Atom a) {
		atomlist.add(a);
	}
	public void deleteAtom(Atom b) {
		if(atomlist.contains(b)) {
			atomlist.remove(b);
		}
	}
	public void setIdToAtom() {
		for(Atom a:atomlist) {
			idToAtom.put(a.id,a);
			
		}
	}
	public void caldipole(float xlo,float xhi,float ylo,float yhi,float zlo,float zhi) {
		float sumdx=(float) 0.0;
		float sumdy=(float) 0.0;
		float sumdz=(float) 0.0;
		
		Atom firstatom = atomlist.get(0);
		
		  for(Atom a:atomlist) { float dx = (float) (a.x-firstatom.x); float dy =
		  (float) (a.y-firstatom.y); float dz = (float) (a.z-firstatom.z); dx = dx -
		  (xhi-xlo)*Math.round(dx/(xhi-xlo)); dy = dy -
		  (yhi-ylo)*Math.round(dy/(yhi-ylo)); dz = dz -
		  (zhi-zlo)*Math.round(dx/(zhi-zlo)); sumdx+=(dx+firstatom.x)*a.charge;
		  sumdy+=(dy+firstatom.y)*a.charge; sumdz+=(dz+firstatom.z)*a.charge; }
		 
		
		/*
		 * for(Atom a:atomlist) {
		 * 
		 * sumdx+=a.x*a.charge; sumdy+=a.y*a.charge; sumdz+=a.z*a.charge; }
		 */
		
		
		float length = (float) Math.sqrt(Math.pow(sumdx, 2)+Math.pow(sumdy, 2)+Math.pow(sumdz, 2));
		this.dipolex=sumdx/length;
		this.dipoley=sumdy/length;
		this.dipolez=sumdz/length;
		
		this.theta = (float) Math.toDegrees(Math.atan(Math.sqrt(this.dipolex*this.dipolex+this.dipoley*this.dipoley)/this.dipolez));
		this.phi = (float) Math.toDegrees(Math.atan(this.dipoley/this.dipolex));
		
		
	}
	public void calCoV() {
		float sumVx=(float) 0.0;
		float sumVy=(float) 0.0;
		float sumVz=(float) 0.0;
		float sumM=0;
		for(Atom a:atomlist) {
			sumVx+=a.vx*a.mass;
			sumVy+=a.vy*a.mass;
			sumVz+=a.vz*a.mass;
			sumM+=a.mass;
		}
		//float length = Math.sqrt(Math.pow(sumdx, 2)+Math.pow(sumdy, 2)+Math.pow(sumdz, 2));
		this.coVx=sumVx/sumM;
		this.coVy=sumVy/sumM;
		this.coVz=sumVz/sumM;

		
	}
	public void calglydipole() {
		float sumdx=(float) 0.0;
		float sumdy=(float) 0.0;
		float sumdz=(float) 0.0;
		float midx=(float) 0.0,midy=(float) 0.0,midz=(float) 0.0;
		float leftx=(float) 0.0,lefty=(float) 0.0,leftz=(float) 0.0;//type 64 id small is left
		float rightx=(float) 0.0,righty=(float) 0.0,rightz=(float) 0.0;//type 64 id large is right
		//glycerol 12 middle carbon, 64 side carbon;
		ArrayList<Atom> carbons =new ArrayList<Atom>();
		for(Atom a:atomlist) {
			if(a.type==12) {
				midx=(float) a.x;
				midy=(float) a.y;
				midz=(float) a.z;
			}else if(a.type==64) {
				carbons.add(a);
			}
		}
		Atom leftC,rightC;
		if(carbons.get(0).id<carbons.get(1).id) {
			leftC=carbons.get(0);
			rightC=carbons.get(1);
		}else {
			leftC=carbons.get(1);
			rightC=carbons.get(0);
		}
		leftx=(float) ((leftC.x+midx)/2);
		lefty=(float) ((leftC.y+midy)/2);
		leftz=(float) ((leftC.z+midz)/2);
		rightx=(float) ((rightC.x+midx)/2);
		righty=(float) ((rightC.y+midy)/2);
		rightz=(float) ((rightC.z+midz)/2);
		sumdx=rightx-leftx;
		sumdy=righty-lefty;
		sumdz=rightz-leftz;
		float length=(float) Math.sqrt(sumdx*sumdx+sumdy*sumdy+sumdz*sumdz);
		this.dipolex=sumdx/length;
		this.dipoley=sumdy/length;
		this.dipolez=sumdz/length;
		
		this.theta = (float) Math.toDegrees(Math.atan(Math.sqrt(this.dipolex*this.dipolex+this.dipoley*this.dipoley)/this.dipolez));
		this.phi = (float) Math.toDegrees(Math.atan(this.dipoley/this.dipolex));
		
		
	}
	
	public void caldipoleV() {
		float sumdx=(float) 0.0;
		float sumdy=(float) 0.0;
		float sumdz=(float) 0.0;
		for(Atom a:atomlist) {
			sumdx+=a.x*a.charge;
			sumdy+=a.y*a.charge;
			sumdz+=a.z*a.charge;
		}
		//float length = Math.sqrt(Math.pow(sumdx, 2)+Math.pow(sumdy, 2)+Math.pow(sumdz, 2));
		this.dipoleVx=sumdx;
		this.dipoleVy=sumdy;
		this.dipoleVz=sumdz;
		this.theta = (float) Math.atan(Math.sqrt(this.dipoleVx*this.dipoleVx+this.dipoleVy*this.dipoleVy)/this.dipoleVz);
		this.phi = (float) Math.atan(this.dipoleVy/this.dipoleVx);
		
	}	
	public ArrayList<Float> calcom(float xlo,float xhi,float ylo,float yhi,float zlo,float zhi) {
		ArrayList<Float> com=new ArrayList<Float>() ;
		float sumM=(float) 0.0;
		float sumx=(float) 0.0;
		float sumy=(float) 0.0;
		float sumz=(float) 0.0;
		Atom firstatom = atomlist.get(0);
		for(Atom a:atomlist) {
			float dx = (float) (a.x-firstatom.x);
			float dy = (float) (a.y-firstatom.y);
			float dz = (float) (a.z-firstatom.z);
			dx = dx - (xhi-xlo)*Math.round(dx/(xhi-xlo));
			dy = dy - (yhi-ylo)*Math.round(dy/(yhi-ylo));
			dz = dz - (zhi-zlo)*Math.round(dx/(zhi-zlo));
			sumx+=(dx+firstatom.x)*a.mass;
			sumy+=(dy+firstatom.y)*a.mass;
			sumz+=(dz+firstatom.z)*a.mass;
			sumM+=a.mass;
		}
		totalM=sumM;
		sumx/=sumM;
		sumy/=sumM;
		sumz/=sumM;
		com.add(sumx);
		com.add(sumy);
		com.add(sumz);
		comx=sumx;
		comy=sumy;
		comz=sumz;
		return com;
	}
	
	public ArrayList<Float> calcoc() {//calculate center of charge
		ArrayList<Float> coc=new ArrayList<Float>() ;
		float sumC=(float) 0.0;
		float sumx=(float) 0.0;
		float sumy=(float) 0.0;
		float sumz=(float) 0.0;
		for(Atom a:atomlist) {
			sumx+=a.x*a.charge;
			sumy+=a.y*a.charge;
			sumz+=a.z*a.charge;
			sumC+=a.charge;
		}
		totalC=sumC;
		sumx/=sumC;
		sumy/=sumC;
		sumz/=sumC;
		coc.add(sumx);
		coc.add(sumy);
		coc.add(sumz);
		cocx=sumx;
		cocy=sumy;
		cocz=sumz;
		return coc;
	}
}
