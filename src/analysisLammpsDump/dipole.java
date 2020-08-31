package analysisLammpsDump;

public class dipole {
	double dipoleX;
	double dipoleY;
	double dipoleZ;
	double theta;
	double phi;
	public dipole(double x,double y,double z) {
		this.dipoleX=x;
		this.dipoleY=y;
		this.dipoleZ=z;
		double tempTheta = Math.toDegrees(Math.atan(Math.sqrt(this.dipoleX*this.dipoleX+this.dipoleY*this.dipoleY)/this.dipoleZ));
		if(tempTheta<0) {
			this.theta=tempTheta+180;
		}else {
			this.theta=tempTheta;
		}
		double tempPhi = Math.toDegrees(Math.atan(this.dipoleY/this.dipoleX));
		if(this.dipoleX>0) {
			if(tempPhi<0) {
				this.phi=tempPhi+360;
			}else {
				this.phi=tempPhi;
			}
		
		}else if(this.dipoleX<0) {
			this.phi=tempPhi+180;
		}

		
		
	}
}
