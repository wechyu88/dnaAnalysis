package analysisLammpsDump;

import java.util.HashMap;

public class timeDipoleList {
	int timestep;
	HashMap<Integer,dipole> moleIdDipole = new HashMap<Integer,dipole> ();
	public timeDipoleList(int timestep,HashMap<Integer,dipole> moleIdDipole) {
		timestep=timestep;
		
	}
}
