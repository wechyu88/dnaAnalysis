package analysisLammpsDump;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

//calculate residence time based on the paper
//Computation of the Mean Residence Time of Water in the
//Hydration Shells of Biomolecules
public class waterResTime {
	int[] AList;//exist in the shell or not arraylist. in the paper it is A array
	int id;//water molecule id;
	long[] BList;//in the paper, AList is converted to BList for calculation.
	// 111, 3,2,1, 1111,4,3,2,1
	int[] EList;
	public waterResTime(int length,int id) {
		AList = new int[length];
		BList = new long[length];
		this.id = id;
		EList = new int[length];		
	}
	public void setAllList() {
		AtoEList();
		EtoBList();
	}
	public void AtoEList() {
		int start=0;
		int end=start;
		int prev = 0;
		//System.out.println("id "+this.id);
		for(int i=0;i<AList.length;i++) {
			if(i==AList.length-1&&AList[i]==1) {
				end=i;
				int templength = end - start +1;
				EList[templength-1]++;
				
			}else if(prev==0&&AList[i]==1) {
				start=i;
				prev=1;
				
			}else if(prev==1&&AList[i]==0) {
				end = i-1;
				int templength = end - start +1;
				EList[templength-1]++;
				prev = 0;
				
			}
			
		}
	}
	public void EtoBList() {
		int[][] M = new int[AList.length][AList.length];
		for(int i=0;i<M.length;i++) {
			for(int j=0;j<M.length;j++) {
				if(j<i) {
					M[i][j]=0;
				}else {
					M[i][j]=j-i+1;
				}
			}
		}
		BList=multiply(M,EList);
		
	}
	public static long[] multiply(int[][] matrix, int[] vector) {
	    return Arrays.stream(matrix)
	                 .mapToLong(row -> 
	                    IntStream.range(0, row.length)
	                             .mapToLong(col -> row[col] * vector[col])
	                             .sum()
	                 ).toArray();
	}
}
