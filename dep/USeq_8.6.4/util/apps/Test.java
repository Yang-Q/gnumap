package util.apps;


import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import util.gen.IO;

public class Test {
	
	public static final void main (String[] args){
		
		String[] pd = {"CL","GV","ICM","MI","MII","MOR","PN","TROPH"};
		File dir = new File("/Users/davidnix/Desktop/RunMe");
		for (int i=0; i< pd.length; i++){
			for (int j=i+1; j< pd.length; j++){
				StringBuilder sb = new StringBuilder();
				sb.append("#e david.nix@hci.utah.edu\n");
				sb.append("#a A1553\n");
				sb.append("pd=/uufs/chpc.utah.edu/common/home/u0028003/Jessie/PointData\n");
				sb.append("mrss=/uufs/chpc.utah.edu/common/home/u0028003/BioApps/USeq/Apps/MultipleReplicaScanSeqs\n");
				sb.append("r=/uufs/chpc.utah.edu/common/home/u0028003/R/R-3.0.1/bin/R\n");
				sb.append("t="+pd[i]+"\n");
				sb.append("c="+pd[j]+"\n");
				sb.append("java -jar -Xmx20G $mrss -r $r -m 20 -p 0 -w 150 -b -s "+pd[i]+"_"+pd[j]+" -t $pd/$t -c $pd/$c\n");
				System.out.println("xxxxxxxxxxxxxxxxx\n"+sb);
				File folder = new File(dir,pd[i]+"_"+pd[j]);
				folder.mkdir();
				File cmd = new File(folder,"cmd.txt");
				IO.writeString(sb.toString(), cmd);
			}
		}
		
		
		
	}
	
	

	
	
	
	
	
	
}
