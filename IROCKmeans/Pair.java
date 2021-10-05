package IROCKmeans;

import java.util.*;
import attributenetwork.vertex;

public class Pair implements Comparable<Pair> {
	
	public double similarity;	
	public Cluster [] member;
	public ArrayList<vertex> structureCN;
	public ArrayList<Integer> featureCN;

	public Pair(){
		member = new Cluster[2];	
		structureCN = new ArrayList<vertex>();
		featureCN = new ArrayList<Integer>();
	}
	
	public int compareTo(Pair p){
		if((member[0] == p.member[0] && member[1] == p.member[1])
			|| (member[0] == p.member[1] && member[1] == p.member[0]))
			return 0;
		else if(similarity < p.similarity)
			return 1;
		else
			return -1;
	}
}
