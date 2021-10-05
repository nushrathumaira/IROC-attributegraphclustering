package DocAGnew;

import java.util.ArrayList;
import attributenetwork.vertex;
import attributenetwork.attributenetwork;

public class Cluster {
	attributenetwork an;
	public ArrayList<vertex> Nodes;
	public ArrayList<Integer> Features;	
	int[][] adjStructure;
	int[][] adjFeature;
//	int[][] FulladjFeature;
	
	public Cluster(attributenetwork an){
		this.an = an;
		Nodes = new ArrayList<vertex>();
		Features = new ArrayList<Integer>();    // Feature order		
//		adjStructure = new int[Nodes.size()][Nodes.size()];
//		adjFeature = new int[Nodes.size()][Features.size()];
	}
	
	public Cluster copy(){
		Cluster newObject = new Cluster(an);
		for(vertex v : Nodes){
//			vertex nv = new vertex(v.id);
			newObject.Nodes.add(v);
		}
		for(Integer f : Features){
			newObject.Features.add(f);
		}	
		return newObject;
	}
	
	public void clear(){
		Nodes.clear();
		Features.clear();
//		adjStructure = new int[Nodes.size()][Nodes.size()];
//		adjFeature = new int[Nodes.size()][Features.size()];
	}
	
	public int[][] FindadjStructure(){
		adjStructure = new int[Nodes.size()][Nodes.size()];
		for(int i = 0; i < Nodes.size(); i ++){
			for(int j = 0; j < Nodes.size(); j ++){
				if(an.exsitEdge(Nodes.get(i).id, Nodes.get(j).id)){
					adjStructure[i][j] = 1;
				}
				else{
					adjStructure[i][j] = 0;
				}
			}
		}
		return adjStructure;
	}
	
	public int[][] FindadjFeature(){
		adjFeature = new int[Nodes.size()][Features.size()];
//		System.out.println("Nodes Size :   " + Nodes.size());
//		System.out.println("Feature Size :   " + Features.size());
//		System.out.println("Adj Size1 :   " + adjFeature.length);
//		System.out.println("Adj Size2 :   " + adjFeature[0].length);
		
		for(int i = 0 ; i < Nodes.size(); i ++){
			for(int j = 0; j < Features.size(); j ++){
//				System.out.println(Nodes.get(i).feature.get(j));
				adjFeature[i][j] = Nodes.get(i).feature.get(Features.get(j));
			}
		}
		return adjFeature;
	}
	
//	public int[][] FindFulladjFeature(){
////		Clustering clu = new Clustering(an);
//		FulladjFeature = new int[Nodes.size()][an.numCategory.length];
//		for(int i = 0 ; i < Nodes.size(); i ++){
//			for(int j = 0; j < an.numCategory.length; j ++){
////				System.out.println(Nodes.get(i).feature.get(j));
//				FulladjFeature[i][j] = Nodes.get(i).feature.get(j);
//			}
//		}
//		return FulladjFeature;
//	}
	
	public boolean equals(Object o){
		Cluster c = (Cluster)o;
		if(c.Nodes.containsAll(Nodes) && Nodes.containsAll(c.Nodes) && c.Features.containsAll(Features) && Features.containsAll(c.Features))
			return true;
		return false;
	}
	
	public boolean equalsNodes(Object o){
		Cluster c = (Cluster)o;
		if(c.Nodes.containsAll(Nodes) && Nodes.containsAll(c.Nodes))
			return true;
		return false;
	}

}
