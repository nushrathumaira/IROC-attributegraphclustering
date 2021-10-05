package AttributeGraph;
import subcat.datagenerator;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import attributenetwork.GraphGenerator;
import attributenetwork.attributenetwork;
import attributenetwork.IO;
import attributenetwork.vertex;
import java.util.*;

public class Main {
	public static void main(String[] args){		
		attributenetwork an = new attributenetwork();
		
//		ArrayList<Cluster> TestList = new ArrayList<Cluster>();
//		ArrayList<vertex> Ver = new ArrayList<vertex>();
//		ArrayList<Integer> Fea = new ArrayList<Integer>();
//		for(int i = 0; i < 10; i ++){
//			vertex v = new vertex(i);
//			Ver.add(v);
//			Fea.add(i);
//		}
//		Cluster C1 = new Cluster(an);		
//		for(int i = 0; i < 3; i ++){						
//			C1.Nodes.add(Ver.get(i));
//			C1.Features.add(Fea.get(i));
//		}
//		TestList.add(C1);
//		Cluster C2 = new Cluster(an);
//		for(int i = 4; i < 7; i ++){						
//			C2.Nodes.add(Ver.get(i));
//			C2.Features.add(Fea.get(i));
//		}
//		TestList.add(C2);
//		Cluster C3 = new Cluster(an);
//		for(int i = 7; i < 9; i ++){			
//			C3.Nodes.add(Ver.get(i));
//			C3.Features.add(Fea.get(i));
//		}
//		TestList.add(C3);
//		Cluster C4 = new Cluster(an);
//		for(int i = 2; i < 7; i ++){		
//			C4.Nodes.add(Ver.get(i));
//			C4.Features.add(Fea.get(i));
//		}
//		TestList.add(C4);
//		Cluster C5 = new Cluster(an);
//		for(int i = 4; i < 9; i ++){					
//			C5.Nodes.add(Ver.get(i));
//			C5.Features.add(Fea.get(i));
//		}
//		TestList.add(C5);
//		for(int i = 0; i < TestList.size(); i ++){
//			for(int j = 0; j < TestList.get(i).Features.size(); j ++){
//				System.out.println(TestList.get(i).Features.get(j));
//			}
//			System.out.println("end");
//		}
//		Clustering clu = new Clustering(an);
//		ArrayList<ArrayList<ArrayList<vertex>>> newNodeID = clu.CombineNodeID(TestList);
//		double Cost = clu.CalculateCC (TestList); 

		
		IO ea = new IO();		
		an = ea.readNetwork("data/syn/test");
		System.out.println("Read Attributed Network");
		CodingCost cc = new CodingCost(an);
		//Initialization
	    Clustering clu = new Clustering(an);	    
	    ArrayList<Cluster> InitialClusList = clu.FindInitialClus();	    
	    
//	    System.out.println("Ego Size:  " + InitialClusList.size());
//	    for(int i = 0; i < InitialClusList.size(); i++){
//	    	for(int j = 0; j < InitialClusList.get(i).Nodes.size(); j ++ ){
//	    		System.out.print(InitialClusList.get(i).Nodes.get(j).id + ",  ");	    		
//	    	}
//	    	System.out.print("\n");
//	    }
	    
	    //Remove Redundency of Cluster List
//    	ArrayList<Cluster> ReClusList = clu.RemoveRedundency(InitialClusList);    	
//    	System.out.println("Removed Redundency:  " +  ReClusList.size());
//  	    for(int i = 0; i < ReClusList.size(); i++){
//  	    	for(int j = 0; j < ReClusList.get(i).Nodes.size(); j ++ ){
//  	    		System.out.print(ReClusList.get(i).Nodes.get(j).id + ",  ");	    		
//  	    	}
//  	    	System.out.print("\n");
//  	    }
//    	System.out.println(ReClusList.size()); 
	    
	    
	    // Find subspaces of all clusters
//	    ArrayList<Cluster> ClusList = clu.Findsubspace(ReClusList);
	    ArrayList<Cluster> ClusList = clu.Findsubspace(InitialClusList);
	    
//	    // Output CLuster List
//	    int id = 0;
//	    for(Cluster c : ClusList){
//	    	System.out.println("Cluster" + id + " : ");
//	    	for(vertex v : c.Nodes)
//	    		System.out.println(" Nodes:  " + v.id);	
//	    	for(Integer f : c.Features)
//	    		System.out.println(" Feaure:  " + f);	    
//	    	id ++;
//	    }
	    	
    	double oldCost = cc.CalculateCC (ClusList);    	 
    	double newCost = 0.0;

    	// Rough Merge
    	HashMap<ArrayList<Cluster>, Double> RoughMerger = new HashMap<ArrayList<Cluster>, Double>(); 
    	
    	ArrayList<Cluster> FeatureClusList = new ArrayList<Cluster>();
    	FeatureClusList.addAll(ClusList);  // ClusList only change features
    	ArrayList<Cluster> StructureClusList = new ArrayList<Cluster>();
    	StructureClusList.addAll(ClusList); // ClusList only change structure          
    	
    	//Initialization of all similarity between any two clusters    
    	TreeSet<Pair> simQueue = new TreeSet<Pair>();
    	for(int i = 0; i < ClusList.size() - 1; i ++){
    		for(int j = i + 1; j < ClusList.size(); j ++){
    			Pair simPair = clu.CalculateSim(ClusList.get(i), ClusList.get(j));
    			simQueue.add(simPair);
    		}
    	}  
    	
    	for(Pair p : simQueue){
    		System.out.println("Similarity: " + p.similarity);
    	}
    	
	    while(simQueue.size() > 1){	    	
	    	Pair mergeCandidate = simQueue.first();
	    	//Output most similar clusters
	    	for(int t = 0; t < mergeCandidate.member.length; t ++){
	    		System.out.println();
	    		for(vertex v : mergeCandidate.member[t].Nodes)
	    			System.out.print(v.id);
	    	}
	    	System.out.println();
	    	
	    	simQueue.remove(mergeCandidate);
	    	// Change Structure Part		
			ArrayList<Cluster> newClusters =  clu.ModifyStructure(StructureClusList, FeatureClusList, mergeCandidate);
			//Output modified clusters
			for(int t = 0; t < newClusters.size(); t ++){
				System.out.println("new Clusters" + t + ":");
				for(vertex v : newClusters.get(t).Nodes){
					System.out.print(v.id);
				}
			}
			System.out.println();
//			ArrayList<Cluster> tempClusListSRe = clu.RemoveRedundency(tempClusListS);
			
			ArrayList<Cluster> tempStructureClusList = new ArrayList<Cluster>();
			tempStructureClusList.addAll(StructureClusList);
			for(int i = 0 ; i < mergeCandidate.member.length; i ++){
				tempStructureClusList.remove(mergeCandidate.member[i]);
			}
			for(Cluster c : newClusters){
				tempStructureClusList.add(c);
			}
			newCost = cc.CalculateCC(tempStructureClusList);
			if(newCost < oldCost){
				StructureClusList.clear();
				StructureClusList.addAll(tempStructureClusList);
				simQueue = clu.ModifyQueue(tempStructureClusList,simQueue,newClusters, mergeCandidate);				
			}			
			System.out.println("cc without subspace:" + newCost);
			
					
			
		    // Find features
			ArrayList<Cluster> tempClusListF = clu.Findsubspace(StructureClusList);
			FeatureClusList.clear();
			FeatureClusList.addAll(StructureClusList);	
			double cctemp = cc.CalculateCC(StructureClusList);
			System.out.println("cc with subspace:" + cctemp);
			RoughMerger.put(StructureClusList, cctemp);
			oldCost = cctemp;			
			System.out.println("Size:  " + FeatureClusList.size());
			
//			
//			Output test
			  int id1 = 0;
			    for(Cluster c : StructureClusList){
			    	System.out.println("Cluster" + id1 + ": ");
			    	System.out.println("Nodes:  ");
			    	for(vertex v : c.Nodes)
			    		System.out.print(v.id+",");
			    	System.out.println("");
			    	System.out.println("Feature:  ");
			    	for(Integer f : c.Features)
			    		System.out.print(f+",");
			    	System.out.println("");
			    	id1 ++;
			    }
		}  
	    
	    double cc_roughmerge = Double.MAX_VALUE;
	    ArrayList<Cluster> RoughClusList = new ArrayList<Cluster>();
	    
	    // Detail Modify
	    
	}
}
