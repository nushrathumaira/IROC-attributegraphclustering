package AttributeGraph;
import subcat.datagenerator;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import attributenetwork.GraphGenerator;
import attributenetwork.attributenetwork;
import attributenetwork.IO;
import attributenetwork.vertex;
import java.util.*; 

public class NewMain {
	public static void main(String[] args){		
		attributenetwork an = new attributenetwork();
		IO ea = new IO();		
		an = ea.readNetwork("data/syn/2clusters400");		
		System.out.println("Read Attributed Network");
		
		long   start   =   System.currentTimeMillis();
		
//		Cluster C1 = new Cluster(an);
//		for(int i = 0; i < 6; i ++){
//			C1.Nodes.add(an.getVertex(i));
//		}
//		for(int i = 0; i < 4; i ++){
//			C1.Features.add(i);
//		}
//		Cluster C2 = new Cluster(an);
//		for(int i = 4; i < 10; i ++){
//			C2.Nodes.add(an.getVertex(i));
//		}
//		for(int i = 2; i < 6; i ++){
//			C2.Features.add(i);
//		}
//		ArrayList<Cluster> TestList = new ArrayList<Cluster>();
//		TestList.add(C1);
//		TestList.add(C2);
		 
//		
//		//Initialization
	   Clustering clu = new Clustering(an);
//	    double testcost = clu.CalculateCC(TestList);
	    
	    double ccentropy = clu.DataEntropy(an);
	    System.out.println(ccentropy);
	    
	    ArrayList<Cluster> InitialClusList = clu.FindInitialClus();
	    //output
//	    System.out.println("Initial Cluster List without subspace");
//	    for(Cluster c : InitialClusList){	    	
//	    	clu.Output(c);
//	    }	    
	    System.out.println("Size of initial cluster List :"  + InitialClusList.size());
	    ArrayList<Cluster> ClusList = clu.RemoveSame(InitialClusList);
	    
	    ArrayList<Cluster> newList = clu.InitialKClusters(ClusList, 15);
//	    ArrayList<Cluster> newList = clu.InitialKClustersAutoMat(ClusList);
	    System.out.println(newList.size());
	    ClusList.clear();
	    ClusList.addAll(newList);
	    
	  //output
//	    System.out.println("Initial Cluster List without subspace");
//	    for(Cluster c : ClusList){	    	
//	    	clu.Output(c);
//	    }
//		ArrayList<Cluster> ClusList = clu.RemoveRedundency(InitialClusList); 
//		System.out.println("Size of cluster List after remove redundency :"  + ClusList.size());
	    //output
//	    System.out.println("remove same:");
//	    for(Cluster c : ClusList){	    	
//	    	clu.Output(c);
//	    }
	   
	    ArrayList<Cluster> clustersnosubspace = new ArrayList<Cluster>(); // All clusters without subspace
	    clustersnosubspace.addAll(ClusList);
	    // Find subspaces of all clusters
	    ClusList = clu.NewFindSubspace(clustersnosubspace, ClusList);
	    
	    
	    CodingCost cc = new CodingCost(an);
	    double oldCost = cc.CalculateCC (ClusList);
//	    ClusList = clu.FindSubspacewithoutOrder(clustersnosubspace, ClusList);
//	    ClusList = clu.Findsubspace(ClusList);
//	    for(Cluster c : ClusList)
//	    	clu.Output(c);
//	    ClusList = clu.FindSubspacewithoutOrder(clustersnosubspace, ClusList, oldCost - cc.ccUnClusterFeature);
//	    System.out.println("Initial Cluster with subspace: ");
//	    for(Cluster cl : ClusList){
//	    	clu.Output(cl);
//	    }
	    oldCost = cc.CalculateCC (ClusList);    	 
    	double newCost = 0.0;
    	
    	TreeSet<Pair> simQueue = new TreeSet<Pair>();
    	for(int i = 0; i < ClusList.size() - 1; i ++){
    		for(int j = i + 1; j < ClusList.size(); j ++){
    			Pair simPair = clu.CalculateSim(ClusList.get(i), ClusList.get(j));
    			simQueue.add(simPair);
    		}
    	}      
    	
	    while(simQueue.size() > 0){	 
//	    	System.out.println("Round:   ");
//	    	for(Pair p : simQueue){
//	    		System.out.println("SimQueue");
//	    		clu.Output(p.member[0]);
//	    		clu.Output(p.member[1]);
//	    	}
//	    	System.out.println(simQueue.size());
	    	
	    	
	 	    Pair mergeCandidate = simQueue.first();
	 	    
	 	    //Output most similar clusters
//	 	    System.out.println("Merge Clusters :   ");
//	 	    for(int t = 0; t < mergeCandidate.member.length; t ++){
//	 	    	System.out.println();
//	 	    	for(vertex v : mergeCandidate.member[t].Nodes)
//	 	    		System.out.print(v.id);
//	 	    	}
//	 	    System.out.println();
	 	    	
	 	    simQueue.remove(mergeCandidate);
	 	    Cluster clus1 = mergeCandidate.member[0];
	 	    Cluster clus2 = mergeCandidate.member[1];
	 	    if(!ClusList.contains(clus1) || !ClusList.contains(clus2))
	 	    	continue;
	 	    
	 	    // Change Structure Part	
	 	   ArrayList<Cluster> modifiedClusters = new ArrayList<Cluster>();	 	   
//	 	   modifiedClusters =  clu.ModifyClusternew(ClusList, mergeCandidate);
//	 	   modifiedClusters =  clu.ModifyClusterwithoutorder(ClusList, mergeCandidate);
	 	   modifiedClusters =  clu.ModifyClusternew(ClusList, mergeCandidate);
//	 	   modifiedClusters =  clu.ModifyClusterAll(ClusList, mergeCandidate);
//	 	   System.out.println("ModifiedCluster Feature:   ");
//	 	   for(Cluster c : modifiedClusters){
//	 		   clu.Output(c);
//	 	   }
//	 		ArrayList<Cluster> newClusters =  clu.ModifyCluster(ClusList, mergeCandidate);
//	 		for(int t = 0; t < newClusters.size(); t ++){
//	 			System.out.println("new Clusters" + t + ":");
//	 			for(vertex v : newClusters.get(t).Nodes){
//	 				System.out.print(v.id);
//	 			}
//	 		}
//	 		System.out.println();
//	 			ArrayList<Cluster> tempClusListSRe = clu.RemoveRedundency(tempClusListS);
	 			
	 		ArrayList<Cluster> newClusList = new ArrayList<Cluster>();
	 		newClusList.addAll(ClusList);
	 		if(modifiedClusters.size() > 0){
	 			for(int i = 0 ; i < mergeCandidate.member.length; i ++){
	 				newClusList.remove(mergeCandidate.member[i]);
	 			}
	 			for(Cluster c : modifiedClusters){
//	 				clu.Output(c);
	 				newClusList.add(c);
	 			}
	 		}
	 			
	 		newCost = cc.CalculateCC(newClusList);
//	 		System.out.println("newCost:  " +  newCost);
	 		if(newCost < oldCost){
//	 			oldCost = newCost;
	 			ClusList.clear();
	 			ClusList.addAll(newClusList);
	 			ClusList = clu.RemoveSame(ClusList);
	 			simQueue = clu.ModifyQueue(newClusList,simQueue,modifiedClusters, mergeCandidate);				
	 		}		 		
//	 		System.out.println(simQueue.size());
	    }
//	    Pair simPair = clu.CalculateSim(ClusList.get(0), ClusList.get(3));
//	    clu.ModifyClusterwithoutorder(ClusList, simPair);
//		Output test
	      System.out.println("Coding cost :  " +cc.CalculateCC(ClusList));
//	      ClusList.remove(1);
//	      ClusList.get(0).Features.add(3);
//	      ClusList.get(0).Features.add(2);
//	      System.out.println("Coding cost :  " +clu.CalculateCC(ClusList));
	      
	      ArrayList<TreeSet<Integer>> IDV = new ArrayList<TreeSet<Integer>>();
	      ArrayList<TreeSet<Integer>> IDF = new ArrayList<TreeSet<Integer>>();
	      for(Cluster c : ClusList){
	    	  TreeSet<Integer> idvc = new TreeSet<Integer>();
	    	  TreeSet<Integer> idfc = new TreeSet<Integer>();
	    	  for(vertex v : c.Nodes){
	    		  idvc.add(v.id);
	    	  }
	    	  IDV.add(idvc);
	    	  for(Integer f : c.Features){
	    		  idfc.add(f);
	    	  }
	    	  IDF.add(idfc);
	      }
	      int numid = 0;
	      for(int i = 0; i < IDV.size(); i ++){
	    	  System.out.println("Cluster" + numid + ": ");
			  System.out.println("Nodes:  "); 			  
			  for(Integer idv : IDV.get(i)){
				  System.out.print(idv+","); 
			  }
			  System.out.println();
			  System.out.println("Feature:  ");
			  for(Integer idf : IDF.get(i)){
				  System.out.print(idf+","); 
			  }
			  System.out.println();
			  numid ++;
	      }		
	      
	      
	      
	      ClusList = clu.PostProcessing(ClusList);
//	      ClusList = clu.PostProcessingnew(ClusList);
	      System.out.println("Coding cost after pp :  " +cc.CalculateCC(ClusList));	      
//	      ClusList.get(0).Nodes.remove(new vertex(8));
//	      ClusList.get(1).Features.add(1);
//	      ClusList.get(1).Features.add(2);
//	      ClusList.get(1).Features.add(3);
//	      System.out.println("Coding cost after pp :  " +cc.CalculateCC(ClusList));
	      ArrayList<TreeSet<Integer>> IDVP = new ArrayList<TreeSet<Integer>>();
	      ArrayList<TreeSet<Integer>> IDFP = new ArrayList<TreeSet<Integer>>();
	      for(Cluster c : ClusList){
	    	  TreeSet<Integer> idvc = new TreeSet<Integer>();
	    	  TreeSet<Integer> idfc = new TreeSet<Integer>();
	    	  for(vertex v : c.Nodes){
	    		  idvc.add(v.id);
	    	  }
	    	  IDVP.add(idvc);
	    	  for(Integer f : c.Features){
	    		  idfc.add(f);
	    	  }
	    	  IDFP.add(idfc);
	      }
	      int numidP = 0;
	      for(int i = 0; i < IDVP.size(); i ++){
	    	  System.out.println("Cluster" + numidP + ": ");
			  System.out.println("Nodes:  "); 			  
			  for(Integer idv : IDVP.get(i)){
				  System.out.print(idv+","); 
			  }
			  System.out.println();
			  System.out.println("Feature:  ");
			  for(Integer idf : IDFP.get(i)){
				  System.out.print(idf+","); 
			  }
			  System.out.println();
			  numidP ++;
	      }	
//	      ClusList = clu.FinalSubspace(ClusList);
	      
	      long   end   =   System.currentTimeMillis(); 
		  System.out.println("All running time:"+Long.toString(end-start)+" ms.");
	}
}
