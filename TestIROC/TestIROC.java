package TestIROC;

import java.util.*;

import DocAG.Cluster;
import DocAG.CodingScheme;
import DocAG.DocAG;
import DocAG.Pair;
import attributenetwork.IO;
import attributenetwork.attributenetwork;
import attributenetwork.vertex;

public class TestIROC {

	attributenetwork an;
	CodingScheme codingScheme;	
	int numVertices;
	int numFeature;
	int[][] AdjNode;
	int[][] AdjFeature;
	ArrayList<Integer> FullFeature;
	ArrayList<Cluster> ClusList;
	
	public TestIROC(attributenetwork an){
		this.an = an;
		this.numVertices = an.numVertices;
		this.numFeature = an.numFeature;
		FullFeature = new ArrayList<Integer>();  // all features
		for (int i = 0; i < numFeature; i++) {
			FullFeature.add(i);
		}	
		AdjNode = new int[numVertices][numVertices];  // adjacency matrix of nodes
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (an.exsitEdge(i, j))
					AdjNode[i][j] = 1;
			}
		}
		AdjFeature = an.featureAdj; // adjacency matrix of attributes
		ClusList = new ArrayList<Cluster>();
		this.codingScheme = new CodingScheme(an, AdjNode, AdjFeature);
	}	
	
	public ArrayList<Cluster> FindSubspace(ArrayList<Cluster> ClusList){
		for(Cluster clus : ClusList){
			FindClusterSubspace(clus, ClusList);
		}
		return ClusList;
	}
	
	public void FindClusterSubspace(Cluster clus, ArrayList<Cluster> ClusList){
		clus.Features.clear();
		ArrayList<ArrayList<Integer>> clusterfeatureorder = FindClusterFeatureOrder(clus); // feature with order
		
		double costmin = Double.MAX_VALUE;
		ArrayList<Integer> featuremin = new ArrayList<Integer>();
		for (int i = 0; i < clusterfeatureorder.size(); i ++) {
			clus.Features.addAll(clusterfeatureorder.get(i));
//			double cost = codingScheme.CalculateCC(ClusList);
			double cost = codingScheme.CCClusterFeature(clus) + codingScheme.CCUnClusterFeature(ClusList);
			if(cost < costmin){
				costmin = cost;
				featuremin.clear();
				featuremin.addAll(clus.Features);
			}
		}
		clus.Features.clear();
		clus.Features.addAll(featuremin);
	}
	
	public ArrayList<Integer> FeatureOrder(HashMap<Integer, Double> featureent) {
		HashMap<Integer, Double> featureentcopy = new HashMap<Integer, Double>();
		featureentcopy.putAll(featureent);
		ArrayList<Integer> Order = new ArrayList<Integer>();
		while (Order.size() < numFeature) {
			double minent = Double.MAX_VALUE;
			int minf = 0;
			for (Integer f : featureentcopy.keySet()) {
				if (featureentcopy.get(f) < minent) {
					minent = featureentcopy.get(f);
					minf = f;
				}
			}
			Order.add(minf);
			featureentcopy.put(minf, Double.MAX_VALUE);
		}
		return Order;
	}
	
	public int [][] FindFullAdjFeature(Cluster c){
		int [][] FulladjFeature = new int[c.Nodes.size()][an.numCategory.length];
		for(int i = 0 ; i < c.Nodes.size(); i ++){
			for(int j = 0; j < an.numCategory.length; j ++){
//				System.out.println(Nodes.get(i).feature.get(j));
				FulladjFeature[i][j] = c.Nodes.get(i).feature.get(j);
			}
		}
		return FulladjFeature;
	}
	
	private double featureentropy(int[] feature, int ID) {
		int sum = feature.length;
		int[] featurecat = new int[an.numCategory[ID]];
		for (int i = 0; i < sum; i++) {
			featurecat[feature[i]]++;
		}
		double ent = 0.0;
		for (int i = 0; i < featurecat.length; i++) {
			if (featurecat[i] != 0)
				ent += featurecat[i]
						* lg2((double) sum / (double) featurecat[i]);
		}
		return ent;
	}
	
	private double lg2(double d) {
		return Math.log(d) / Math.log(2.0);
	}
	
	public ArrayList<ArrayList<Integer>> FindClusterFeatureOrder(Cluster c){		
		int[][] adjFeature = FindFullAdjFeature(c);
		HashMap<Integer, Double> featureEnt = new HashMap<Integer, Double>();
		for (Integer i : FullFeature) {			
			int[] featurei = new int[c.Nodes.size()];
			for (int j = 0; j < c.Nodes.size(); j++) {
				featurei[j] = adjFeature[j][i];
			}
			double entvalue = featureentropy(featurei, i);
			featureEnt.put(i, entvalue);
		}
		ArrayList<Integer> Order = FeatureOrder(featureEnt);			
		int[] featureLabel = new int[Order.size()]; // Feature Label
		for (int i = 0; i < featureLabel.length; i++) {
			featureLabel[i] = -1;
		}
		ArrayList<Integer> Ordercopy = new ArrayList<Integer>();
		Ordercopy.addAll(Order);
		int num = 0;
		while (Ordercopy.size() > 0) {
			featureLabel[Ordercopy.get(0)] = num;
			for (Integer o : Ordercopy) {
				if (featureEnt.get(Ordercopy.get(0)).equals(featureEnt.get(o))) {
					featureLabel[o] = num;
				}
			}
			for (int i = 0; i < featureLabel.length; i++) {
				if (featureLabel[i] == num) {
					Ordercopy.remove((Object) i);
				}
			}
			num++;
		}
		int labelnum = featureLabel[Order.get(Order.size() - 1)];			
		ArrayList<ArrayList<Integer>> featureGroup = new ArrayList<ArrayList<Integer>>();			
		for(int i = 0; i < labelnum + 1; i ++){
			ArrayList<Integer> featurewithSameEnt = new ArrayList<Integer>();
			for (int j = 0; j < featureLabel.length; j++) {
				if (featureLabel[j] == i) {
					featurewithSameEnt.add(j);
				}
			}
			featureGroup.add(featurewithSameEnt);
		}	
		return featureGroup;
	}
	
	public ArrayList<Cluster> FindInitialClus() {
		ArrayList<Cluster> InitialClus = new ArrayList<Cluster>();
		for (int i = 0; i < numVertices; i++) {
			ArrayList<vertex> Ego = new ArrayList<vertex>();
			Ego.add(an.getVertex(i));
			Ego.addAll(an.getVertex(i).getNeighbor());
			Cluster Clus = new Cluster(an);
			Clus.Nodes.addAll(Ego);
			InitialClus.add(Clus);
		}
		return InitialClus;
	}	
	
	public ArrayList<Cluster> RemoveSame(ArrayList<Cluster> ClusList){
		int [] label = new int [ClusList.size()];
		for(int i = 0; i < ClusList.size()-1; i ++){
			if(label[i] == 0){
				for(int j = i + 1; j < ClusList.size(); j ++){
					if(i != j && ClusList.get(i).equals(ClusList.get(j))){
						label[j] = 1;
					}
				}
			}
		}
		ArrayList<Cluster> newClusList = new ArrayList<Cluster>();
		for(int i=0; i<label.length; i++)
			if(label[i] == 0)
				newClusList.add(ClusList.get(i));

		return newClusList;
	}
	
	public ArrayList<Cluster> InitialKClustersAutoMatNew(ArrayList<Cluster> ClusList){
		int[] Label = new int[numVertices];
		ArrayList<vertex> selectVertex = new ArrayList<vertex>();	
		ArrayList<Cluster> kCluster = new ArrayList<Cluster>();
		ArrayList<Cluster> selectList = new ArrayList<Cluster>();
		while(true){
			double ccGraph = Double.MAX_VALUE;
			Cluster selectCluster = new Cluster(an);
			for(Cluster cl : ClusList){
//				Output(cl);
				selectList.clear();
				selectList.addAll(kCluster);
				selectList.add(cl);
				double ccselectList = codingScheme.CCStructure(selectList);
				if(ccselectList < ccGraph){
					ccGraph = ccselectList;
					selectCluster = cl;
//					Output(selectCluster);
				}					
			}
//			Output(selectCluster);
//			selectList.add(selectCluster);
			Cluster chosenCluster = new Cluster(an);
			chosenCluster = selectCluster.copy();
			kCluster.add(chosenCluster);
			ClusList.remove(chosenCluster);
			/*Remove redundant cluster*/			
			for(vertex v : chosenCluster.Nodes){
				if(! selectVertex.contains(v))
					selectVertex.add(v);
			}				
			
			ArrayList<Cluster> redundantClus = new ArrayList<Cluster>();
			for(Cluster cl : ClusList){
				if(selectVertex.containsAll(cl.Nodes))
					redundantClus.add(cl);
			}
			ClusList.removeAll(redundantClus);
//			for(Cluster c : ClusList){
//				Output(c);
//			}
			for(vertex v : chosenCluster.Nodes){
				Label[v.id] = 1;
			}
			int sum = 0;
			for(int i = 0; i < Label.length; i ++){
				sum += Label[i];
			}
			if(sum == numVertices)
				break;
			System.out.println(sum+" "+numVertices);
		}
//		System.out.print(kCluster.size());
		return kCluster;		
	}
	
	public ArrayList<vertex> NodeCommonNeighbor(ArrayList<vertex> a,
			ArrayList<vertex> b) {
		ArrayList<vertex> CN = new ArrayList<vertex>();
		for (vertex v : a) {
			if (b.contains(v)) {
				CN.add(v);
			}
		}
		return CN;
	}
	
	public Pair CalculateSimStructure(Cluster a, Cluster b){
		Pair simPair = new Pair();
		simPair.member[0] = a;
		simPair.member[1] = b;
		simPair.structureCN = NodeCommonNeighbor(a.Nodes, b.Nodes);
		simPair.similarity = (double)simPair.structureCN.size()/(double)(a.Nodes.size()+ b.Nodes.size());
		return simPair;
	}
	
	public ArrayList<ArrayList<Cluster>>  FindBig(Cluster MergeC){
		ArrayList<ArrayList<Cluster>> candiBig = new ArrayList<ArrayList<Cluster>>();
		Cluster candiMerge = new Cluster(an);
		candiMerge = MergeC.copy();
		ArrayList<Integer> selectFeature = new ArrayList<Integer>();
		while(candiBig.size() < numFeature){
			ArrayList<Cluster> splitClus = new ArrayList<Cluster>();  //split cluster			
			int[][] adjMerge = new int[candiMerge.Nodes.size()][numFeature];  // feature matrix
			for(int i = 0; i < candiMerge.Nodes.size(); i ++){
				for(int j =0 ; j < numFeature; j ++){
					adjMerge[i][j] = an.featureAdj[candiMerge.Nodes.get(i).id][j];
				}
			}
			HashMap<Integer,int[]> cateOrder  = new HashMap<Integer,int[]>();
			for(int i = 0; i < numFeature; i ++){
				int[] cateNum = new int[an.numCategory[i]];
				for(vertex v : candiMerge.Nodes){
					cateNum[an.featureAdj[v.id][i]]++;
				}
				cateOrder.put(i, cateNum);
			}	
			int maxValue = Integer.MIN_VALUE;
			int[] maxID = new int[2];
		
			for(Integer f : cateOrder.keySet()){
				for(int i = 0; i <cateOrder.get(f).length; i ++){
					if(cateOrder.get(f)[i] > maxValue){
						maxValue = cateOrder.get(f)[i];
						maxID[0] = f;
						maxID[1] = i;
					}
				}
			}
			selectFeature.add(maxID[0]);
			Cluster selectClus = new Cluster(an);
			Cluster remainClus = new Cluster(an);
			for(vertex v : candiMerge.Nodes){
				if(an.featureAdj[v.id][maxID[0]]== maxID[1])
					selectClus.Nodes.add(v);				
				else
					remainClus.Nodes.add(v);
			}
			selectClus.Features.addAll(selectFeature);
			splitClus.add(selectClus);
			splitClus.add(remainClus);
			candiMerge.clear();
			candiMerge = selectClus.copy();
			candiBig.add(splitClus);
		}		
		return candiBig;
	}
	
	public ArrayList<Cluster> ModifyMerge(Pair mergeCandidate){
		double ccorg = codingScheme.CalculateCC(ClusList);		
	
		ArrayList<Cluster> remainClusList = new ArrayList<Cluster>();
		remainClusList.addAll(ClusList);
		// Remove merging clusters
		for (int m = 0; m < mergeCandidate.member.length; m ++) {
			remainClusList.remove(mergeCandidate.member[m]);
		}// cluster list remove merging clusters
		
		// Merge two clusters into one
		Cluster MergeC = new Cluster(an);
		ArrayList<vertex> NodesMergeC = new ArrayList<vertex>();	
		NodesMergeC.addAll(mergeCandidate.member[0].Nodes);
		for(vertex v : mergeCandidate.member[1].Nodes)
			if(!NodesMergeC.contains(v))
				NodesMergeC.add(v);
		MergeC.Nodes.addAll(NodesMergeC);
		ArrayList<Cluster> mergeClusList = new ArrayList<Cluster>();
		mergeClusList.addAll(remainClusList);
		mergeClusList.add(0, MergeC);		
		double ccold =  codingScheme.CalculateCC(mergeClusList);
		
		ArrayList<Cluster> modifiedClusters= new ArrayList<Cluster>();
		if(ccold > ccorg){
			modifiedClusters.add(mergeCandidate.member[0]);
			modifiedClusters.add(mergeCandidate.member[1]);
		}
		else{			
			ArrayList<ArrayList<Cluster>> candiClus =  FindBig(MergeC);
			
			for(int i = 0; i < candiClus.size(); i ++){
				
				ArrayList<Cluster> candiClusList = new ArrayList<Cluster>();
				candiClusList.addAll(remainClusList);
				candiClusList.addAll(candiClus.get(i));
				double cccandi = codingScheme.CalculateCC(candiClusList);
				if(cccandi < ccold){
					ccold = cccandi;
					modifiedClusters.clear();
					modifiedClusters.addAll(candiClus.get(i));
				}
			}			
		}
		return modifiedClusters;
	}
	
	public TreeSet<Pair> ModifyQueue(ArrayList<Cluster> StructureClusList, TreeSet<Pair> simQueue, ArrayList<Cluster> newClusters, Pair mergeCandidate){
		for(int i = 0; i < newClusters.size(); i ++){
    		for(int j = 0; j < StructureClusList.size(); j ++){
    			if(! StructureClusList.get(j).equals(newClusters.get(i))){
//    				Pair simPair = CalculateSim(newClusters.get(i), StructureClusList.get(j));
    				Pair simPair = CalculateSimStructure(newClusters.get(i), StructureClusList.get(j));
    				simQueue.add(simPair);
    			}
    		}
    	}
		return simQueue; 
	}
	
	
	public void run(){
		/** Find initial clusters, which are Egonet of each vertices */
		ArrayList<Cluster> InitialClusList = FindInitialClus();		
		InitialClusList = RemoveSame(InitialClusList);		
		/** Greedily Find top K clusters that cover all the vertices 
		 * with minimum coding cost on structure part of graph */
		ClusList = InitialKClustersAutoMatNew(InitialClusList);
		System.out.println("There are "+ClusList.size()+" initial clusters!");		
		double oldCost = codingScheme.CalculateCC(ClusList);
		System.out.println("old coding cost :  " + oldCost);

		
		/** For each pair of clusters, calculate their similarities and rank */
    	TreeSet<Pair> simQueue = new TreeSet<Pair>();
    	for(int i = 0; i < ClusList.size() - 1; i ++){
    		for(int j = i + 1; j < ClusList.size(); j ++){
//    			Pair simPair = CalculateSim(ClusList.get(i), ClusList.get(j));
    			Pair simPair = CalculateSimStructure(ClusList.get(i), ClusList.get(j));
    			simQueue.add(simPair);
    		}
    	} 
    	
    	/** Merge or Split pairs of clusters, if coding cost is decreased */
	    while(simQueue.size() > 0){
	    	/** The candidate pair of clusters is that with biggest similarity */
	    	Pair mergeCandidate = simQueue.first();
	    	
	 	    simQueue.remove(mergeCandidate);
	 	    
	 	    /** Get two clusters, continue if they are not in the clusters list any more */
	 	    Cluster clus1 = mergeCandidate.member[0];
	 	    Cluster clus2 = mergeCandidate.member[1];
	 	    if(!ClusList.contains(clus1) || !ClusList.contains(clus2))
	 	    	continue;
	 	    
	 	    /** Merge or Split the candidate pair of cluster, when coding cost is comparably small */
	 	    ArrayList<Cluster> modifiedClusters = ModifyMerge(mergeCandidate);
	 	    
	 	    if(!modifiedClusters.isEmpty()){
		 		/** Make the clusters list */
		 	    ArrayList<Cluster> newClusList = new ArrayList<Cluster>(ClusList);
		 	    newClusList.remove(mergeCandidate.member[0]);
		 	    newClusList.remove(mergeCandidate.member[1]);
	 			for(Cluster c : modifiedClusters){
	 				newClusList.add(c);
	// 				Output(c);
	 			}
	 			
	 			/** Calculate the new coding cost */
	 			newClusList = RemoveSame(newClusList);
	 			double newCost = codingScheme.CalculateCC(newClusList);
		 		if(newCost < oldCost){
		 			System.out.println(newCost);
		 			oldCost = newCost;
		 			ClusList = newClusList;
		 			simQueue = ModifyQueue(newClusList, simQueue, modifiedClusters, mergeCandidate);				
		 		}
	 	    }
	    }
	    
	    /** Postprocessing, reassign each object to each cluster, if coding cost decreases */
 		postProcessing();
		ArrayList<Cluster> emptyClus = new ArrayList<Cluster>();
		for(Cluster clus : ClusList){
			if(clus.Nodes.isEmpty())
				emptyClus.add(clus);
		}
		ClusList.removeAll(emptyClus);
 		
//		ClusList.get(2).Nodes.remove(0);
//		System.out.println(codingScheme.CalculateCC(ClusList));
	}
	
	public void postProcessing(){
		double cc = codingScheme.CalculateCC(ClusList);         // total coding cost		
		
		for(int i = 0; i < ClusList.size(); i ++){	
			Cluster modiclus = ClusList.get(i);				
//			Output(modiclus);		
//			Output(keepmodiclus);					
			for(int j = 0; j < numVertices; j ++){
				vertex v = an.getVertex(j);		
				ArrayList<Integer> modiclusFeature = new ArrayList<Integer>(modiclus.Features);	
				if(modiclus.Nodes.contains(v)){
					modiclus.Nodes.remove(v);
					modiclus.Features.clear();
//					Output(modiclus);
					ClusList = FindSubspace(ClusList);
//					Output(modiclus);					
					double ccnew = codingScheme.CalculateCC(ClusList);				
					if(ccnew >= cc){	
						modiclus.Nodes.add(v);
						modiclus.Features = modiclusFeature;
					}
					else{
						cc = ccnew;
						System.out.println(ccnew);
					}
				} else {
					modiclus.Nodes.add(v);
					modiclus.Features.clear();
//					Output(modiclus);
					ClusList = FindSubspace(ClusList);		
					double ccnew = codingScheme.CalculateCC(ClusList);
					if(ccnew >= cc){
						modiclus.Nodes.remove(v);
						modiclus.Features = modiclusFeature;
					}
					else{
						cc = ccnew;
						System.out.println(ccnew);
					}
				}				
			}			
		}
		
		ArrayList<Cluster> emptyClus = new ArrayList<Cluster>();
		for(Cluster clus : ClusList){
			if(clus.Nodes.isEmpty())
				emptyClus.add(clus);
		}
		ClusList.removeAll(emptyClus);
	}
	
	public static void main(String[] args){
		attributenetwork an = new attributenetwork();
		IO ea = new IO();	
//		ea.readFacebookLabel("data/facebook/0");
		an = ea.readNetwork("data/syn/2clusterGraph");
//		an = ea.readNetwork("data/real/1912");
//		ea.convertToGraphml(an, "data/syn/2clusterGraph.graphml");
//		ea.writeAllNodes("data/real/0", an);
//		ea.writeAllEdges("data/real/0", an);
//		an = ea.readNetwork("data/real/0");
//		an = ea.readGraphml("data/real/dblp_top.graphml");
//		Scanner sc = new Scanner(System.in);
//		System.out.println("Press any key to continue ...");
//		sc.next();

		System.out.println("Read Attributed Network");
		TestIROC ti = new TestIROC(an);
		long   start   =   System.currentTimeMillis();
		ti.run();
		
	      ArrayList<TreeSet<Integer>> IDV = new ArrayList<TreeSet<Integer>>();
	      ArrayList<TreeSet<Integer>> IDF = new ArrayList<TreeSet<Integer>>();
	      for(Cluster c : ti.ClusList){
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
	    	  System.out.println("Cluster" + numid + ": "+IDV.get(i).size());
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
			long end   =   System.currentTimeMillis(); 
		  	System.out.println("All running time:"+Long.toString(end-start)+" ms.");
		  	ea.writeLabel("0", IDV);
		  	ea.writesubspace("0", IDF);
//		ArrayList<Cluster> ClusList = ea.readClusterFromStructure("data/syn/3clusterGraph_0_40_30_40_60_40", an);
//		for(Cluster clus : ClusList)
//			docag.Output(clus);
//		ClusList = docag.FindSubspace(ClusList);
//		for(Cluster clus : ClusList)
//			docag.Output(clus);
	}

}
