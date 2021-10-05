package DocAG;

import java.util.*;
import attributenetwork.IO;
import attributenetwork.attributenetwork;
import attributenetwork.vertex;

public class DocAG {
	attributenetwork an;
	CodingScheme codingScheme;
	
	int numVertices;
	int numFeature;
	int[][] AdjNode;
	int[][] AdjFeature;
	ArrayList<Integer> FullFeature;
	ArrayList<Cluster> ClusList;
	
	public DocAG(attributenetwork an){
		this.an = an;
		this.numVertices = an.numVertices;
		this.numFeature = an.numFeature;
		FullFeature = new ArrayList<Integer>();
		for (int i = 0; i < numFeature; i++) {
			FullFeature.add(i);
		}
		AdjNode = new int[numVertices][numVertices];
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (an.exsitEdge(i, j))
					AdjNode[i][j] = 1;
			}
		}
		AdjFeature = new int[numVertices][numFeature];
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numFeature; j++) {
				AdjFeature[i][j] = an.getVertex(i).feature.get(j);
			}
		}
		ClusList = new ArrayList<Cluster>();
		this.codingScheme = new CodingScheme(an, AdjNode, AdjFeature);
	}
	
	public void Output(Cluster c){
		System.out.print("Nodes: ");
		TreeSet<Integer> nodes = new TreeSet<Integer>();
		for(vertex v : c.Nodes){
			nodes.add(v.id);
		}
		System.out.println(nodes);
		System.out.print("Feature: ");
		TreeSet<Integer> features = new TreeSet<Integer>();
		features.addAll(c.Features);
		System.out.println(features);
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
	
//	public double DataEntropy(){
//		double cctotal = 0.0;
//		double ccstructure = 0.0;
//		int numones = 0;
//		for(int i = 0; i < numVertices; i ++){
//			for(int j = 0; j < numVertices; j ++){
//				if(i != j && an.exsitEdge(i, j))
//					numones ++;
//			}
//		}
//		int numzeros = numVertices * numVertices - numones;
//		if(numones != 0)
//			ccstructure += numones * lg2((double)numVertices * numVertices / (double)numones);
//		if(numzeros != 0)
//			ccstructure += numzeros * lg2((double)numVertices * numVertices / (double)numzeros);
//		double ccstructurepara = 0.5 * lg2(numVertices * numVertices);
//		
//		double ccfeature = 0.0;
//		double ccfeaturepara = 0.0;
//		for(int i = 0; i < numFeature; i ++){
//			int[] numcate = new int[an.numCategory[i]];
//			for(int j = 0; j < numVertices; j ++){
//				numcate[an.vertexMap.get(j).feature.get(i)] ++;
//			}
//			int numnonzero = 0;
//			for(int k = 0; k < numcate.length; k ++){
//				if(numcate[k] != 0)
//					ccfeature += numcate[k] * lg2((double)numVertices /(double)numcate[k]);	
//				if(numcate[k] == 0)
//					numnonzero ++;
//			}
//			if(numnonzero == 1)
//				ccfeaturepara += 0.5*lg2(numcate.length);
//			else
//				ccfeaturepara += 0.5*(numnonzero -1) * lg2(numcate.length);		
//		}
//		cctotal = ccstructure + ccstructurepara + ccfeature + ccfeaturepara;
//		return cctotal;
//	}
	
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
	
	public ArrayList<Cluster> InitialKClusters(ArrayList<Cluster> ClusList, int K ){
		ArrayList<Cluster> kCluster = new ArrayList<Cluster>();
		ArrayList<Cluster> selectList = new ArrayList<Cluster>();
		while(kCluster.size() < K){
			double ccGraph = Double.MAX_VALUE;
			Cluster selectCluster = new Cluster(an);
			for(Cluster cl : ClusList){
//				Output(cl);
				selectList.clear();
				selectList.add(cl);
				double ccselectList = codingScheme.CCStructure(selectList);
				if(ccselectList < ccGraph){
					ccGraph = ccselectList;
					selectCluster = cl;
//					Output(selectCluster);
				}
			}
			selectList.add(selectCluster);
			kCluster.add(selectCluster);
			ClusList.remove(selectCluster);
		}
//		System.out.print(kCluster.size());
		return kCluster;		
	}
	
	public ArrayList<Cluster> InitialKClustersAutoMat(ArrayList<Cluster> ClusList){
		int[] Label = new int[numVertices];
			
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
	
	public ArrayList<Cluster> InitialRandom(ArrayList<Cluster> ClusList){
		int[] Label = new int[numVertices];		
		ArrayList<Cluster> kCluster = new ArrayList<Cluster>();
		ArrayList<Integer> kindex = new ArrayList<Integer>();
		while(true){
			Random r = new Random();
			int index = r.nextInt(ClusList.size());
			if(! kindex.contains(index)){
				kindex.add(index);
				kCluster.add(ClusList.get(index));
				for(vertex v : ClusList.get(index).Nodes){
					Label[v.id] = 1;
				}
			}			
			int sum = 0;
			for(int i = 0; i < Label.length; i ++){
				sum += Label[i];
			}
			if(sum == numVertices)
				break;
		}
		return kCluster;
	}
	
	public ArrayList<Cluster> InitialBigRange(ArrayList<Cluster> ClusList){
		int[] Label = new int[numVertices];			
		int clusSize = 0;	
		Cluster firstClus = new Cluster(an);
		for(Cluster clu : ClusList){
			if(clu.Nodes.size() > clusSize){
				clusSize = clu.Nodes.size();	
				firstClus.clear();
				firstClus = clu.copy();
			}				
		}	
		ArrayList<Cluster> selectClus = new ArrayList<Cluster>();
		ClusList.remove(firstClus);
		selectClus.add(firstClus);	
		
		while(true){
			ArrayList<Integer> range = new ArrayList<Integer>();
			for(Cluster cluS : selectClus){
				for(vertex v : cluS.Nodes){
					Label[v.id] = 1;
				}
			}			
			ArrayList<Integer> noCover = new ArrayList<Integer>();
			for(int i = 0; i < Label.length; i ++){
				if(Label[i] == 0)
					noCover.add(i);
			}
			if(noCover.size() == 0)
				break;
			for(Cluster clu : ClusList){
				int numCN = CommonNodes(clu,noCover);
				range.add(numCN);
			}
			int maxNum = 0;
			int maxID = 0;
			for(int i = 0; i < range.size(); i ++){
				if(range.get(i) > maxNum){
					maxNum = range.get(i);
					maxID = i;
				}
			}
			selectClus.add(ClusList.get(maxID));
			ClusList.remove(maxID);
			
		}
		return selectClus;
	}
	
	public int CommonNodes(Cluster A, ArrayList<Integer> B){
		int numCN = 0;
		for(vertex v1 : A.Nodes){
			for(Integer i : B){
				if(v1.id == i)
					numCN ++;
			}
		}
		return numCN;
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

	public ArrayList<Integer> FeatureCommonNeighbor(ArrayList<Integer> a,
			ArrayList<Integer> b) {
		ArrayList<Integer> CN = new ArrayList<Integer>();
		for (Integer i : a) {
			if (b.contains(i)) {
				CN.add(i);
			}
		}
		return CN;
	}
	
	public Pair CalculateSimEntropy(Cluster a, Cluster b){
		Pair simPair = new Pair();
		simPair.member[0] = a;
		simPair.member[1] = b;
		simPair.structureCN = NodeCommonNeighbor(a.Nodes, b.Nodes);
		simPair.featureCN = FeatureCommonNeighbor(a.Features, b.Features);
	
		int[][] adjCN = new int [simPair.structureCN.size()][simPair.structureCN.size()];	
		for(int i = 0; i < simPair.structureCN.size(); i ++){
			for(int j = 0 ; j < simPair.structureCN.size(); j ++){
				if(simPair.structureCN.get(i).existsEdge(simPair.structureCN.get(j)))
					adjCN[i][j] = 1;
				else
					adjCN[i][j] = 0;
			}			
		}		
		double structuresimilarity = 0.0;
		if(simPair.structureCN.size() > 0)
			structuresimilarity = AdjEntropy(adjCN);

		double featuresimilarity = 0.0;
		int[][] featurematrixCN = new int[simPair.structureCN.size()][simPair.featureCN.size()];
		for (int i =0 ; i < simPair.structureCN.size(); i ++ ) {
			for(int j = 0; j < simPair.featureCN.size(); j ++){
				featurematrixCN[i][j] = simPair.structureCN.get(i).feature.get(simPair.featureCN.get(j));
			}
		}	
		if(simPair.structureCN.size() > 0 && simPair.featureCN.size() > 0)
			featuresimilarity = (double)FeatureEntropy(featurematrixCN,simPair.featureCN);
		simPair.similarity = structuresimilarity * featuresimilarity; 
		if(structuresimilarity == 0 && featuresimilarity == 0)
			simPair.similarity = -1;
//		System.out.println(structuresimilarity+" "+featuresimilarity);
		return simPair;
	}
	
	public double FeatureEntropy(int[][] featurematrix, ArrayList<Integer> feature){
		double ent = 0.0;
		for(int j = 0; j < feature.size(); j ++){
			int f = feature.get(j);
			int[] cate = new int[an.numCategory[f]];
			for(int i = 0; i < featurematrix.length; i ++){
				cate[featurematrix[i][j]] ++;
			}
			for(int i = 0 ;i < cate.length; i++)
				if((double)cate[i]/(double)featurematrix.length != 0)
					ent += -(double)cate[i]/(double)featurematrix.length * lg2((double)cate[i]/(double)featurematrix.length);
		}
		return ent;
	}
	
	public double AdjEntropy(int[][] adj){
		int numone = 0;
		int numzero = 0;
		for(int i = 0; i < adj.length; i ++){
			for(int j = 0; j < adj.length; j ++){
				if(adj[i][j] == 1)
					numone ++;
				else
					numzero ++;
			}
		}
		int sum = numone + numzero;
		double p1 = (double)numone/(double)sum;
		double p0 = (double)numzero/(double)sum;
		double ent = 0.0;
		if(p1 != 0 && p0 != 0)
			ent = - p1 * lg2(p1) - p0 * lg2(p0);
		return ent;
	}
	
	public Pair CalculateSim(Cluster a, Cluster b){
		Pair simPair = new Pair();
		simPair.member[0] = a;
		simPair.member[1] = b;
		simPair.structureCN = NodeCommonNeighbor(a.Nodes, b.Nodes);
		simPair.featureCN = FeatureCommonNeighbor(a.Features, b.Features);
		
		HashSet<vertex> combinedNodes = new HashSet<vertex>();
		combinedNodes.addAll(a.Nodes);
		combinedNodes.addAll(b.Nodes);
		double structuresimilarity = (double) simPair.structureCN.size()/ (double) (a.Nodes.size() + b.Nodes.size());

		double featuresimilarity = 0.0;
		for (Integer f : simPair.featureCN) {
			int[] featurecat = new int[an.numCategory[f]];
			for(vertex v : combinedNodes)
//			for(vertex v : simPair.structureCN)
				featurecat[v.feature.get(f)] ++;
			int maxcat = 0;
			for(int jj = 0; jj < featurecat.length; jj ++){
				if(featurecat[jj] > maxcat)
					maxcat = featurecat[jj];
			}
			featuresimilarity += (double)maxcat/combinedNodes.size();
//			featuresimilarity += (double)maxcat/simPair.structureCN.size();
		}
		if(simPair.featureCN.size() != 0)
			featuresimilarity = (double) featuresimilarity/(double) simPair.featureCN.size();		
		simPair.similarity = structuresimilarity * featuresimilarity; 
//		System.out.println(simPair.similarity);
		return simPair;
	}
	
	public TreeSet<Pair> ModifyQueue(ArrayList<Cluster> StructureClusList, TreeSet<Pair> simQueue, ArrayList<Cluster> newClusters, Pair mergeCandidate){
		for(int i = 0; i < newClusters.size(); i ++){
    		for(int j = 0; j < StructureClusList.size(); j ++){
    			if(! StructureClusList.get(j).equals(newClusters.get(i))){
//    				Pair simPair = CalculateSim(newClusters.get(i), StructureClusList.get(j));
    				Pair simPair = CalculateSimEntropy(newClusters.get(i), StructureClusList.get(j));
    				simQueue.add(simPair);
    			}
    		}
    	}
		return simQueue; 
	}
	
	public ArrayList<Cluster> ModifyClusterwithoutPP(Pair mergeCandidate){
		double ccorg = codingScheme.CalculateCC(ClusList);
		
		ArrayList<Cluster> tempClusList = new ArrayList<Cluster>();
		ArrayList<Cluster> remainClusList = new ArrayList<Cluster>();
		tempClusList.addAll(ClusList);
		// Remove merging clusters
		for (int m = 0; m < mergeCandidate.member.length; m ++) {
			tempClusList.remove(mergeCandidate.member[m]);
		}
		remainClusList.addAll(tempClusList); // cluster list remove merging clusters
		
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
//		mergeClusList = FindSubspace(mergeClusList);
		FindClusterSubspace(MergeC, mergeClusList);
		double ccinitial = codingScheme.ccClusterStructure + codingScheme.ccClusterFeature 
		- codingScheme.CCCluster(mergeCandidate.member[0]) - codingScheme.CCCluster(mergeCandidate.member[1]);
		double ccold =  ccinitial + codingScheme.CCCluster(MergeC) + codingScheme.CCUnCluster(mergeClusList);
//		double ccold = codingScheme.CalculateCC(mergeClusList);
//		System.out.println(ccold+" "+ccold1);
		
		//Split
//		ArrayList<vertex> NodesCommon = new ArrayList<vertex>();
//		ArrayList<vertex> NodesSplit = new ArrayList<vertex>();
//		for(vertex v : NodesMergeC){
//			if(mergeCandidate.member[0].Nodes.contains(v)&& mergeCandidate.member[1].Nodes.contains(v))
//				NodesCommon.add(v);
//			else
//				NodesSplit.add(v);
//		}		
		ArrayList<vertex> moveNodes = new ArrayList<vertex>();
		moveNodes.addAll(MergeC.Nodes);
		Cluster RemainC = new Cluster(an);
//		if(!cnNodes.isEmpty()){
		mergeClusList.add(0, RemainC);
//		double ccnew = 0.0;
		boolean stop = true;
		while(stop){
			stop = false;				
			for(vertex v : moveNodes){		
				double ccnew = ccinitial;
				ArrayList<Integer> MergeCFeature = new ArrayList<Integer>(MergeC.Features);	
				ArrayList<Integer> RemainCFeature = new ArrayList<Integer>(RemainC.Features);	
				MergeC.Nodes.remove(v);
				RemainC.Nodes.add(v);			
				FindClusterSubspace(MergeC, mergeClusList);
				FindClusterSubspace(RemainC, mergeClusList);
//				mergeClusList = FindSubspace(mergeClusList);
				
				ccnew += codingScheme.CCCluster(MergeC);
				ccnew += codingScheme.CCCluster(RemainC);
				ccnew += codingScheme.CCUnCluster(mergeClusList);
//				ccnew = codingScheme.CalculateCC(mergeClusList);
				if(ccnew < ccold){
					ccold = ccnew;					
					stop = true;			
				}
				else{
					MergeC.Nodes.add(v);
					RemainC.Nodes.remove(v);
					MergeC.Features.clear();
					MergeC.Features.addAll(MergeCFeature);
					RemainC.Features.clear();
					RemainC.Features.addAll(RemainCFeature);					
				}
			}
			for(vertex v : RemainC.Nodes){
				if(moveNodes.contains(v))
					moveNodes.remove(v);
			}
//			mergeClusList = FindSubspace(mergeClusList);
		}
//		}
		
		ArrayList<Cluster> modifiedClusters= new ArrayList<Cluster>();
		if(ccold > ccorg)
			return modifiedClusters;
		
		if(!MergeC.Nodes.isEmpty())
			modifiedClusters.add(MergeC);
		if(!RemainC.Nodes.isEmpty())
			modifiedClusters.add(RemainC);
//		Output(MergeC);
//		Output(RemainC);
		return modifiedClusters;
	}
	
	public ArrayList<Cluster> ModifyClusterPair(Pair mergeCandidate){
		HashMap<ArrayList<Cluster>, Double> Cases = new HashMap<ArrayList<Cluster>, Double>();
		ArrayList<Cluster> tempClusList = new ArrayList<Cluster>();
		ArrayList<Cluster> remainClusList = new ArrayList<Cluster>();
		tempClusList.addAll(ClusList);		
				
		// Remove merging clusters
		for (int m = 0; m < mergeCandidate.member.length; m ++) {
			tempClusList.remove(mergeCandidate.member[m]);
		}
		remainClusList.addAll(tempClusList); // cluster list remove merging clusters
		
		
		// Case1	
		ArrayList<Cluster> case1ClusList = new ArrayList<Cluster>();
		case1ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters1 = new ArrayList<Cluster>(); // Produced new clusters
		Cluster C1 = new Cluster(an);
		ArrayList<vertex> NodesC1 = new ArrayList<vertex>();
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			Cluster Clus = mergeCandidate.member[i];
			for (vertex v : Clus.Nodes) {
				if (!NodesC1.contains(v)) {
					NodesC1.add(v);
				}
			}
		}
		C1.Nodes.addAll(NodesC1);	
		newClusters1.add(C1);
		if(C1.Nodes.size()>0){
			case1ClusList.addAll(newClusters1);
			case1ClusList = FindSubspace(case1ClusList);		
			double cccase1 = codingScheme.CalculateCC(case1ClusList);
			Cases.put(newClusters1, cccase1);
		}
		
		//Case 2
		ArrayList<Cluster> case2ClusList = new ArrayList<Cluster>();
		case2ClusList.addAll(remainClusList);
		ArrayList<vertex> CNNodes = mergeCandidate.structureCN;			
		ArrayList<Cluster> newClusters2 = new ArrayList<Cluster>(); // Produced new clusters	
		Cluster C21 = new Cluster(an);
		Cluster C22 = new Cluster(an);
		ArrayList<vertex> NodesC21 = new ArrayList<vertex>();
		NodesC21.addAll(mergeCandidate.member[0].Nodes);
		ArrayList<vertex> NodesC22 = new ArrayList<vertex>();		
		for(vertex v : mergeCandidate.member[1].Nodes){
			if(! NodesC22.contains(v) && ! CNNodes.contains(v))
				NodesC22.add(v);
		}
		C21.Nodes.addAll(NodesC21);
		C22.Nodes.addAll(NodesC22);	
		newClusters2.add(C21);
		newClusters2.add(C22);	
		if(C21.Nodes.size() > 0 && C22.Nodes.size() > 0){
			case2ClusList.addAll(newClusters2);
			case2ClusList = FindSubspace(case2ClusList);		
			double cccase2 = codingScheme.CalculateCC(case2ClusList);
			Cases.put(newClusters2, cccase2);
		}
		
		//Case3
		ArrayList<Cluster> case3ClusList = new ArrayList<Cluster>();
		case3ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters3 = new ArrayList<Cluster>(); // Produced new clusters	
		Cluster C31 = new Cluster(an);
		Cluster C32 = new Cluster(an);
		ArrayList<vertex> NodesC31 = new ArrayList<vertex>();
		NodesC31.addAll(mergeCandidate.member[1].Nodes);
		ArrayList<vertex> NodesC32 = new ArrayList<vertex>();		
		for(vertex v : mergeCandidate.member[0].Nodes){
			if(! NodesC32.contains(v) && ! CNNodes.contains(v))
				NodesC32.add(v);
		}
		C31.Nodes.addAll(NodesC31);
		C32.Nodes.addAll(NodesC32);	
		newClusters3.add(C31);
		newClusters3.add(C32);
		if(C31.Nodes.size() > 0 && C32.Nodes.size() > 0){
			case3ClusList.addAll(newClusters3);
			case3ClusList = FindSubspace(case3ClusList);		
			double cccase3 = codingScheme.CalculateCC(case3ClusList);
			Cases.put(newClusters3, cccase3);	
		}
		
		// Case4: Both clusters are redundant
		ArrayList<Cluster> case4ClusList = new ArrayList<Cluster>();
		case4ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters4 = new ArrayList<Cluster>(); // Produced new clusters
		Cluster C41 = new Cluster(an);
		Cluster C42 = new Cluster(an);
		Cluster C43 = new Cluster(an);
		ArrayList<vertex> NodesC41 = new ArrayList<vertex>();
		for(vertex v : mergeCandidate.member[1].Nodes){
			if(! NodesC41.contains(v) && ! CNNodes.contains(v))
				NodesC41.add(v);
		}
		ArrayList<vertex> NodesC42 = new ArrayList<vertex>();		
		for(vertex v : mergeCandidate.member[0].Nodes){
			if(! NodesC42.contains(v) && ! CNNodes.contains(v))
				NodesC42.add(v);
		}		
		C41.Nodes.addAll(NodesC41);		
		C42.Nodes.addAll(NodesC42);
		C43.Nodes.addAll(CNNodes);
		newClusters4.add(C41);
		newClusters4.add(C42);
		newClusters4.add(C43);
		if(C41.Nodes.size() > 0 && C42.Nodes.size() > 0 && C43.Nodes.size() > 0){
			case4ClusList.addAll(newClusters4);
			case4ClusList = FindSubspace(case4ClusList);		
			double cccase4 = codingScheme.CalculateCC(case4ClusList);
			Cases.put(newClusters4, cccase4);	
		}	
		
		double ccfinalmin = Double.MAX_VALUE;
		ArrayList<Cluster> minCluster = new ArrayList<Cluster>();
		int num = 0;
		for(ArrayList<Cluster> n : Cases.keySet()){
//			for(Cluster c : n){
//				Output(c);
//			}
			if(Cases.get(n) < ccfinalmin){
				ccfinalmin = Cases.get(n);
				minCluster.clear();
				minCluster.addAll(n);
				num ++;
			}
		}	
		System.out.println("Choose Cases:" + num);
		
		return minCluster;	
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
	
	public void run(){
		/** Find initial clusters, which are Egonet of each vertices */
		ArrayList<Cluster> InitialClusList = FindInitialClus();		
		InitialClusList = RemoveSame(InitialClusList);
//		for(Cluster c : InitialClusList)
//			Output(c);
		
		/** Greedily Find top K clusters that cover all the vertices 
		 * with minimum coding cost on structure part of graph */
//		ClusList = InitialKClusters(InitialClusList, 40);
//		ClusList = InitialKClustersAutoMat(InitialClusList);
		ClusList = InitialKClustersAutoMatNew(InitialClusList);
//		ClusL.ist = InitialRandom(InitialClusList);
//		ClusList = InitialBigRange(InitialClusList);
//		ClusList = InitialClusList;
		System.out.println("There are "+ClusList.size()+" initial clusters!");
		
		/** Find subspaces for these K clusters */
		ClusList = FindSubspace(ClusList);
		double oldCost = codingScheme.CalculateCC(ClusList);
		System.out.println(oldCost);
		System.out.println("test!");
		
		/** For each pair of clusters, calculate their similarities and rank */
    	TreeSet<Pair> simQueue = new TreeSet<Pair>();
    	for(int i = 0; i < ClusList.size() - 1; i ++){
    		for(int j = i + 1; j < ClusList.size(); j ++){
//    			Pair simPair = CalculateSim(ClusList.get(i), ClusList.get(j));
    			Pair simPair = CalculateSimEntropy(ClusList.get(i), ClusList.get(j));
    			simQueue.add(simPair);
    		}
    	} 
    	
    	/** Merge or Split pairs of clusters, if coding cost is decreased */
	    while(simQueue.size() > 0){
	    	/** The candidate pair of clusters is that with biggest similarity */
	    	Pair mergeCandidate = simQueue.first();
	 	    simQueue.remove(mergeCandidate);
//	 	    Output(mergeCandidate.member[0]);
//	 	    Output(mergeCandidate.member[1]);
	 	    
	 	    /** Get two clusters, continue if they are not in the clusters list any more */
	 	    Cluster clus1 = mergeCandidate.member[0];
	 	    Cluster clus2 = mergeCandidate.member[1];
	 	    if(!ClusList.contains(clus1) || !ClusList.contains(clus2))
	 	    	continue;
	 	    
	 	    /** Merge or Split the candidate pair of cluster, when coding cost is comparably small */
//	 	    ArrayList<Cluster> modifiedClusters = ModifyClusterPair(mergeCandidate);
	 	    ArrayList<Cluster> modifiedClusters = ModifyClusterwithoutPP(mergeCandidate);
//	 	    for(Cluster c : modifiedClusters)
//	 	    	Output(c);
	 	    
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
// 		postProcessing();
		ArrayList<Cluster> emptyClus = new ArrayList<Cluster>();
		for(Cluster clus : ClusList){
			if(clus.Nodes.isEmpty())
				emptyClus.add(clus);
		}
		ClusList.removeAll(emptyClus);
 		
//		ClusList.get(2).Nodes.remove(0);
//		System.out.println(codingScheme.CalculateCC(ClusList));
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
		DocAG docag = new DocAG(an);
		long   start   =   System.currentTimeMillis();
		docag.run();
		
	      ArrayList<TreeSet<Integer>> IDV = new ArrayList<TreeSet<Integer>>();
	      ArrayList<TreeSet<Integer>> IDF = new ArrayList<TreeSet<Integer>>();
	      for(Cluster c : docag.ClusList){
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
