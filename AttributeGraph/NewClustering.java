package AttributeGraph;

import java.util.*;
import attributenetwork.edge;
import attributenetwork.attributenetwork;
import attributenetwork.vertex;


public class NewClustering {
	attributenetwork an;
	int numVertices;
	int numFeature;
	ArrayList<Integer> FullFeature;
	int[][] AdjNode;
	int[][] AdjFeature;
	
	NewCodingCost cc;
	static int minVertexInClus = 2;

	public NewClustering(attributenetwork an) {
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
		
		cc = new NewCodingCost(an);
	}
	
	public double DataEntropy(attributenetwork an){
		double cctotal = 0.0;
		double ccstructure = 0.0;
		int numones = 0;
		for(int i = 0; i < numVertices; i ++){
			for(int j = 0; j < numVertices; j ++){
				if(i != j && an.exsitEdge(i, j))
					numones ++;
			}
		}
		int numzeros = numVertices * numVertices - numones;
		if(numones != 0)
			ccstructure += numones * lg2((double)numVertices * numVertices / (double)numones);
		if(numzeros != 0)
			ccstructure += numzeros * lg2((double)numVertices * numVertices / (double)numzeros);
		double ccstructurepara = 0.5 * lg2(numVertices * numVertices);
		
		double ccfeature = 0.0;
		double ccfeaturepara = 0.0;
		for(int i = 0; i < numFeature; i ++){
			int[] numcate = new int[an.numCategory[i]];
			for(int j = 0; j < numVertices; j ++){
				numcate[an.vertexMap.get(j).feature.get(i)] ++;
			}
			int numnonzero = 0;
			for(int k = 0; k < numcate.length; k ++){
				if(numcate[k] != 0)
					ccfeature += numcate[k] * lg2((double)numVertices /(double)numcate[k]);	
				if(numcate[k] == 0)
					numnonzero ++;
			}
			if(numnonzero == 1)
				ccfeaturepara += 0.5*lg2(numcate.length);
			else
				ccfeaturepara += 0.5*(numnonzero -1) * lg2(numcate.length);		
		}
		cctotal = ccstructure + ccstructurepara + ccfeature + ccfeaturepara;
		return cctotal;
	}

	public ArrayList<Cluster> FindInitialClus() {
		ArrayList<Cluster> InitialClus = new ArrayList<Cluster>();
		for (int i = 0; i < numVertices; i++) {
			ArrayList<vertex> Ego = new ArrayList<vertex>();
			Ego.add(an.getVertex(i));
			Ego.addAll(an.getVertex(i).getNeighbor());

			// for(int t = 0; t < Ego.size(); t ++){
			// System.out.println(Ego.get(t).id);
			// }

//			ArrayList<Integer> Features = new ArrayList<Integer>();
//			for (int j = 0; j < numFeature; j++) {
//				Features.add(j);
//			}
			Cluster Clus = new Cluster(an);
			Clus.Nodes.addAll(Ego);
//			Clus.Features.addAll(Features);
			InitialClus.add(Clus);
		}
		return InitialClus;
	}		

	public ArrayList<Cluster> RemoveRedundency(ArrayList<Cluster> ClusList) {
		ArrayList<Cluster> removeCluster = new ArrayList<Cluster>();
		// ArrayList<Integer> removeID = new ArrayList<Integer>();
		for (int i = 0; i < ClusList.size() - 1; i++) {
			for (int j = i + 1; j < ClusList.size(); j++) {
				Cluster Clus1 = ClusList.get(i);
				Cluster Clus2 = ClusList.get(j);
				if (Clus1.Nodes.containsAll(Clus2.Nodes)) {
					removeCluster.add(Clus2);
					// removeID.add(j);
				}
				if (Clus2.Nodes.containsAll(Clus1.Nodes)) {
					removeCluster.add(Clus1);
					// removeID.add(i);
				}
			}
		}
		for (Cluster clu : removeCluster) {
			ClusList.remove(clu);
		}
		// for(Integer ID : removeID){
		// System.out.println("ID:  " + ID);
		// }
		return ClusList;
	}

	public ArrayList<Cluster> InitialKClusters(ArrayList<Cluster> ClusList, int K ){
		CodingCost cc = new CodingCost(an);		
		ArrayList<Cluster> kCluster = new ArrayList<Cluster>();
		ArrayList<Cluster> selectList = new ArrayList<Cluster>();
		while(kCluster.size() < K){
			double ccGraph = Double.MAX_VALUE;
			Cluster selectCluster = new Cluster(an);
			for(Cluster cl : ClusList){
//				Output(cl);
				selectList.clear();
				selectList.add(cl);
				double ccselectList = cc.CCStructure(selectList);
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

	public ArrayList<Pair> SimClusterPair(ArrayList<Cluster> ClusList) {
		ArrayList<Pair> SimList = new ArrayList<Pair>();
		for (int i = 0; i < ClusList.size() - 1; i++) {
			for (int j = i + 1; j < ClusList.size(); j++) {
				Pair SimPair = new Pair();
				if (i != j) {
					SimPair.member[0] = ClusList.get(i);
					SimPair.member[1] = ClusList.get(j);
					SimPair.structureCN = NodeCommonNeighbor(
							ClusList.get(i).Nodes, ClusList.get(j).Nodes);
					SimPair.featureCN = FeatureCommonNeighbor(
							ClusList.get(i).Features, ClusList.get(j).Features);
					ArrayList<vertex> CombinedNodes = new ArrayList<vertex>();
					CombinedNodes.addAll(ClusList.get(i).Nodes);
					for (vertex v : ClusList.get(j).Nodes) {
						if (!CombinedNodes.contains(v))
							CombinedNodes.add(v);
					}
					double structuresimilarity = (double) SimPair.structureCN.size()/ (double) (ClusList.get(i).Nodes.size() + ClusList.get(j).Nodes.size());

					double featuresimilarity = 0.0;
					for (Integer f : SimPair.featureCN) {
						int[] CommonFeature = new int[CombinedNodes.size()];
						for (int m = 0; m < CombinedNodes.size(); m++) {
							CommonFeature[m] = CombinedNodes.get(m).feature.get(f);
						}						
						int sum = CommonFeature.length;
						int[] featurecat = new int[an.numCategory[f]];
						for (int ii = 0; ii < sum; ii++) {
							featurecat[CommonFeature[ii]]++;
						}
						int maxcat = 0;
						for(int jj = 0; jj < featurecat.length; jj ++){
							if(featurecat[jj] > maxcat)
								maxcat = featurecat[jj];
						}
						featuresimilarity += (double)maxcat/(double)sum;
					}
					if(SimPair.featureCN.size() != 0)
						featuresimilarity = (double) featuresimilarity/(double) SimPair.featureCN.size();					
					SimPair.similarity = structuresimilarity * featuresimilarity; // The smaller the better
					// System.out.println("Similarity" + SimEgoPair.similarity);
				}
				SimList.add(SimPair);
			}
		}
		TreeSet<Pair> tmp = new TreeSet<Pair>();
		tmp.addAll(SimList);
		for(Pair p : tmp){
			System.out.println(p.similarity);
		}
		return SimList;
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
//		for (vertex v : b.Nodes) {
//			if (!combinedNodes.contains(v))
//				combinedNodes.add(v);
//		}
		double structuresimilarity = (double) simPair.structureCN.size()/ (double) (a.Nodes.size() + b.Nodes.size());

		double featuresimilarity = 0.0;
		for (Integer f : simPair.featureCN) {
			int[] commonFeature = new int[combinedNodes.size()];
//			for (int m = 0; m < combinedNodes.size(); m++) {
//				commonFeature[m] = combinedNodes.get(m).feature.get(f);
//			}		
			int m = 0;
			for(vertex v : combinedNodes)
				commonFeature[m++] = v.feature.get(f);
			int sum = commonFeature.length;
			int[] featurecat = new int[an.numCategory[f]];
			for (int ii = 0; ii < sum; ii++) {
				featurecat[commonFeature[ii]]++;
			}
			int maxcat = 0;
			for(int jj = 0; jj < featurecat.length; jj ++){
				if(featurecat[jj] > maxcat)
					maxcat = featurecat[jj];
			}
			featuresimilarity += (double)maxcat/(double)sum;
		}
		if(simPair.featureCN.size() != 0)
			featuresimilarity = (double) featuresimilarity/(double) simPair.featureCN.size();		
		simPair.similarity = structuresimilarity * featuresimilarity; 
		return simPair;
	}

	public Pair FindMostSimPair(ArrayList<Cluster> ClusList) {
		ArrayList<Pair> SimilarityEgo = SimClusterPair(ClusList);
		double minSim = Double.MAX_VALUE;
		int minid = 0;
		for (int i = 0; i < SimilarityEgo.size(); i++) {
			if (SimilarityEgo.get(i).similarity < minSim) {
				minSim = SimilarityEgo.get(i).similarity;
				minid = i;
			}
		}
		Pair MostSimPair = SimilarityEgo.get(minid);
		return MostSimPair;
	}
	
	public ArrayList<Cluster> ModifyClusternew(ArrayList<Cluster> ClusList, Pair mergeCandidate){
		CodingCost cc = new CodingCost(an);
		ArrayList<Cluster> newCluster = new ArrayList<Cluster>();	
		double cc_compare = cc.CalculateCC(ClusList);    // Compared coding cost
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
			case1ClusList = NewFindSubspace(newClusters1, case1ClusList);		
			double cccase1 = cc.CalculateCC(case1ClusList);
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
			case2ClusList = NewFindSubspace(newClusters2, case2ClusList);		
			double cccase2 = cc.CalculateCC(case2ClusList);
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
			case3ClusList = NewFindSubspace(newClusters3, case3ClusList);		
			double cccase3 = cc.CalculateCC(case3ClusList);
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
			case4ClusList = NewFindSubspace(newClusters4, case4ClusList);		
			double cccase4 = cc.CalculateCC(case4ClusList);
			Cases.put(newClusters4, cccase4);	
		}	
		
		double ccfinalmin = Double.MAX_VALUE;
		ArrayList<Cluster> minCluster = new ArrayList<Cluster>();
		for(ArrayList<Cluster> n : Cases.keySet()){
//			for(Cluster c : n){
//				Output(c);
//			}
			if(Cases.get(n) < ccfinalmin){
				ccfinalmin = Cases.get(n);
				minCluster.clear();
				minCluster.addAll(n);
			}
		}		
		if(ccfinalmin < cc_compare)			
			return minCluster;	
		else
			return newCluster;
	}
	
	public ArrayList<Cluster> ModifyCluster(ArrayList<Cluster> ClusList, Pair mergeCandidate){		
		ArrayList<Cluster> newCluster = new ArrayList<Cluster>();		
		HashMap<ArrayList<Cluster>, Double> Cases = new HashMap<ArrayList<Cluster>, Double>();	// Four cases for merging and splitting	
		ArrayList<Cluster> tempClusList = new ArrayList<Cluster>();
		ArrayList<Cluster> remainClusList = new ArrayList<Cluster>();
		tempClusList.addAll(ClusList);		
		double cc_compare = cc.CalculateCC(ClusList);    // Compared coding cost		
		
		// Remove merging clusters
		for (int m = 0; m < mergeCandidate.member.length; m ++) {
			tempClusList.remove(mergeCandidate.member[m]);
		}
		remainClusList.addAll(tempClusList); // cluster list remove merging clusters	
	
		
		// Case1	
		ArrayList<Cluster> newClusters1 = new ArrayList<Cluster>(); // Produced new clusters
		Cluster C1 = new Cluster(an);
		HashSet<vertex> NodesC1 = new HashSet<vertex>();
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			Cluster Clus = mergeCandidate.member[i];
			for (vertex v : Clus.Nodes) {
				NodesC1.add(v);
			}
		}
		C1.Nodes.addAll(NodesC1);		
		ArrayList<ArrayList<Integer>> featureOrderC1  = FindClusterFeatureOrder(C1);
		ArrayList<Integer> C1Feature = new ArrayList<Integer>();
		HashMap<Cluster, Double> candidateC1 = new HashMap<Cluster, Double>();
		for(int i = 0; i < featureOrderC1.size(); i ++){
			Cluster C1S = new Cluster(an);
			C1S.Nodes.addAll(C1.Nodes);
			C1Feature.addAll(featureOrderC1.get(i));			
			C1S.Features.addAll(C1Feature);
			tempClusList.add(C1S);
			double cctemp = cc.CalculateCC(tempClusList);
			candidateC1.put(C1S, cctemp);
			tempClusList.remove(C1S);
		}
		double mincc = Double.MAX_VALUE;
		Cluster minC = new Cluster(an);
		for(Cluster c : candidateC1.keySet()){
			System.out.println(c.Nodes.size());
			if(candidateC1.get(c) < mincc){
				mincc = candidateC1.get(c);
				minC.clear();
				minC.Nodes.addAll(c.Nodes);
				minC.Features.addAll(c.Features);
				System.out.println(minC.Nodes.size());
			}
		}
		System.out.println(minC.Nodes.size());
		newClusters1.add(minC);		
		double cccase1 = mincc;
		Cases.put(newClusters1, cccase1);		
		for(ArrayList<Cluster> nc : Cases.keySet()){
			System.out.println(nc.get(0).Nodes.size());
			System.out.println(Cases.get(nc));
		}
		
		
		//Case 2
		ArrayList<vertex> CNNodes = mergeCandidate.structureCN;			
		ArrayList<Cluster> newClusters2 = new ArrayList<Cluster>(); // Produced new clusters
		tempClusList.clear();
		tempClusList.addAll(remainClusList);
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
		if(C21.Nodes.size() != 0 && C22.Nodes.size() != 0){
			ArrayList<ArrayList<Integer>> featureOrderC21  = FindClusterFeatureOrder(C21);
			ArrayList<Integer> C21Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC21 = new HashMap<Cluster, Double>();
			tempClusList.add(C22);
			for(int i = 0; i < featureOrderC21.size(); i ++){
				Cluster C21S = new Cluster(an);
				C21Feature.addAll(featureOrderC21.get(i));				
				C21S.Features.addAll(C21Feature);
				tempClusList.add(C21S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC21.put(C21S, cctemp);
				tempClusList.remove(C21S);
			}
			double mincc21 = Double.MAX_VALUE;
			Cluster minC21 = new Cluster(an);
			for(Cluster c : candidateC21.keySet()){
				if(candidateC21.get(c) < mincc21){
					mincc21 = candidateC21.get(c);
					minC21.Nodes.addAll(c.Nodes);
					minC21.Features.addAll(c.Features);
				}
			}
			newClusters2.add(minC21);
			tempClusList.clear();
			tempClusList.addAll(remainClusList);
			tempClusList.add(minC21);
			ArrayList<ArrayList<Integer>> featureOrderC22  = FindClusterFeatureOrder(C22);
			ArrayList<Integer> C22Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC22 = new HashMap<Cluster, Double>();		
			for(int i = 0; i < featureOrderC22.size(); i ++){
				Cluster C22S = new Cluster(an);
				C22Feature.addAll(featureOrderC22.get(i));			
				C22S.Features.addAll(C22Feature);
				tempClusList.add(C22S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC22.put(C22S, cctemp);
				tempClusList.remove(C22S);
			}
			double mincc22 = Double.MAX_VALUE;
			Cluster minC22 = new Cluster(an);
			for(Cluster c : candidateC22.keySet()){
				if(candidateC22.get(c) < mincc22){
					mincc22 = candidateC22.get(c);
					minC22.Nodes.addAll(c.Nodes);
					minC22.Features.addAll(c.Features);
				}
			}
			newClusters2.add(minC22);		
			double cccase2 = mincc22;
			Cases.put(newClusters2, cccase2);
		}
		
		//Case3
		ArrayList<Cluster> newClusters3 = new ArrayList<Cluster>(); // Produced new clusters
		tempClusList.clear();
		tempClusList.addAll(remainClusList);
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
		if(C21.Nodes.size() != 0 && C22.Nodes.size() != 0){
			ArrayList<ArrayList<Integer>> featureOrderC31  = FindClusterFeatureOrder(C31);
			ArrayList<Integer> C31Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC31 = new HashMap<Cluster, Double>();
			tempClusList.add(C32);
			for(int i = 0; i < featureOrderC31.size(); i ++){
				Cluster C31S = new Cluster(an);
				C31Feature.addAll(featureOrderC31.get(i));			
				C31S.Features.addAll(C31Feature);
				tempClusList.add(C31S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC31.put(C31S, cctemp);
				tempClusList.remove(C31S);
			}
			double mincc31 = Double.MAX_VALUE;
			Cluster minC31 = new Cluster(an);
			for(Cluster c : candidateC31.keySet()){
				if(candidateC31.get(c) < mincc31){
					mincc31 = candidateC31.get(c);
					minC31.Nodes.addAll(c.Nodes);
					minC31.Features.addAll(c.Features);
				}
			}
			newClusters3.add(minC31);
			tempClusList.clear();
			tempClusList.addAll(remainClusList);
			tempClusList.add(minC31);
			ArrayList<ArrayList<Integer>> featureOrderC32  = FindClusterFeatureOrder(C32);
			ArrayList<Integer> C32Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC32 = new HashMap<Cluster, Double>();		
			for(int i = 0; i < featureOrderC32.size(); i ++){
				Cluster C32S = new Cluster(an);
				C32Feature.addAll(featureOrderC32.get(i));				
				C32S.Features.addAll(C32Feature);
				tempClusList.add(C32S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC32.put(C32S, cctemp);
				tempClusList.remove(C32S);
			}
			double mincc32 = Double.MAX_VALUE;
			Cluster minC32 = new Cluster(an);
			for(Cluster c : candidateC32.keySet()){
				if(candidateC32.get(c) < mincc32){
					mincc32 = candidateC32.get(c);
					minC32.Nodes.addAll(c.Nodes);
					minC32.Features.addAll(c.Features);
				}
			}
			newClusters3.add(minC32);		
			double cccase3 = mincc32;		
			Cases.put(newClusters3, cccase3);
		}
		// Case4: Both clusters are redundant
		ArrayList<Cluster> newClusters4 = new ArrayList<Cluster>(); // Produced new clusters
		tempClusList.clear();
		tempClusList.addAll(remainClusList);
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
		if(C41.Nodes.size() != 0 && C42.Nodes.size() != 0 && C43.Nodes.size() != 0){
			ArrayList<ArrayList<Integer>> featureOrderC41  = FindClusterFeatureOrder(C41);
			ArrayList<Integer> C41Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC41 = new HashMap<Cluster, Double>();
			tempClusList.add(C42);
			tempClusList.add(C43);
			for(int i = 0; i < featureOrderC41.size(); i ++){
				Cluster C41S = new Cluster(an);
				C41Feature.addAll(featureOrderC41.get(i));			
				C41S.Features.addAll(C41Feature);
				tempClusList.add(C41S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC41.put(C41S, cctemp);
				tempClusList.remove(C41S);
			}
			double mincc41 = Double.MAX_VALUE;
			Cluster minC41 = new Cluster(an);
			for(Cluster c : candidateC41.keySet()){
				if(candidateC41.get(c) < mincc41){
					mincc41 = candidateC41.get(c);
					minC41.Nodes.addAll(c.Nodes);
					minC41.Features.addAll(c.Features);
				}
			}
			newClusters4.add(minC41);
			tempClusList.clear();
			tempClusList.addAll(remainClusList);
			tempClusList.add(minC41);
			tempClusList.add(C43);
			ArrayList<ArrayList<Integer>> featureOrderC42  = FindClusterFeatureOrder(C42);
			ArrayList<Integer> C42Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC42 = new HashMap<Cluster, Double>();		
			for(int i = 0; i < featureOrderC42.size(); i ++){
				Cluster C42S = new Cluster(an);
				C42Feature.addAll(featureOrderC42.get(i));		
				C42.Features.addAll(C42Feature);
				tempClusList.add(C42S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC42.put(C42S, cctemp);
				tempClusList.remove(C42S);
			}
			double mincc42 = Double.MAX_VALUE;
			Cluster minC42 = new Cluster(an);
			for(Cluster c : candidateC42.keySet()){
				if(candidateC42.get(c) < mincc42){
					mincc42 = candidateC42.get(c);
					minC42.Nodes.addAll(c.Nodes);
					minC42.Features.addAll(c.Features);
				}
			}
			newClusters4.add(minC42);
			tempClusList.clear();
			tempClusList.addAll(remainClusList);
			tempClusList.add(minC41);
			tempClusList.add(minC42);
			ArrayList<ArrayList<Integer>> featureOrderC43  = FindClusterFeatureOrder(C43);
			ArrayList<Integer> C43Feature = new ArrayList<Integer>();
			HashMap<Cluster, Double> candidateC43 = new HashMap<Cluster, Double>();		
			for(int i = 0; i < featureOrderC43.size(); i ++){
				Cluster C43S = new Cluster(an);
				C43Feature.addAll(featureOrderC43.get(i));			
				C43S.Features.addAll(C43Feature);
				tempClusList.add(C43S);
				double cctemp = cc.CalculateCC(tempClusList);
				candidateC43.put(C43S, cctemp);
				tempClusList.remove(C43S);
			}
			double mincc43 = Double.MAX_VALUE;
			Cluster minC43 = new Cluster(an);
			for(Cluster c : candidateC42.keySet()){
				if(candidateC43.get(c) < mincc43){
					mincc43 = candidateC43.get(c);
					minC43.Nodes.addAll(c.Nodes);
					minC43.Features.addAll(c.Features);
				}
			}
			newClusters4.add(minC43);		
			double cccase4 = mincc43;
			Cases.put(newClusters4, cccase4);
		}
			
		double ccfinalmin = Double.MAX_VALUE;
		ArrayList<Cluster> selectnewCluster = new ArrayList<Cluster>();
		for(ArrayList<Cluster> n : Cases.keySet()){
			System.out.println(Cases.get(n));		
			if(Cases.get(n) < ccfinalmin){
				ccfinalmin = Cases.get(n);
				selectnewCluster.clear();
				selectnewCluster.addAll(n);
			}
		}		
		if(ccfinalmin < cc_compare)
			newCluster.addAll(selectnewCluster);
			return newCluster;		
	}	
	
	public ArrayList<Cluster> ModifyClusterwithoutorder(ArrayList<Cluster> ClusList, Pair mergeCandidate){
		ArrayList<Cluster> newCluster = new ArrayList<Cluster>();	
		double cc_compare = cc.CalculateCC(ClusList);    // Compared coding cost
		
		HashMap<ArrayList<Cluster>, Double> Cases = new HashMap<ArrayList<Cluster>, Double>();
		ArrayList<Cluster> tempClusList = new ArrayList<Cluster>();
		ArrayList<Cluster> remainClusList = new ArrayList<Cluster>();
		tempClusList.addAll(ClusList);		
				
		// Remove merging clusters
		for (int m = 0; m < mergeCandidate.member.length; m ++) {
			tempClusList.remove(mergeCandidate.member[m]);
		}
		remainClusList.addAll(tempClusList); // cluster list remove merging clusters
		
		double ccCluster = cc.ccClusterStructure + cc.ccClusterFeature;
		ccCluster -= cc.CCClusterStructure(mergeCandidate.member[0]) + cc.CCClusterFeature(mergeCandidate.member[0]);
		ccCluster -= cc.CCClusterStructure(mergeCandidate.member[1]) + cc.CCClusterFeature(mergeCandidate.member[1]);
		
		
		// Case1	
		ArrayList<Cluster> case1ClusList = new ArrayList<Cluster>();
		case1ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters1 = new ArrayList<Cluster>(); // Produced new clusters
		Cluster C1 = new Cluster(an);
		HashSet<vertex> NodesC1 = new HashSet<vertex>();
		NodesC1.addAll(mergeCandidate.member[0].Nodes);
		NodesC1.addAll(mergeCandidate.member[1].Nodes);
		C1.Nodes.addAll(NodesC1);
		newClusters1.add(C1);
		case1ClusList.addAll(newClusters1);
		double initialCost1 = ccCluster+cc.CCClusterStructure(C1)+cc.CCUnClusterStructure(case1ClusList);
		case1ClusList = FindSubspacewithoutOrder(newClusters1, case1ClusList, initialCost1);		
//		double cccase1 = cc.CalculateCC(case1ClusList);
		double cccase1 = initialCost1 + cc.CCClusterFeature(C1)+cc.CCUnClusterFeature(case1ClusList);
		Cases.put(newClusters1, cccase1);
		
		//Case 2
		ArrayList<Cluster> case2ClusList = new ArrayList<Cluster>();
		case2ClusList.addAll(remainClusList);
		HashSet<vertex> CNNodes = new HashSet<vertex>(mergeCandidate.structureCN);			
		ArrayList<Cluster> newClusters2 = new ArrayList<Cluster>(); // Produced new clusters	
		Cluster C21 = new Cluster(an);
		Cluster C22 = new Cluster(an);
		HashSet<vertex> NodesC21 = new HashSet<vertex>();
		NodesC21.addAll(mergeCandidate.member[0].Nodes);
		HashSet<vertex> NodesC22 = new HashSet<vertex>();		
		for(vertex v : mergeCandidate.member[1].Nodes){
			if(! CNNodes.contains(v))
				NodesC22.add(v);
		}
		C21.Nodes.addAll(NodesC21);
		C22.Nodes.addAll(NodesC22);
		double initialCost2 = ccCluster;
		if(C21.Nodes.size() >= minVertexInClus){
			newClusters2.add(C21);
			initialCost2 += cc.CCClusterStructure(C21);
		}
		if(C22.Nodes.size() >= minVertexInClus){
			newClusters2.add(C22);
			initialCost2 += cc.CCClusterStructure(C22);
		}
		case2ClusList.addAll(newClusters2);
		initialCost2 += cc.CCUnClusterStructure(case2ClusList);
		case2ClusList = FindSubspacewithoutOrder(newClusters2, case2ClusList, initialCost2);		
//		double cccase2 = cc.CalculateCC(case2ClusList);
		double cccase2 = initialCost2+cc.CCUnClusterFeature(case2ClusList);
		if(C21.Nodes.size() >= minVertexInClus)
			cccase2 += cc.CCClusterFeature(C21);
		if(C22.Nodes.size() >= minVertexInClus)
			cccase2 += cc.CCClusterFeature(C22);
		Cases.put(newClusters2, cccase2);
		
		//Case3
		ArrayList<Cluster> case3ClusList = new ArrayList<Cluster>();
		case3ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters3 = new ArrayList<Cluster>(); // Produced new clusters	
		Cluster C31 = new Cluster(an);
		Cluster C32 = new Cluster(an);
		HashSet<vertex> NodesC31 = new HashSet<vertex>();
		NodesC31.addAll(mergeCandidate.member[1].Nodes);
		HashSet<vertex> NodesC32 = new HashSet<vertex>();		
		for(vertex v : mergeCandidate.member[0].Nodes){
			if(! CNNodes.contains(v))
				NodesC32.add(v);
		}
		C31.Nodes.addAll(NodesC31);
		C32.Nodes.addAll(NodesC32);	
		double initialCost3 = ccCluster;
		if(C31.Nodes.size() >= minVertexInClus){
			newClusters3.add(C31);
			initialCost3 += cc.CCClusterStructure(C31);
		}
		if(C32.Nodes.size() >= minVertexInClus){
			newClusters3.add(C32);
			initialCost3 += cc.CCClusterStructure(C32);
		}
		case3ClusList.addAll(newClusters3);
		initialCost3 += cc.CCUnClusterStructure(case3ClusList);
		case3ClusList = FindSubspacewithoutOrder(newClusters3, case3ClusList, initialCost3);	
		double cccase3 = initialCost3+cc.CCUnClusterFeature(case3ClusList);
		if(C31.Nodes.size() >= minVertexInClus)
			cccase3 += cc.CCClusterFeature(C31);
		if(C32.Nodes.size() >= minVertexInClus)
			cccase3 += cc.CCClusterFeature(C32);
//		double cccase3 = cc.CalculateCC(case3ClusList);
		Cases.put(newClusters3, cccase3);	
		
		// Case4: Both clusters are redundant
		ArrayList<Cluster> case4ClusList = new ArrayList<Cluster>();
		case4ClusList.addAll(remainClusList);
		ArrayList<Cluster> newClusters4 = new ArrayList<Cluster>(); // Produced new clusters
		Cluster C41 = new Cluster(an);
		Cluster C42 = new Cluster(an);
		Cluster C43 = new Cluster(an);
		C41.Nodes.addAll(NodesC22);		
		C42.Nodes.addAll(NodesC32);
		C43.Nodes.addAll(CNNodes);
		double initialCost4 = ccCluster;
		if(C41.Nodes.size() >= minVertexInClus){
			newClusters4.add(C41);
			initialCost4 += cc.CCClusterStructure(C41);
		}
		if(C42.Nodes.size() >= minVertexInClus){
			newClusters4.add(C42);
			initialCost4 += cc.CCClusterStructure(C42);
		}
		if(C43.Nodes.size() >= minVertexInClus){
			newClusters4.add(C43);
			initialCost4 += cc.CCClusterStructure(C43);
		}
		case4ClusList.addAll(newClusters4);
		initialCost4 += cc.CCUnClusterStructure(case4ClusList);
		case4ClusList = FindSubspacewithoutOrder(newClusters4, case4ClusList, initialCost4);		
//		double cccase4 = cc.CalculateCC(case4ClusList);
		double cccase4 = initialCost4+cc.CCUnClusterFeature(case4ClusList);
		if(C41.Nodes.size() >= minVertexInClus)
			cccase4 += cc.CCClusterFeature(C41);
		if(C42.Nodes.size() >= minVertexInClus)
			cccase4 += cc.CCClusterFeature(C42);
		if(C43.Nodes.size() >= minVertexInClus)
			cccase4 += cc.CCClusterFeature(C43);
		Cases.put(newClusters4, cccase4);	
		
		
		double ccfinalmin = Double.MAX_VALUE;
		ArrayList<Cluster> minCluster = new ArrayList<Cluster>();
		for(ArrayList<Cluster> n : Cases.keySet()){
//			for(Cluster c : n){
//				Output(c);
//			}
//			System.out.print(Cases.get(n)+" ");
			if(Cases.get(n) < ccfinalmin){
				ccfinalmin = Cases.get(n);
				minCluster.clear();
				minCluster.addAll(n);
			}
		}		
//		System.out.println();
		if(ccfinalmin < cc_compare)
			if(Math.abs(ccfinalmin-cc_compare) > 0.000000001)
				return minCluster;
			else
				return newCluster;
		else
			return newCluster;
	}
	
	
	public ArrayList<ArrayList<Integer>> FindClusterFeatureOrder(Cluster c){		
		int[][] adjFeature = c.FindFulladjFeature();
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
	

	public ArrayList<Cluster> ModifyStructure(ArrayList<Cluster> StructureClusList, ArrayList<Cluster> FeatureClusList, Pair mergeCandidate) {
		CodingCost cc = new CodingCost(an);
		// The coding cost of feature do not changed
		double ccfeature = cc.CCFeature(FeatureClusList);	
		
		// Merged Cluster Common Feature
		ArrayList<Integer> mergefeature = mergeCandidate.featureCN;
		
//		// Output most similar CLuster pair
//		for(int i = 0; i < mergeCandidate.member.length; i ++){
//			Cluster Clu = mergeCandidate.member[i];
//			System.out.print("The most similar pair:  ");
//			for(vertex v : mergeCandidate.member[i].Nodes){
//				System.out.print(v.id);
//			}
//			System.out.print("\n");
//		}
		
		ArrayList<Cluster> newClusters = new ArrayList<Cluster>(); // Produced new clusters
		
		
		ArrayList<Cluster> tempClusListS = new ArrayList<Cluster>(); // Temporary cluster list
		tempClusListS.addAll(StructureClusList);

		
		HashMap<ArrayList<Cluster>, Double> MergeCases = new HashMap<ArrayList<Cluster>, Double>();  // all possible of merge

		// Case1 : Merger two clusters into one
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
		C1.Features.addAll(mergefeature);		
		
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			tempClusListS.remove(mergeCandidate.member[i]);
		}
		tempClusListS.add(C1);
		newClusters.add(C1);
		double ccstructureCase1 = cc.CCStructure(tempClusListS);
		double cccase1 = ccstructureCase1 + ccfeature + cc.CCNoisenew(tempClusListS, FeatureClusList);
		MergeCases.put(newClusters, cccase1);

		
		
		ArrayList<vertex> CNNodes = mergeCandidate.structureCN;		
		// Case2: The first cluster is original, the second is redundant
		newClusters.clear();
		tempClusListS.clear();
		tempClusListS.addAll(StructureClusList);
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
		C21.Features.addAll(mergefeature);
		C22.Nodes.addAll(NodesC22);
		C22.Features.addAll(mergefeature);
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			tempClusListS.remove(mergeCandidate.member[i]);
		}
		if(C21.Nodes.size() > 0){
			tempClusListS.add(C21);
			newClusters.add(C21);
		}
		if(C22.Nodes.size() > 0){
			tempClusListS.add(C22);
			newClusters.add(C22);
		}
		double ccstructureCase2 = cc.CCStructure(tempClusListS);
		double cccase2 = ccstructureCase2 + ccfeature + cc.CCNoisenew(tempClusListS, FeatureClusList);
		MergeCases.put(newClusters, cccase2);
		
		// Case3: The second cluster is original, the first is redundant
		newClusters.clear();
		tempClusListS.clear();
		tempClusListS.addAll(StructureClusList);
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
		C31.Features.addAll(mergefeature);
		C32.Nodes.addAll(NodesC32);
		C32.Features.addAll(mergefeature);
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			tempClusListS.remove(mergeCandidate.member[i]);
		}
		if(C31.Nodes.size() > 0){
			tempClusListS.add(C31);
			newClusters.add(C31);
		}
		if(C32.Nodes.size() > 0){
			tempClusListS.add(C32);
			newClusters.add(C32);
		}		
		double ccstructureCase3 = cc.CCStructure(tempClusListS);
		double cccase3 = ccstructureCase3 + ccfeature + cc.CCNoisenew(tempClusListS, FeatureClusList);
		MergeCases.put(newClusters, cccase3);
		
		// Case4: Both clusters are redundant
		newClusters.clear();
		tempClusListS.clear();
		tempClusListS.addAll(StructureClusList);
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
		C41.Features.addAll(mergefeature);
		C42.Nodes.addAll(NodesC42);
		C42.Features.addAll(mergefeature);
		C43.Nodes.addAll(CNNodes);
		C43.Features.addAll(mergefeature);
		for (int i = 0; i < mergeCandidate.member.length; i ++) {
			tempClusListS.remove(mergeCandidate.member[i]);
		}
		if(C41.Nodes.size() > 0){
			tempClusListS.add(C41);
			newClusters.add(C41);
		}
		if(C42.Nodes.size() > 0){
			tempClusListS.add(C42);
			newClusters.add(C42);
		}
		if(C43.Nodes.size() > 0){
			tempClusListS.add(C43);
			newClusters.add(C43);
		}		
		double ccstructureCase4 = cc.CCStructure(tempClusListS);
		double cccase4 = ccstructureCase4 + ccfeature + cc.CCNoisenew(tempClusListS, FeatureClusList);
		MergeCases.put(newClusters, cccase4);

		// Find minimum cost
		double mincc = Double.MAX_VALUE;
		ArrayList<Cluster> MinnewCluster = new ArrayList<Cluster>();
		for(ArrayList<Cluster> cl : MergeCases.keySet()){
			if(MergeCases.get(cl) < mincc){
				mincc = MergeCases.get(cl);
				MinnewCluster.addAll(cl);
			}
		}
		
		double ccold = cc.CalculateCC(StructureClusList);	
		
		return MinnewCluster;
	}
	
	public TreeSet<Pair> ModifyQueue(ArrayList<Cluster> StructureClusList, TreeSet<Pair> simQueue, ArrayList<Cluster> newClusters, Pair mergeCandidate){
//		ArrayList<Pair> removePair = new ArrayList<Pair>();
//		for(Pair p : simQueue){
//			if(p.member[0].equals(mergeCandidate.member[0])||p.member[0].equals(mergeCandidate.member[1])
//			||p.member[1].equals(mergeCandidate.member[0])||p.member[1].equals(mergeCandidate.member[1]))
//				removePair.add(p);
//		}
//		for(Pair p : removePair){
//			simQueue.remove(p);
//		}
		for(int i = 0; i < newClusters.size(); i ++){
    		for(int j = 0; j < StructureClusList.size(); j ++){
    			if(! StructureClusList.get(j).equals(newClusters.get(i))){
    				Pair simPair = CalculateSim(newClusters.get(i), StructureClusList.get(j));
    				simQueue.add(simPair);
    			}
    		}
    	}
		return simQueue; 
	}

	public HashMap<vertex, Double> Similarity(ArrayList<vertex> CandidateNode,
			Cluster Core) {
		HashMap<vertex, Double> Sim = new HashMap<vertex, Double>();
		for (vertex v : CandidateNode) {
			double simstructure = 0.0;
			double simfeature = 0.0;
			for (vertex c : Core.Nodes) {
				if (an.exsitEdge(v.id, c.id))
					simstructure++;
				for (Integer f : Core.Features) {
					if (v.feature.contains(c.feature.get(f)))
						simfeature++;
				}
			}
			simstructure = simstructure / Core.Nodes.size();
			simfeature = simfeature
					/ (Core.Features.size() * Core.Nodes.size());
			Sim.put(v, simstructure + simfeature);
		}
		return Sim;
	}

	public ArrayList<Cluster> Findsubspace(ArrayList<Cluster> ClusList) {
		CodingCost cc = new CodingCost(an);
		double InitialCost = cc.CalculateCC(ClusList);
		for (Cluster c : ClusList) {

			// for(vertex v : c.Nodes){
			// System.out.println("Clus c Nodes:  " + v.id);
			// for(Integer fv : v.feature){
			// System.out.println("Nodes v feature:  " + fv);
			// }
			// }
			// for(Integer f : c.Features){
			// System.out.println("Clus c Feature:  " + f);
			// }

			HashMap<Integer, Double> featureent = new HashMap<Integer, Double>();
			for (Integer i : FullFeature) {
				int[][] adjFeature = c.FindFulladjFeature();
				int[] featurei = new int[c.Nodes.size()];
				for (int j = 0; j < c.Nodes.size(); j++) {
					featurei[j] = adjFeature[j][i];
				}
				double entvalue = featureentropy(featurei, i);
				featureent.put(i, entvalue);
			}

//			//
//			for (Integer fe : featureent.keySet()) {
//				System.out.println("Feature Entropy:   " + featureent.get(fe));
//			}
//			//
//
			ArrayList<Integer> Order = FeatureOrder(featureent);
//			//
//			for (Integer o : Order) {
//				System.out.println("Order :   " + o);
//			}
//			//

			// Feature Label
			int[] FeatureLabel = new int[Order.size()];
			for (int i = 0; i < FeatureLabel.length; i++) {
				FeatureLabel[i] = -1;
			}
			ArrayList<Integer> Ordercopy = new ArrayList<Integer>();
			Ordercopy.addAll(Order);
			int num = 0;
			while (Ordercopy.size() > 0) {
				FeatureLabel[Ordercopy.get(0)] = num;
				for (Integer o : Ordercopy) {
					if (featureent.get(Ordercopy.get(0)).equals(
							featureent.get(o))) {
						FeatureLabel[o] = num;
					}
				}
				for (int i = 0; i < FeatureLabel.length; i++) {
					if (FeatureLabel[i] == num) {
						Ordercopy.remove((Object) i);
					}
				}
				num++;
			}

			int labelnum = FeatureLabel[Order.get(Order.size() - 1)];
			c.Features.clear();
			for (int j = 0; j < FeatureLabel.length; j++) {
				if (FeatureLabel[j] == 0) {
					c.Features.add(j);
				}
			}
			double Cost = cc.CalculateCC(ClusList);

			for (int i = 1; i < labelnum; i++) {
				ArrayList<Integer> group = new ArrayList<Integer>();
				for (int j = 0; j < FeatureLabel.length; j++) {
					if (FeatureLabel[j] == i) {
						group.add(j);
					}
				}
				c.Features.addAll(group);
				double NewCost = cc.CalculateCC(ClusList);
				if (NewCost < Cost)
					Cost = NewCost;
				else {
					c.Features.removeAll(group);
					break;
				}
			}
		}
		return ClusList;
	}
	
	public ArrayList<Cluster> NewFindSubspace(ArrayList<Cluster> Clusters, ArrayList<Cluster> ClusList){	
		ArrayList<Cluster> tempClusters = new ArrayList<Cluster>();
		tempClusters.addAll(Clusters);
		
		while(tempClusters.size() > 0){
//			System.out.println("Size:" + tempClusters.size());
			for(Cluster cl : ClusList){
				if(tempClusters.contains(cl))
					cl.Features.clear();			
			}   // All feature of cluster candidate is null
			Cluster cluster = FindClusterSubspace(tempClusters, ClusList);	
			for(Cluster cl : ClusList){
				if(cl.equalsNodes(cluster)){
					cl.Features.clear();
					cl.Features.addAll(cluster.Features);
//					Output(cl);
				}
			}
			for(Cluster c : tempClusters){
				if(c.equalsNodes(cluster))
					cluster = c;
			}	
			tempClusters.remove(cluster);		
		}
		return ClusList;
	}
	
	public void Output(Cluster c){
		System.out.print("Nodes: ");
		for(vertex v : c.Nodes){
			System.out.print(v.id + ",");
		}
		System.out.println();
		System.out.print("Feature: ");
		for(Integer f : c.Features){
			System.out.print(f);
		}
		System.out.println();
	}
	
	public ArrayList<Cluster> FindSubspacewithoutOrder(ArrayList<Cluster> Clusters, ArrayList<Cluster> ClusList, double initialCost){
		for (Cluster c : Clusters) {			
			ArrayList<ArrayList<Integer>> clusterfeatureorder = FindClusterFeatureOrder(c); // feature with order
			double costmin = Double.MAX_VALUE;			
			for(Cluster cl : ClusList){
				if(cl.equals(c)){
					ArrayList<Integer> featuremin = new ArrayList<Integer>();					
					for (int i = 0; i < clusterfeatureorder.size(); i ++) {						
						cl.Features.addAll(clusterfeatureorder.get(i));
						double ccUnClusterFeature = cc.CCUnClusterFeature(ClusList);
						double ccClusterFeaturecl = cc.CCClusterFeature(cl);
						double cost = initialCost + ccUnClusterFeature + ccClusterFeaturecl;
//						double cost = cc.CalculateCC(ClusList);
						if(cost < costmin){
							costmin = cost;
							featuremin.clear();
							featuremin.addAll(cl.Features);
//							initialCost += ccClusterFeaturecl;
//							cc.ccUnClusterFeature = ccUnClusterFeature;
						}
					}
					cl.Features.clear();
					cl.Features.addAll(featuremin);	
					break;
				}	
			}	
//			TreeSet<Integer> points = new TreeSet<Integer>();
//			for(vertex v : c.Nodes)
//				points.add(v.id);
//			System.out.println(points);
//			System.out.println(c.Features);
		}		
		return ClusList;
	}
	
	public ArrayList<Cluster> FindSubspacewithoutOrder(ArrayList<Cluster> Clusters, ArrayList<Cluster> ClusList){
		for (Cluster c : Clusters) {			
			ArrayList<ArrayList<Integer>> clusterfeatureorder = FindClusterFeatureOrder(c); // feature with order
			double costmin = Double.MAX_VALUE;			
			for(Cluster cl : ClusList){
				if(cl.equals(c)){
					ArrayList<Integer> featuremin = new ArrayList<Integer>();					
					for (int i = 0; i < clusterfeatureorder.size(); i ++) {						
						cl.Features.addAll(clusterfeatureorder.get(i));
						double cost = cc.CalculateCC(ClusList);
						if(cost < costmin){
							costmin = cost;
							featuremin.clear();
							featuremin.addAll(cl.Features);
						}
					}
					cl.Features.clear();
					cl.Features.addAll(featuremin);	
					break;
				}	
			}			
		}		
		return ClusList;
	}
	
	public Cluster FindClusterSubspace(ArrayList<Cluster> Clusters, ArrayList<Cluster> ClusList){		
		CodingCost cc = new CodingCost(an);
		ArrayList<Cluster> subspaceClusList = new ArrayList<Cluster>(); // Subspace fixed clusters			
		// select each cluster and calculate subspace and subspace of other clusters set to null
		Cluster scluster = new Cluster(an);
		double minCostallCluster = Double.MAX_VALUE;
		for (Cluster c : Clusters) {
			subspaceClusList.clear();
			subspaceClusList.addAll(ClusList);
			ArrayList<ArrayList<Integer>> clusterfeatureorder = FindClusterFeatureOrder(c); // feature with order
			double costmin = Double.MAX_VALUE;
			Cluster clusterwithfeature = new Cluster(an);
			for(Cluster cl : subspaceClusList){
				if(cl.equals(c)){
					ArrayList<Integer> featuremin = new ArrayList<Integer>();					
					for (int i = 0; i < clusterfeatureorder.size(); i ++) {						
						cl.Features.addAll(clusterfeatureorder.get(i));
//						System.out.println(cl.Features);
						double cost = cc.CalculateCC(subspaceClusList);
						if(cost < costmin){
							costmin = cost;
							featuremin.clear();
							featuremin.addAll(cl.Features);
						}
					}
					cl.Features.clear();
					cl.Features.addAll(featuremin);	
					clusterwithfeature = cl.copy();
					cl.Features.clear();
//					Output(clusterwithfeature);
					break;
				}	
			}
			if(costmin < minCostallCluster){
				minCostallCluster = costmin;
				scluster.clear();
				scluster = clusterwithfeature.copy();
			}
		}
//		Output(scluster);
		
		return scluster;
	}


	public ArrayList<Cluster> RemoveSame(ArrayList<Cluster> ClusList){
		ArrayList<Cluster> remove = new ArrayList<Cluster>();
		for(int i = 0; i < ClusList.size()-1; i ++){
			for(int j = i + 1; j < ClusList.size(); j ++){
				if(i != j && ClusList.get(i).equals(ClusList.get(j)))
					remove.add(ClusList.get(j));
			}
		}
		for(Cluster r : remove){
			ClusList.remove(r);
		}
		return ClusList;
	}
	
	public ArrayList<Cluster> PostProcessing(ArrayList<Cluster> ClusList){		
		NewCodingCost cd = new NewCodingCost(an);			
		double cc_cluster = 0.0;
		double cc_noise = 0.0;
		double cc = 0.0;
		double ccremain = 0.0;
		for(int i = 0; i < ClusList.size(); i ++){
			Cluster modiclus = ClusList.get(i);	
			Output(modiclus);
			Cluster keepmodiclus = modiclus.copy(); 
			Output(keepmodiclus);
//			cc_noise = cd.CCUnClusterStructure(ClusList) + cd.CCUnClusterFeature(ClusList); // coding cost of un-clusterarea
//			cc_cluster = cd.CCClusterStructure(modiclus) + cd.CCClusterFeature(modiclus); // coding cost of cluster
			cc = cd.CalculateCC(ClusList);         // total coding cost
//			ccremain = cc - cc_cluster - cc_noise;  // coding cost of other clusters 
			ccremain = cd.ccClusterFeature + cd.ccClusterStructure - cd.CCClusterFeature(modiclus) - cd.CCClusterStructure(modiclus);
			for(int j = 0; j < numVertices; j ++){
				vertex v = an.getVertex(j);
				ArrayList<Cluster> modi = new ArrayList<Cluster>(); 
				
				if(modiclus.Nodes.contains(v)){
					modiclus.Nodes.remove(v);
					Output(modiclus);
					modi.add(modiclus);
					double cctemp = cd.CCClusterStructure(modiclus) + ccremain + cd.CCUnClusterStructure(ClusList);
					modi = FindSubspacewithoutOrder(modi, ClusList, cctemp);
					Output(modiclus);
					cc_cluster = cd.CCClusterStructure(modiclus) + cd.CCClusterFeature(modiclus);
					cc_noise = cd.CCUnClusterStructure(ClusList) + cd.CCUnClusterFeature(ClusList);
					double ccnew = cc_cluster + cc_noise + ccremain;
					double tmp = cd.CalculateCC(ClusList);
					if(ccnew > cc)						
						modiclus = keepmodiclus.copy() ;				
					else{
						keepmodiclus = modiclus.copy();
						cc = ccnew;
					}
					continue;
				}
				if(!modiclus.Nodes.contains(v)){
					modiclus.Nodes.add(v);
					Output(modiclus);
					modi.add(modiclus);
					double cctemp = cd.CCClusterStructure(modiclus) + ccremain + cd.CCUnClusterStructure(ClusList);
					modi = FindSubspacewithoutOrder(modi, ClusList, cctemp);
					cc_cluster = cd.CCClusterStructure(modiclus) + cd.CCClusterFeature(modiclus);
					cc_noise = cd.CCUnClusterStructure(ClusList) + cd.CCUnClusterFeature(ClusList);
					double ccnew = cc_cluster + cc_noise + ccremain;
					if(ccnew > cc)
						modiclus = keepmodiclus.copy();
					else
						keepmodiclus = modiclus.copy();
				}				
			}			
		}
		return ClusList;
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

	public ArrayList<vertex> NodeCommonNeighbor(ArrayList<vertex> a,
			ArrayList<vertex> b) {
		HashSet<vertex> hb = new HashSet<vertex>(b);
		ArrayList<vertex> CN = new ArrayList<vertex>();
		for (int i = 0; i < a.size(); i++) {
			if (hb.contains(a.get(i))) {
				CN.add(a.get(i));
			}
		}
		return CN;
	}

	public ArrayList<Integer> FeatureCommonNeighbor(ArrayList<Integer> a,
			ArrayList<Integer> b) {
		HashSet<Integer> hb = new HashSet<Integer>(b);
		ArrayList<Integer> CN = new ArrayList<Integer>();
		for (int i = 0; i < a.size(); i++) {
			if (hb.contains(a.get(i))) {
				CN.add(a.get(i));
			}
		}
		return CN;
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

}
