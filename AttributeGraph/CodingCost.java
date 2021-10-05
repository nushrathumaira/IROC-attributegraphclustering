package AttributeGraph;

import java.util.*;

import attributenetwork.attributenetwork;
import attributenetwork.vertex;

public class CodingCost {
	attributenetwork an;
	int numVertices;
	int numFeature;
	ArrayList<Integer> FullFeature;
	int[][] AdjNode;
	int[][] AdjFeature;
	
	public CodingCost(attributenetwork an){
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
	}
	
	public double CalculateCC(ArrayList<Cluster> ClusList) {
//		double cc_structure = CCStructure(ClusList);
		double cc_structure = CCStructurenoCombine(ClusList);
//		double cc_feature = CCFeature(ClusList);
		double cc_feature = CCFeaturenoCombine(ClusList);
//		double cc_noise = CCNoise(ClusList);
		double cc_noise = CCNoiseP(ClusList);
		double cc = cc_structure + cc_feature + cc_noise;
		if(Double.isNaN(cc))
			System.out.print("");
		return cc;
	}	
	
	public double CCNoise(ArrayList<Cluster> ClusList) {
		double cc_noise = 0.0;
		double cc_strucutenoise_para = 0.0;
		double cc_featurenoise_para = 0.0;
		int[][] NodeLabel = new int[numVertices][numVertices];
		int[][] FeatureLabel = new int[numVertices][numFeature];
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);			
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					NodeLabel[Clus.Nodes.get(j).id][Clus.Nodes.get(k).id] = 1;
				}
			}
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Features.size(); k++) {
					FeatureLabel[Clus.Nodes.get(j).id][Clus.Features.get(k)] = 1;
				}
			}

		}
		int numallnoisenode = 0;
		int numnoiseone = 0;
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (NodeLabel[i][j] == 0) {
					numallnoisenode++;
					if (AdjNode[i][j] == 1)
						numnoiseone++;
				}
			}
		}
		if(numallnoisenode != 0){
			if(numallnoisenode == numnoiseone || numnoiseone == 0)
				cc_strucutenoise_para += 1;
			else
				cc_strucutenoise_para += 0.5 * lg2(numallnoisenode);
		}			
		double cc_noise_feature = 0.0;
		for (int i = 0; i < numFeature; i++) {
			int[] categoryinonefeature = new int[an.numCategory[i]];
			int numnoiseinonefeature = 0;
			for (int j = 0; j < numVertices; j++) {
				if (FeatureLabel[j][i] == 0) {
					numnoiseinonefeature++;
					categoryinonefeature[AdjFeature[j][i]]++;
				}
			}
			int catenum = 0;
			for(int n = 0; n < categoryinonefeature.length; n ++){
				if(categoryinonefeature[n] != 0)
					catenum ++;
			}
			if(numnoiseinonefeature != 0){
				if(catenum == 1)
					cc_featurenoise_para += 0.5 * lg2(numnoiseinonefeature);
				else
					cc_featurenoise_para += 0.5* (catenum -1)*lg2(numnoiseinonefeature);
			}
				
			for (int k = 0; k < categoryinonefeature.length; k++) {
				if (categoryinonefeature[k] != 0) {
					cc_noise_feature += categoryinonefeature[k]
							* lg2((double) numnoiseinonefeature
									/ (double) categoryinonefeature[k]);
				}
			}
		}
		double cc_noise_structure_one = 0.0;
		if (numnoiseone != 0) {
			cc_noise_structure_one = numnoiseone
					* lg2((double) numallnoisenode / (double) numnoiseone);
		}
		double cc_noise_structure_zero = 0.0;
		if (numallnoisenode - numnoiseone != 0) {
			cc_noise_structure_zero = (numallnoisenode - numnoiseone)
					* lg2((double) numallnoisenode
							/ (double) (numallnoisenode - numnoiseone));
		}
		double cc_noise_structure = cc_noise_structure_one
				+ cc_noise_structure_zero;
		cc_noise = cc_noise_structure + cc_noise_feature + cc_strucutenoise_para + cc_featurenoise_para;
		return cc_noise;
	}
	
	public double CCNoiseP(ArrayList<Cluster> ClusList){
		double cc_noise = 0.0;
		double cc_strucutenoise_para = 0.0;
		double cc_featurenoise_para = 0.0;
		int[][] NodeLabel = new int[numVertices][numVertices];
		int[][] FeatureLabel = new int[numVertices][numFeature];
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);			
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					NodeLabel[Clus.Nodes.get(j).id][Clus.Nodes.get(k).id] = 1;
				}
			}
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Features.size(); k++) {
					FeatureLabel[Clus.Nodes.get(j).id][Clus.Features.get(k)] = 1;
				}
			}

		}
		int numallnoisenode = 0;
		int numnoiseone = 0;
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (NodeLabel[i][j] == 0) {
					numallnoisenode++;
					if (AdjNode[i][j] == 1)
						numnoiseone++;
				}
			}
		}
		if(numallnoisenode != 0){
			if(numallnoisenode == numnoiseone || numnoiseone == 0)
				cc_strucutenoise_para += 1;
			else
				cc_strucutenoise_para += 0.5 * lg2(numallnoisenode);
		}			
		double cc_noise_feature = 0.0;
		for (int i = 0; i < numFeature; i++) {			
			int[] categoryinonefeature = new int[an.numCategory[i]];
			int numnoiseinonefeature = 0;
			for (int j = 0; j < numVertices; j++) {
				if (FeatureLabel[j][i] == 0) {
					numnoiseinonefeature++;
					categoryinonefeature[AdjFeature[j][i]]++;
				}
			}
			int catenum = 0;
			for(int n = 0; n < categoryinonefeature.length; n ++){
				if(categoryinonefeature[n] != 0)
					catenum ++;
			}
			if(numnoiseinonefeature != 0){
				if(catenum == 1)
					cc_featurenoise_para += 0.5 * lg2(numnoiseinonefeature);
				else
					cc_featurenoise_para += 0.5* (catenum -1)*lg2(numnoiseinonefeature);
			}
				
			for (int k = 0; k < categoryinonefeature.length; k++) {
				if (categoryinonefeature[k] != 0) {
					double cateprob = an.pcategory.get(i).get(k);
					cc_noise_feature += categoryinonefeature[k]* cateprob * lg2(cateprob);
				}
			}
		}
		double cc_noise_structure_one = 0.0;
		if (numnoiseone != 0) {
			cc_noise_structure_one = numnoiseone * lg2((double) numallnoisenode / (double) numnoiseone);
		}
		double cc_noise_structure_zero = 0.0;
		if (numallnoisenode - numnoiseone != 0) {
			cc_noise_structure_zero = (numallnoisenode - numnoiseone)
					* lg2((double) numallnoisenode
							/ (double) (numallnoisenode - numnoiseone));
		}
		double cc_noise_structure = cc_noise_structure_one
				+ cc_noise_structure_zero;
		cc_noise = cc_noise_structure + cc_noise_feature + cc_strucutenoise_para + cc_featurenoise_para;
		return cc_noise;
	}

	public double CCStructure(ArrayList<Cluster> ClusList) {
		double cc_structure = 0.0;
		double cc_structure_para = 0.0;
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);
			int ClusSize = Clus.Nodes.size();
			int ClusAdjSize = Clus.Nodes.size() * Clus.Nodes.size();
			// structure information
			int numOnes = 0;
			int numZeros = 0;
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					if (AdjNode[Clus.Nodes.get(j).id][Clus.Nodes.get(k).id] == 1)
						numOnes++;
					else
						numZeros++;
				}
			}
			if(numOnes != 0)
				cc_structure += numOnes * lg2((double) ClusAdjSize / (double) numOnes);
			if(numZeros != 0)
				cc_structure += numZeros * lg2((double) ClusAdjSize / (double) numZeros);
			if(numOnes == ClusAdjSize || numZeros == ClusAdjSize)
				cc_structure_para += 1;
			else
				cc_structure_para += 0.5 * lg2(ClusAdjSize);			
		}
		ArrayList<ArrayList<ArrayList<vertex>>> newNodeID = CombineNodeID(ClusList);
		double cc_nodeid = 0.0;
		for (int i = 0; i < newNodeID.size(); i++) {
			ArrayList<ArrayList<vertex>> nodetable = newNodeID.get(i);
			double cc_nodetable = 0.0;
			int remainnode = 0; // Whether all the nodes seize all the table
			int sumnode = 0;
			for (int j = 0; j < nodetable.size(); j++) {
				sumnode += nodetable.get(j).size();
			}
			remainnode = numVertices - sumnode;
			if (remainnode == 0) {
				for (int j = 0; j < nodetable.size(); j++) {
					if(nodetable.get(j).size() != 0)
					cc_nodetable += nodetable.get(j).size()
							* lg2((double) numVertices
									/ (double) nodetable.get(j).size());
				}
			}
			if (remainnode > 0) {
				for (int j = 0; j < nodetable.size(); j++) {
					if(nodetable.get(j).size() != 0)
					cc_nodetable += nodetable.get(j).size()
							* lg2((double) numVertices
									/ (double) nodetable.get(j).size());
				}
				cc_nodetable += remainnode
						* lg2((double) numVertices / (double) remainnode);
			}
			cc_nodeid += cc_nodetable;
		}

		double CCstructure = cc_structure + cc_nodeid + cc_structure_para;
		
		if(Double.isNaN(CCstructure))
			System.out.print("");

		return CCstructure;
	}
	
	public double CCStructurenoCombine(ArrayList<Cluster> ClusList) {
		double cc_structure = 0.0;
		double cc_structure_para = 0.0;
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);
			int ClusSize = Clus.Nodes.size();
			int ClusAdjSize = Clus.Nodes.size() * Clus.Nodes.size();
			// structure information
			int numOnes = 0;
			int numZeros = 0;
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					if (AdjNode[Clus.Nodes.get(j).id][Clus.Nodes.get(k).id] == 1)
						numOnes++;
					else
						numZeros++;
				}
			}
			if(numOnes != 0)
				cc_structure += numOnes * lg2((double) ClusAdjSize / (double) numOnes);
			if(numZeros != 0)
				cc_structure += numZeros * lg2((double) ClusAdjSize / (double) numZeros);
			if(numOnes == ClusAdjSize || numZeros == ClusAdjSize)
				cc_structure_para += 1;
			else
				cc_structure_para += 0.5 * lg2(ClusAdjSize);			
		}		
		double cc_nodeid = 0.0;
		for(Cluster c : ClusList){
			if(c.Nodes.size() != 0 && c.Nodes.size() != numVertices)
				cc_nodeid += c.Nodes.size() * lg2((double)numVertices/(double)c.Nodes.size()) + (numVertices - c.Nodes.size()) * lg2((double)numVertices/(double)(numVertices - c.Nodes.size()));
		}

		double CCstructure = cc_structure + cc_nodeid + cc_structure_para;
		
		if(Double.isNaN(CCstructure))
			System.out.print("");

		return CCstructure;
	}
	
	public double CCFeaturenoCombine(ArrayList<Cluster> ClusList) {
		double cc_feature = 0.0;
		double cc_feature_para = 0.0;
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);
			int ClusSize = Clus.Nodes.size();
			for (int j = 0; j < Clus.Features.size(); j++) {
				int featureid = Clus.Features.get(j);
				int[] Categorynum = new int[an.numCategory[featureid]];
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					Categorynum[Clus.Nodes.get(k).feature.get(featureid)]++;
				}
				for (int m = 0; m < Categorynum.length; m++) {
					if (Categorynum[m] != 0)
						cc_feature += Categorynum[m] * lg2((double) ClusSize / (double) Categorynum[m]);
				}
				int numpara = 0;
				for(int n = 0; n < Categorynum.length; n ++){
					if(Categorynum[n] != 0)
						numpara ++;
				}
				if(numpara == 1)
					cc_feature_para +=  0.5 * lg2(ClusSize);
				else
					cc_feature_para += 0.5* (numpara-1) * lg2(ClusSize);
				// System.out.println(cc_feature);
			}
		}	
		double cc_featureid = 0.0;
		for(Cluster c : ClusList){				
			if(c.Features.size() != 0 && c.Features.size() != numFeature)
				cc_featureid += c.Features.size() * lg2((double)numFeature/(double)c.Features.size()) + (numFeature-c.Features.size())*lg2((double)numFeature/((double)(numFeature-c.Features.size())));
		}
		double CCfeature = cc_feature + cc_featureid + cc_feature_para;
		if(Double.isNaN(CCfeature))
			System.out.print("");

		return CCfeature;
	}


	public double CCFeature(ArrayList<Cluster> ClusList) {
		double cc_feature = 0.0;
		double cc_feature_para = 0.0;
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);
			int ClusSize = Clus.Nodes.size();
			for (int j = 0; j < Clus.Features.size(); j++) {
				int featureid = Clus.Features.get(j);
				int[] Categorynum = new int[an.numCategory[featureid]];
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					Categorynum[Clus.Nodes.get(k).feature.get(featureid)]++;
				}
				for (int m = 0; m < Categorynum.length; m++) {
					if (Categorynum[m] != 0)
						cc_feature += Categorynum[m] * lg2((double) ClusSize / (double) Categorynum[m]);
				}
				int numpara = 0;
				for(int n = 0; n < Categorynum.length; n ++){
					if(Categorynum[n] != 0)
						numpara ++;
				}
				if(numpara == 1)
					cc_feature_para +=  0.5 * lg2(ClusSize);
				else
					cc_feature_para += 0.5* (numpara-1) * lg2(ClusSize);
				// System.out.println(cc_feature);
			}
		}
		ArrayList<ArrayList<ArrayList<Integer>>> newFeatureID = CombineFeatureID(ClusList);
		double cc_featureid = 0.0;
		for (int i = 0; i < newFeatureID.size(); i++) {
			ArrayList<ArrayList<Integer>> featuretable = newFeatureID.get(i);
			double cc_featuretable = 0.0;
			int remainfeature = 0;
			int sumfeature = 0;
			for (int j = 0; j < featuretable.size(); j++) {
				sumfeature += featuretable.get(j).size();
			}
			remainfeature = numFeature - sumfeature;
			if (remainfeature == 0) {
				for (int j = 0; j < featuretable.size(); j++) {
					if(featuretable.get(j).size() != 0)
						cc_featuretable += featuretable.get(j).size()* lg2(numFeature / featuretable.get(j).size());
				}
			}
			if (remainfeature > 0) {
				for (int j = 0; j < featuretable.size(); j++) {
//					System.out.println(featuretable.get(j).size());
					if(featuretable.get(j).size() != 0)
						cc_featuretable += (double)(featuretable.get(j).size())* lg2((double)numFeature / (double)(featuretable.get(j).size()));
				}
				cc_featuretable += remainfeature * lg2(numFeature / remainfeature);
			}
			cc_featureid += cc_featuretable;
		}
		double CCfeature = cc_feature + cc_featureid + cc_feature_para;

		return CCfeature;
	}

	public ArrayList<ArrayList<ArrayList<vertex>>> CombineNodeID(ArrayList<Cluster> ClusList) {
		ArrayList<ArrayList<ArrayList<vertex>>> newNodeID = new ArrayList<ArrayList<ArrayList<vertex>>>();
		ArrayList<ArrayList<vertex>> NodeID = new ArrayList<ArrayList<vertex>>();
		for (int i = 0; i < ClusList.size(); i++) {
			NodeID.add(ClusList.get(i).Nodes);
		}
		while (NodeID.size() > 0) {			
			ArrayList<vertex> top = NodeID.get(0);
			NodeID.remove(top);
			ArrayList<ArrayList<vertex>> combine = new ArrayList<ArrayList<vertex>>();
			combine.add(top);

			if (NodeID.size() == 0) {
				newNodeID.add(combine);
				break;
			}

			ArrayList<ArrayList<vertex>> removeClus = new ArrayList<ArrayList<vertex>>();
			for (int i = 0; i < NodeID.size(); i++) {
				if (!haveOverlapNode(combine, NodeID.get(i))) {
					combine.add(NodeID.get(i));
					removeClus.add(NodeID.get(i));
				}
			}
			for (int i = 0; i < removeClus.size(); i++) {
				NodeID.remove(removeClus.get(i));
			}
			newNodeID.add(combine);
		}
		return newNodeID;
	}

	public ArrayList<ArrayList<ArrayList<Integer>>> CombineFeatureID(
			ArrayList<Cluster> ClusList) {
		ArrayList<ArrayList<ArrayList<Integer>>> newFeatureID = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ArrayList<ArrayList<Integer>> FeatureID = new ArrayList<ArrayList<Integer>>();
		for (int i = 0; i < ClusList.size(); i++) {
			FeatureID.add(ClusList.get(i).Features);
		}
		while (FeatureID.size() > 0) {
			ArrayList<Integer> top = FeatureID.get(0);
			FeatureID.remove(top);
			ArrayList<ArrayList<Integer>> combine = new ArrayList<ArrayList<Integer>>();
			combine.add(top);

			if (FeatureID.size() == 0) {
				newFeatureID.add(combine);
				break;
			}

			ArrayList<ArrayList<Integer>> removeClus = new ArrayList<ArrayList<Integer>>();
			for (int i = 0; i < FeatureID.size(); i++) {
				if (!haveOverlapFeature(combine, FeatureID.get(i))) {
					combine.add(FeatureID.get(i));
					removeClus.add(FeatureID.get(i));
				}
			}
			for (int i = 0; i < removeClus.size(); i++) {
				FeatureID.remove(removeClus.get(i));
			}
			newFeatureID.add(combine);
		}
		return newFeatureID;
	}

	public boolean haveOverlapNode(ArrayList<ArrayList<vertex>> a,
			ArrayList<vertex> b) {
		ArrayList<vertex> CN = new ArrayList<vertex>();
		for (int i = 0; i < b.size(); i++) {
			for (int j = 0; j < a.size(); j++) {
				if (a.get(j).contains(b.get(i))) {
					CN.add(b.get(i));
					return true;
				}
			}
		}		
		return false;
	}

	public boolean haveOverlapFeature(ArrayList<ArrayList<Integer>> a, ArrayList<Integer> b) {
		ArrayList<Integer> CN = new ArrayList<Integer>();
		for (int i = 0; i < b.size(); i++) {
			for (int j = 0; j < a.size(); j++) {
				if (a.get(j).contains(b.get(i))){
					CN.add(b.get(i));
					return true;
				}
			}
		}
		return false;
	}
	
	private double lg2(double d) {
		return Math.log(d) / Math.log(2.0);
	}
	
	
	
	public double CCNoisenew(ArrayList<Cluster> StructureClusList, ArrayList<Cluster> FeatureClusList){
		double cc_noise = 0.0;
		int[][] NodeLabel = new int[numVertices][numVertices];
		int[][] FeatureLabel = new int[numVertices][numFeature];
		// Structure Part
		for (int i = 0; i < StructureClusList.size(); i++) {
			Cluster ClusN = StructureClusList.get(i);
			for (int j = 0; j < ClusN.Nodes.size(); j++) {
				for (int k = 0; k < ClusN.Nodes.size(); k++) {
					NodeLabel[ClusN.Nodes.get(j).id][ClusN.Nodes.get(k).id] = 1;
				}
			}
		}
		int numallnoisenode = 0;
		int numnoiseone = 0;
		for (int i = 0; i < numVertices; i++) {
			for (int j = 0; j < numVertices; j++) {
				if (NodeLabel[i][j] == 0) {
					numallnoisenode++;
					if (AdjNode[i][j] == 1)
						numnoiseone++;
				}
			}
		}
		
		//Feature Part
		for(int i = 0; i < FeatureClusList.size(); i ++){
			Cluster ClusF = FeatureClusList.get(i);
			for (int j = 0; j < ClusF.Nodes.size(); j++) {
				for (int k = 0; k < ClusF.Features.size(); k++) {
					FeatureLabel[ClusF.Nodes.get(j).id][ClusF.Features.get(k)] = 1;
				}
			}
		}
		
		double cc_noise_feature = 0.0;
		for (int i = 0; i < numFeature; i++) {
			int[] categoryinonefeature = new int[an.numCategory[i]];
			int numnoiseinonefeature = 0;
			for (int j = 0; j < numVertices; j++) {
				if (FeatureLabel[j][i] == 0) {
					numnoiseinonefeature++;
					categoryinonefeature[AdjFeature[j][i]]++;
				}
			}
			for (int k = 0; k < categoryinonefeature.length; k++) {
				if (categoryinonefeature[k] != 0) {
					cc_noise_feature += categoryinonefeature[k]
							* lg2((double) numnoiseinonefeature
									/ (double) categoryinonefeature[k]);
				}
			}
		}
		double cc_noise_structure_one = 0.0;
		if (numnoiseone != 0) {
			cc_noise_structure_one = numnoiseone
					* lg2((double) numallnoisenode / (double) numnoiseone);
		}
		double cc_noise_structure_zero = 0.0;
		if (numallnoisenode - numnoiseone != 0) {
			cc_noise_structure_zero = (numallnoisenode - numnoiseone)
					* lg2((double) numallnoisenode
							/ (double) (numallnoisenode - numnoiseone));
		}
		double cc_noise_structure = cc_noise_structure_one
				+ cc_noise_structure_zero;
		cc_noise = cc_noise_structure + cc_noise_feature;
		return cc_noise;		
	}
	
	public double CalculateCCSeperate(ArrayList<Cluster> StructureClusList, ArrayList<Cluster> FeatureClusList){
		double cc_structure = CCStructure(StructureClusList);
		double cc_feature = CCFeature(FeatureClusList);
		double cc_noise = CCNoisenew(StructureClusList, FeatureClusList);
		double cc = cc_structure + cc_feature + cc_noise;
		return cc;
	}

}
