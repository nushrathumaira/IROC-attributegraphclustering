package AttributeGraph;

import java.util.*;

import attributenetwork.attributenetwork;

public class NewCodingCost {
	attributenetwork an;
	int numVertices;
	int numFeature;
	ArrayList<Integer> FullFeature;
	int[][] AdjNode;
	int[][] AdjFeature;
	
	double ccClusterStructure;
	double ccClusterFeature;
	double ccUnClusterStructure;
	double ccUnClusterFeature;
	
	public NewCodingCost(attributenetwork an){
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
	
	public double CalculateCC(ArrayList<Cluster> ClusList){
		double ccClusterStructure = 0.0;
		double ccClusterFeature = 0.0;
		double ccUnClusterStructure = 0.0;
		double ccUnClusterFeature = 0.0;
		
		for(Cluster Clus : ClusList){
			ccClusterStructure += CCClusterStructure(Clus);
			ccClusterFeature += CCClusterFeature(Clus);
		}
		
		ccUnClusterStructure += CCUnClusterStructure(ClusList);
		ccUnClusterFeature += CCUnClusterFeature(ClusList);
		
		this.ccClusterStructure = ccClusterStructure;
		this.ccClusterFeature = ccClusterFeature;
		this.ccUnClusterStructure = ccUnClusterStructure;
		this.ccUnClusterFeature = ccUnClusterFeature;
		
		return ccClusterStructure+ccUnClusterStructure+ccClusterFeature+ccUnClusterFeature;
	}
	
	public double CCClusterStructure(Cluster Clus){
		double cc_structure = 0.0;
		double cc_structure_para = 0.0;
		double cc_nodeid = 0.0;
		
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
	
		if(Clus.Nodes.size() != 0 && Clus.Nodes.size() != numVertices)
			cc_nodeid += Clus.Nodes.size() * lg2((double)numVertices/(double)Clus.Nodes.size()) + (numVertices - Clus.Nodes.size()) * lg2((double)numVertices/(double)(numVertices - Clus.Nodes.size()));
		
		double CCstructure = cc_structure + cc_nodeid + cc_structure_para;
		
		return CCstructure;
	}
	
	public double CCClusterFeature(Cluster Clus){		
		double cc_feature = 0.0;
		double cc_feature_para = 0.0;
		double cc_featureid = 0.0;
		
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
		
		if(Clus.Features.size() != 0 && Clus.Features.size() != numFeature)
			cc_featureid += Clus.Features.size() * lg2((double)numFeature/(double)Clus.Features.size()) + (numFeature-Clus.Features.size())*lg2((double)numFeature/((double)(numFeature-Clus.Features.size())));
		
		double CCfeature = cc_feature + cc_featureid + cc_feature_para;
		
		return CCfeature;
	}
	
	public double CCUnClusterStructure(ArrayList<Cluster> ClusList){
		double cc_strucutenoise_para = 0.0;
		int[][] NodeLabel = new int[numVertices][numVertices];
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);			
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Nodes.size(); k++) {
					NodeLabel[Clus.Nodes.get(j).id][Clus.Nodes.get(k).id] = 1;
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
		return cc_noise_structure + cc_strucutenoise_para;
	}
	
	public double CCUnClusterFeature(ArrayList<Cluster> ClusList){
		double cc_featurenoise_para = 0.0;
		int[][] FeatureLabel = new int[numVertices][numFeature];
		for (int i = 0; i < ClusList.size(); i++) {
			Cluster Clus = ClusList.get(i);			
			for (int j = 0; j < Clus.Nodes.size(); j++) {
				for (int k = 0; k < Clus.Features.size(); k++) {
					FeatureLabel[Clus.Nodes.get(j).id][Clus.Features.get(k)] = 1;
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

		return cc_noise_feature + cc_featurenoise_para;
	}
	
	private double lg2(double d) {
		return Math.log(d) / Math.log(2.0);
	}

}
