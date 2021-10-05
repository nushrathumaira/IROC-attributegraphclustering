package dbcsc.base;

import dbcsc.graph.Node;

import java.util.ArrayList;
import java.util.HashSet;

public class DenseClusterCandidates extends DenseCluster{
	protected ArrayList<DenseCluster> parents; //Cluster auf dem Pfad
	protected ArrayList<DenseCluster> better_parents; //Cluster auf dem Pfad, zu dem Candidate redundant ist
	
	public DenseClusterCandidates(HashSet<Node> nodes, Subspace subspace, double q_max, ArrayList<DenseCluster> parents, ArrayList<DenseCluster> better_parents) {
		//Quality ist hier die maximal mögliche Qualität, Subspace der Subspace von nodes
		super(nodes,subspace,q_max);
		this.parents = parents;
		this.better_parents = better_parents;
	}
	
	public DenseClusterCandidates(HashSet<Node> nodes, HashSet<Node> borders, Subspace subspace, double q_max, ArrayList<DenseCluster> parents, ArrayList<DenseCluster> better_parents) {
		//Quality ist hier die maximal mögliche Qualität, Subspace der Subspace von nodes
		super(nodes,borders,subspace,q_max);
		this.parents = parents;
		this.better_parents = better_parents;
	}

	public ArrayList<DenseCluster> getParents() {
		return parents;
	}

	public ArrayList<DenseCluster> getBetterParents() {
		return better_parents;
	}
}
