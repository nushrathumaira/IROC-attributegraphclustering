package dbcsc.main;

import dbcsc.base.DenseCluster;
import dbcsc.base.DenseClusterCandidates;
import dbcsc.base.DenseClusterComparator;
import dbcsc.base.Log;
import dbcsc.base.Parameter;
import dbcsc.base.Subspace;
import dbcsc.base.Timer;
import dbcsc.graph.Graph;
import dbcsc.graph.Node;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.PriorityQueue;
import java.util.Properties;

public class DBCSC {
	
	public static void main(String[] args) {
        if (args.length < 1) {
            System.out.print("Required argument: Properties file");
            System.exit(0);
        }
        String config_file = args[0];
        String filename = "";

        Properties properties;
        try {
            properties = new Properties();
            FileInputStream stream = new FileInputStream(config_file);
            properties.load(stream);
            stream.close();

            Parameter.k_min = Integer.parseInt(properties.getProperty("k_min"));
            Parameter.epsilon = Double.parseDouble(properties.getProperty("epsilon"));
            Parameter.min_pts = Integer.parseInt(properties.getProperty("min_pts"));
            Parameter.r_dim = Double.parseDouble(properties.getProperty("r_dim"));
            Parameter.r_obj = Double.parseDouble(properties.getProperty("r_obj"));
            Parameter.s_min = Integer.parseInt(properties.getProperty("s_min"));

            filename = properties.getProperty("2clusters400");
        } catch (FileNotFoundException e) {
            System.out.println("Properties file not found");
            e.printStackTrace();
            System.exit(0);
        } catch (IOException e) {
            System.out.println("Properties file could not be read");
            e.printStackTrace();
            System.exit(0);
        }

        //Einlesen des Graphen
        Graph myGraph = null;
        try {
            myGraph =dbcsc.graphReader.GraphReader.loadGraphFromGraphMLFile(filename);
            Parameter.numberOfAtts = myGraph.getNumberOfAtts();

        } catch (Exception e) {
            e.printStackTrace();
            System.out.print("Graph konnte nicht eingelesen werden.");
            System.exit(0);
        }

        dbcsc(config_file,myGraph);

    }
	
	
	
	public static void dbcsc(String config_file,Graph myGraph) {
		
		//Log-Datei
		String foundfile = config_file.replace(".properties", "_dbcsc.found");
		Log found = new Log(foundfile,true,false,true);
    	
    	Timer timer= new Timer();
		timer.start();
		
		//Initiales Pruning fuer density-based clustering
    	myGraph.densityPruning();
		
    	//Finden der Cluster
    	DenseClusterComparator comp = new DenseClusterComparator();
    	PriorityQueue<DenseCluster> clustering = new PriorityQueue<DenseCluster>(10,comp);
    	
    	//Result set
    	ArrayList<DenseCluster> non_red_clusters = new ArrayList<DenseCluster>();
    	
    	Subspace empty_sub = new Subspace();
    	for (int dim=0;dim<Parameter.numberOfAtts;dim++){
    		empty_sub.removeDimension(dim);
    	}
    	
    	for (int dim =0;dim < Parameter.numberOfAtts;dim++){
    		Subspace one_dim = empty_sub.copy();
    		one_dim.setDimension(dim, Double.NaN, Double.NaN);
    		
    		dfs_traversal(one_dim,myGraph.getNodes(),new ArrayList<DenseCluster>(), clustering, myGraph,false);
    	}
    	
    	while(!clustering.isEmpty()) {
    		DenseCluster cluster1 = clustering.poll();
    		if(cluster1 instanceof DenseClusterCandidates) {
    			//Falls keiner der Better_Parents ausgegeben wurde, werte die Kandidaten aus und fuege sie in die Queue ein, sonst verwerde Kandidaten
    			
    			boolean better_cluster_in_output =false;
    			
    			for (DenseCluster better_cluster : ((DenseClusterCandidates) cluster1).getBetterParents() ){
    				if(non_red_clusters.contains(better_cluster)){
    					better_cluster_in_output = true;
    					break;
    				}
    			}
    			
    			if(!better_cluster_in_output) {
    				//Teilbaum auswerten
    				dfs_traversal(cluster1.getSubspace(), cluster1.getNodes(),((DenseClusterCandidates) cluster1).getParents(), clustering, myGraph, false);
    			} 
    			//else: Kandidat verwerfen, nichts tun
    			
    		} else {	
	    		double quality1= cluster1.getQuality();
	    		//Vergleiche Cluster1 mit allen nicht-redundanten Clusters, sie haben alle eine Qualitaet >=quality1
	    		boolean redundant = false; 
	    		for(DenseCluster cluster2 : non_red_clusters) {
	    			double quality2=cluster2.getQuality();
	    			if(quality1 < quality2) {
	    				//Ueberlappung der Knoten ueberpruefen
	    				ArrayList<Node> intersection = new ArrayList<Node>(cluster1.getNodes());
	    				intersection.retainAll(cluster2.getNodes());
	    				if(intersection.size()>=Parameter.r_obj*cluster1.getNodes().size()) {
	    					//Ueberlappung der Dimensionen ueberpruefen
	    					boolean[]dim1 = cluster1.getSubspace().getDimensions();
	    					boolean[]dim2 = cluster2.getSubspace().getDimensions();
	    					int intersectDims = 0;
	    					for(int i =0;i<dim1.length;i++) {
	    						if(dim1[i] && dim2[i]){
	    							intersectDims++;
	    						}
	    					}
	    					if(intersectDims >= Parameter.r_dim*cluster1.getSubspace().size()) {
		    					//group1 ist redundant gegen group2
		    					redundant=true;
		    					break;
	    					}
	    				}
		    		} else {
		    			break;
		    		}
	    		}
	    		if(!redundant){
					non_red_clusters.add(cluster1);
					found.log(cluster1.toString()+"\n");
				}
    		}
    	}
    	String runtime = timer.toString();
    	found.log(runtime);
    }
		
		
	public static void dfs_traversal(Subspace sub, Collection<Node> cands, Collection<DenseCluster> parents, PriorityQueue<DenseCluster> queue, Graph orig_graph, boolean output_enriched){
		
		//System.out.println("DFS traversal: Subspace "+sub.toString());
		
		ArrayList<DenseCluster> found_clusters = new ArrayList<DenseCluster>();
		ArrayList<Collection<Node>> prelim_clusters = new ArrayList<Collection<Node>>();
		prelim_clusters.add(cands);
		
		while(!prelim_clusters.isEmpty()){
			HashSet<Node> current_orig = new HashSet<Node>(prelim_clusters.get(0));
			prelim_clusters.remove(0);
			
			HashSet<Node> current = new HashSet<Node>();
			//Subgraph muss Knoten aus current_orig, aber deren Kanten aus dem Originalgraphen!!! enthalten
			for(Node node : current_orig){
				Node new_node = node.copyWithoutNeighbours();
				current.add(new_node);
			}
			
			for(Node node : orig_graph.getNodes()){
				
				for(Node current_node: current){
					if(current_node.getID() == node.getID()){
						for(Node neighbor : node.getNeighbors()){
							for (Node current_neighbor : current){
								if(current_neighbor.getID() == neighbor.getID()){
									current_node.addNeighbor(current_neighbor);
									current_neighbor.addNeighbor(current_node);
									break;
								}
							}
						}
						break;
					}
				}
			}
			
			
			//Enriched subgraph berechnen
			//System.out.println("Compute Enriched Subgraph");
			Graph enriched = getEnrichedSubgraph(new Graph(current,Parameter.numberOfAtts), sub);
			//System.out.println("Finished");
			
			//Darin (MinPts-1)-cores finden
			ArrayList<HashSet<Node>> cores = enriched.detect_cores();
			
			if(cores.size()==1 && cores.get(0).size() == current.size()){
				//Cluster muss Knoten aus Original-Graphen enthalten, sonst geht Redundanz-Berechnung nicht 
				HashSet<Node> node_set = new HashSet<Node>();
				for(Node node : current){
					node_set.add(orig_graph.getNodeHavingID(node.getID()));
				}
				
				DenseCluster cluster = new DenseCluster(node_set, sub, DenseCluster.quality(current.size(), sub.size()));
				found_clusters.add(cluster);
				
			} else {
				prelim_clusters.addAll(cores);
			}
			
		}
		if(sub.size()>=Parameter.s_min) {
			queue.addAll(found_clusters);
		}
		
		for(DenseCluster cluster: found_clusters){
			//Set Enumeration Tree f�r Subspaces -> nur "dahinter" liegende Dimensionen k�nnen aufgenommen werden
			for (int dim=sub.get_max_d()+1;dim<Parameter.numberOfAtts;dim++){
				
				Subspace new_sub = sub.copy();
				new_sub.setDimension(dim, Double.NaN, Double.NaN);
				int addable_dims = Parameter.numberOfAtts - new_sub.get_max_d() - 1; //Wie viele Dimensionen koennen noch aufgenommen werden?
				
				//Parents-Menge f�r Subtree unter cluster, cluster wird hinzugef�gt
				ArrayList<DenseCluster> new_parents = new ArrayList<DenseCluster>(parents);
				if(sub.size()>=Parameter.s_min){
					new_parents.add(cluster);
				}
				
				//Parents finden, zu denen Subtree auf jeden Fall redundant ist
				ArrayList<DenseCluster> new_better_parents =new ArrayList<DenseCluster>();
				double q_max=cluster.getNodes().size() * (new_sub.size() + addable_dims); //Max. Qualitaet, die in Subtree auftreten kann
				
				for(DenseCluster parent : new_parents){
					//testen, ob subtree zu parent redundant ist
					if(parent.getQuality()> q_max && parent.getSubspace().size() >= Parameter.r_dim * (new_sub.size() + addable_dims)){
						new_better_parents.add(parent);
					}
				}
				
				if (new_better_parents.isEmpty()){
					dfs_traversal(new_sub, new HashSet<Node>(cluster.getNodes()), new_parents , queue, orig_graph, output_enriched);
				} else {
					//Cluster-Candidates einf�gen, Subtree erstmal nicht durchsuchen
					DenseClusterCandidates subtree = new DenseClusterCandidates(new HashSet<Node>(cluster.getNodes()), new_sub, q_max, new_parents, new_better_parents);
					queue.add(subtree);
				}
			}
		}
	}
	
	
	
	//Hilfsfunktionen
	
	 public static Graph getEnrichedSubgraph(Graph baseGraph, Subspace sub) {
		 Graph result = getEnrichedSubgraph(baseGraph, Parameter.k_min);
		 
		 //F�r alle k-Nachbarn �berpr�fen, ob sie auch auch �hnlich sind
		
		 for (Node node :  result.getNodes()){
			 for (Iterator<Node> it= node.getNeighbors().iterator();it.hasNext();){
				 Node neighbour = it.next();
				 for(int i=0;i<Parameter.numberOfAtts;i++){
					 if(sub.hasDimension(i)){
						 if(Math.abs(node.getAttribute(i) - neighbour.getAttribute(i))>Parameter.epsilon){
							 it.remove();
							 neighbour.getNeighbors().remove(node);
							 break;
						 }
					 }
				 }
			 }
		 }
		 
		 return result;
	 }

	
	 private static Graph getEnrichedSubgraph(Graph baseGraph, int k_min) {
	      
	        Graph enrichedGraph = baseGraph.copy();
	        while (k_min > 1) {
	            Graph tempGraph = enrichedGraph.copy();
	            for (Node node : enrichedGraph.getNodes()) {
	                for (Node neighbour : node.getNeighbors()) {
	                    for (Node extendedNeighbour : neighbour.getNeighbors()) {
	                        // check 1. if the new proposed edge is already an existing edge
	                        if (node.getNeighbors().contains(extendedNeighbour) == false && node.getID() != extendedNeighbour.getID()) {
	                            Node n1 = tempGraph.getNodeHavingID(node.getID());
	                            Node n2 = tempGraph.getNodeHavingID(extendedNeighbour.getID());
	                            tempGraph.addEdge(n1, n2);
	                        }
	                    }
	                }
	            }
	            k_min--;
	            enrichedGraph = tempGraph;
	        }
	        return enrichedGraph;
	    }


}
