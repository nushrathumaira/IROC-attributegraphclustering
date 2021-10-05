package attributenetwork;

import edu.uci.ics.jung.algorithms.generators.random.KleinbergSmallWorldGenerator;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

import java.util.Collection;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.collections15.Factory;
import subcat.datagenerator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.ArrayList;
import subcat.subcluster;

public class GraphGenerator {
    int numVertices; 
    int numCluster;
    double[][] pconnection;
    int[][] clusteraasign; 
    int edgeCounter = 0;
	
	public GraphGenerator(String fn) {
		readParam(fn);		
	}	
	
	public Graph<Integer, Integer> generateRandomGraph2Cluster(int startNode1, int numNodes1, double edgeProb1, int startNode2, int numNodes2, double edgeProb2, double betweenProb){
	    Graph<Integer, Integer> g = new UndirectedSparseGraph<Integer, Integer>();
	    int numEdgesC1 = (int) Math.abs(edgeProb1 * (numNodes1 * (numNodes1 - 1)) / 2.0);
	    int numEdgesC2 = (int) Math.abs(edgeProb2 * (numNodes2 * (numNodes2 - 1)) / 2.0);
	    int numEdgesBetween = (int) Math.abs(betweenProb * (numNodes1 * numNodes2)); 
	    int numEdge = numEdgesC1 + numEdgesC2 + numEdgesBetween;
	    for(int i = startNode1; i < startNode1 + numNodes1; i ++){
	    	if(!g.containsVertex(i))
	    		g.addVertex(i);	    		
	    }
	    for(int i = startNode2; i < startNode2 + numNodes2; i ++){
	    	if(!g.containsVertex(i))
	    		g.addVertex(i);	    		
	    }	    
	    int numVertices = g.getVertexCount();    
	  
	    int edgeCounter = 0;
	    Random r = new Random();
	    //add random Edges in Cluster 1
	    while (edgeCounter < numEdgesC1) {
	    	int i = startNode1 + r.nextInt(numNodes1);
	        int j = startNode1 + r.nextInt(numNodes1);
	        if ((!g.isNeighbor(i, j))&&(i != j)) {
	        	g.addEdge(edgeCounter, i, j);
	            edgeCounter++;
	        }
	     }
	    //add random Edges in Cluster 2      
	    while (edgeCounter < numEdgesC1+numEdgesC2) {
	    	int i = startNode2 + r.nextInt(numNodes2);
	        int j = startNode2 + r.nextInt(numNodes2);
	        if ((!g.isNeighbor(i, j))&&(i != j)) {
	        	g.addEdge(edgeCounter, i, j);
	            edgeCounter++;
	        }
	    }
	    //add random Edges between Clusters     
	    while (edgeCounter < numEdgesC1+numEdgesC2+numEdgesBetween) {
	    	int i = startNode1 + r.nextInt(numNodes1);
	        int j = startNode2 + r.nextInt(numNodes2);
	        if ((!g.isNeighbor(i, j))&&(i != j)) {
	        	g.addEdge(edgeCounter, i, j);
	        	edgeCounter++;
	        }
	    }    	
	    return g;
     }	
	
	public Graph<Integer, Integer> generateRandomGraphCluster(Graph<Integer,Integer> g,int startNode, int numNodes, double edgeProb){
		int numEdgesC = (int) Math.abs(edgeProb * (numNodes * (numNodes - 1)) / 2.0);	  	    
	    Random r = new Random();
	    //add random Edges in Cluster 1
	    int edgeCounterC = 0;
	    while (edgeCounterC < numEdgesC) {
	    	int i = startNode + r.nextInt(numNodes);
	        int j = startNode + r.nextInt(numNodes);
	        if ((!g.isNeighbor(i, j))&&(i != j)) {
	        	g.addEdge(edgeCounter, i, j);
	        	edgeCounter ++;
	        	edgeCounterC++;
	        }
	     }
	    return g;
	}
	
	public Graph<Integer, Integer> generateRandomGraphBetweenCluster(Graph<Integer,Integer> g,int startNode1, int numNodes1, int startNodes2, int numNodes2, double betweenProb){
		int edgeCounterBetween = 0;
		Random r = new Random();
		int numEdgesBetween = (int) Math.abs(betweenProb * (numNodes1 * numNodes2));
		while (edgeCounterBetween < numEdgesBetween) {
			int i = startNode1 + r.nextInt(numNodes1);
		    int j = startNodes2 + r.nextInt(numNodes2);
		    if ((!g.isNeighbor(i, j))&&(i != j)) {
		    	g.addEdge(edgeCounter, i, j);
		        edgeCounter++;
		        edgeCounterBetween ++;
		    }
		}
		return g;
	}
	
	public void readParam(String fn){
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
			String line;
			int tmp = 0;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				String[] str = line.split("\\,");
				if(tmp == 0){
					numVertices = Integer.parseInt(str[0]);
					numCluster = Integer.parseInt(str[1]);					
				} else if(tmp == 1){
					pconnection = new double[numCluster][numCluster];
					for(int i=0; i < numCluster; i++){
						for(int j = 0; j < numCluster; j++){
							pconnection[i][j] = Double.parseDouble(str[i*numCluster+j]);
						}
					}						
				} else {
					clusteraasign = new int[numCluster][2];	
					for(int i =0; i < numCluster; i ++){
						clusteraasign[i][0] =  Integer.parseInt(str[i*2]);
						clusteraasign[i][1] =  Integer.parseInt(str[i*2+1]);
					}					
				}
				tmp ++;
			}
		} catch (IOException ex) {
			ex.printStackTrace();
		}		
	}
	
	public Graph<Integer,Integer> generateGraph(){
		Graph<Integer,Integer> g = new UndirectedSparseGraph<Integer, Integer>();		
		for(int i = 0; i < numVertices; i ++){
			if(!g.containsVertex(i))
				g.addVertex(i);	    		
		}
		for(int i = 0; i < numCluster; i ++){
			for(int j = 0; j < numCluster; j ++){
				if(i == j){
					double edgeProb = pconnection[i][j];
					int startNode = clusteraasign[i][0];
					int numNodes = clusteraasign[i][1];
					g = generateRandomGraphCluster(g,startNode, numNodes, edgeProb);
				}
				if(i != j){
					double betweenProb = pconnection[i][j];
					int startNode1 = clusteraasign[i][0];
					int numNodes1 = clusteraasign[i][1];
					int startNode2 = clusteraasign[j][0];
					int numNodes2 = clusteraasign[j][1];
					g = generateRandomGraphBetweenCluster(g,startNode1, numNodes1, startNode2, numNodes2, betweenProb);
				}
			}
		}
		return g;
	}
}
