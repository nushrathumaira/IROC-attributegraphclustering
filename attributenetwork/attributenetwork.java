package attributenetwork;

import java.util.HashSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;

import attributenetwork.edge;
import attributenetwork.attributenetwork;

import attributenetwork.vertex;

public class attributenetwork {
	 public int numVertices;
	 public int numFeature;	
	 public int[] numCategory;
	 public int numEdges;
     public HashMap<Integer, vertex> vertexMap;
     public ArrayList<HashMap<Integer,Double>> pcategory;
     public int[][] featureAdj;
     
     public attributenetwork(){ 
    	 numCategory = new int[numFeature];
    	 vertexMap = new HashMap<Integer, vertex>();
    	 pcategory = new ArrayList<HashMap<Integer,Double>>();
    	 featureAdj = new int[numVertices][numFeature];
 	 }   

	/**
	 * Obtain vertex id
	 * @param id
	 * @return id of vertex
	 */
	public vertex getVertex(int id){
		assert(vertexMap.containsKey(id));
		return vertexMap.get(id);
	}
	
	/**
	 * Add a vertex
	 * @param v
	 */
	public void addVertex(vertex v){
		vertexMap.put(v.id, v);
		numVertices ++;	
	}
	
	/**
	 * Add vertex with id
	 * @param id
	 */
	public void addVertex(int id){
		vertex v = new vertex(id);
		vertexMap.put(id, v);
		numVertices ++;				
	}
	
	/**
	 * judge whether the edge is exist
	 * @param i
	 * @param j
	 * @return
	 */
	public boolean exsitEdge(int i, int j){
		vertex v1 = vertexMap.get(i);
		vertex v2 = vertexMap.get(j);
		return v1.existsEdge(v2);
	}
	
	/**
	 * Add an edge
	 * @param i
	 * @param j
	 */	
	public void addEdge(int i, int j){
		vertex v1 = vertexMap.get(i);
		vertex v2 = vertexMap.get(j);
		v1.addEdge(v2);
		v2.addEdge(v1);
		numEdges ++;
	}
	
	/**
	 * Remove an edge
	 * @param i
	 * @param j
	 */
	public void removeEdge(int i, int j){
		vertex v1 = vertexMap.get(i);
		vertex v2 = vertexMap.get(j);
		v1.removeEdge(new edge(v1,v2));
		v2.removeEdge(new edge(v2,v1));
		numEdges --;
	}
	
	//
	public void removeVertex(vertex v){
		for(edge e : v.edges){
			vertex neighbor = vertexMap.get(e.end.id);
			neighbor.removeEdge(new edge(neighbor,v));
			numEdges --;
		}
		
		vertexMap.remove(v.id);
		numVertices --;	
	}
	
//	public void mergeSameV(vertex v, HashSet<vertex> cluster){
//		for(Iterator<vertex> it=cluster.iterator(); it.hasNext(); ){
//			vertex merged = it.next();
//			v.vertices.add(merged.id);
//			removeVertex(merged);
//		}
//	}
	
	public ArrayList<vertex> commonNeighbor(vertex v1, vertex v2){
		ArrayList<vertex> cn = new ArrayList<vertex>();
		
		ArrayList<vertex> nv1 = v1.getNeighbor();
		ArrayList<vertex> nv2 = v2.getNeighbor();
		for(vertex v : nv1){
			if(nv2.contains(v))
				cn.add(v);
		}
		
		return cn;
	}
	
	public double featuresimilarity(vertex v1, vertex v2){
		double fs = 0;
		for (int i=0; i<numFeature; i++){
			if(v1.feature.get(i) == v2.feature.get(i))
				fs ++;	
		}
		fs = fs / numFeature;
	    return fs;		
	}

}
