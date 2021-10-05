package attributenetwork;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Collections;

import attributenetwork.edge;
import attributenetwork.vertex;

public class vertex {
	public int id;
	public ArrayList<edge> edges;
	public ArrayList<Integer> feature;
	
	public vertex(int id){
		this.id = id;		
		this.edges = new ArrayList<edge>();
		this.feature = new ArrayList<Integer>();		
	}
	
//	public vertex clone(){
//		vertex v = new vertex(this.id);
//		
//	}
    
	public boolean existsEdge(vertex v){
		edge toAdd = new edge(this, v);
		int index = Collections.binarySearch(edges, toAdd);
		return(index >= 0);
	}
	
    public void addEdge(vertex v){
    	edge toAdd = new edge(this, v);
    	int index = Collections.binarySearch(edges, toAdd);
    	assert(index < 0);
    	if(index < 0)
    		edges.add(-index-1, toAdd);
    }
    
    public void removeEdge(edge e){
    	int index = Collections.binarySearch(edges, e);
    	assert(index >= 0);
    	if(index >= 0)
    		edges.remove(index);
    }
    
    public ArrayList<vertex> getNeighbor(){
        ArrayList<vertex> neighbors = new ArrayList<vertex>();
        for(edge e : this.edges)
        	neighbors.add(e.end);        
        return neighbors;
    }
    
    public boolean equals(Object o)
    {
    	return (this.id == ((vertex) o).id);
    }

}
