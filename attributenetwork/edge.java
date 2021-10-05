package attributenetwork;

import attributenetwork.edge;
import attributenetwork.vertex;

public class edge implements Comparable<edge>{
//	public int weight;
	public vertex start, end;
	
    public edge(vertex start, vertex end)
    {
    	this.start = start;
    	this.end = end;
//    	this.weight = weight;
    }
	
    public boolean equals(Object o)
    {
        edge e = (edge) o;
        if((e.start.id == start.id) && (e.end.id == end.id))
        	return true;
        if((e.start.id == end.id) && (e.end.id == start.id))
        	return true;
        return false;
    }
    
    public int compareTo(edge e) 
    {
        //Edge e = (Edge) o;
		if(e.end == null)	
			System.out.println();
        if(end.id > e.end.id)
            return 1;
        else if(end.id < e.end.id)
            return -1;
        return 0;
    }
}
