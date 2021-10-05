package dbcsc.graph;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;


public class Node implements Comparable<Node>{

    private int ID;
    private double[] attributes;
    private HashSet<Node> neighbors = new HashSet<Node>();

    public int getID() {
        return ID;
    }

    public void setID(int ID) {
        this.ID = ID;
    }

    public Node copyWithoutNeighbours() {
        Node neu = new Node(ID, attributes);
        neu.neighbors = new HashSet<Node>();
        return neu;
    }

    public void setAttribute(int key, double value) {
        attributes[key] = value;
    }

    public void setAttributes(double[] values) {
        attributes = values;
    }

    public double getAttribute(int key) {
        return attributes[key];
    }

    public double[] getAttributes() {
        return attributes;
    }

    public HashSet<Node> getNeighbors() {
        return neighbors;
    }

    public void addNeighbor(Node x) {
        neighbors.add(x);
    }

    @Override
	public String toString() {
        return Integer.toString(ID);
    }

    public Node(int myID, double[] atts) {
        ID = myID;
        attributes = atts;
    }

    public Node() {
        ID = -1;
    }

    public Node(int ID) {
        this.ID = ID;
    }



  

    //k-neighborhood ohne Subspace-Pruning
    public HashSet<Node> k_neighborhood(Collection<Node> cands, int k) {
        HashSet<Node> k_neighbors = new HashSet<Node>();
        k_neighborhood_rek(new HashSet<Node>(cands), k, k_neighbors);
        k_neighbors.remove(this);
        return k_neighbors;
    }
    
    //k-neighborhood ohne Subspace-Pruning (cands wird verändert!!!)

    public void k_neighborhood_rek(HashSet<Node> cands, int k, HashSet<Node> k_neighbors) {
        //Tiefensuche
        if (cands.contains(this)) {
            k_neighbors.add(this);
            cands.remove(this);
        }
        if (k <= 0 || cands.isEmpty()) {
            return;
        }
        //Es werden ALLE Nachbarn betrachtet, der Pfad zu einem Knoten aus cands muss nicht in cands liegen!
        for (Node node : this.neighbors) {
            node.k_neighborhood_rek(cands, k - 1, k_neighbors);
        }
    }

  
    public HashSet<Node> connected_component(Collection<Node> cands) {
        if (cands.isEmpty()) {
            return new HashSet<Node>();
        }
        //Breitensuche
        Node firstnode = cands.iterator().next();
        HashSet<Node> visited = new HashSet<Node>();
        LinkedList<Node> queue = new LinkedList<Node>();
        queue.add(firstnode);
        while (!queue.isEmpty()) {
            Node node = queue.poll();
            visited.add(node);
            for (Node neighbor : node.getNeighbors()) {
                if (cands.contains(neighbor) && !visited.contains(neighbor) && !queue.contains(neighbor)) {
                    queue.add(neighbor);
                }
            }
        }
        return visited;
    }


    public int getDegree(){
        return getNeighbors().size();
    }

    @Override
	public int compareTo(Node o) {
        return this.getID() - o.getID();
    }
    
}
