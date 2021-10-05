package attributenetwork;

import subcat.datagenerator;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;


public class GenerateNetwork {
	public static void main(String[] args){
		GraphGenerator g = new GraphGenerator("data/syn/5clusterGraph_structure.txt");
		Graph<Integer,Integer> G = g.generateGraph();
		IO ea = new IO();
		ea.writeGraph(G, "data/syn/5clusterGraph_network.txt");
		ea.writeSparseGraph(G, "data/syn/5clusterGraph_sparse.txt");
		datagenerator dg = new datagenerator("data/syn/5clusterGraph_parameter.txt");
		dg.writeIntData("data/syn/5clusterGraph_feature.txt");		
		System.out.println("Network Ready");		
	}
}
