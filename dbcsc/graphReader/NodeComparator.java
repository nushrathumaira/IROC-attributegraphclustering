package dbcsc.graphReader;

import java.util.Comparator;
import dbcsc.graph.*;

//Vergleicht 2 Knoten nach ihrer ID
public class NodeComparator implements Comparator<Node> {
	
	@Override
	public int compare(Node node1,Node node2) {
		if(node1 ==null && node2 ==null) {
			return 0;
		}
		if(node1 ==null && node2 !=null) {
			return -1;
		}
		if(node1 !=null && node2 ==null) {
			return 1;
		}
		if(node1.getID() == node2.getID()) {
			return 0;
		}
		if(node1.getID() < node2.getID()){
			return -1;
		}
		return 1;
		
	}

}
