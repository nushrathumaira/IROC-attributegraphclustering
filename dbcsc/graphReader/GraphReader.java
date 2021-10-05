package dbcsc.graphReader;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;
import edu.uci.ics.jung.io.GraphMLReader;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.MapIterator;
import org.apache.commons.collections15.Transformer; 

import java.util.Iterator;
import java.util.ArrayList;
import java.util.Collections;

public class GraphReader {
 
       
    @SuppressWarnings("unchecked")
	public static dbcsc.graph.Graph loadGraphFromGraphMLFile(String filename) {

    	dbcsc.graph.Graph myGraph = null;
        //Leeren JUNG-Graphen erzeugen
        UndirectedSparseGraph<dbcsc.graph.Node, ALink> g = new UndirectedSparseGraph<dbcsc.graph.Node, ALink>();

        try {

            //Knoten- und Kantenfactories erzeugen
            Factory<dbcsc.graph.Node> vertexFactory = new Factory<dbcsc.graph.Node>() {

                @Override
				public dbcsc.graph.Node create() {
                	dbcsc.graph.Node node = new dbcsc.graph.Node();

                    return node;
                }
            };

            Factory<ALink> edgeFactory = new Factory<ALink>() {

                @Override
				public ALink create() {
                    ALink link = new ALink();

                    return link;
                }
            };

            //GraphML-File lesen
            GraphMLReader<Graph<dbcsc.graph.Node, ALink>, dbcsc.graph.Node, ALink> reader = new GraphMLReader<Graph<dbcsc.graph.Node, ALink>, dbcsc.graph.Node, ALink>(vertexFactory, edgeFactory);
            reader.load(filename, g);
            
            myGraph = new dbcsc.graph.Graph(reader.getVertexMetadata().size());
            //Die Menge der Attribute
            ArrayList<String> attributeSet = new ArrayList<String>(reader.getVertexMetadata().keySet());
            Collections.sort(attributeSet);
            
            //Mit XML-Tranformern die Attribute der Knoten auslesen
            Transformer<dbcsc.graph.Node,String>[] transformers = new Transformer[attributeSet.size()];
            int count = 0;
            for(Iterator<String> it = attributeSet.iterator();it.hasNext();) {
            	String key = it.next();
            	transformers[count] = reader.getVertexMetadata().get(key).transformer;
            	count++;
            }
            
            
            //Iteriere über Knoten und lese Attribute
            MapIterator<dbcsc.graph.Node, String> iVertices = reader.getVertexIDs().mapIterator();

            while (iVertices.hasNext()) {
            	dbcsc.graph.Node node = iVertices.next();
                
                node.setID(Integer.parseInt(iVertices.getValue()));
                //alle Attribute einlesen und zum Knoten hinzufuegen
                double [] atts = new double[attributeSet.size()];
                for(int i = 0;i<attributeSet.size();i++) {
                	if(transformers[i].transform(node)==null){
                		//Nicht vorhandene Attribute werden auf NaN gesetzt
                		atts[i]=Double.NaN;
                	} else {
                		atts[i]=Double.parseDouble(transformers[i].transform(node));
                	}
                }       
                node.setAttributes(atts);
                myGraph.addNode(node);
            }


            //Kanten einlesen
            for (ALink current : g.getEdges()) {
                Pair<dbcsc.graph.Node> nodes = g.getEndpoints(current);

                myGraph.addEdge((nodes.getFirst()), (nodes.getSecond()));

                
            }


        } catch (Exception e) {
            if (e.getMessage() == null)
                System.out.print("Error while parsing GraphML-File: Some of the edges and/or vertices do not have appropriate attributes defined. This case is unhandled by this parse for now.\n");
            else
                System.out.print("Error while parsing GraphML-File: " + e.getMessage() + "\n");
        }
        
        //Knoten in richtiger Reihenfolge ausgeben (nach ID sortiert)
        NodeComparator nodeComparator = new NodeComparator();
        Collections.sort(myGraph.getNodes(), nodeComparator);
        
        return myGraph;
    }
    
}
