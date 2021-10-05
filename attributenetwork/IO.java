package attributenetwork;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import AttributeGraph.Cluster;

import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLDouble;

import edu.uci.ics.jung.graph.Graph;

public class IO {	
	
	public attributenetwork readNetwork(String file){
		String networkFile = file + "_network.txt";
//		String featureFile = file + "_newfeature.txt";		
		String featureFile = file + "_feature.txt";
		
		System.out.println("Read network from file: "+featureFile);	
		
		attributenetwork an = new attributenetwork();		
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(featureFile)));
			String line;
			int num = 0;
			while((line=br.readLine())!=null){
				line = line.trim();
				String[] str = line.split(",");
//				String[] str = line.split("\\s+");
				ArrayList<Integer> feature = new ArrayList<Integer>();
				an.numFeature = str.length;
				for(int i=0; i<str.length; i++)
					feature.add(Integer.parseInt(str[i]));				
				an.addVertex(num);
				an.getVertex(num).feature = feature;
				num++;
			}
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }
		
		System.out.println("Read network from file: "+networkFile);	

	    try {
	    	BufferedReader br = new BufferedReader(new FileReader(new File(networkFile)));
	    	String line;
	    	int i=0;
	    	int numEdges = 0;
	    	while((line=br.readLine())!=null){
	    		line = line.trim();        		
	    		String[] str = line.split(",");
	    		an.numVertices = str.length;
	    		for(int j=0; j<str.length; j++){
	    			double data = Double.parseDouble(str[j]);
	    			if((data == 1)&&(!an.exsitEdge(i, j))){
	    				an.addEdge(i, j);	
	    				numEdges ++;
	    			}
	    		}    		
	    		i++;
	    	}
	    	an.numEdges = numEdges;
	    } catch (IOException ex) {
	        ex.printStackTrace();
	    }	
	    
	    for(Integer i : an.vertexMap.keySet()){
	    	vertex v = an.vertexMap.get(i);
	    	v.getNeighbor();
	    }
	    
	    int[][] Feature = new int[an.numVertices][an.numFeature];
	    for (int i = 0; i < an.numVertices; i ++){
	    	for (int j =0 ; j < an.numFeature; j ++){
	    		Feature[i][j] = an.getVertex(i).feature.get(j);
	    	}
	    }
	    an.featureAdj = Feature;
	    int[] numCategory = new int[an.numFeature];
	    for(int i = 0; i < an.numFeature; i ++){
	    	ArrayList<Integer> Category = new ArrayList<Integer>();
	    	for(int j = 0; j < an.numVertices; j ++){
	    		if(! Category.contains(Feature[j][i])){
	    			Category.add(Feature[j][i]);
	    		}
	    	}
	    	numCategory[i] = Category.size();
	    	
	    	HashMap<Integer, Double> ponecategory = new HashMap<Integer, Double>();	    	
	    	for(Integer k : Category){
	    		int catenum = 0;
	    		for(int j = 0; j < an.numVertices; j ++){	    			
	    			if(Feature[j][i] == k){
	    				catenum ++;
	    			}	    				
	    		}
	    		double pcate = (double)catenum/(double)an.numVertices;
	    		ponecategory.put(k, pcate);
	    	}
	    	an.pcategory.add(ponecategory);
	    }
	    an.numCategory = numCategory;	   
		return an;
	}
	
	public attributenetwork readGraphml(String fn){
		attributenetwork an = new attributenetwork();
		ArrayList<String> nodeId = new ArrayList<String>();
		ArrayList<ArrayList<String>> nodeFeature = new ArrayList<ArrayList<String>>();
		ArrayList<String> allFeature = new ArrayList<String>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn)));
			String line;
			int num = 0;
			while((line=br.readLine())!=null){
				line = line.trim();
				if(line.startsWith("<node")){
					String[] str = line.split("\\s+");
					int start = str[1].indexOf("\"");
					int end = str[1].lastIndexOf("\"");
					nodeId.add(str[1].substring(start+1, end));
					an.addVertex(num++);
					ArrayList<String> feature = new ArrayList<String>();
					for(int i=2; i<str.length-1; i++){
						start = 0;
						end = str[i].indexOf("=");
						feature.add(str[i].substring(start, end));
						if(!allFeature.contains(str[i].substring(start, end)))
							allFeature.add(str[i].substring(start, end));
					}
					nodeFeature.add(feature);
				}
				if(line.startsWith("<edge")){
					String[] str = line.split("\\s+");
					int start1 = str[1].indexOf("\"");
					int end1 = str[1].lastIndexOf("\"");
					int start2 = str[2].indexOf("\"");
					int end2 = str[2].lastIndexOf("\"");
					an.addEdge(nodeId.indexOf(str[1].substring(start1+1, end1)), nodeId.indexOf(str[2].substring(start2+1, end2)));
				}
			}
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }
		for(int i=0; i<nodeId.size(); i++){
			ArrayList<Integer> feat = new ArrayList<Integer>();
			ArrayList<String> feature = nodeFeature.get(i);
			for(int j=0; j<allFeature.size(); j++){
				if(feature.contains(allFeature.get(j)))
					feat.add(1);
				else
					feat.add(0);
			}
			an.getVertex(i).feature = feat;
		}
		an.numVertices = nodeId.size();
		an.numFeature = allFeature.size();
		
		int[][] Feature = new int[an.numVertices][an.numFeature];
		    for (int i = 0; i < an.numVertices; i ++){
		    	for (int j =0 ; j < an.numFeature; j ++){
		    		Feature[i][j] = an.getVertex(i).feature.get(j);
		    	}
		    }	  
		    int[] numCategory = new int[an.numFeature];
		    for(int i = 0; i < an.numFeature; i ++){
		    	ArrayList<Integer> Category = new ArrayList<Integer>();
		    	for(int j = 0; j < an.numVertices; j ++){
		    		if(! Category.contains(Feature[j][i])){
		    			Category.add(Feature[j][i]);
		    		}
		    	}
		    	numCategory[i] = Category.size();
		    }
		    an.numCategory = numCategory;
		return an;
	}
	
	public void convertToGraphml(attributenetwork an, String fn){
		 try{
			 FileOutputStream fout = new FileOutputStream(new File(fn));	
			 fout.write(("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n").getBytes());
			 fout.write(("<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\r\n").getBytes());
			 fout.write(("xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \r\n").getBytes());
			 fout.write(("xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns \r\n").getBytes());
			 fout.write(("http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\"> \r\n").getBytes());
			 fout.write(("<graph id=\"G\" edgedefault=\"undirected\">\r\n").getBytes());
			 for(int id = 0; id<an.numVertices; id++){
				 vertex v = an.getVertex(id);
				 fout.write(("\t<node id=\""+id+"\" ").getBytes());
				 for(int i=0; i<v.feature.size(); i++)
					 fout.write(("att"+i+"="+"\""+v.feature.get(i)+"\" ").getBytes());
				 fout.write(("/>\r\n").getBytes());
			 }
			 for(int id1=0; id1<an.numVertices; id1++){
				 for(int id2=id1+1; id2 < an.numVertices; id2++){
					 if(an.exsitEdge(id1, id2)){
						 fout.write(("\t<edge source=\""+id1+"\" ").getBytes());
						 fout.write(("target=\""+id2+"\" />\r\n").getBytes());
					 }
				 }
			 }
			 fout.write(("</graph> \r\n").getBytes());
			 fout.write(("</graphml> \r\n").getBytes());
			 } catch (IOException ex) {
				 ex.printStackTrace();
			 } 
	}
	
	
	 public void writeGraph(Graph<Integer, Integer> g, String fn) {		 
		 int numObj = g.getVertexCount();
		 int[][] data = new int[numObj][numObj];
		 try{
				FileOutputStream fout = new FileOutputStream(new File(fn));
				for(int i=0;i<numObj;i++){
					for(int j=0;j<numObj;j++){
						if (g.isNeighbor(i, j)) {
		                    data[i][j] = 1;
		                } else {
		                    data[i][j] = 0;
		                }
						fout.write((data[i][j]+",").getBytes());
					}
					fout.write(("\r\n").getBytes());
				}			
			} catch (IOException ex) {
			    ex.printStackTrace();
			} 
	 }
	 
	 public void writeSparseGraph(Graph<Integer, Integer> g, String fn){
		 int numObj = g.getVertexCount();		
		 try{
			 FileOutputStream fout = new FileOutputStream(new File(fn));
			 for(int i=0;i<numObj;i++){
				for(int j=0;j<numObj;j++){
					if (g.isNeighbor(i, j)) {					
						fout.write((i +"," + j + "," + 1).getBytes());
						fout.write(("\n").getBytes());
					}
				}				
			 }			
		 } catch (IOException ex) {
			 ex.printStackTrace();
		 } 
	 }
	 
	 public void writeDemonNode(int numVertices, String fn){			
		 try{
			 FileOutputStream fout = new FileOutputStream(new File(fn));
			 for(int i = 0; i < numVertices; i ++){									
				fout.write((i+"").getBytes());
				fout.write(("\n").getBytes());						
			 }			
		 } catch (IOException ex) {
			 ex.printStackTrace();
		 } 
	 }
	  
	 public void writeDemonEdge(attributenetwork an, String fn){			
		 try{
			 FileOutputStream fout = new FileOutputStream(new File(fn));
			 for(int i = 0; i < an.numVertices; i ++){	
				 for(int j = 0; j < an.numVertices; j ++){
					 if(an.exsitEdge(i, j)){
						 fout.write((i +"," + j + "," + 1).getBytes());
						 fout.write(("\n").getBytes());
					 }				
				 }								
			 }			
		 } catch (IOException ex) {
			 ex.printStackTrace();
		 } 
	 }
	 
	 public void convertGraphgmlresult(String file){
		 String resultFile = file + "_gamer.found";			
			attributenetwork an = new attributenetwork();		
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(resultFile)));
				String line;
				int num = 0;
				while((line=br.readLine())!=null){
					line = line.trim();
					String[] str = line.split(",");
					ArrayList<Integer> feature = new ArrayList<Integer>();
					an.numFeature = str.length;
					for(int i=0; i<str.length; i++)
						feature.add(Integer.parseInt(str[i]));				
					an.addVertex(num);
					an.getVertex(num).feature = feature;
					num++;
				}
			} catch (IOException ex) {
		        ex.printStackTrace();
		    } 
	 }

	 
	 public void writeNetwork(String file, ArrayList<edge> edges){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_network.txt"));
			for(edge e : edges){
				fout.write((e.start.id+"\t").getBytes());
				fout.write((e.end.id+"").getBytes());
				fout.write(("\r\n").getBytes());
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}	

	
	public void writeResult(String file, int [] label){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_result.txt"));
			for(int i=0; i<label.length ; i++){
				fout.write((label[i]+"\t").getBytes());
				fout.write(("\r\n").getBytes());
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}
	
	public void writeLabel(String file, ArrayList<TreeSet<Integer>> IDV){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_label.txt"));
			for(int i = 0; i < IDV.size(); i ++){
				TreeSet<Integer> clusid = IDV.get(i);
				for(Integer id : clusid)				
					fout.write((id+",").getBytes());
				fout.write(("\r\n").getBytes());				
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}	
	
	public void writesubspace(String file, ArrayList<TreeSet<Integer>> IDF){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_featurelabel.txt"));
			for(int i = 0; i < IDF.size(); i ++){
				TreeSet<Integer> featureid = IDF.get(i);
				for(Integer id : featureid)				
					fout.write((id+",").getBytes());
				fout.write(("\r\n").getBytes());				
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}
	
	public ArrayList<Cluster> readClusterFromStructure(String structureFile, attributenetwork an){
		ArrayList<Cluster> ClusList = new ArrayList<Cluster>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(structureFile+"_structure.txt")));
			String line;
			int numLine = 0;
			while((line=br.readLine())!=null){
				if(numLine == 2){
					line = line.trim();
					String [] str = line.split(",");
					for(int i=0; i<str.length; i+=2){
						Cluster clus = new Cluster(an);
						int start = Integer.parseInt(str[i]);
						int numObj = Integer.parseInt(str[i+1]);
						for(int j=start; j<start+numObj; j++)
							clus.Nodes.add(an.getVertex(j));
						ClusList.add(clus);
					}
				}
				numLine ++;
			}
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }
		
		return ClusList;
	}
	

	
	public ArrayList<ArrayList<Integer>> readFacebookLabel(String fn){
		ArrayList<ArrayList<Integer>> facebookLabel = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+".circles")));
			String line;		
			while((line=br.readLine())!=null){	
				line = line.trim();
				String [] str = line.split("\\s+");
				ArrayList<Integer> circle = new ArrayList<Integer>();
				for(int i=1; i<str.length; i++){			
					circle.add(Integer.parseInt(str[i]));
				}			
				facebookLabel.add(circle);
			}
			
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+".feat")));
			String line;
			ArrayList<Integer> trueNodes = new ArrayList<Integer>();
			while((line=br.readLine())!=null){	
				line = line.trim();
				String [] str = line.split("\\s+");				
				trueNodes.add(Integer.parseInt(str[0]));				
			}
			for(int i = 0; i < trueNodes.size(); i ++){
				for(int j = 0; j < facebookLabel.size(); j ++){
					ArrayList<Integer> flabel = facebookLabel.get(j);
					for(Integer fl : flabel){
						if(fl == trueNodes.get(i))
							fl = i;
					}
				}
			}
			
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }
		return facebookLabel;
	}
	
	public ArrayList<ArrayList<Integer>> readPICSlabel(String fn){
		ArrayList<ArrayList<Integer>> label = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+"_label.txt")));
			String line;		
			while((line=br.readLine())!=null){	
				line = line.trim();
				String [] str = line.split("\\s+");
				ArrayList<Integer> clusterlabel = new ArrayList<Integer>();
				for(int i=0; i<str.length; i++){	
					if(! clusterlabel.contains(Integer.parseInt(str[i])))
						clusterlabel.add(Integer.parseInt(str[i]));					
				}	
				for(Integer l : clusterlabel){
					ArrayList<Integer> clus = new ArrayList<Integer>();
					for(int i=1; i<str.length; i++){
						if(Integer.parseInt(str[i])== l){
							clus.add(i);
						}
					}
					label.add(clus);
				}				
			}
			
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }		
		return label;
	}
	
	public ArrayList<ArrayList<Integer>> readBAGClabel(String fn){
		ArrayList<ArrayList<Integer>> label = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+"_label.txt")));
			String line;	
			ArrayList<Integer> alllabel = new ArrayList<Integer>();
			while((line=br.readLine())!=null){	
				line = line.trim();
				alllabel.add(Integer.parseInt(line));
			}
			ArrayList<Integer> clusterlabel = new ArrayList<Integer>();
			for(int i=1; i<alllabel.size(); i++){	
				if(! clusterlabel.contains(alllabel.get(i)))
					clusterlabel.add(alllabel.get(i));					
			}	
			for(Integer l : clusterlabel){
					ArrayList<Integer> clus = new ArrayList<Integer>();
					for(int i=1; i<alllabel.size(); i++){
						if(alllabel.get(i)== l){
							clus.add(i);
						}
					}
					label.add(clus);
			}				
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }		
		return label;
	}
	
	public ArrayList<ArrayList<Integer>> readDBCSClabel(String fn, int numFeatures){
		ArrayList<ArrayList<Integer>> label = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+"_dbcsc.found")));
			String line;		
			while((line=br.readLine())!=null){	
				line = line.trim();
				String [] str = line.split("\\s+");
				ArrayList<Integer> clus = new ArrayList<Integer>();
				for(int i= numFeatures +1; i<str.length; i++){	
					clus.add(Integer.parseInt(str[i]));					
				}		
				label.add(clus);						
			}
			
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }		
		return label;
	}
	
	public ArrayList<ArrayList<Integer>> readsubspace(String fn, int numFeatures){
		ArrayList<ArrayList<Integer>> label = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+"_dbcsc.found")));
			String line;		
			while((line=br.readLine())!=null){					
				line = line.trim();
				String [] str = line.split("\\s+");
				if(str[0] == "time")
					break;
				ArrayList<Integer> clus = new ArrayList<Integer>();
				for(int i= 0; i< numFeatures; i++){	
					if(Integer.parseInt(str[i]) == 1)
						clus.add(i);					
				}		
				label.add(clus);						
			}
			
		} catch (IOException ex) {
	        ex.printStackTrace();
	    }		
		return label;
	}
	
	public void writeAllNodes(String file, attributenetwork an){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_nodes.txt"));
			for(int i = 0; i < an.numVertices; i ++){
				fout.write((i+"").getBytes());
				fout.write(("\r\n").getBytes());
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}
	
	public void writeAllEdges(String file, attributenetwork an){
		try{
			FileOutputStream fout = new FileOutputStream(new File(file+"_edges.txt"));
			for(int i = 0; i < an.numVertices; i ++){
				for(int j = 0; j < an.numVertices; j ++){
					if(i!=j && an.exsitEdge(i, j)){
						fout.write((i+" ").getBytes());
						fout.write((j+"").getBytes());
						fout.write(("\r\n").getBytes());
					}
				}			
			}
		} catch (IOException ex) {
		    ex.printStackTrace();
		} 
	}
	
	public static void main(String[] args){
//		attributenetwork an = new attributenetwork();
//		IO ea = new IO();		
//		an = ea.readNetwork("3clusterGraph");
//		ea.convertToGraphml(an, "3clusterGraph.graphml");
		
		IO ea = new IO();
		attributenetwork an = ea.readGraphml("data/real/dblp_top.graphml");
	}
}

