package AttributeGraph;

	import java.io.BufferedReader;
	import java.io.File;
	import java.io.FileReader;
	import java.io.IOException;
	import java.util.ArrayList;

import attributenetwork.IO;
import attributenetwork.attributenetwork;

	import subcat.subcluster;
import subcat.utilities;

	public class Evaluation {
		attributenetwork an;
		int numCluster;
		int numVertics;
		int numFeatures;
		public ArrayList<ArrayList<Integer>> trueclus;
		public ArrayList<ArrayList<Integer>> findclus;
		public ArrayList<ArrayList<Integer>> truesubspace;
		public ArrayList<ArrayList<Integer>> findsubspace;
		
		public Evaluation(){
			trueclus = new ArrayList<ArrayList<Integer>>();
			findclus = new ArrayList<ArrayList<Integer>>();
			truesubspace = new ArrayList<ArrayList<Integer>>();
			findsubspace = new ArrayList<ArrayList<Integer>>();
		}
		
		// Read acquired cluster label
		public void readLabel(String fn){	
			String Flabel = fn + "_label.txt";
//			String Flabel = fn;
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(Flabel)));
				String line;
				while ((line = br.readLine()) != null) {
					line = line.trim();
//					String[] str = line.split("\\s+");		
					String[] str = line.split(",");	
					ArrayList<Integer> cluselement = new ArrayList<Integer>();
					for(int i=0; i<str.length; i++)
						cluselement.add(Integer.parseInt(str[i]));				
					findclus.add(cluselement);
				}
			} catch (IOException ex) {
				ex.printStackTrace();
			}			
		}
		
		//Read true labels
		public void readTure(String fn){
			String FTlabel = fn + "_truelabel.txt";
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(FTlabel)));
				String line;
				int tmp = 0;				
				while ((line = br.readLine()) != null) {
					line = line.trim();
					String[] str = line.split("\\s+");
//					String[] str = line.split(",");
					if(tmp == 0){
						numVertics = Integer.parseInt(str[0]);
						numCluster = Integer.parseInt(str[1]);
					}
					if(tmp > 0){
						if(str.length < 3){
							int startP = Integer.parseInt(str[0]);
							int numP = Integer.parseInt(str[1]);					
							ArrayList<Integer> tclus = new ArrayList<Integer>();
							for(int i = startP; i < startP + numP; i ++)
								tclus.add(i);
							trueclus.add(tclus);
						}
						if(str.length > 3){
							ArrayList<Integer> tclus = new ArrayList<Integer>();
							for(int i = 0; i < str.length; i +=2){
								int startP = Integer.parseInt(str[i]);
								int numP = Integer.parseInt(str[i+1]);								
								for(int j = startP; j < startP + numP; j ++)
									tclus.add(j);
							}
							trueclus.add(tclus);
						}
					}
					tmp ++;
				}
				
					
			} catch (IOException ex) {
				ex.printStackTrace();
			}	
		}
		
		

		//Read true labels
		public void readTureSubspace(String fn){
			String FTlabel = fn + "_featuretruelabel.txt";
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(FTlabel)));
				String line;
				int tmp = 0;				
				while ((line = br.readLine()) != null) {
					line = line.trim();
					String[] str = line.split("\\s+");
					if(tmp == 0){
						numFeatures = Integer.parseInt(str[0]);						
					}
					if(tmp > 0){
						int startF = Integer.parseInt(str[0]);
						int numF = Integer.parseInt(str[1]);					
						ArrayList<Integer> tclus = new ArrayList<Integer>();
						for(int i = startF; i < startF + numF; i ++)
							tclus.add(i);
						truesubspace.add(tclus);
					}
					tmp ++;
				}
				
					
			} catch (IOException ex) {
				ex.printStackTrace();
			}	
		}
		
		public void FmeasureDBCSCsubspace(){
			double Fm = 0.0;
			IO ea = new IO();			
			ArrayList<ArrayList<Integer>> dbsubspace = ea.readsubspace("D:/AttributedGraph/src/DB/newnew2clusterGraph", numFeatures);
			ArrayList<ArrayList<Integer>> dbcluster = ea.readDBCSClabel("D:/AttributedGraph/src/DB/newnew2clusterGraph", numFeatures);
			for(int i = 0; i < dbcluster.size(); i ++){
				ArrayList<Integer> clus = dbcluster.get(i);
				int[] label = new int[trueclus.size()];
				for(Integer j : clus){
					for(int k = 0; k < trueclus.size(); k ++){
						if(trueclus.get(k).contains(j))
							label[k]++;
					}
				}
				int maxlabel = 0;
				int maxvalue = 0;
				for(int l = 0; l < label.length; l ++){
					if(label[l] > maxvalue){
						maxvalue = label[l];
						maxlabel = l;
					}
				}
				double fmeasurec = 0.0;
				ArrayList<Integer> truelist = truesubspace.get(maxlabel);
				ArrayList<Integer> findlist = dbsubspace.get(i);
				double precision = 0.0;
				double recall = 0.0;
				double fmeasure = 0.0;
				
				int TP = 0;
				int TN = 0;
				int FP = 0;
				int FN = 0;
				for(int a = 0; a < numFeatures; a ++){
					for(int b = a + 1; b < numFeatures; b ++){
						boolean flagF  = false;						
						if(findlist.contains(a) && findlist.contains(b))
							flagF = true;						
						boolean flagT = false;						
						if(truelist.contains(a) && truelist.contains(b))
							flagT = true;						
						if(flagF && flagT)
							TP ++;
						if(!flagF && !flagT)
							TN ++;
						if(!flagF && flagT)
							FN ++;
						if(flagF && !flagT)
							FP ++;
					}
				}
				if((TP + FP) > 0)
					precision = (double)TP/(double)(TP + FP);
				if((TP + FN) > 0)
					recall = (double)TP/(double)(TP + FN);
				double beta = 1;
				if(precision > 0 &&  recall > 0)
				fmeasure = (double)(1+beta*beta)*precision*recall/(double)(beta*beta*precision + recall);
				Fm += fmeasure;
			}
			System.out.println("Fmeasure:  " + Fm/(double)dbsubspace.size());
			System.out.println("Cluster Size:  " + dbsubspace.size());
		}
		
		public void Fmeasure(){
			double precision = 0.0;
			double recall = 0.0;
			double fmeasure = 0.0;
			
			int TP = 0;
			int TN = 0;
			int FP = 0;
			int FN = 0;
			for(int i = 0; i < numVertics; i ++){
				for(int j = i + 1; j < numVertics; j ++){
					boolean flagF  = false;
					for(int m = 0; m < findclus.size(); m ++){
						if(findclus.get(m).contains(i) && findclus.get(m).contains(j))
							flagF = true;
					}
					boolean flagT = false;
					for(int m = 0; m < trueclus.size(); m ++){
						if(trueclus.get(m).contains(i) && trueclus.get(m).contains(j))
							flagT = true;
					}
					if(flagF && flagT)
						TP ++;
					if(!flagF && !flagT)
						TN ++;
					if(!flagF && flagT)
						FN ++;
					if(flagF && !flagT)
						FP ++;
				}
			}
			
			precision = (double)TP/(double)(TP + FP);
			recall = (double)TP/(double)(TP + FN);
			double beta = 1;
			fmeasure = (double)(1+beta*beta)*precision*recall/(double)(beta*beta*precision + recall);
			System.out.println("Precision:  " + precision);
			System.out.println("Recall:  " + recall);
			System.out.println("Fmeasure:  " + fmeasure);
		}
		
		public void Modularity(){
			attributenetwork an = new attributenetwork();
			IO ea = new IO();
			an = ea.readNetwork("D:/research/AttributedGraph/data/real/0");
			Evaluation ev = new Evaluation();
			ev.readLabel("D:/research/AttributedGraph/data/real/0_result.txt"); 
			double modu = 0.0;
			for(int i = 0; i < findclus.size(); i ++){
				double fkk = 0.0;
				double flk = 0.0;
				ArrayList<Integer> clusi = findclus.get(i);
				int edgein = 0;
				int edgebetween = 0;
				for(Integer v1 : clusi){
					for(Integer v2 : clusi){
						if(v1 != v2 && an.exsitEdge(v1, v2))
							edgein ++;
					}
					for(int j = 0 ; j < findclus.size(); j ++){
						if(i!=j){
							for(Integer v3 : findclus.get(j)){
								if(v1 != v3 && an.exsitEdge(v1, v3)){
									edgebetween ++;
								}
							}
						}
						flk += (double)edgebetween/(2*(double)an.numEdges);
					}
				}
				fkk = (double)edgein/(double)an.numEdges;
				
				modu += fkk - Math.pow(flk, 2);
			}
			System.out.println("Modularity:  " + modu);
		}
		
		public static void main(String[] args){	
			Evaluation ev = new Evaluation();
			ev.readTure("D:/AttributedGraph/data/syn/ResultNew/newnew3clusterGraph");
			
//			ev.readLabel("D:/AttributedGraph/107");
//			ev.readTureSubspace("D:/AttributedGraph/data/syn/ResultNew/newnew2clusterGraph");
			IO ea = new IO();
			ev.findclus = ea.readBNMTFlabel("D:/AttributedGraph/data/syn/BNMTF");     
//			ev.findclus = ea.readBAGClabel("D:/AttributedGraph/data/real/BAGC_107");
//			ev.findclus = ea.readPICSlabel("D:/AttributedGraph/data/realgp/PICSgp7");
//			ev.findclus = ea.readCongaLabel("D:/AttributedGraph/data/syn/ResultNew/new3clusterGraph");
//			ev.findclus = ea.readDBCSClabel("src/DB/gp53",7);
//			ev.readLabel(fn);
//			ev.readTure("data/syn/Result/5clusterGraph");
//			ev.trueclus = ea.readFacebookLabel("D:/AttributedGraph/data/facebook/107");
//			ev.trueclus = ea.readGoogleplusLabel("D:/AttributedGraph/data/googleplus/gplus/gp53");
//			ev.numVertics = 1084;
			ev.Fmeasure();
//			ev.FmeasureDBCSCsubspace();
//			ev.Modularity();
		}
	}