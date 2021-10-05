package DocAGnew;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class readFacebook {
	String fn;
	ArrayList<Integer> nodes;
	
	public readFacebook(String fn){
		this.fn = fn;
		this.nodes = new ArrayList<Integer>();
		
		readFeat();
		readEdges();
	}
	
	public void readFeat(){
		ArrayList<ArrayList<Integer>> featList = new ArrayList<ArrayList<Integer>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+".feat")));
			String line;
			while((line=br.readLine())!=null){
				line = line.trim();
				String[] str = line.split("\\s+");
				ArrayList<Integer> feature = new ArrayList<Integer>();
				nodes.add(Integer.parseInt(str[0]));
				for(int i=1; i<str.length; i++)
					feature.add(Integer.parseInt(str[i]));
				featList.add(feature);
			}
		}  catch (IOException ex) {
	        ex.printStackTrace();
	    }
			
		try{
			FileOutputStream fout = new FileOutputStream(new File(fn+"_feature.txt"));
			for(int i=0; i<featList.size(); i++){
				for(int j=0; j<featList.get(i).size(); j++){
					fout.write((featList.get(i).get(j)+",").getBytes());
				}
				fout.write(("\r\n").getBytes());
			}			
		 } catch (IOException ex) {
			 ex.printStackTrace();
		 } 
	}
	
	public void readEdges(){
		if(nodes.isEmpty()){
			System.out.println("Read feature file first!");
			return;
		}
		
		int [][] adj = new int [nodes.size()][nodes.size()];
			
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(fn+".edges")));
			String line;
			while((line=br.readLine())!=null){
				line = line.trim();
				String[] str = line.split("\\s+");
				if(str.length > 2){
					System.out.println("Check the edge file!");
					return;
				}
				int id0 = Integer.parseInt(str[0]);
				int id1 = Integer.parseInt(str[1]);
				if(nodes.contains(id0) && nodes.contains(id1)){
					adj[nodes.indexOf(id0)][nodes.indexOf(id1)] = 1;
					adj[nodes.indexOf(id1)][nodes.indexOf(id0)] = 1;
				} else {
					System.out.println(id0+" "+id1);
				}
			}
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		
		try{
			FileOutputStream fout = new FileOutputStream(new File(fn+"_network.txt"));
			for(int i=0; i<adj.length; i++){
				for(int j=0; j<adj.length; j++){
					fout.write((adj[i][j]+",").getBytes());
				}
				fout.write(("\r\n").getBytes());
			}
		 } catch (IOException ex) {
			 ex.printStackTrace();
		 } 
	}
	
	public static void main(String[] args){
		String fn = "3980";
		readFacebook rf = new readFacebook(fn);
	}
}
