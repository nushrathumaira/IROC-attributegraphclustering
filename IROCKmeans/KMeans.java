package IROCKmeans;

import java.util.HashSet;
import java.util.Random;
import java.util.ArrayList;

public class KMeans {
	private int k;
	private int num;
	private int [] label;
	private int [] preLabel;
	private double [][] dist;
	
	public KMeans(int k, double [][] dist){
		this.k = k;
		this.num = dist.length;
		this.label = new int [num];
		this.preLabel = label.clone();
		this.dist = dist;
		
		/** Initialization */
		ArrayList<Integer> center = randomInitCenter();
		for(int i=0; i<num; i++){
			double min = Double.MAX_VALUE;
			for(int j=0; j<k; j++){
				if(dist[i][center.get(j)] < min){
					min = dist[i][center.get(j)];
					label[i] = j;
				}
			}
		}
		
		int iteration = 0;
		boolean converge = false;
		while(!converge && iteration < 1000){
			preLabel = label.clone();
			
			for(int i=0; i<num; i++){
				double min = Double.MAX_VALUE;
				for(int j=0; j<k; j++){
					double distij = averageDist(i, j);
					if(distij < min){
						min = distij;
						label[i] = j;
					}
				}
			}
			
			if(Evaluation.NMI(label, preLabel) > 0.99999999)
				converge = true;
			
			iteration ++;
		}
	}
	
	public int [] getLabel(){
		return label;
	}
	
	public double averageDist(int obj, int clus){
		int numClus = 0;
		double avgDist = 0.0;
		for(int i=0; i<num; i++){
			if((i != obj) && (preLabel[i] == clus)){
				numClus ++;
				avgDist += dist[obj][i];
			}
		}
		
		return avgDist/numClus;
	}
	
	public ArrayList<Integer> randomInitCenter(){
		HashSet<Integer> center = new HashSet<Integer>();
		Random ran = new Random();
		while(center.size() < k){
			center.add(ran.nextInt(num));
		}
		return new ArrayList<Integer>(center);
	}
	
	public static void main(String[] args) {	
//		/** Read Musk Data */
////		readMusk r = new readMusk("Musk1.data", 476, 166);
//		readMusk r = new readMusk("Musk2.data", 6598, 166);
////		readFace r = new readFace("imageData5_5_gray.txt", 640*25, 5);
//		MultiInsData [] data = r.getData();
//		int [] gold = r.getGold();
//		
//		int max = 0;
//		for(int i=0; i<data.length; i++){
//			if(max < data[i].size())
//				max = data[i].size();
//		}
//		DistanceFunctions.QWX = max;
//		
//		DistanceFunction<MultiInstance> metric = null;
//			
////		metric = DistanceFunctions.INFHAMMING;
////		double [][] dist = new double [data.length][data.length];
////		for(int i=0; i<data.length; i++){
////			for(int j=i+1; j<data.length; j++){
////				dist[i][j] = metric.calculate(data[i], data[j]);
////				dist[j][i] = dist[i][j];
////			}
////		}
//
//		
//		int nn = 1;
//		ArrayList<double [][]> distList = new ArrayList<double [][]>();
//		for(int m=0; m<nn; m++){
//			System.out.println(m);
//			metric = DistanceFunctions.QUANTILE;
//			double [][] dist = new double [data.length][data.length];
//			for(int i=0; i<data.length; i++){
//				for(int j=i+1; j<data.length; j++){
//					dist[i][j] = metric.calculate(data[i], data[j]);
//					dist[j][i] = dist[i][j];
//				}
//			}
//			distList.add(dist);
//		}
//		
//		int itr = 100;
//		double [] nmi = new double [itr];
//		double [] purity = new double [itr];
//		double [] fmeasure = new double [itr];
//		for(int m=0; m<itr; m++){			
//			KMeans km = new KMeans(2, distList.get(m/100));
//			purity[m] = Evaluation.Purity(gold, km.getLabel());
//			nmi[m] = Evaluation.NMI(gold, km.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, km.getLabel());
//		}
//		System.out.println("K = 2;");
//		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
//		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
//		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
//		System.out.println();
//		
//		for(int m=0; m<itr; m++){			
//			KMeans km = new KMeans(3, distList.get(m/100));
//			purity[m] = Evaluation.Purity(gold, km.getLabel());
//			nmi[m] = Evaluation.NMI(gold, km.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, km.getLabel());
//		}
//		System.out.println("K = 3;");
//		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
//		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
//		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
//		System.out.println();
//		
//		for(int m=0; m<itr; m++){
//			KMeans km = new KMeans(4, distList.get(m/100));
//			purity[m] = Evaluation.Purity(gold, km.getLabel());
//			nmi[m] = Evaluation.NMI(gold, km.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, km.getLabel());
//		}
//		System.out.println("K = 4;");
//		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
//		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
//		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
//		System.out.println();
//		
//		for(int m=0; m<itr; m++){			
//			KMeans km = new KMeans(5, distList.get(m/100));
//			purity[m] = Evaluation.Purity(gold, km.getLabel());
//			nmi[m] = Evaluation.NMI(gold, km.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, km.getLabel());
//		}
//		System.out.println("K = 5;");
//		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
//		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
//		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
//		System.out.println();
//		
//		for(int m=0; m<itr; m++){			
//			KMeans km = new KMeans(6, distList.get(m/100));
//			purity[m] = Evaluation.Purity(gold, km.getLabel());
//			nmi[m] = Evaluation.NMI(gold, km.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, km.getLabel());
//		}
//		System.out.println("K = 6;");
//		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
//		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
//		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
//		System.out.println();
	}
}
