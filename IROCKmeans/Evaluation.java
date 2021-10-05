package IROCKmeans;

import java.util.*;

public class Evaluation {
	
	public static double Mean(double [] nmi){
		double mean = 0.0;
		for(int i=0; i<nmi.length; i++){
			mean += nmi[i];
		}
		return mean/nmi.length;
	}
	
	public static double Var(double [] nmi){
		double var = 0.0;
		double mean = Mean(nmi);
		for(int i=0; i<nmi.length; i++){
			var += Math.pow(nmi[i]-mean, 2);
		}
		
		return var/nmi.length;
	}
	
	public static double Stdv(double [] nmi){
		return Math.sqrt(Var(nmi));
	}

	public static double NMI(int [] g, int [] t){
		int num = g.length;
		
		int [] gold = g.clone();
		int [] test = t.clone();
		
		/** Re-order the label array in case there are missing labels */
		ArrayList<Integer> gold_element = new ArrayList<Integer>();
		ArrayList<Integer> test_element = new ArrayList<Integer>();
		for(int i=0; i<num; i++){
			if(!gold_element.contains(gold[i])) gold_element.add(gold[i]);
			if(!test_element.contains(test[i])) test_element.add(test[i]);
		}
		
		/** Get the number of true classes in each array */
		int K1 = gold_element.size();
		int K2 = test_element.size();
		
		/** Count the number of entries in each class and finish re-ordering */
		int [] NK1 = new int [K1];
		int [] NK2 = new int [K2];		
		for(int i=0; i<num; i++){
			gold[i] = gold_element.indexOf(gold[i]);
			NK1[gold[i]] ++;
			test[i] = test_element.indexOf(test[i]);
			NK2[test[i]] ++;
		}
		
		/** Calculate the entropy */
		double HK1 = 0.0;
		double HK2 = 0.0;
		for(int i=0; i<K1; i++)
			HK1 += - ((double)NK1[i]/num) * Math.log((double)NK1[i]/num);
		for(int i=0; i<K2; i++)
			HK2 += - ((double)NK2[i]/num) * Math.log((double)NK2[i]/num);
		
		/** Build contingency table */
		int [][] T = new int [K1][K2];
		for(int i=0; i<num; i++)
			T[gold[i]][test[i]] ++;
		
		/** Calculate the MI */
		double MI = 0.0;
		for(int i=0; i<K1; i++){
			for(int j=0; j<K2; j++){
				if(T[i][j] > 0) MI += T[i][j] * Math.log((double)T[i][j] * num / (NK1[i] * NK2[j]));
			}
		}
		MI /= num;
		
		/** Calculate the NMI */
		double NMI = MI / ((HK1+HK2)/2);
		
		return NMI;
	}
	
	public static double Purity(int [] g, int [] t){
		int num = g.length;
		
		int [] gold = g.clone();
		int [] test = t.clone();
		
		/** Re-order the label array in case there are missing labels */
		ArrayList<Integer> gold_element = new ArrayList<Integer>();
		ArrayList<Integer> test_element = new ArrayList<Integer>();
		for(int i=0; i<num; i++){
			if(!gold_element.contains(gold[i])) gold_element.add(gold[i]);
			if(!test_element.contains(test[i])) test_element.add(test[i]);
		}
		
		/** Get the number of true classes in each array */
		int K1 = gold_element.size();
		int K2 = test_element.size();
		
		/** Finish re-ordering */	
		for(int i=0; i<num; i++){
			gold[i] = gold_element.indexOf(gold[i]);
			test[i] = test_element.indexOf(test[i]);
		}
		
		/** Build contingency table */
		int [][] T = new int [K1][K2];
		for(int i=0; i<num; i++)
			T[gold[i]][test[i]] ++;
		
		/** Count the number of correctly assigned points */
		double count = 0;
		for(int i=0; i<K2; i++){
			int max = 0;
			for(int j=0; j<K1; j++){
				if(T[j][i] > max)
					max = T[j][i];
			}
			count += max;
		}
		
		/** Calculate the Purity */
		double purity = count / num;
		
		return purity;
	}
	
	public static double FMeasure(int [] g, int [] t){
		int num = g.length;
		
		int [] gold = g.clone();
		int [] test = t.clone();
		
		/** Re-order the label array in case there are missing labels */
		ArrayList<Integer> gold_element = new ArrayList<Integer>();
		ArrayList<Integer> test_element = new ArrayList<Integer>();
		for(int i=0; i<num; i++){
			if(!gold_element.contains(gold[i])) gold_element.add(gold[i]);
			if(!test_element.contains(test[i])) test_element.add(test[i]);
		}
		
		/** Finish re-ordering */	
		for(int i=0; i<num; i++){
			gold[i] = gold_element.indexOf(gold[i]);
			test[i] = test_element.indexOf(test[i]);
		}
		
		/** Calculate true positive, false positive, false negative and true negative */
		double TP = 0;
		double FP = 0;
		double FN = 0;
		double TN = 0;
		for(int i=0; i<num; i++){
			for(int j=i+1; j<num; j++){
				if(test[i] == test[j]){
					if(gold[i] == gold[j])
						TP ++;
					else
						FP ++;
				} else {
					if(gold[i] == gold[j])
						FN ++;
					else
						TN ++;
				}
			}
		}
		
		/** Calculate the FMeasure */
		double precision = TP / (TP + FP);
		double recall = TP / (TP + FN);
		double fmeasure = 2 * precision * recall / (precision + recall);
		
		return fmeasure;
	}
}
