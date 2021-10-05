package clustering;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import util.Evaluation;
import util.readMusk;
import data.MultiInsData;
import distance.Metrics;

/**
 * Partitioning Around Medoids. Cluster data elements by aiming to minimize the average
 * dissimilarity of objects to their closest selected element (medoid).
 * Data elements are indexable.
 * Independent of Cytoscape.
 * NB   This class is an implementation of the PAM algorithm described in the course notes
 *      at <www.cs.umb.edu/cs738/pam1.pdf> (accessed 2012-07-31).
 *      There appears to be different variants of the PAM algorithm varying in details
 *      within the build and swap phases.
 *      The original algorithm is described in chapter 2 of Kaufman and Rousseeuw (1990),
 *      whose original Fortran code was translated and augmented as part of the cluster R package.
 *      The cluster results from current implementation can differ from the implementation 
 *      in R's cluster::pam. (See PAMTest for details.)
 * @author djh.shih
 * @comment
 *
 */
public class PAM {	
	
	// clustering cost, to be minimized
	double cost;
	
	// distance between element and closest medoid
	protected double[] nearestDistances;
	// distance between element and second closest medoid
	double[] nextNearestDistances;
	
	// nearest medoid of each element
	int[] nearestMedoids;
	// next-nearest medoid of each element 
	int[] nextNearestMedoids;
	
	// set of medoids
	HashSet<Integer> medoids;
	
	// set of non-meoids (maintain for finding swap candidates)
	HashSet<Integer> nonmedoids;
	
	// set of all indexed elements
	// required since Java's HashSet cannot use native types
	Integer[] elements;
	
	int maxSwaps = 1000;
	
	double [][] dist;
	
	int nClusters;
	int nObjects;
	
	int [] label;
	
	public PAM(int k, double [][] dist){
		this.dist = dist;
		this.nClusters = k;
		this.nObjects = dist.length;
		
		initialize();
//		buildPhase();
		bulidPhaseRandom();
		swapPhase();
		
		int count = 0;
		int [] center = new int [k];
		for(Integer i : medoids)
			center[count++] = i;
		label = new int [nObjects];
		for(int i=0; i<nObjects; i++){
			double min = Double.MAX_VALUE;
			for(int j=0; j<k; j++){
				if(dist[i][center[j]] < min){
					min = dist[i][center[j]];
					label[i] = j;
				}
			}
		}
	}
	
	public int [] getLabel(){
		return label;
	}
	
//	/**
//	 * Calculate the clustering cost: sum of distances to cluster medoids.
//	 * @return cost
//	 */
//	private double getCost() {
//		double c = 0;
//		for (int i = 0; i < nearestDistances.length; ++i) {
//			c += nearestDistances[i];
//		}
//		return c;
//	}
	
	private void initialize() {
		int m = nObjects;
		nearestDistances = new double[m];
		nextNearestDistances = new double[m];
		nearestMedoids = new int[m];
		nextNearestMedoids = new int[m];
		
		
		elements = new Integer[m];
		medoids = new HashSet<Integer>();
		nonmedoids = new HashSet<Integer>();
		
		for (int ii = 0; ii < m; ++ii) {
			// initialize distances to infinity
			nearestDistances[ii] = nextNearestDistances[ii] = Double.POSITIVE_INFINITY;
			// initialize medoids to non-valid indices, s.t. unexpected bugs trigger indexing error
			nearestMedoids[ii] = nextNearestMedoids[ii] = -1;
			
			elements[ii] = new Integer(ii);
			
			// all (indexed) data elements are initially non-medoids
			nonmedoids.add( elements[ii] );
		}
	}
	
	private void bulidPhaseRandom() {
		HashSet<Integer> center = new HashSet<Integer>();
		Random ran = new Random();
		while(center.size() < nClusters){
			center.add(ran.nextInt(nObjects));
		}
		for(Integer i : center)
			addMedoid(i);
	}
	
//	/**
//	 * BUILD phase. Select a initial set of k medoids.
//	 */
//	private void buildPhase() {
//		int m = nObjects;
//		
//		// select first medoid
//		
//		// find element with minimum total distance to all other elements
//		double[] totalDistances = new double[m];
//		for (int ii = 0; ii < m; ++ii) {
//			// sum distances to all other elements
//			// assume distance to itself is 0
//			double d = 0;
//			for (int jj = 0; jj < m; ++jj) {
//				d += dist[ii][jj];
//			}
//			totalDistances[ii] = d;
//		}
//		double minDistance = totalDistances[0];
//		int minIndex = 0;
//		for (int ii = 0; ii < m; ++ii) {
//			if (totalDistances[ii] < minDistance) {
//				minDistance = totalDistances[ii];
//				minIndex = ii;
//			}
//		}
//		// add element to medoid set
//		addMedoid(minIndex);
//		
//		
//		// select remaining k - 1 medoids
//		
//		double[] gains = new double[m];
//		
//		for (int kk = 1; kk < nClusters; ++kk) {
//		
//			// consider each i as medoid candidate
//			for (int ii = 0; ii < m; ++ii) {
//				// if ii is already a medoid, it has negative gain to prevent it from being selected again
//				if (medoids.contains(elements[ii])) {
//					gains[ii] = -1.0;
//				} else {
//					double gain = 0;
//					// for each non-medoid j != i, calculate the gain
//					for (int jj = 0; jj < m; ++jj) {
//						if (jj == ii || medoids.contains(elements[jj]) ) continue;
//						if (nearestDistances[jj] > dist[ii][jj]) {
//							// add i will improve j's nearest distances
//							// (if selected, i will be the new nearest neighbour of j)
//							gain += nearestDistances[jj] - dist[ii][jj];
//						}
//					}
//					gains[ii] = gain;
//				}
//			}
//			// select candidate with maximum gain
//			double maxGain = Double.NEGATIVE_INFINITY;
//			int maxIndex = -1;
//			for (int ii = 0; ii < m; ++ii) {
//				if (gains[ii] > maxGain) {
//					maxGain = gains[ii];
//					maxIndex = ii;
//				}
//			}
//			// add element to medoid set
//			addMedoid(maxIndex);
//			
//		}
//		
//		// check that the number of medoids match the expected
//		if (nClusters != medoids.size()) {
//			throw new RuntimeException("Expected error in BUILD phase: Number of medoids does not match parameter k.");
//		}
//		
//	}
	
	/**
	 * SWAP phase. Attempt to improve clustering quality by exchanging medoids with non-medoids.
	 */
	private void swapPhase() {
		boolean notConverged = true;
		boolean continueLoop = true;
		int nSwaps = 0;
		
		while (notConverged && continueLoop) {
			notConverged = false;
			continueLoop = false;
	
			Iterator<Integer> medIt = medoids.iterator();
			while (medIt.hasNext() && continueLoop) {
				int ii = medIt.next().intValue();
				
				Iterator<Integer> nonmedIt = nonmedoids.iterator();
				while (nonmedIt.hasNext()) {
					int hh = nonmedIt.next().intValue();
				
					// Consider swapping medoid i and nonmedoid h
					// by calculating gains by all other elements
					
					// Calculate cumulative change to distance to nearest medoid for all nonmedoids j != h
					double change = 0;
					Iterator<Integer> nonmedIt2 = nonmedoids.iterator();
					while (nonmedIt2.hasNext()) {
						int jj = nonmedIt2.next().intValue();
						if (jj == hh) continue;
					
						double d = nearestDistances[jj];
						if (dist[ii][jj] > d) {
							// if removed, i will have no impact
							if (dist[ii][hh] < d) {
								// if selected, h will improve nearest distance for j
								change += dist[jj][hh] - d;
							}
						} else {
							// i cannot be closer than the nearest neighbour for j;
							// therefore, distances[i][j] == d
							// and i is currently the nearest neighbour for j
							double e = nextNearestDistances[jj];
							if (dist[jj][hh] < e) {
								// if i and h are swapped, h will become the nearest neighbour
								// nearest distance for j may improve or worsen
								change += dist[jj][hh] - d;
							} else {
								// if i is removed, the current next-nearest of j will be promoted to nearest
								change += e - d;
							}
						}
					}
					
					if (change < 0) {
						// distance to nearest medoid summed over all nonmedoids is improved: swap
						swap(hh, ii);
						//System.out.print("Swap " + hh + " and " + ii + " for change = " + change + "\n");
						
						// non-convergence if any swap occurs, up to a maximum number of swaps (to guard against swap cycles)
						if (nSwaps++ < maxSwaps) {
							notConverged = true;
						} else {
							continueLoop = false;
						}
						
						// reset iterator
						medIt = medoids.iterator();
						// break out of inner loop to consider next medoid
						break;
					}
					
				}
				
			}
			
		}
	}
	
	private void addMedoid(int add) {
		medoids.add( elements[add] );
		nonmedoids.remove( elements[add] );
		updateNearest(add, -1);
	}
	
	private void swap(int add, int remove) {
		medoids.add( elements[add] );
		nonmedoids.remove( elements[add] );
		medoids.remove( elements[remove] );
		nonmedoids.add( elements[remove] );
		updateNearest(add, remove);
	}
	
	/**
	 * Update nearest and next-nearest distances.
	 * Does not check whether {@code added} or {@ removed} have been added to or removed from the medoid set.
	 * @param added Index of element added to medoid set (-1 for none)
	 * @param removed Index of element removed from medoid set (-1 for none)
	 */
	private void updateNearest(int added, int removed) {
		int m = nObjects;
		if (added >= 0) {
			// added index is valid
			
			// check if any nearest distance improves
			for (int ii = 0; ii < m; ++ii) {
				double d = dist[ii][added];
				if (d < nearestDistances[ii]) {
					// element i is nearer to added medoid than previous nearest: update
					double oldDistance = nearestDistances[ii];
					int oldMedoid = nearestMedoids[ii];
					nearestMedoids[ii] = added;
					nearestDistances[ii] = d;
					// pump nearest distance to next-nearest distance
					nextNearestMedoids[ii] = oldMedoid;
					nextNearestDistances[ii] = oldDistance;
				} else if (d < nextNearestDistances[ii]) {
					// element i is nearer to added medoid than previous next-nearest: update
					nextNearestMedoids[ii] = added;
					nextNearestDistances[ii] = d;
				}
			}
			
		}
		
		if (removed >= 0) {
			// removed index is valid
			
			// check if the removed medoid is the nearest or next-nearest of any element
			for (int ii = 0; ii < m; ++ii) {
				if (nearestMedoids[ii] == removed) {
					// promote next-nearest to nearest
					nearestMedoids[ii] = nextNearestMedoids[ii];
					nearestDistances[ii] = nextNearestDistances[ii];
					// find new next-nearest
					updateNextNearest(ii);
				} else if (nextNearestMedoids[ii] == removed) {
					// find new next-nearest
					updateNextNearest(ii);
				}
			}
			
		}
	}
	
	/**
	 * Update next nearest for element i.
	 * Assume nearest medoid is already set.
	 * @param ii element index to be updated
	 */
	private void updateNextNearest(int ii) {
		int nearestMedoid = nearestMedoids[ii];
		
		// find the next-nearest
		Iterator<Integer> it = medoids.iterator();
		double minDistance = Double.POSITIVE_INFINITY;
		int nextNearestMedoid = -1;
		while (it.hasNext()) {
			int jj = it.next().intValue();
			// ignore if j is the nearestMedoid, since we are interested in the next-nearest
			if (jj == nearestMedoid) continue;
			if (dist[ii][jj] < minDistance) {
				minDistance = dist[ii][jj];
				nextNearestMedoid = jj;
			}
		}
	
		// update
		nextNearestDistances[ii] = minDistance;
		nextNearestMedoids[ii] = nextNearestMedoid;
	}

	public static void main(String[] args) {	
		/** Read Musk Data */
		readMusk rm = new readMusk("Musk1.data", 476, 166);
//		readMusk rm = new readMusk("Musk2.data", 6598, 166);
		MultiInsData [] data = rm.getData();
		int [] gold = rm.getGold();
		double [][] dist = new double [data.length][data.length];
		for(int i=0; i<data.length; i++){
			System.out.println(i);
			for(int j=i+1; j<data.length; j++){
				dist[i][j] = Metrics.InfHamming(data[i], data[j]);
				dist[j][i] = dist[i][j];
			}
		}
		
		int itr = 1000;
		double [] nmi = new double [itr];
		double [] purity = new double [itr];
		double [] fmeasure = new double [itr];
		for(int m=0; m<itr; m++){
			PAM pam = new PAM(2, dist);
			purity[m] = Evaluation.Purity(gold, pam.getLabel());
			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
		}
		System.out.println("K = 2;");
		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
		System.out.println();
		
		for(int m=0; m<itr; m++){
			PAM pam = new PAM(3, dist);
			purity[m] = Evaluation.Purity(gold, pam.getLabel());
			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
		}
		System.out.println("K = 3;");
		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
		System.out.println();
		
		for(int m=0; m<itr; m++){
			PAM pam = new PAM(4, dist);
			purity[m] = Evaluation.Purity(gold, pam.getLabel());
			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
		}
		System.out.println("K = 4;");
		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
		System.out.println();
		
		for(int m=0; m<itr; m++){
			PAM pam = new PAM(5, dist);
			purity[m] = Evaluation.Purity(gold, pam.getLabel());
			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
		}
		System.out.println("K = 5;");
		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
		System.out.println();
		
		for(int m=0; m<itr; m++){
			PAM pam = new PAM(6, dist);
			purity[m] = Evaluation.Purity(gold, pam.getLabel());
			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
		}
		System.out.println("K = 6;");
		System.out.println("Purity, Mean: "+Evaluation.Mean(purity)+", Stdv: "+Evaluation.Stdv(purity));
		System.out.println("NMI, Mean: "+Evaluation.Mean(nmi)+", Stdv: "+Evaluation.Stdv(nmi));
		System.out.println("FMeasure, Mean: "+Evaluation.Mean(fmeasure)+", Stdv: "+Evaluation.Stdv(fmeasure));
		System.out.println();
		
//		for(int m=0; m<itr; m++){
//			PAM pam = new PAM(10, dist);
//			purity[m] = Evaluation.Purity(gold, pam.getLabel());
//			nmi[m] = Evaluation.NMI(gold, pam.getLabel());
//			fmeasure[m] = Evaluation.FMeasure(gold, pam.getLabel());
//		}
//		System.out.println("Mean: "+Evaluation.Mean(nmi));
//		System.out.println("Stdv: "+Evaluation.Stdv(nmi));
	}
}

