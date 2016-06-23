/**
 * @author Tao Chen
 *
 * Implementation of minhash and locality sensitive hashing (lsh) to find similar objects.
 *
 * The LSH should first construct a signature matrix. Based on this, LSH is performed resulting in a mapping of band ids to hash tables (stored in bandToBuckets).
 * From this bandsToBuckets mapping, the most similar items should then be retrieved.
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

public class LSH extends SimilaritySearcher{
	List<Map<String, Set<Integer>>> bandToBuckets;

	/**
	 * Construct an LSH similarity searcher.
	 * 
	 * @param objectMapping objects and their set representations of which similarity should be searched
	 * @param numHashes number of hashes to use to construct the signature matrix
	 * @param numBands number of bands to use during locality sensitive hashing
	 * @param numValues the number of unique values that occur in the objects' set representations (i.e. the number of rows of the original characteristic matrix)
	 * @param rand should be used to generate any random numbers needed
	 */
	public LSH(Map<Integer, Set<Integer>> objectMapping, int numHashes, int numBands, int numValues, Random rand){
		super(objectMapping);
		
		int prime  = Primes.findLeastPrimeNumber(numValues);
		int[][] hashValues = LSH.constructHashTable(numHashes, numValues, prime, rand);
		int[][] signatureMatrix = LSH.constructSignatureMatrix(objectMapping, hashValues);
		bandToBuckets = LSH.lsh(signatureMatrix, numBands);
	}

	/**
	 * Returns the band to buckets mapping.
	 * @return
	 */
	public List<Map<String, Set<Integer>>> getBandToBuckets(){
		return bandToBuckets;
	}

	/**
	 * Construct the table of hash values needed to construct the signature matrix.
	 * Position (i,j) contains the result of applying function j to row number i.
	 * @param numHashes number of hashes that will be used in the signature matrix
	 * @param numValues number of unique values that occur in the object set representations (i.e. number of rows of the characteristic matrix)
	 * @param prime prime needed for universal hashing
	 * @param rand object to generate random numbers
	 * @return the (numValues x numHashes) matrix of hash values
	 */
	public static int[][] constructHashTable(int numHashes, int numValues, int prime, Random rand) {
		int[][] hashes = new int[numValues][numHashes];

		/* Fill in here */
		System.out.println("Constructing hash table...");
		for (int j=0; j<numHashes; j++) {
			for (int i=0; i<numValues; i++) {
			/* Universal hashing: h_a,b(x) = ((a*x + b) mod p) mod N */
				int a = rand.nextInt(numValues);
				if (a==0) a=1;
				int b = rand.nextInt(numValues);
				hashes[i][j] = (int) (((long) a*i + b) % prime) % numValues;
			}
		}
		System.out.println("done");
		System.out.println("--------------");
		
		return hashes;
	}

	/**
	 * Constructing the signature matrix.
	 * 
	 * @param objectMapping objects and their set representations for which the signature matrix should be constructed
	 * @param hashValues (numValues x numHashes) matrix of hash values
	 * @return the (numHashes x numObjects) signature matrix
	 */
	public static int[][] constructSignatureMatrix(Map<Integer, Set<Integer>> objectMapping, int[][] hashValues) {
		int numHashes = hashValues[0].length;
		int numObjects = objectMapping.size();

		// initialize it to max int values
		System.out.println("Constructing Signature Matrix...");
		int[][] signatureMatrix = new int[numHashes][numObjects];
		for (int i=0; i<numHashes; i++) {
			for (int j=0; j<numObjects; j++) {
				signatureMatrix[i][j] = Integer.MAX_VALUE;
			}
		}

		/* Fill in here */
		int numValues = hashValues.length;
		System.out.println("Number of feature values: " + numValues);
		for (int i=0; i<numValues; i++) {
			//if (i%1000==0) System.out.println(i + " feature values done");
			if (i % (numValues/10) == 0) System.out.println((i*100 / numValues) + "% feature values done");
			for (int j=0; j<numObjects; j++) {
				if (objectMapping.get(j).contains(i))
					for (int k=0; k<numHashes; k++)
						//if (signatureMatrix[k][j] > hashValues[i][k]) signatureMatrix[k][j] = hashValues[i][k];
						signatureMatrix[k][j] = Math.min(signatureMatrix[k][j], hashValues[i][k]);
			}
		}
		System.out.println("done");//
		System.out.println("--------------");

		return signatureMatrix;

	}

	/**
	 * Perform locality sensitive hashing.
	 * 
	 * @param signatureMatrix previously constructed signature matrix
	 * @param numBands the number of bands to use
	 * @return mapping of bands to corresponding hash maps, these hash maps in turn map keys (being parts of the signature matrix) to sets of ids of potentially similar objects
	 */
	public static List<Map<String, Set<Integer>>> lsh(int[][] signatureMatrix, int numBands) {
		List<Map<String, Set<Integer>>> bandToBuckets =
		    new ArrayList<Map<String, Set<Integer>>>();

		/* Fill in here */
		System.out.println("Performing locality sensitive hashing...");//
		int numHashes = signatureMatrix.length;
		int numObjects = signatureMatrix[0].length;
		int prime = Primes.findLeastPrimeNumber(numObjects);
		//int prime = Primes.findLeastPrimeNumber(10000000);
		int Rows = (int) numHashes / numBands;	// rows per band
		int b = numObjects * 5; // number of buckets per band
		System.out.println("Number of hashes: " + numHashes);
		System.out.println("Number of objects: " + numObjects);
		System.out.println("Number of rows per band: " + Rows);
		
		for (int i=0; i<numBands; i++) {
			int[] hv = new int[numObjects];
			for (int k=0; k<numObjects; k++) hv[k] = 0;
			//String[] keys = new String[numObjects];//
			//for (int k=0; k<numObjects; k++) keys[k] = "";//
			
			for (int j=0; j<Rows; j++) {
				for (int k=0; k<numObjects; k++) {
					hv[k] = (int) ((hv[k] + (long) signatureMatrix[i*Rows+j][k]*prime) % b);
					//keys[k] = keys[k] + Integer.toString(signatureMatrix[i*Rows+j][k]);//
				}
			}

			Map<String, Set<Integer>> btb = new HashMap<String, Set<Integer>>();
			for (int k=0; k<numObjects; k++) {
				if (btb.containsKey(Integer.toString(hv[k]))) {
					btb.get(Integer.toString(hv[k])).add(k);
				}else{
					Set<Integer> s = new HashSet<Integer>();
					s.add(k);
					btb.put(Integer.toString(hv[k]), s);
				}
			}
			bandToBuckets.add(btb);
		}
		System.out.println("Length of bandToBuckets: " + bandToBuckets.size());
		System.out.println("done");//
		System.out.println("--------------");

		return bandToBuckets;

	}

	/**
	 * Returns the pairs with similarity above threshold (approximate).
	 */
	@Override
	public Set<SimilarPair> getSimilarPairsAboveThreshold(double threshold) {
		List<SimilarPair> cands = new ArrayList<SimilarPair>();
		
		/* Fill in here */
		System.out.println("Computing similar pairs...");
		int fp = 0;
		for (int i=0; i<bandToBuckets.size(); i++) {
			System.out.println(i + "th band");//
			
			Map<String, Set<Integer>> btb = bandToBuckets.get(i);
			List<Set<Integer>> bList = new ArrayList<Set<Integer>>(btb.values());
			for (int j=0; j<bList.size(); j++) {
				List<Integer> candsBucket = new ArrayList<Integer>(bList.get(j));
				for (int k=0; k<candsBucket.size(); k++) {
					for (int l=k+1; l<candsBucket.size(); l++) {
						/* compute similarity of candsBucket.get(k) and candsBucket.get(l) */
						double sim = this.jaccardSimilarity(this.objectMapping.get(candsBucket.get(k)), this.objectMapping.get(candsBucket.get(l)));
						if (sim > threshold) {
							//cands.add(new SimilarPair(Math.min(candsBucket.get(k), candsBucket.get(l)), Math.max(candsBucket.get(k), candsBucket.get(l)), sim));
							
							int index = 0;
							while (index<cands.size()) {
								if (cands.get(index).getId1()==Math.min(candsBucket.get(k), candsBucket.get(l)) && cands.get(index).getId2()==Math.max(candsBucket.get(k), candsBucket.get(l)))
									break;
								index ++;
							}
							if (index == cands.size())
								cands.add(new SimilarPair(Math.min(candsBucket.get(k), candsBucket.get(l)), Math.max(candsBucket.get(k), candsBucket.get(l)), sim));
							
						}else{ fp ++;}
					}
				}
			}
		}
		/* delete duplicates */
		/* Alternative method which improves the speed, but still has bugs. 
		List<SimilarPair> simPairs = new ArrayList<SimilarPair>(cands);
		Collections.sort(simPairs, Collections.reverseOrder());
		List<Integer> index = new ArrayList<Integer>();
		int i=0;
		while (i<(simPairs.size()-1)) {
			index.add(i);;
			int j=i+1;
			while (j<simPairs.size()) {
				if (simPairs.get(i).getId1() == simPairs.get(j).getId1()) {
					if (simPairs.get(i).getId2() == simPairs.get(j).getId2()) {
						j++;
						if (j>=simPairs.size()) i=j;
					}else{
						i = j;
						break;
					}
				}else{
					i = j;
					break;
				}
			}
		}
		List<SimilarPair> simPairs2 = new ArrayList<SimilarPair>();
		for (int k=0; k<index.size(); k++) {
			simPairs2.add(simPairs.get(index.get(k)));
		}
		
		//System.out.println("Number of candidates: " + simPairs2.size());
		*/
		
		System.out.println("Number of candidates: " + cands.size());
		
		System.out.println("Number of false positives: " + fp);
		System.out.println("--------------");

		return new HashSet<SimilarPair>(cands);
		//return new HashSet<SimilarPair>(simPairs2);

	}

	/**
	 * Get the objects that have a similarity above threshold thr to the object identified by the given object id objID
	 * @param objID the object of which we want to search neighbors
	 * @param thr the similarity threshold
	 * @return the objects with similarity above thr
	 */
	@Override
	public Set<Neighbor> getNeighborsAboveThreshold(int internalID, double thr) {
		Set<Neighbor> candidateNeighbors = new HashSet<Neighbor>();
		
		/* Fill in here (only required for movie rating predictions */ 
		
		Set<Integer> added = new HashSet<Integer>(); // for checking neighbors that have been added
		
		for (Map<String, Set<Integer>> btb: bandToBuckets)
			for (Set<Integer> s: btb.values())
				if (s.contains(internalID))
					for (Integer id: s) {
						double sim = this.jaccardSimilarity(this.objectMapping.get(internalID), this.objectMapping.get(id));
						if (sim > thr && !added.contains(id)){
							candidateNeighbors.add(new Neighbor(id, sim));
							added.add(id);
						}
					}
		
		return candidateNeighbors;
	}


}
