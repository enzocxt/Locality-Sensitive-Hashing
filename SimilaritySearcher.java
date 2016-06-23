import java.util.HashSet;
import java.util.Set;
import java.util.Map;


/**
 * Searching similar objects. Objects should be represented as a mapping from an object identifier to a set containing the associated values.
 * 
 * @author Toon Van Craenendonck
 *
 */
public abstract class SimilaritySearcher {

	Map<Integer,Set<Integer>> objectMapping;

	public SimilaritySearcher(Map<Integer,Set<Integer>> objectMapping){
		this.objectMapping = objectMapping;
	}

	// /**
	//  * Returns the k most similar pairs in the given objectMapping.
	//  * @param k
	//  * @return a set containing the most similar pairs
	//  */
	// abstract public Set<SimilarPair> getMostSimilarPairs(int k);
	
	/**
	 * Returns the pairs of the objectMapping that have a similarity coefficient exceeding threshold
	 * @param threshold the similarity threshold
	 * @return the pairs with similarity above the threshold
	 */
	abstract public Set<SimilarPair> getSimilarPairsAboveThreshold(double threshold);
	
	/**
	 * Get the objects that have a similarity above threshold thr to the object identified by the given object id objID
	 * @param objID the object of which we want to search neighbors
	 * @param thr the similarity threshold
 	 * @return the objects with similarity above threshold
	 */
	abstract public Set<Neighbor> getNeighborsAboveThreshold(int objId, double threshold);

	/**
	 * Jaccard similarity between two sets.
	 * @param set1
	 * @param set2
	 * @return the similarity
	 */
	public <T> double jaccardSimilarity(Set<T> set1, Set<T> set2) {
		Set<T> union = new HashSet<T>(set1);
		union.addAll(set2);

		Set<T> intersection = new HashSet<T>(set1);
		intersection.retainAll(set2);

		if (union.size() == 0){
			return 0;
		}
		return (double) intersection.size() / union.size();
	}

}
