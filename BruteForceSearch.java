import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

/**
 * Brute force implementation of the similarity searcher. The Jaccard similarity is computed for all pairs and the most similar ones are selected.
 * 
 * @author Toon Van Craenendonck
 *
 */
public class BruteForceSearch extends SimilaritySearcher{

	public BruteForceSearch(Map<Integer,Set<Integer>> objectMapping){
		super(objectMapping);
	}
	
	/**
	 * Get pairs of objects with similarity above threshold.
	 * @param threshold the similarity threshold
	 * @return the pairs
	 */
	@Override
	public Set<SimilarPair> getSimilarPairsAboveThreshold(double threshold) {
		Set<SimilarPair> cands = new HashSet<SimilarPair>();
		for (Integer obj1 : objectMapping.keySet()){
			for (Integer obj2 : objectMapping.keySet()){
				if (obj1 < obj2){
					double sim = jaccardSimilarity(objectMapping.get(obj1),objectMapping.get(obj2));
					if (sim > threshold){
						cands.add(new SimilarPair(obj1, obj2, sim));
					}
	   			}
			}
		}
		return cands;
	}

	/**
	 * Get the objects that have a similarity above threshold thr to the object identified by the given object id objID
	 * @param objID the object of which we want to search neighbors
	 * @param thr the similarity threshold
  	 * @return the objects with similarity above threshold
	 */
	@Override
	public Set<Neighbor> getNeighborsAboveThreshold(int objID, double thr) {
		Set<Neighbor> candidateNeighbors = new HashSet<Neighbor>();
		for (int otherObj = 0; otherObj < objectMapping.size(); otherObj++){
			double sim = jaccardSimilarity(objectMapping.get(objID),objectMapping.get(otherObj));
			if (sim > thr){
				candidateNeighbors.add(new Neighbor(otherObj,sim));
			}
		}
		return candidateNeighbors;
	}
	

}
