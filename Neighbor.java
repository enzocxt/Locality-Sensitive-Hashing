/**
 * Simple class introduced for convenience, a neighbor has an ID and a Jaccard similarity to a certain user.
 * 
 * @author Toon Van Craenendonck
 *
 */
public class Neighbor implements Comparable<Neighbor>{

	int id;
	double similarity;
	
	public Neighbor(int id, double similarity) {
		this.id = id;
		this.similarity = similarity;
	}

	public int compareTo(Neighbor nb) {
		if (similarity < nb.getSimilarity()){
			return -1; 
		}else if (similarity == nb.getSimilarity()){
			return 0;
		}else{
			return 1;
		}
	}

	public double getSimilarity() {
		return similarity;
	}
	
	public int getUserID() {
		return id;
	}
	

}
