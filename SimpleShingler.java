import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A SimpleShingler constructs the shingle representations of documents.
 * It takes all substrings of length k of the document, and maps these substrings to an integer value that is inserted into the documents shingle set.
 * 
 * @author Toon Van Craenendonck
 *
 */
public class SimpleShingler {
	
	Map<String, Integer> shingleMap = new HashMap<String, Integer>(); // maps the k-shingles to integers
	int k;
	
	/**
	 * Construct a shingler.
	 * @param k number of characters in one shingle
	 */
	public SimpleShingler(int k){
		this.k = k;
	}

	/**
	 * Hash a k-shingle to an integer.
	 * @param shingle shingle to hash
	 * @return integer that the shingle maps to
	 */
	private int hashShingle(String shingle){
		if (shingleMap.containsKey(shingle)){
			return shingleMap.get(shingle);
		}else{
			shingleMap.put(shingle, shingleMap.size());
			return shingleMap.size() - 1;
		}
	}
	
	/**
	 * Get the shingle set representation of a document.
	 * @param fn filename of the document that should be shingled
	 * @return set of integers being the hash maps of the shingles
	 */
	public Set<Integer> shingle(String fn){
		Set<Integer> shingled = new HashSet<Integer>();

		String completeDocument = "";
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(fn));
			String line;
			while ((line = br.readLine()) != null) {
				completeDocument += line + " ";
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		for (int i = 0; i < completeDocument.length() - k; i++){
				String toHash = Character.toString(completeDocument.charAt(i));
				for (int j = 0; j < k; j++){
					toHash += Character.toString(completeDocument.charAt(i+j+1));
				}
				shingled.add(hashShingle(toHash));
		}
		return shingled;
	}

	/**
	 * Get the number of unique shingles this shingler has processed.
	 * @return number of unique shingles
	 */
	public int getNumShingles() {
		return shingleMap.size();
	}
}
