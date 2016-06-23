import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Reads a directory of documents and constructs shingle representations for these documents.
 * 
 * @author Toon Van Craenendonck
 *
 */
public class DocumentHandler {

	Map<Integer, Set<Integer>> fileToShingle = new HashMap<Integer, Set<Integer>>(); // maps a fileID (which in this implementation simply corresponds to the filename)
																				     // to its shingle set representation
	SimpleShingler shingler; 
	int numDocuments;

	/**
	 * Constructing a DocumentHandler causes it to read maxFiles documents from directory fileDir and construct shingle set representations for them.
	 * @param fileDir file directory
	 * @param maxFiles number of files to read
	 */
	public DocumentHandler(String fileDir, int maxFiles, int shingleLength){
		this.numDocuments = maxFiles;
		shingler = new SimpleShingler(shingleLength);
		
		for (int fileID = 0; fileID < maxFiles; fileID++){
			Set<Integer> shingle = shingler.shingle(fileDir + "/" + fileID);
			fileToShingle.put(fileID,shingle);
		}
	}
	
	/**
	 * Get the mapping of filename (which in this case is equal to the object id) to its set representation.
	 * @return the mapping
	 */
	public Map<Integer, Set<Integer>> getObjectMapping() {
		return fileToShingle;
	}

	/**
	 * Get the number of unique shingles that were processed.
	 * @return the number of unique shingles
	 */
	public int getNumShingles(){
		return shingler.getNumShingles();
	}

	/**
	 * Get the number of documents that were processed.
	 * @return the number of documents.
	 */
	public int getNumDocuments(){
		return numDocuments;
	}
}
