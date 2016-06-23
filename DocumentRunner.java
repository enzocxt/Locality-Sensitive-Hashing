import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * The DocumentRunner can be ran from the commandline to find the most similar pairs of documents in a directory.
 * Example command to run with brute force similarity search:
 * 				java DocumentRunner -threshold 0.5 -method bf -maxFiles 100 -dir articles -shingleLength 5 
 * Example command to run with LSH similarity search:
 * 				java DocumentRunner -threshold 0.5 -method lsh -maxFiles 100 -dir articles -shingleLength 5 -numHashes 100 -numBands 20
 *
 * @author Toon Van Craenendonck
 */

public class DocumentRunner {

	public static void main(String[] args) {	

		SimilaritySearcher searcher = constructSimilaritySearcher(args); // can be a brute force (bf) searcher or an LSH (lsh) searcher
		Set<SimilarPair> similarItems = new HashSet<SimilarPair>();
		
		if(args[0].equals("-threshold")){
		    similarItems = searcher.getSimilarPairsAboveThreshold(Double.parseDouble(args[1]));
		}
		
		printPairs(similarItems);

	}
	
	/**
	 * Constructs a similarity searcher.
	 */
	public static SimilaritySearcher constructSimilaritySearcher(String[] args){
		
		String fileDir = "";
		String method = "";
		int numHashes = -1;
		int numBands = -1;
		int seed = -1;
		int maxFiles = -1;
		int shingleLength = -1;
		
		int i = 0;
		while (i < args.length && args[i].startsWith("-")) {
			String arg = args[i];
			if (arg.equals("-method")) {
				if (!args[i+1].equals("bf") && !args[i+1].equals("lsh")){
					System.err.println("The search method should either be brute force (bf) or minhash and locality sensitive hashing (lsh)");
				}
                method = args[i+1];
            }else if(arg.equals("-numHashes")){
            	numHashes = Integer.parseInt(args[i+1]);
            }else if(arg.equals("-numBands")){
            	numBands = Integer.parseInt(args[i+1]);
            }else if(arg.equals("-seed")){
            	seed = Integer.parseInt(args[i+1]);
            }else if(arg.equals("-dir")){
            	fileDir = args[i+1];
            }else if(arg.equals("-maxFiles")){
            	maxFiles = Integer.parseInt(args[i+1]);
            }else if(arg.equals("-shingleLength")){
            	shingleLength = Integer.parseInt(args[i+1]);
            }
			i += 2;
		}

		DocumentHandler dh = new DocumentHandler(fileDir, maxFiles, shingleLength);
		Random rand = new Random(seed);
		
		if (method.equals("bf")){
			return new BruteForceSearch(dh.getObjectMapping());
		}else if(method.equals("lsh")){
		    if(numHashes == -1 || numBands == -1){
			throw new Error("Both -numHashes and -numBands are mandatory arguments for the LSH method"); 
		    }
		    return new LSH(dh.getObjectMapping(), numHashes, numBands, dh.getNumShingles(), rand);
		}
		return null;

	}
	
	/**
	 * Prints pairs and their similarity.
	 * @param similarItems
	 */
	public static void printPairs(Set<SimilarPair> similarItems){
		System.out.println("Number of candidate neighnors: " + similarItems.size());//
		List<SimilarPair> sim = new ArrayList<SimilarPair>(similarItems);
		Collections.sort(sim, Collections.reverseOrder());
		for(SimilarPair p : sim){
			System.out.println(p.getId1() + "," + p.getId2() + "," + p.getSimilarity());
		}
	}
	
}
