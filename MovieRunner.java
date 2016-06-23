import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * The MovieRunner can be ran from the commandline to perform minhash and LSH, and to afterwards predict user ratings.
 * Example command to run with LSH:
 * 			java MovieRunner -method lsh -numHashes 100 -numBands 20 -trainingFile ../r1.train -testFile ../r1.test -minRatingCount 3 -threshold 0.1
 * 
 * @author Toon Van Craenendonck
 */

public class MovieRunner {
	
	static MovieHandler ratings;
	static SimilaritySearcher searcher;
	static double threshold;
	static int minRatingCount;
	static String testFile;
	
	public static void main(String[] args) {	

		searcher = constructSimilaritySearcher(args);
		evaluate(testFile);
		
	}
	
	/**
	 * Constructs a similarity searcher.
	 */
	public static SimilaritySearcher constructSimilaritySearcher(String[] args){
		
		String method = "";
		String trainingFile = "";
		int numHashes = -1;
		int numBands = -1;
		int seed = -1;
		
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
            }else if(arg.equals("-trainingFile")){
            	trainingFile = args[i+1];
            }else if(arg.equals("-testFile")){
            	testFile = args[i+1];
            }else if(arg.equals("-threshold")){
            	threshold = Double.parseDouble(args[i+1]);
            }else if(arg.equals("-minRatingCount")){
            	minRatingCount = Integer.parseInt(args[i+1]);
            }
			i += 2;
		}
		
		ratings = new MovieHandler(trainingFile);	
		Random rand = new Random(seed);
		
		if (method.equals("bf")){
			if(numHashes == -1 || numBands == -1){
				throw new Error("Both -numHashes and -numBands are mandatory arguments for the LSH method"); 
			}
			return new BruteForceSearch(ratings.getObjectMapping());
		}else if(method.equals("lsh")){
			return new LSH(ratings.getObjectMapping(), numHashes, numBands, ratings.getNumValues(), rand);//
		}
		return null;

	}
	
	/**
	 * Predict the rating of user with external id externUserID for movie with id movieID.
	 * @param externUserID external id of user whose rating should be predict
	 * @param movieID movie for which the rating should be predicted
	 * @return the predicted rating
	 */
	public static double predictRating(int externUserID, int movieID){


	    /* Update this method */
		int internUserID = ratings.getUserIDs().indexOf(externUserID);
		Set<Neighbor> neighbors = searcher.getNeighborsAboveThreshold(internUserID, threshold);
		double avgRating = ratings.getAverageRating(externUserID);
		double predict = 0.0;
		double sum = 0.0;
		double a = 0.0;
		
		for (Neighbor nb: neighbors) {
			int nbID = ratings.getUserIDs().get(nb.getUserID()); // get external user id
			double avgNbRating = ratings.getAverageRating(nbID); // use external user id
			double sim = nb.getSimilarity();
			a += Math.abs(sim);
			double rating = 0.0;
			List<MovieRating> ratingList = ratings.getUsersToRatings().get(nbID); // use external user id
			for (MovieRating mr: ratingList)
				if (mr.getMovieID() == movieID) {
					rating = mr.getRating();
				}
			if (rating==0.0) rating = avgNbRating;
			sum += sim * (rating - avgNbRating);
		}
		predict = avgRating + (1.0/a) * sum;

	    //return MovieHandler.DEFAULT_RATING;
		return predict;

	}
	
	/**
	 * For each user/movie combination in the test set, predict the users' rating for the movie and compare to the true rating. 
	 * Prints the current mean absolute error (MAE) after every 50 users. 
	 * @param testFile path to file containing test set 
	 */
	public static void evaluate(String testFile){

		double summedErrorRecommenderSq = 0;
		double summedErrorAvgSq = 0;

		int avg_used = 0;
		int est_used = 0;
		int ctr = 0;
		
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(testFile));
			String line;
			while ((line = br.readLine()) != null) {
				String[] tokens = line.split("::|\t");

				int userID = Integer.parseInt(tokens[0]);
				int movieID = Integer.parseInt(tokens[1]);
				double rating = Double.parseDouble(tokens[2]);

				double avgRating = ratings.getMovieAverageRating(movieID);
				double estimate = predictRating(userID, movieID);
				System.out.println("predict: " + estimate);//

				summedErrorRecommenderSq += Math.pow(rating - estimate,2);
				summedErrorAvgSq += Math.pow(rating - avgRating, 2);
				ctr++;

				if (avgRating == estimate){
					avg_used++;
				}else{
					est_used++;
				}
				if ((ctr % 50) == 0){
					System.out.println("RMSE (default): " + Math.sqrt(summedErrorAvgSq/ctr) + " RMSE (recommender): " + Math.sqrt(summedErrorRecommenderSq/ctr));
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}
	
}
