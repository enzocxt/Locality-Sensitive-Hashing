/**
 * @author Tao Chen
 *
 * The MovieHandler reads the MovieLens data and constructs several mappings:
 * 		- usersToRatings: maps the true user ID to a list of movie ratings
 *
 * 		- userIDs: maps the internal user ID to the true user ID
 * 		- movieIDs: maps the internal movie ID to the true movie ID
 *
 * 		- usersToSets: maps the internal user ID to a set containing the users' likes and dislikes
 *
 * The internal IDs are introduced to make sure that the IDs used to perform minhashing and LSH nicely go from 0 to num_users or num_movies.
 * Also, creating this mapping here ensures that the produced signature matrix and LSH tables are equal for equal inputs.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class MovieHandler{

	String ratingFile;

	Map<Integer, List<MovieRating>> usersToRatings = new HashMap<Integer, List<MovieRating>>(); 
	Map<Integer, Set<Integer>> usersToSets = new HashMap<Integer, Set<Integer>>();

	ArrayList<Integer> movieIDs;
	ArrayList<Integer> userIDs;

	Map<Integer, Double> movieAverageRatings; 

	static double DEFAULT_RATING = 2.5; 

	/**
	 * Create new MovieHandler. 
	 * Creating a MovieHandler object results in reading the rating and movie files, and constructing the above described structures.
	 * 
	 * @param fileName name of file containing the ratings
	 * @param movieFile name of file containing the additional movie information, we will only process the title
	 */
	public MovieHandler(String fileName) {
		this.ratingFile = fileName;

		long startTime = System.currentTimeMillis();
		System.out.println("Reading data.. ");
		this.readData();
		System.out.println("done, took " +  (System.currentTimeMillis() - startTime)/1000.0 + "seconds.");
		System.out.println("--------------");

		startTime = System.currentTimeMillis();
		System.out.println("Converting to set representation.. ");		
		this.convertToSetRepresentation();
		System.out.println("done, took " +  (System.currentTimeMillis() - startTime)/1000.0 + "seconds.");
		System.out.println("--------------");

	}

	/**
	 * Returns internal ID to true ID mapping.
	 * @return userIDs the mapping
	 */
	public ArrayList<Integer> getUserIDs(){
		return userIDs;
	}

	/**
	 * Reads the MovieLens data into a map, mapping user IDs to lists of movie ratings.
	 * Also creates internal to true ID mappings for users and movies. 
	 */
	private void readData() {
		Set<Integer> movieSet = new HashSet<Integer>();

		BufferedReader br;
		//int flag = 0;//
		try {
			br = new BufferedReader(new FileReader(ratingFile));
			String line;
			while ((line = br.readLine()) != null) {
				String[] tokens = line.split("::|\t");

				int userID = Integer.parseInt(tokens[0]);
				int movieID = Integer.parseInt(tokens[1]);
				double rating = Double.parseDouble(tokens[2]);

				movieSet.add(movieID);
				if (!usersToRatings.containsKey(userID)){
					List<MovieRating> ratingList = new ArrayList<MovieRating>();
					ratingList.add(new MovieRating(movieID,rating));
					usersToRatings.put(userID, ratingList);
				}else{
					usersToRatings.get(userID).add(new MovieRating(movieID, rating));
				}
			}
			br.close();

		} catch (IOException e) {
			e.printStackTrace();
		}

		// store the user ids as a sorted list (just to make sure that we have a unique ordering)
		userIDs = new ArrayList<Integer>(usersToRatings.keySet());
		Collections.sort(userIDs);

		// same for movie ids
		movieIDs = new ArrayList<Integer>(movieSet);
		Collections.sort(movieIDs);

		// precompute average ratings for each movie
		computeMovieAverageRatings(); 
	}

	/**
	 * Converts the user to MovieRating mapping to a user to set mapping. 
	 * The set contains the movies that the user either liked or disliked. 
	 * A user likes a movie if he rated it >= 3, and dislikes a movie if he rated it < 3. 
	 * If a user likes a movie, the element 2*id is added to the users' set.
	 * If a users dislikes a movie, the element 2*id + 1 is added to the users' set.
	 */
	private void convertToSetRepresentation(){
		// constructing these mappings once avoids an indexOf lookup for every movie/user
		Map<Integer, Integer> trueToInternalMovieIDs = new HashMap<Integer, Integer>();
		for (int i = 0; i < movieIDs.size(); i++){
			trueToInternalMovieIDs.put(movieIDs.get(i),i);
		}
		Map<Integer, Integer> trueToInternalUserIDs = new HashMap<Integer, Integer>();
		for (int i = 0; i < userIDs.size(); i++){
			trueToInternalUserIDs.put(userIDs.get(i),i);
		}
		
		for (Integer userID : usersToRatings.keySet()) {
			Set<Integer> newSet = new HashSet<Integer>();
			for (MovieRating rating : usersToRatings.get(userID)){
				if (rating.getRating() >= this.getAverageRating(userID)){
					newSet.add(2*trueToInternalMovieIDs.get(rating.getMovieID())); 
				}else{
					newSet.add(2*trueToInternalMovieIDs.get(rating.getMovieID()) + 1);
				}
			}
			usersToSets.put(trueToInternalUserIDs.get(userID), newSet);  
		}
		//System.out.println("usersToSets: " + usersToSets.get(1));//
	}

	/**
	 * Computes the average rating for a user.
	 * @param userID external userID 
	 * @return the users' average rating
	 */
	public double getAverageRating(Integer userID) {
		double avg = 0;
		for (MovieRating r : usersToRatings.get(userID)){
			avg += r.getRating();
		}
		return avg / usersToRatings.get(userID).size();
	}

	/** 
	 * Fetch the average movie rating from the cache (see computeMovieAverageRatings())
	 */ 
	public double getMovieAverageRating(Integer movieID){
		if(movieAverageRatings.containsKey(movieID))
			return movieAverageRatings.get(movieID);
		else
			return DEFAULT_RATING;
	}


	/**
	 * Returns the user to set mapping. Internal IDs are used here.
	 * @return the mapping
	 */
	public Map<Integer, Set<Integer>> getObjectMapping() {
		return usersToSets;
	}

	/**
	 * Returns the user to movie mapping. External IDs are used here.
	 * @return the mappings
	 */
	public Map<Integer, List<MovieRating>> getUsersToRatings() {
		return usersToRatings;
	}

	/**
	 * Returns the number of movies that were read.
	 * @return the number of movies
	 */
	public int getNumMovies() {
		return movieIDs.size();
	}

	/**
	 * Returns the number of unique values that can be found in the object mapping sets.
	 * @return the number of values
	 */
	public int getNumValues() {
		return 2*movieIDs.size();
	}

	/**
	 * Computes and caches movies average ratings. 
	 */
	private void computeMovieAverageRatings() {
		movieAverageRatings = new HashMap<Integer, Double>(); 

		Map<Integer, Double> ratingSum = new HashMap<Integer, Double>();
		Map<Integer, Integer> ratingCount = new HashMap<Integer, Integer>();

		for(Integer userID : userIDs){
			for (MovieRating r : usersToRatings.get(userID)){
				if(! ratingSum.containsKey(r.getMovieID())){
					ratingSum.put(r.getMovieID(), (double)r.getRating()); 
					ratingCount.put(r.getMovieID(), 1); 
				}
				else{
					ratingSum.put(r.getMovieID(), (ratingSum.get(r.getMovieID()) + r.getRating())); 
					ratingCount.put(r.getMovieID(), ratingCount.get(r.getMovieID()) + 1);
				}
			}
		}

		for (Map.Entry<Integer, Double> entry : ratingSum.entrySet()) {
			int movieID = entry.getKey(); 
			double avg = entry.getValue() / ratingCount.get(movieID);

			movieAverageRatings.put(movieID, avg);
		}

	}


}
