/**
 * @author Yuki Kitayama
 * */
package graph;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;
import java.util.Iterator;
import javafx.util.Pair;
import util.GraphLoader;

public class FBGraph implements Graph {
	
	/** 
	 * Member variables 
	 * */
    private HashMap<Integer, HashSet<Integer>> network;
    private HashMap<Integer, FBGraph> mapSCC;
    private Integer root;
    private List<Graph> SCCs;
    private HashSet<Integer> setA;
    private HashSet<Integer> setB;
    private HashSet<Integer> initialAdopters;
    
    /** 
     * Constructor 
     * */
    public FBGraph() {
    	network = new HashMap<Integer, HashSet<Integer>>();
    	mapSCC = new HashMap<Integer, FBGraph>();
    	SCCs = new ArrayList<Graph>();
    }
    
    /** 
     * Add a new vertex in graph, used in GraphLoader loadGraph method 
     * to import vertices to this network.
     * 
     * @param num vertex
     * */
    @Override
    public void addVertex(int num) {
    	if (!network.containsKey(num)) {
    		HashSet<Integer> edge = new HashSet<Integer>();
    		network.put(num, edge);
    	}
    }
    
    /** 
     * Add a new edge to a given vertex, used in GraphLoade loadGraph
     * method to import edges to this network.
     * 
     * @param from this vertex
     * @param to one neighbor
     * */
    @Override
    public void addEdge(int from, int to) {
    	if (network.containsKey(from)) {
    		HashSet<Integer> edge = network.get(from);
    		if (!edge.contains(to)) {
    			edge.add(to);
    		}
    	}
    }
    
    /** 
     * Method to get Egonet.
     * 
     * @param center central vertex of an Egonet.
     * @return Graph class which this class implements 
     * */
    @Override
    public Graph getEgonet(int center) {
    	// initialize graph to be returned
    	Graph graph = new FBGraph();
    	// if the center is not in the graph
    	if (!network.containsKey(center)) {
    		return graph;
    	}
    	// if the center IS in the graph
    	graph.addVertex(center);
    	HashSet<Integer> friends = network.get(center);
    	for (Integer friend : friends) {
    		// add edge between center and friend
    		graph.addEdge(center, friend);
    		// connections between friends of friends
    		HashSet<Integer> friendsOfFriend = network.get(friend);
    		for (Integer friendOfFriend : friendsOfFriend) {
    			if (friends.contains(friendOfFriend) || center == friendOfFriend) {
    				// add friend to graph
    				graph.addVertex(friend);
    				// add edge between friend and friend of friend
    				graph.addEdge(friend, friendOfFriend);
    			}
    		}
    	}
    	return graph;
    }
    
    /** 
     * Method to get Strongly Connected Components.
     * 
     * @return list of Graph classes
     * */
    @Override
    public List<Graph> getSCCs() {
    	// initialize vertices
    	Stack<Integer> vertices = new Stack<Integer>();
    	// add every vertex in the current network to the stack
    	for (Integer vertex : network.keySet()) {
    		// push just put element to stack, but add not only puts but also return boolean
    		vertices.push(vertex);
    	}
    	// get Stack from the first DFS
    	Stack<Integer> finishedOne = DFS(this, vertices);
    	// transpose network
    	FBGraph transposedGraph = transposeGraph();
    	// get Stack from the second DFS
    	Stack<Integer> finishedSecond = StackFromSCC(transposedGraph, finishedOne);
    	// finally get subgraph SCCs and add it to List
    	HashMap<Integer, FBGraph> map = transposedGraph.getMapSCC();
    	for (FBGraph subgraph : map.values()) {
    		SCCs.add(subgraph);
    	}
    	return SCCs;
    }
    
    /** 
     * When GraphLoader loadGraph is called, addVertex and addEdge methods
     * are called and add them to the network member variable in this
     * class. This method just gets the variable.
     * 
     * @return hash map with key vertices and value edges
     * */
    @Override
    public HashMap<Integer, HashSet<Integer>> exportGraph() {
    	return network;
    }
    
    /** 
     * Helper method to get Strongly Connected Components. Do Depth First
     * Search.
     * 
     * @param graph Graph class from UCSD
     * @param vertices stack of vertices to trace SCC during exploration
     * @return Stack of finished vertices 
     * */
    private Stack<Integer> DFS(Graph graph, Stack<Integer> vertices) {
    	HashSet<Integer> visited = new HashSet<Integer>();
    	Stack<Integer> finished = new Stack<Integer>();
    	while (!vertices.empty()) {
    		// pop removes element at the top of the stack
    		Integer vertex = vertices.pop();
    		if (!visited.contains(vertex)) {
    			DFSVisit(graph, vertex, visited, finished);
    		}
    	}
    	return finished;
    }
    
    /** 
     * Helper method for Depth First Search method
     * 
     * @param graph Graph class from UCSD
     * @param vertex a vertex available from DFS method
     * @param visited hashset of vertices to let DFS happen
     * @param finished stack of vertices to trace SCC method
     * */
    private void DFSVisit(
    		Graph graph, 
    		Integer vertex, 
    		HashSet<Integer> visited, 
    		Stack<Integer> finished) {
    	visited.add(vertex);
    	HashMap<Integer, HashSet<Integer>> currNetwork = graph.exportGraph();
    	HashSet<Integer> friends = currNetwork.get(vertex);
    	for (Integer friend : friends) {
    		if (!visited.contains(friend)) {
    			DFSVisit(graph, friend, visited, finished);
    		}
    	}
    	finished.push(vertex);
    }
    
    /** 
     * Helper method for Depth First Search. This is used in the second 
     * DFS from the Stack from the first DFS.
     * 
     * @param graph this graph
     * @param vertex a vertex available from DFS method
     * @param visited hashset of vertices to let DFS happen
     * @param finished stack of vertices to trace SCC method 
     * */
    private void DFSVisitSCC(
    		FBGraph graph,
    		Integer vertex,
    		HashSet<Integer> visited,
    		Stack<Integer> finished) {
    	visited.add(vertex);
    	Integer root = graph.getRoot();
    	HashMap<Integer, FBGraph> mapSCC = graph.getMapSCC();
    	FBGraph subgraph = mapSCC.get(root);
    	// get transposed network
    	HashMap<Integer, HashSet<Integer>> transposedNetwork = graph.exportGraph();
    	HashSet<Integer> friends = transposedNetwork.get(vertex);
    	// add vertex and edges to subgraph
    	for (Integer friend : friends) {
    		if (!finished.contains(friend)) {
    			subgraph.addVertex(vertex);
    			subgraph.addEdge(vertex, friend);
    			// connections between friends of friends
    			HashSet<Integer> friendsOfFriend = transposedNetwork.get(friend);
    			for (Integer friendOfFriend : friendsOfFriend) {
    				if (friends.contains(friendOfFriend) || vertex == friendOfFriend) {
    					subgraph.addVertex(friend);
    					subgraph.addEdge(friend, friendOfFriend);
    				}
    			}
    		}
    		if (!visited.contains(friend)) {
    			DFSVisitSCC(graph, friend, visited, finished);
    		}
    	}
    	finished.push(vertex);
    }
    
    /** 
     * Helper method to return transposed graph of this network
     * 
     * @return this graph
     * */
    private FBGraph transposeGraph() {
    	// initialize
    	FBGraph transposedGraph = new FBGraph();
    	// add vertices and edges
    	for (Integer friend : network.keySet()) {
    		HashSet<Integer> friendsOfFriend = network.get(friend);
    		for (Integer friendOfFriend : friendsOfFriend) {
    			transposedGraph.addVertex(friendOfFriend);
    			transposedGraph.addEdge(friendOfFriend, friend);
    		}
    	}
    	// add missing vertices from the original network
    	// but leave edge empty
    	for (Integer friend : network.keySet()) {
    		if (!transposedGraph.exportGraph().keySet().contains(friend)) {
    			transposedGraph.addVertex(friend);
    		}
    	}
    	return transposedGraph;
    }
    
    /** 
     * Helper method in getting Strongly Connected Components.
     * 
     * @param graph This graph class implementing Graph
     * @param vertices stack of integer available from SCC method
     * @return Stack of vertices produced during getting SCC (will not be used) 
     * */
    private Stack<Integer> StackFromSCC(FBGraph graph, Stack<Integer> vertices) {
    	// initialize
    	HashSet<Integer> visited = new HashSet<Integer>();
    	Stack<Integer> finished = new Stack<Integer>();
    	// traverse
    	while (!vertices.empty()) {
    		Integer vertex = vertices.pop();
    		if (!visited.contains(vertex)) {
    			FBGraph subgraph = new FBGraph();
    			subgraph.addVertex(vertex);
    			// update SCC
    			graph.setMapSCC(vertex, subgraph);
    			graph.setRoot(vertex);
    			// recursive
    			DFSVisitSCC(graph, vertex, visited, finished);
    		}
    	}
    	return finished;
    }
    
    /** 
     * Method to suggest friends for easy question in the project. Pick one
     * user from the network and seek pair of people who are user's friend 
     * but are not friends with each other. When the network is big, returned
     * list is long. So only randomly pick 5 pairs as suggestion.
     * 
     * @param user one vertex from the network
     * @return a list of pairs of unconnected potential friends
     * */
    public ArrayList<Pair<Integer, Integer>> suggestFriend(int user) {
    	// initialize returned list
    	ArrayList<Pair<Integer, Integer>> suggestion = new ArrayList<Pair<Integer, Integer>>();
    	// there's duplicated pair, for example 1=2 and 2=1, so if reverse is the same, should omit it
    	HashSet<Integer> keySet = new HashSet<Integer>();
    	// get a list of friends of a user
    	HashSet<Integer> friends = network.get(user);
    	for (Integer x : friends) {
    		for (Integer y : friends) {
    			if (x != y && !network.get(x).contains(y) && !keySet.contains(y)) {
    				Pair<Integer, Integer> pair = new Pair<Integer, Integer>(x, y);
    				suggestion.add(pair);
    				keySet.add(x);
    			}
    		}
    	}
    	// if no suggestion, show that there is not friend suggestion
    	if (suggestion.size() == 0) {
    		System.out.println("*** There is no friends suggestion ***");
    		return suggestion;
    	}
    	else if (suggestion.size() <= 5) {
    		return suggestion;
    	}
    	// with real data, returned ArrayList<Pair> becomes long
    	// randomly sample only 5 pairs of friends as suggestion
    	// select random items without repetitions list element
    	// https://www.geeksforgeeks.org/randomly-select-items-from-a-list-in-java/
    	else {
	    	Random rand = new Random();
	    	long s = 123;
	    	rand.setSeed(s);
	    	int totalItems = 5;
	    	ArrayList<Pair<Integer, Integer>> sampledList = new ArrayList<Pair<Integer, Integer>>();
	    	for (int i = 0; i < totalItems; i++) {
	    		int randomIndex = rand.nextInt(suggestion.size());
	    		sampledList.add(suggestion.get(randomIndex));
	    		// remove selected element from original list
	    		suggestion.remove(randomIndex);
	    	}
	    	return sampledList;
    	}
    }
        
    /**
     * Iteratively cascade a new technology A around the network.
     * Before running this method, initial adopters of A must be initialized by setInitialAdopters method
     * During cascade, people converted from B to A are added to setA member variable.
     * After cascade, people who did not convert are stored in setB as reference.
     * 
     * @param payoffA relative payoff to B, e.g. 2 means twice as good as B, 
     *                should experiment various numbers
     * @param payoffB relative payoff to A, but usually fixed at 1.
     * */
    public void cascade(double payoffA, double payoffB) {
    	// initialize threshold to compare with p
    	double threshold = payoffB / (payoffA + payoffB);
    	// count how many iterations it takes to equilibrium
    	int iteration = 0;
    	// boolean to keep running iteration until equilibrium
    	boolean switching = true;
    	// store cascade result
    	setA = initialAdopters;
    	// start iteration of cascade
    	while (switching) {
    		iteration += 1;
	        // initialize empty set which will contain candidate to join new tech group
	    	HashSet<Integer> candidateSetA = new HashSet<Integer>();
    		// count how many new people join a new tech group in each iteration
    		int counterAddition = 0;
    		Iterator<Integer> SetAIter = setA.iterator();
    		while (SetAIter.hasNext()) {
	    		HashSet<Integer> neighbors = network.get(SetAIter.next());
	        	for (Integer neighbor : neighbors) {
	        		int counterA = 0;
	        		HashSet<Integer> neighborsOfNeighbor = network.get(neighbor);
	        		int d = neighborsOfNeighbor.size();
	        		for (Integer neighborOfNeighbor : neighborsOfNeighbor) {
	        			if (setA.contains(neighborOfNeighbor)) {
	        				counterA += 1;
	        			}
	        		}
	            	// p is fraction of neighbors who use a new technology A
	        		double p = counterA / (double) d;
	        		if (p > threshold && !setA.contains(neighbor) && !candidateSetA.contains(neighbor)) {
	        			candidateSetA.add(neighbor);
	        			counterAddition += 1;
	        		}
	    	    }
	    	}
	    	// after for loop add newCandidates to newSet
	    	setA.addAll(candidateSetA);
	    	if (counterAddition == 0) {
	    		switching = false;
	    	}
    	}
    	// set a set of vertices which was not affected by cascade
    	setB = new HashSet<Integer>(network.keySet());
    	setB.removeAll(setA);
    	// reference info
    	System.out.println("*** cascade info ***");
    	System.out.printf("threshold: %.3f\n", threshold);
    	System.out.printf("iterations to get equilibrium: %,d\n", iteration);
    	System.out.println("********************");
    }
    
    /**
     * Check whether a complete cascade occurred.
     * Run this only after cascade method.
     * 
     * @return Boolean of whether or not complete cascade occurred.
     * */
    public boolean isCompleteCascade() {
    	int numVertexNetwork = getNumVertex();
    	int numVertexSetA = setA.size();
    	if (numVertexNetwork == numVertexSetA) {
    		return true;
    	}
    	else {
    		return false;
    	}
    }
    
    /**
     * Define cascade ratio as number of vertices using A after cascade 
     * divided by total number of vertices. Depending on the network 
     * structure, a complete cascade sometimes does not happen. Checking
     * the ratio can give you sense of how cascade diffused the innovation.
     * 
     *  @return cascade ratio
     * */
    public double getCascadeRatio() {
    	int numVertex = getNumVertex();
    	int numVertexSetA = setA.size();
    	double ratio = numVertexSetA / (double) numVertex;
    	return ratio;
    }
    
    /**
     * Generates initial adopters and set it to initialAdopters member variable.
     * Randomly pick vertices from the set of integers loaded from GraphLoader.
     * 
     * @param size the number of initial adopters before cascade
     * @param seed to observe experiments
     * */
    public void setInitialAdopters(int size, long seed) {
    	// initialize result
    	initialAdopters = new HashSet<Integer>();
    	// convert set to list
    	List<Integer> list = new ArrayList<Integer>(network.keySet());
    	// initialize random
    	Random random = new Random(seed);
    	for (int i = 0; i < size; i++ ) {
    	    int vertex = list.get(random.nextInt(list.size()));
    	    if (!initialAdopters.contains(vertex)) {
    	    	initialAdopters.add(vertex);
    	    }
    	}
    }
    
    /**
     * Find isolated networks within the network. To let a complete cascade
     * happen, all the vertices need to belong to the main network. But
     * some network like UCSD Facebook has a small number of vertices which
     * only has one neighbor, and the neighbor also has one neighbor, which
     * is the first vertex. This method finds those special vertices. 
     * 
     * @return ArrayList of Pair of integers, no duplicates
     * */
    public ArrayList<Pair<Integer, Integer>> findOneToOne() {
    	// initialize result
    	ArrayList<Pair<Integer, Integer>> listOneToOne = new ArrayList<Pair<Integer, Integer>>();
    	// initialize set to avoid duplicated pair
    	HashSet<Integer> keySet = new HashSet<Integer>();
    	for (Integer vertex : network.keySet()) {
    		HashSet<Integer> thisEdge = network.get(vertex);
    		// check this friend size
    		if (thisEdge.size() == 1) {
    			Integer neighbor = thisEdge.toArray(new Integer[1])[0];
    			// check that friend size
    			HashSet<Integer> thatEdge = network.get(neighbor);
    			// second condition avoid duplicated pairs
    			if (thatEdge.size() == 1 && !keySet.contains(neighbor)) {
    				// this condition can add them to pair
    				Pair<Integer, Integer> pair = new Pair<Integer, Integer>(vertex, neighbor);
    				listOneToOne.add(pair);
    				keySet.add(vertex);
    			}
    		}
    	}
    	return listOneToOne;
    }
    
    /**
     * Setters
     * */
    public void setMapSCC(Integer vertex, FBGraph subgraph) {
    	mapSCC.put(vertex, subgraph);
    }
    
    /** helper method to update root
     * */
    public void setRoot(Integer vertex) {
    	root = vertex;
    }    

    /**
     * Getters
     * */
    public HashMap<Integer, FBGraph> getMapSCC() {
    	return mapSCC;
    }
    
    public Integer getRoot() {
    	return root;
    }     
    
    public HashSet<Integer> getInitialAdopters() {
    	return initialAdopters;
    }
    
    public HashSet<Integer> getSetA() {
    	return setA;
    }
    
    public HashSet<Integer> getSetB() {
    	return setB;
    }
    
    /**
     * get the number of vertices in the graph. 
     * */
    public int getNumVertex() {
    	return network.size();
    }
    
    /**
     * get the average number of edges among vertices.
     * */
    public double getAvgEdge() {
    	int n = 0;
    	int sum = 0;
    	double avg;
    	for (HashSet<Integer> edges : network.values()) {
    		sum += edges.size();
    		n += 1;
    	}
    	avg = sum / (double) n;
    	return avg;
    }
    
    /**
     * get the average number of edges in set B*/
    public double getAvgEdgeSetB() {
    	int n = 0;
    	int sum = 0;
    	double avg;
    	for (Integer vertex : setB) {
    		HashSet<Integer> neighbors = network.get(vertex);
    		sum += neighbors.size();
    		n += 1;
    	}
    	avg = sum / (double) n;
    	return avg;
    }
        
    /** test
     * */
    public static void main(String[] args) {
    	// test getSCCs
//    	System.out.println("Test getSCCs started...");
//    	FBGraph graphSCC = new FBGraph();
//    	GraphLoader.loadGraph(graphSCC, "data/test_scc.txt");
//    	List<Graph> SCCs = graphSCC.getSCCs();
//    	for (Graph subgraph : SCCs) {
//    		System.out.println("subgraph.exportGraph() " + subgraph.exportGraph());
//    	}
//    	System.out.println("Test getSCCs ended...");
//    	System.out.println("**********");
    	
    	// test suggestFriend
//    	System.out.println("Test suggestFriend method started...");
//    	FBGraph graphSuggest = new FBGraph();
//    	GraphLoader.loadGraph(graphSuggest, "data/test_easy_1.txt");
//    	int user = 0;
//    	ArrayList<Pair<Integer, Integer>> suggestion = graphSuggest.suggestFriend(user);
//    	for (Pair<Integer, Integer> pair : suggestion) {
//    		System.out.println("pair " + pair.toString());
//    	}
//    	System.out.println("Test suggestFriend method ended...");
//    	System.out.println("**********");
    	
    	// test suggestFriend with a bit more complicated test cast
//    	System.out.println("Test suggestFriend method with difficult case started...");
//    	FBGraph graphSuggest2 = new FBGraph();
//    	GraphLoader.loadGraph(graphSuggest2, "data/test_easy_2.txt");
    	// one of the users case
//    	int user2 = 0;
//    	System.out.println("Friend suggestion that user " + user2 + " can do is,");
//    	ArrayList<Pair<Integer, Integer>> suggestion2 = graphSuggest2.suggestFriend(user2);
//    	for (Pair<Integer, Integer> pair : suggestion2) {
//    		System.out.println("pair " + pair.toString());
//    	}
    	// another user case
//    	int user3 = 4;
//    	System.out.println("Friend suggestion that user " + user3 + " can do is,");
//    	ArrayList<Pair<Integer, Integer>> suggestion3 = graphSuggest2.suggestFriend(user3);
//    	for (Pair<Integer, Integer> pair : suggestion3) {
//    		System.out.println("pair " + pair.toString());
//    	}
    	// one more user case
//    	int user4 = 2;
//    	System.out.println("Friend suggestion that user " + user4 + " can do is,");
//    	ArrayList<Pair<Integer, Integer>> suggestion4 = graphSuggest2.suggestFriend(user4);
//    	for (Pair<Integer, Integer> pair : suggestion4) {
//    		System.out.println("pair " + pair.toString());
//    	}
//    	System.out.println("Test suggestFriend method with difficult case ended...");
//    	System.out.println("**********");
    	
    	// test suggestFriend with real data, facebook_ucsd.txt
//    	System.out.println("Apply suggestFriend method to facebook_ucsd.txt started...");
//    	FBGraph graphUCSD = new FBGraph();
//    	GraphLoader.loadGraph(graphUCSD, "data/facebook_ucsd.txt");
//    	int userUCSD1 = 0;
//    	System.out.println("Friend suggestion of UCSD student " + userUCSD1 + " can do is,");
//    	ArrayList<Pair<Integer, Integer>> suggestionUCSD1 = graphUCSD.suggestFriend(userUCSD1);
//    	for (Pair<Integer, Integer> pair : suggestionUCSD1) {
//    		System.out.println("pair " + pair.toString());
//    	}
//    	System.out.println("Apply suggestFriend method to facebook_ucsd.txt ended...");
//    	System.out.println("**********");
    	
    	// test informationFlow with sample graph
//    	System.out.println("Test informationFlow with sample graph started...");
//    	FBGraph graphIF = new FBGraph();
//    	GraphLoader.loadGraph(graphIF, "data/test_information_flow.txt");
//    	graphIF.informationFlow();
//    	System.out.println("**********");
//    	System.out.println("Test overloaded method");
//    	int rewardA = 3;
//    	int rewardB = 1;
//    	HashSet<Integer> newSet = new HashSet<Integer>();
//    	newSet.add(25);
//    	newSet.add(18);
//    	graphIF.informationFlow(newSet, rewardA, rewardB);
//    	System.out.println("Test informationFlow with sample graph ended...");
//    	System.out.println("**********");
    	
    	// test informationFlow with real data facebook_ucsd.txt
//    	System.out.println("Test informationFlow with UCSD facebook data started...");
//    	FBGraph graphIfUcsd = new FBGraph();
//    	GraphLoader.loadGraph(graphIfUcsd, "data/facebook_ucsd.txt");
//    	HashSet<Integer> newSetUcsd = new HashSet<Integer>();
//    	newSetUcsd.add(0);
//    	int rewardAUcsd = 1;
//    	int rewardBUcsd = 1;
//    	graphIfUcsd.informationFlow(newSetUcsd, rewardAUcsd, rewardBUcsd);
//    	System.out.println("Test informationFlow with UCSD facebook data ended...");
//    	System.out.println("**********");
    	
    	// test utility methods
//    	System.out.println("Test utility methods started...");
//    	FBGraph graphTest = new FBGraph();
//    	GraphLoader.loadGraph(graphTest, "data/facebook_ucsd.txt");
//    	int numVertex = graphTest.getNumVertex();
//    	double avgEdge = graphTest.getAvgEdge();
//    	System.out.printf("numVertex: %,d\n", numVertex);
//    	System.out.printf("avgEdge: %.1f\n", avgEdge);
//    	System.out.println("Test utility methods ended...");
//    	System.out.println("**********");
    	
    	// test cascade
    	System.out.println("Test cascade started...");
    	FBGraph graphTest = new FBGraph();
    	// load data
    	String data = "data/facebook_ucsd.txt";
    	// String data = "data/test_cascade.txt";
    	 GraphLoader.loadGraph(graphTest, data);
    	// basic information of the graph
    	int numVertex = graphTest.getNumVertex();
    	double avgEdge = graphTest.getAvgEdge();
    	ArrayList<Pair<Integer, Integer>> listOneToOne = graphTest.findOneToOne();
    	System.out.println("*** basic info ***");
    	System.out.printf("numVertex: %,d\n", numVertex);
    	System.out.printf("avgEdge: %.1f\n", avgEdge);
    	System.out.println("listOneToOne: " + listOneToOne);
    	System.out.println("******************");
    	// initialize initial adopters, setInitialAdopters(size, seed)
    	graphTest.setInitialAdopters(1000, 1);
    	// run cascade, cascade(payoffA, payoffB)
    	graphTest.cascade(10, 1);
    	// check result of cascade
    	boolean cascadeCompleteBool = graphTest.isCompleteCascade();
    	double cascadeRatio = graphTest.getCascadeRatio();
    	HashSet<Integer> setB = graphTest.getSetB();
    	double avgEdgeSetB = graphTest.getAvgEdgeSetB();
    	System.out.println("*** cascade result ***");
    	System.out.println("Is there complete cascade?: " + cascadeCompleteBool);
    	System.out.printf("Cascade ratio: %.3f\n", cascadeRatio);
    	System.out.println("Set B: " + setB.toString());
    	System.out.printf("Average number of edges in set B: %.3f\n", avgEdgeSetB);
    	System.out.println("***********************");
    	// HashSet<Integer> initialAdopters = graphTest.getInitialAdopters();
    	// System.out.println("initialAdopters: " + initialAdopters.toString());    	
    	System.out.println("Test cascade ended...");
    	System.out.println("**********");
    }
}
