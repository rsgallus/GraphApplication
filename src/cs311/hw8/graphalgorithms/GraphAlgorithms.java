
package cs311.hw8.graphalgorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Stack;
import cs311.hw8.graph.Graph;
import cs311.hw8.graph.IGraph;
import cs311.hw8.graph.IGraph.Edge;
import cs311.hw8.graph.IGraph.Vertex;
import cs311.hw8.graph.IWeight;

/**
 * 
 * This class contains various algorithms to perform operations on graphs,
 * including topological sort, Krsucal's algorithm for minimum spanning tree,
 * and Dijkstra's algorithm for shortest path.
 * 
 * @author Ryan Gallus
 *
 */
public class GraphAlgorithms {

	/**
	 * This method performs a topological sort on a graph, returning one
	 * possible sort.
	 * 
	 * @param g
	 *            the graph to perform a topological sort on
	 * @return a list of vertices representing a topological sort of the graph
	 */
	public static <V, E> List<Vertex<V>> TopologicalSort(IGraph<V, E> g) {

		// check if directed graph
		if (!g.isDirectedGraph()) {
			throw new ImproperGraphException();
		}

		// empty stack to contain sorted vertices
		Stack<Vertex<V>> stack = new Stack<Vertex<V>>();

		// map to track visits
		HashMap<Vertex<V>, String> visited = new HashMap<Vertex<V>, String>();

		// initialize the map
		List<Vertex<V>> vertexList = g.getVertices();
		for (int i = 0; i < vertexList.size(); i++) {
			visited.put(vertexList.get(i), "Not Visited");
		}

		// visit all the vertices
		for (int i = 0; i < vertexList.size(); i++) {
			if (visited.get(vertexList.get(i)) != "Visited") {
				topoHelper(visited, vertexList.get(i), stack, g);
			}
		}

		// convert stack to list
		List<Vertex<V>> sortedGraph = new ArrayList<Vertex<V>>();
		while (stack.size() > 0) {
			sortedGraph.add(stack.pop());
		}
		return sortedGraph;
	}

	// helper method for topological sort
	private static <V, E> void topoHelper(HashMap<Vertex<V>, String> visited, Vertex<V> vertex, Stack<Vertex<V>> stack,
			IGraph<V, E> graph) {

		// check for cycles
		if (visited.get(vertex) == "Visiting") {
			throw new ImproperGraphException();
		} else if (visited.get(vertex) != "Visited") {

			// mark current vertex
			visited.put(vertex, "Visiting");

			// visit all neighbors
			List<Vertex<V>> neighborList = graph.getNeighbors(vertex.getVertexName());
			for (int i = 0; i < neighborList.size(); i++) {
				topoHelper(visited, neighborList.get(i), stack, graph);
			}

			// mark current vertex
			visited.put(vertex, "Visited");
			stack.push(vertex);
		}
	}

	/**
	 * This method uses Kruscal's algorithm to return the minimum spanning of a
	 * graph.
	 * 
	 * @param g
	 *            a graph to perform Kruscal's algorithm for minimum spanning
	 *            tree on
	 * @return a graph representing a minimum spanning tree of the input graph
	 */
	public static <V, E extends IWeight> IGraph<V, E> Kruscal(IGraph<V, E> g) {

		// make list of vertex names
		ArrayList<String> stringList = new ArrayList<String>();
		for (Vertex<V> vertex : g.getVertices()) {
			stringList.add(vertex.getVertexName());
		}

		// create subsets
		Subset subset = new Subset(stringList);

		// create new priority queue with edge comparator
		Comparator<Edge<E>> comparator = new EdgeComparator<E>();
		PriorityQueue<Edge<E>> queue = new PriorityQueue<Edge<E>>(g.getEdges().size(), comparator);
		Graph<V, E> result = new Graph<V, E>();
		int count = 0;

		// sort edges by edge weight using quicksort
		ArrayList<Edge<E>> sortedList = (ArrayList<Edge<E>>) g.getEdges();
		quickSort(sortedList, 0, sortedList.size() - 1);

		// build priority queue
		for (Edge<E> edge : sortedList) {
			queue.add(edge);
		}

		// add vertices to result graph
		for (Vertex<V> vertex : g.getVertices()) {
			result.addVertex(vertex.getVertexName());
		}

		// go through vertices
		while (count < g.getVertices().size() - 1 || !(queue.isEmpty())) {
			Edge<E> minEdge = queue.poll();

			// if vertices aren't in the same subset, add to result
			if (!subset.find(minEdge.getVertexName1()).equals(subset.find(minEdge.getVertexName2()))) {
				result.addEdge(minEdge.getVertexName1(), minEdge.getVertexName2(), minEdge.getEdgeData());
				count++;
				subset.union(minEdge.getVertexName1(), minEdge.getVertexName2());
			}
		}

		return (IGraph<V, E>) result;
	}

	// quicksort algorithm to sort edges by IWeight
	private static <V, E extends IWeight> ArrayList<Edge<E>> quickSort(ArrayList<Edge<E>> list, int lowerIndex, int upperIndex) {
		int i = lowerIndex;
		int j = upperIndex;

		Edge<E> pivot = list.get(lowerIndex + (upperIndex - lowerIndex) / 2);
		while (i <= j) {
			while (list.get(i).getEdgeData().getWeight() < pivot.getEdgeData().getWeight()) {
				i++;
			}
			while (list.get(j).getEdgeData().getWeight() > pivot.getEdgeData().getWeight()) {
				j--;
			}
			if (i <= j) {
				Edge<E> temp = list.get(i);
				list.set(i, list.get(j));
				list.set(j, temp);
				i++;
				j--;
			}
		}

		if (lowerIndex < j) {
			quickSort(list, lowerIndex, j);
		}
		if (i < upperIndex) {
			quickSort(list, i, upperIndex);
		}

		return null;
	}

	/**
	 * This method uses Dijkstra's algorithm to return a list of edges
	 * representing the shortest path between two vertices in a graph.
	 * 
	 * @param g
	 *            a graph to find the shortest path in
	 * @param vertexStart
	 *            a vertex to begin at
	 * @param vertexEnd
	 *            a vertex to end at
	 * @return a list of edges representing a path from the starting vertex to
	 *         the ending vertex
	 */
	public static <V, E extends IWeight> List<Edge<E>> ShortestPath(IGraph<V, E> g, String vertexStart, String vertexEnd) {

		// check for null inputs
		if (g == null || vertexStart == null || vertexEnd == null) {
			throw new NullPointerException();
		}

		// instantiate maps for distance and predecessor
		HashMap<String, Double> dist = new HashMap<String, Double>();
		HashMap<String, String> pred = new HashMap<String, String>();

		// instantiate lists for open and closed sets
		ArrayList<String> open = new ArrayList<String>();
		ArrayList<String> closed = new ArrayList<String>();

		// add all vertices as max distance and null predecessor
		for (Vertex<V> vertex : g.getVertices()) {
			dist.put(vertex.getVertexName(), Double.MAX_VALUE);
			pred.put(vertex.getVertexName(), null);
		}
		dist.put(vertexStart, 0.0);
		open.add(vertexStart);

		// while open set is not empty
		while (!open.isEmpty()) {

			// a = x such that x is in open and x has minimum cost
			String a = getMin(open, dist);
			closed.add(a);
			open.remove(a);

			// for each neighbor of a
			for (Vertex<V> vertex : g.getNeighbors(a)) {
				if (!closed.contains(vertex)) {
					Double alt = dist.get(a) + g.getEdge(a, vertex.getVertexName()).getEdgeData().getWeight();
					if (alt < dist.get(vertex.getVertexName())) {
						dist.put(vertex.getVertexName(), alt);
						open.add(vertex.getVertexName());
						pred.put(vertex.getVertexName(), a);
					}
				}
			}
		}

		// use predecessor list to build path
		ArrayList<Edge<E>> path = new ArrayList<Edge<E>>();
		String currentVertex = vertexEnd;
		while (pred.get(currentVertex) != null) {
			path.add(g.getEdge(pred.get(currentVertex), currentVertex));
			currentVertex = pred.get(currentVertex);
		}

		// reverse path and return
		Collections.reverse(path);
		return path;
	}

	// helper method to find vertex in open set with minimum cost
	private static String getMin(ArrayList<String> open, HashMap<String, Double> dist) {
		Double minDistance = Double.MAX_VALUE;
		String minVertex = null;

		// check each vertex in open set
		for (String vertex : open) {
			if (dist.get(vertex) < minDistance) {
				minVertex = vertex;
				minDistance = dist.get(vertex);
			}
		}

		if (minVertex == null) {
			minVertex = open.get(0);
		}

		return minVertex;
	}

	// inner class to represent a set
	public static class Subset {

		private HashMap<String, Node> nodeMap;

		// create set from list of vertices
		public Subset(List<String> vertices) {
			this.nodeMap = new HashMap<String, Node>();
			for (String vertex : vertices) {
				nodeMap.put(vertex, new Node(vertex));
			}
		}

		// union of two sets
		public void union(String set1, String set2) {
			Node x = find(this.nodeMap.get(set1));
			Node y = find(this.nodeMap.get(set2));
			x.parent = y;
		}

		// find node in set
		public String find(String s) {
			Node n = nodeMap.get(s);
			return this.find(n).data;
		}

		// helper method to find node set
		private Node find(Node n) {
			while (n != n.parent) {
				n = n.parent;
			}
			return n;
		}

		// node class used by subsets
		private class Node {
			protected Node parent;
			protected String data;

			public Node(String data) {
				this.parent = this;
				this.data = data;
			}
		}
	}

	// edge comparator to sort by edge weight
	public static class EdgeComparator<E extends IWeight> implements Comparator<Edge<E>> {

		@Override
		public int compare(Edge<E> x, Edge<E> y) {

			if (x.getEdgeData().getWeight() < y.getEdgeData().getWeight()) {
				return -1;
			}
			if (x.getEdgeData().getWeight() > y.getEdgeData().getWeight()) {
				return 1;
			}
			return 0;
		}
	}

	// exception for improper graph
	public final static class ImproperGraphException extends RuntimeException {
		private static final long serialVersionUID = 1L;
	}
}
