package cs311.hw8.graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

public class Graph<V, E> implements IGraph<V, E> {

	private boolean directed = false;
	private final HashMap<String, Vertex<V>> vmap = new HashMap<>();
	private final LinkedHashMap<Vertex<V>, LinkedHashMap<Vertex<V>, E>> edges = new LinkedHashMap<>();

	@Override
	public void setDirectedGraph() {
		directed = true;
	}

	@Override
	public void setUndirectedGraph() {
		if (directed) {
			directed = false;
			List<Edge<E>> edgelist = getEdges();
			Iterator<Edge<E>> iter = edgelist.iterator();
			while (iter.hasNext()) {
				Edge<E> ed = iter.next();
				setEdgeData(ed.getVertexName1(), ed.getVertexName2(), null);
				if (!isEdge(ed.getVertexName2(), ed.getVertexName1())) {
					addEdge(ed.getVertexName1(), ed.getVertexName2());
				}
			}
		}
	}

	@Override
	public boolean isDirectedGraph() {
		return directed;
	}

	@Override
	public void addVertex(String vertexName) throws DuplicateVertexException {
		if (vmap.containsKey(vertexName))
			throw new DuplicateVertexException();
		vmap.put(vertexName, new Vertex<>(vertexName, null));
		edges.put(vmap.get(vertexName), new LinkedHashMap<>());
	}

	@Override
	public void addVertex(String vertexName, V vertexData) throws DuplicateVertexException {
		if (vmap.containsKey(vertexName))
			throw new DuplicateVertexException();
		vmap.put(vertexName, new Vertex<>(vertexName, vertexData));
		edges.put(vmap.get(vertexName), new LinkedHashMap<>());
	}

	@Override
	public void addEdge(String vertex1, String vertex2) throws DuplicateEdgeException, NoSuchVertexException {
		Vertex<V> v1 = vmap.get(vertex1);
		Vertex<V> v2 = vmap.get(vertex2);
		if (v1 == null || v2 == null)
			throw new NoSuchVertexException();
		if (edges.get(v1).containsKey(v2))
			throw new DuplicateEdgeException();
		edges.get(v1).put(v2, null);
	}

	@Override
	public void addEdge(String vertex1, String vertex2, E edgeData) throws DuplicateEdgeException, NoSuchVertexException {
		Vertex<V> v1 = vmap.get(vertex1);
		Vertex<V> v2 = vmap.get(vertex2);
		if (v1 == null || v2 == null)
			throw new NoSuchVertexException();
		if (edges.get(v1).containsKey(v2))
			throw new DuplicateEdgeException();
		edges.get(v1).put(v2, edgeData);
	}

	@Override
	public V getVertexData(String vertexName) throws NoSuchVertexException {
		if (!vmap.containsKey(vertexName))
			throw new NoSuchVertexException();
		return vmap.get(vertexName).getVertexData();
	}

	@Override
	public void setVertexData(String vertexName, V vertexData) throws NoSuchVertexException {
		if (!vmap.containsKey(vertexName))
			throw new NoSuchVertexException();
		vmap.remove(vertexName);
		vmap.put(vertexName, new Vertex<>(vertexName, vertexData));
	}

	@Override
	public E getEdgeData(String vertex1, String vertex2) throws NoSuchVertexException, NoSuchEdgeException {
		Vertex<V> v1 = vmap.get(vertex1);
		Vertex<V> v2 = vmap.get(vertex2);
		if (v1 == null || v2 == null)
			throw new NoSuchVertexException();
		if (!edges.get(v1).containsKey(v2))
			throw new NoSuchEdgeException();
		return edges.get(v1).get(v2);
	}

	@Override
	public void setEdgeData(String vertex1, String vertex2, E edgeData) throws NoSuchVertexException, NoSuchEdgeException {
		Vertex<V> v1 = vmap.get(vertex1);
		Vertex<V> v2 = vmap.get(vertex2);
		if (v1 == null || v2 == null)
			throw new NoSuchVertexException();
		if (!edges.get(v1).containsKey(v2))
			throw new NoSuchEdgeException();
		edges.get(v1).put(v2, edgeData);
	}

	@Override
	public Vertex<V> getVertex(String vertexName) throws NoSuchVertexException {
		if (!vmap.containsKey(vertexName))
			throw new NoSuchVertexException();
		return vmap.get(vertexName);
	}

	@Override
	public Edge<E> getEdge(String vertexName1, String vertexName2) {
		Vertex<V> v1 = vmap.get(vertexName1);
		Vertex<V> v2 = vmap.get(vertexName2);
		if (v1 == null || v2 == null)
			throw new NoSuchVertexException();
		if (!edges.get(v1).containsKey(v2))
			throw new NoSuchEdgeException();
		return new Edge<E>(vertexName1, vertexName2, edges.get(v1).get(v2));
	}

	@Override
	public List<Vertex<V>> getVertices() {
		return new ArrayList<>(vmap.values());
	}

	@Override
	public List<Edge<E>> getEdges() {
		ArrayList<Edge<E>> retval = new ArrayList<>();
		vmap.values().forEach(v -> edges.get(v).forEach((nv, ed) -> retval.add(new Edge<>(v.getVertexName(), nv.getVertexName(), ed))));
		return retval;
	}

	@Override
	public List<Vertex<V>> getNeighbors(String vertex) {
		Vertex<V> v1 = vmap.get(vertex);
		if (v1 == null)
			throw new NoSuchVertexException();
		ArrayList<Vertex<V>> retval = new ArrayList<>();
		edges.get(v1).forEach((nv, ed) -> retval.add(nv));
		return retval;
	}

	public boolean isEdge(String v1, String v2) {
		Vertex<V> v1v = vmap.get(v1);
		Vertex<V> v2v = vmap.get(v2);
		return isEdge(v1v, v2v);
	}

	public boolean isEdge(Vertex<?> v1, Vertex<?> v2) {
		return edges.get(v1).containsKey(v2);
	}

	public boolean isEdge(Edge<E> e) {
		return isEdge(e.getVertexName1(), e.getVertexName2());
	}

}