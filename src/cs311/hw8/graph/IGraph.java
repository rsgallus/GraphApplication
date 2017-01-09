package cs311.hw8.graph;

import java.util.List;
import java.util.Objects;

public interface IGraph<V, E>
{
    /**
     * Set the graph to be a directed graph.  Edge (x, y) is different than edge (y, x)
     */
    public void setDirectedGraph();
    
    
    /**
     * Set the graph to be an undirected graph.  Edge (x, y) is in the graph
     * if and only if edge (y, x) is in the graph.  Note that when implementing this
     * and there are already edges defined in the graph, care must be taken to 
     * resolve conflicts and inconsistencies in the overall implementation.
     */
    public void setUndirectedGraph();
    
    
    /**
     * 
     * @return true if the graph is directed.
     */
    public boolean isDirectedGraph();
    
    /**
     * Adds a vertex to the graph with name given by the vertexName.  vertexNames,
     * must be unique in the graph.
     * 
     * @param vertexName The unique name of the vertex.
     * 
     * @throws cs311.hw6.graph.IGraph.DuplicateVertexException 
     */
    public void addVertex( String vertexName ) throws DuplicateVertexException;
    
    
    /**
     * Adds a vertex to the graph with name given by the vertexName.  vertexNames,
     * must be unique in the graph.  The vertexData of generic type is associated with
     * this vertex.
     * 
     * @param vertexName
     * @param vertexData
     * @throws cs311.hw6.graph.IGraph.DuplicateVertexException 
     */
    public void addVertex( String vertexName, V vertexData ) throws DuplicateVertexException;
    
    
    /**
     * Adds an edge to the graph by specifying the two vertices that comprise the 
     * edge.  If the graph is undirected then edge (x, y) or edge (y, x) may be used
     * to add the single edge.  If the graph is undirected and edge (x,y) is added
     * followed by a subsequent edge (y, x), the later add would throw a 
     * DuplicateEdgeException.
     * 
     * @param vertex1 The first vertex in the edge.
     * @param vertex2 The second vertex in the edge.
     * 
     * @throws cs311.hw6.graph.IGraph.DuplicateEdgeException
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException 
     */
    public void addEdge( String vertex1, String vertex2 ) throws DuplicateEdgeException, NoSuchVertexException;
    
    
    /**
     * Adds an edge to the graph by specifying the two vertices that comprise the 
     * edge.  If the graph is undirected then edge (x, y) or edge (y, x) may be used
     * to add the single edge.  If the graph is undirected and edge (x,y) is added
     * followed by a subsequent edge (y, x), the later add would throw a 
     * DuplicateEdgeException.  The edgeDaa parameter is used to associate generic
     * edge data with the edge.
     * 
     * @param vertex1 The first vertex in the edge.
     * @param vertex2 The second vertex in the edge.
     * @param edgeData Thegeneric edge data.
     * 
     * @throws cs311.hw6.graph.IGraph.DuplicateEdgeException
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException 
     */
    public void addEdge( String vertex1, String vertex2, E edgeData ) throws DuplicateEdgeException, NoSuchVertexException;
    
    
    /**
     * Returns the generic vertex data associated with the specified vertex.  If no
     * vertex data is associated with the vertex, then null is returned.
     * 
     * @param vertexName  Name of vertex to get data for
     * 
     * @return The generic vertex data
     * 
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException 
     */
    public V getVertexData( String vertexName ) throws NoSuchVertexException;
    
    
    /**
     * Sets the generic vertex data of the specified vertex.
     * 
     * @param vertexName The name of the vertex.
     * 
     * @param vertexData The generic vertex data.
     * 
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException 
     */
    public void setVertexData( String vertexName, V vertexData) throws NoSuchVertexException;
    
    
    /**
     * Returns the generic edge data associated with the specified edge.  If no
     * edge data is associated with the edge, then null is returned.
     * 
     * @param vertex1 Vertex one of the edge.
     * @param vertex2 Vertex two of the edge.
     * 
     * @return The generic edge data.
     * 
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException
     * @throws cs311.hw6.graph.IGraph.NoSuchEdgeException 
     */
    public E getEdgeData( String vertex1, String vertex2 ) throws NoSuchVertexException, NoSuchEdgeException;
    
    
    /**
     * Sets the generic edge data of the specified edge.
     * 
     * @param vertex1  Vertex one of the edge.
     * @param vertex2  Vertex two of the edge.
     * 
     * @param edgeData The generic edge data.
     * 
     * @throws cs311.hw6.graph.IGraph.NoSuchVertexException
     * @throws cs311.hw6.graph.IGraph.NoSuchEdgeException 
     */
    public void setEdgeData( String vertex1, String vertex2, E edgeData ) throws NoSuchVertexException, NoSuchEdgeException;

    /**
     * Returns an encapsulated Vertex data type based on the vertex name
     * 
     * @param VertexName The name of the vertex.
     * 
     * @return The encapsulated vertex. 
     */
    public Vertex<V> getVertex( String VertexName );
    
    
    /**
     * Returns an encapsulated Edge data type based on the specified edge.
     * 
     * @param vertexName1 Vertex one of edge.
     * @param vertexName2 Vertex two of edge.
     * 
     * @return Encapsulated edge.
     */
    public Edge<E> getEdge( String vertexName1, String vertexName2);
    
    /**
     * Returns a list of all the vertices in the graph.
     * 
     * @return The List<Vertex> of vertices.
     */
    public List<Vertex<V>> getVertices();
    
    /**
     * Returns all the edges in the graph.
     * 
     * @return The List<Edge<E>> of edges. 
     */
    public List<Edge<E>> getEdges();
    
    /**
     * Returns all the neighbors of a specified vertex.
     * 
     * @param vertex The vertex to return neighbors for.
     * 
     * @return The list of vertices that are the neighbors of the specified vertex.
     */
    public List<Vertex<V>> getNeighbors(String vertex);
    
    
    /**
     * This class represents a vertex.  Do not change this class or this interface.
     * @param <V> 
     */
    public static final class Vertex<V>
    {
        private final String vertex;
        private final V vertexData;
        
        protected Vertex(String v, V d)
        {
            vertex = v;
            vertexData = d;
        }
        
        public String getVertexName()
        {
            return vertex;
        }
        
        public V getVertexData()
        {
            return vertexData;
        }
        
        @Override
        public int hashCode()
        {
            int hash = 3;
            hash = 23 * hash + Objects.hashCode(this.vertex);
            return hash;
        }

        @Override
        public boolean equals(Object obj)
        {
            if (this == obj)
            {
                return true;
            }
            if (obj == null)
            {
                return false;
            }
            if (getClass() != obj.getClass())
            {
                return false;
            }
            final Vertex<?> other = (Vertex<?>) obj;
            return Objects.equals(this.vertex, other.vertex);
        }
    }
    
    
    /**
     * This class represents an edge.  Do not change this class or interface.
     * @param <E> 
     */
    public static final class Edge<E>
    {
        private final String vertex1;
        private final String vertex2;
        private final E edgeData;
        
        protected Edge(String v1, String v2, E ed)
        {
            vertex1 = v1;
            vertex2 = v2;
            edgeData = ed;
        }
        
        public String getVertexName1()
        {
            return vertex1;
        }
        
        public String getVertexName2()
        {
            return vertex2;
        }
        
        public E getEdgeData()
        {
            return edgeData;
        }

        
        @Override
        public int hashCode()
        {
            int hash = 7;
            hash = 59 * hash + Objects.hashCode(this.vertex1);
            hash = 59 * hash + Objects.hashCode(this.vertex2);
            return hash;
        }

        @Override
        public boolean equals(Object obj)
        {
            if (this == obj)
            {
                return true;
            }
            if (obj == null)
            {
                return false;
            }
            if (getClass() != obj.getClass())
            {
                return false;
            }
            final Edge<?> other = (Edge<?>) obj;
            if (!Objects.equals(this.vertex1, other.vertex1))
            {
                return false;
            }
            return Objects.equals(this.vertex2, other.vertex2);
        }
    }
    
    // Exceptions used i nthe interface.
    
    public final static class DuplicateVertexException extends RuntimeException
    {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
    }

    public final static class DuplicateEdgeException extends RuntimeException
    {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
    }

    public final static class NoSuchVertexException extends RuntimeException
    {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
    }

    public final static class NoSuchEdgeException extends RuntimeException
    {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
    }
}
