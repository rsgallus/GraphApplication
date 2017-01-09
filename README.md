# GraphApplication

### Overview
This Java application parses Open Street Map xml data and finds the quickest route between two latitude/longitude locations.

### Graph Algorithms
Included in the application is an algorithm for topological sort, an implementation of Kruscal's algorithm for minimum spanning tree, and Djikstra's algorithm for shortest path.

### Usage
The first command line argument points to an Open Street Map xml data file which will be parsed to create a graph. The second command line argument point to a text file with a list of latitude/longitude locations. The application will find the closest street to each location and output a list of street names which represents the shortest street route between each given location.
