# Social Network Data Analysis

This repository uses social network data and researches the methods to analyze the data structure.

## Code overview

* Graph.java
  * Java Interface includes methods to make and deal with Graph structure. These methods are implemented in Java FBGraph Class.

* GraphLoader.java
  * Imports social network data. The addVertex and addEdge methods in FBGraph Class are used in loadGraph method in this Class. It needs to be calleed before analyzing social network data.

* FBGraph.java
  * Implements all the methods necessary for analyzing social network. It contains a cascade method which experiments how a new technology diffuses within a certain network. 
