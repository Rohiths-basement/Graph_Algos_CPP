# Graph_Algos_CPP
The **Primary objective** of this was to develop a generic graph library (File) in C++ that would provide support for various Graph Algorithms and types. 
The **Secondary objective** (and more interesting one) was to a) type out the code for each algorithm manually and b) get ChatGPT (GPT-3) to generate the code individually for each algo as well, and then, compare the two to better understand the differences in methodology & approach, and to combine them both to give a more refined code. 

In this, algorithms are written in terms of data types that are to be specified later which are then instantiated when needed for specific types provided as parameters. 

The things included in this are: 
* Graph Creation
* Support for Directed and Undirected graphs
* Support for Node and Edge Creation
* Algorithm Execution 
* Detection of cycles
* Node colouring
* Edge colouring
* Detect and complete edges
* Connected components Detection 
* Graph centrality (Katz's)
* Prim's algorithm
* Kruskal's algorithm
* Iterative Depth First Search (DFS)
* Uniform Cost Search (UCS)
* A* Search


This project can be further worked on to make it callable from python. This can be done by providing python bindings (Using Pybind11 or a similar library) 

The input graph has been hardcoded and below is the Graph. 
![Graph](https://github.com/Rohiths-basement/Graph_Algos_CPP/assets/83353135/42d3dc77-1abe-4e6c-a965-d9f580789813)


To run this, just run it as you would a normal CPP file ie, - Compile it and then run it. On running the code, the algorithms get executed on their own and their respective results are displayed. 
 
