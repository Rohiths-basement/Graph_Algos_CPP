#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <queue>
#include <climits>
#include <algorithm>
#include <stack>

//#include <pybind11> // include the pybind11 header

// A generic graph class template
template <typename T>
class Graph {
    // A struct to represent an edge
    struct Edge {
        T src; // source node
        T dest; // destination node
        int weight; // weight of the edge
        Edge(T s, T d, int w) : src(s), dest(d), weight(w) {} // constructor
    };

    // A map to store the adjacency list of each node
    std::map<T, std::list<Edge>> adj;

    // A boolean flag to indicate if the graph is directed or not
    bool directed;

public:
    // Constructor
    Graph(bool dir = false) : directed(dir) {}

    // A method to add a node to the graph
    void addNode(T node) {
        // If the node is not already in the graph, add it with an empty list
        if (adj.find(node) == adj.end()) {
            adj[node] = std::list<Edge>();
        }
    }

    // A method to add an edge to the graph
    void addEdge(T src, T dest, int weight = 1) {
        // Add the source and destination nodes if they are not already in the graph
        addNode(src);
        addNode(dest);

        // Create an edge object and add it to the adjacency list of the source node
        Edge e(src, dest, weight);
        adj[src].push_back(e);

        // If the graph is undirected, also add the reverse edge
        if (!directed) {
            Edge e2(dest, src, weight);
            adj[dest].push_back(e2);
        }
    }

    // A helper method to perform depth-first search and check for cycles
    bool dfs(T node, std::map<T, bool>& visited, std::map<T, bool>& recStack) {
        // Mark the node as visited and add it to the recursion stack
        visited[node] = true;
        recStack[node] = true;

        // Iterate over the adjacency list of the node
        for (auto e : adj[node]) {
            // If the destination node is not visited, recursively call dfs on it
            if (!visited[e.dest]) {
                if (dfs(e.dest, visited, recStack)) {
                    return true; // cycle found
                }
            }
            // If the destination node is already in the recursion stack, cycle found
            else if (recStack[e.dest]) {
                return true;
            }
        }

        // Remove the node from the recursion stack and return false
        recStack[node] = false;
        return false;
    }

    // A method to check if the graph contains a cycle
    bool hasCycle() {
        // Create two maps to store the visited and recursion stack status of each node
        std::map<T, bool> visited;
        std::map<T, bool> recStack;

        // Initialize both maps with false values for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            visited[it->first] = false;
            recStack[it->first] = false;
        }

        // Iterate over each node and call dfs on it if it is not visited
        for (auto it = adj.begin(); it != adj.end(); it++) {
            if (!visited[it->first]) {
                if (dfs(it->first, visited, recStack)) {
                    std::cout<<"A Cycle Found"<<std::endl;
                    return true; // cycle found
                }
            }
        }

        // No cycle found
        std::cout<<"No Cycle Found"<<std::endl;
        return false;
        
    }

    // A method to assign colours to the nodes of the graph
    void nodeColouring() {
        // Create a map to store the colour of each node
        std::map<T, int> colour;

        // Initialize the colour map with -1 for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            colour[it->first] = -1;
        }

        // Create a vector to store the available colours for each node
        std::vector<bool> available(adj.size(), true);

        // Iterate over each node and assign a colour to it
        for (auto it = adj.begin(); it != adj.end(); it++) {
            // Mark the colours of the adjacent nodes as unavailable
            for (auto e : it->second) {
                if (colour[e.dest] != -1) {
                    available[colour[e.dest]] = false;
                }
            }

            // Find the first available colour and assign it to the current node
            int c;
            for (c = 0; c < adj.size(); c++) {
                if (available[c]) {
                    break;
                }
            }
            colour[it->first] = c;

            // Reset the available vector for the next node
            for (auto e : it->second) {
                if (colour[e.dest] != -1) {
                    available[colour[e.dest]] = true;
                }
            }
        }

        // Print the colour of each node
        for (auto it = colour.begin(); it != colour.end(); it++) {
            std::cout << "Node " << it->first << " has colour " << it->second << "\n";
        }
    }

    // A method to assign colours to the edges of the graph
    void edgeColouring() {
        // Create a map to store the colour of each edge
        std::map<std::pair<T, T>, int> colour;

        // Initialize the colour map with -1 for each edge
        for (auto it = adj.begin(); it != adj.end(); it++) {
            for (auto e : it->second) {
                colour[{e.src, e.dest}] = -1;
            }
        }

        // Create a vector to store the available colours for each edge
        std::vector<bool> available(adj.size(), true);

        // Iterate over each node and assign colours to its edges
        for (auto it = adj.begin(); it != adj.end(); it++) {
            // Mark the colours of the incident edges as unavailable
            for (auto e : it->second) {
                if (colour[{e.src, e.dest}] != -1) {
                    available[colour[{e.src, e.dest}]] = false;
                }
            }

            // Find the first available colour and assign it to each edge
            int c;
            for (c = 0; c < adj.size(); c++) {
                if (available[c]) {
                    break;
                }
            }
            for (auto e : it->second) {
                if (colour[{e.src, e.dest}] == -1) {
                    colour[{e.src, e.dest}] = c;
                    c++;
                }
            }

            // Reset the available vector for the next node
            for (auto e : it->second) {
                if (colour[{e.src, e.dest}] != -1) {
                    available[colour[{e.src, e.dest}]] = true;
                }
            }
        }

        // Print the colour of each edge
        for (auto it = colour.begin(); it != colour.end(); it++) {
            std::cout << "Edge (" << it->first.first << ", " << it->first.second << ") has colour " << it->second << "\n";
        }
    }

    // A method to print the graph
    void print() {
        // Iterate over the map and print each node and its adjacency list
        for (auto it = adj.begin(); it != adj.end(); it++) {
            std::cout << it->first << ": ";
            for (auto e : it->second) {
                std::cout << "(" << e.dest << ", " << e.weight << ") ";
            }
            std::cout << "\n";
        }
    }

    // A helper method to perform breadth-first search and mark visited nodes
    void bfs(T node, std::map<T, bool>& visited) {
        // Create a queue to store the nodes to be visited
        std::queue<T> q;

        // Mark the node as visited and enqueue it
        visited[node] = true;
        q.push(node);

        // While the queue is not empty, dequeue a node and visit its adjacent nodes
        while (!q.empty()) {
            T u = q.front();
            q.pop();

            for (auto e : adj[u]) {
                // If the destination node is not visited, mark it as visited and enqueue it
                if (!visited[e.dest]) {
                    visited[e.dest] = true;
                    q.push(e.dest);
                }
            }
        }
    }

    // A method to detect and complete edges in the graph
    void detectAndCompleteEdges() {
        // Create a map to store the visited status of each node
        std::map<T, bool> visited;

        // Initialize the visited map with false for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            visited[it->first] = false;
        }

        // Iterate over each node and call bfs on it if it is not visited
        for (auto it = adj.begin(); it != adj.end(); it++) {
            if (!visited[it->first]) {
                bfs(it->first, visited);
            }
        }

        // Iterate over each node and check if there are any missing edges
        for (auto it = adj.begin(); it != adj.end(); it++) {
            for (auto jt = adj.begin(); jt != adj.end(); jt++) {
                // If the nodes are different and both are visited, check if there is an edge between them
                if (it->first != jt->first && visited[it->first] && visited[jt->first]) {
                    bool found = false;
                    for (auto e : it->second) {
                        if (e.dest == jt->first) {
                            found = true;
                            break;
                        }
                    }
                    // If there is no edge between them, add one with weight 0
                    if (!found) {
                        addEdge(it->first, jt->first, 0);
                    }
                }
            }
        }

        // Print the updated graph
        print();
    }

    // A method to find the number of connected components in the graph
    int connectedComponents() {
        // Create a map to store the visited status of each node
        std::map<T, bool> visited;

        // Initialize the visited map with false for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            visited[it->first] = false;
        }

        // Initialize a counter for the number of connected components
        int count = 0;

        // Iterate over each node and call bfs on it if it is not visited
        for (auto it = adj.begin(); it != adj.end(); it++) {
            if (!visited[it->first]) {
                bfs(it->first, visited);
                count++; // increment the counter for each connected component
            }
        }

        // Return the number of connected components
        return count;
    }

    // A method to find the graph centrality (Katz's) of each node
    void katzCentrality() {
        // Create a map to store the centrality of each node
        std::map<T, double> centrality;

        // Initialize the centrality map with 1 for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            centrality[it->first] = 1.0;
        }

        // Choose a constant alpha that is less than the inverse of the largest eigenvalue of the adjacency matrix
        double alpha = 0.1;

        // Choose a constant beta that is any positive value
        double beta = 1.0;

        // Choose a tolerance value for convergence
        double tol = 0.0001;

        // Create a boolean flag to indicate if the centrality values have converged
        bool converged = false;

        // Repeat until convergence
        while (!converged) {
            // Create a map to store the updated centrality values
            std::map<T, double> new_centrality;

            // Initialize the new centrality map with beta for each node
            for (auto it = adj.begin(); it != adj.end(); it++) {
                new_centrality[it->first] = beta;
            }

            // Iterate over each node and update its centrality value based on its adjacent nodes
            for (auto it = adj.begin(); it != adj.end(); it++) {
                for (auto e : it->second) {
                    new_centrality[e.dest] += alpha * centrality[e.src];
                }
            }

            // Check if the difference between the old and new centrality values is less than the tolerance for each node
            converged = true;
            for (auto it = adj.begin(); it != adj.end(); it++) {
                if (std::abs(new_centrality[it->first] - centrality[it->first]) > tol) {
                    converged = false;
                    break;
                }
            }

            // Copy the new centrality values to the old ones
            for (auto it = adj.begin(); it != adj.end(); it++) {
                centrality[it->first] = new_centrality[it->first];
            }
        }

        // Print the centrality of each node
        for (auto it = centrality.begin(); it != centrality.end(); it++) {
            std::cout << "Node " << it->first << " has Katz's centrality " << it->second << "\n";
        }
    }

    // A method to find the minimum spanning tree of the graph using Prim's algorithm
    void primMST() {
        // Create a map to store the parent of each node in the MST
        std::map<T, T> parent;

        // Create a map to store the key value of each node
        std::map<T, int> key;

        // Create a map to store the visited status of each node
        std::map<T, bool> visited;

        // Initialize the parent map with -1 for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            parent[it->first] = -1;
        }

        // Initialize the key map with a large value for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            key[it->first] = INT_MAX;
        }

        // Initialize the visited map with false for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            visited[it->first] = false;
        }

        // Choose any node as the starting node and set its key value to 0
        T start = adj.begin()->first;
        key[start] = 0;

        // Repeat until all nodes are visited
        while (true) {
            // Find the node with the minimum key value that is not visited
            T u = -1;
            int min_key = INT_MAX;
            for (auto it = adj.begin(); it != adj.end(); it++) {
                if (!visited[it->first] && key[it->first] < min_key) {
                    u = it->first;
                    min_key = key[it->first];
                }
            }

            // If no such node is found, break the loop
            if (u == -1) {
                break;
            }

            // Mark the node as visited
            visited[u] = true;

            // Iterate over the adjacent nodes and update their key and parent values if needed
            for (auto e : adj[u]) {
                if (!visited[e.dest] && e.weight < key[e.dest]) {
                    key[e.dest] = e.weight;
                    parent[e.dest] = u;
                }
            }
        }

        // Print the edges and weights of the MST
        int total_weight = 0;
        for (auto it = parent.begin(); it != parent.end(); it++) {
            if (it->second != -1) {
                std::cout << "Edge (" << it->second << ", " << it->first << ") has weight " << key[it->first] << "\n";
                total_weight += key[it->first];
            }
        }
        std::cout << "The total weight of the MST is " << total_weight << "\n";
    }

    // A helper method to find the representative of a node in a disjoint set
    T find(std::map<T, T>& parent, T node) {
        // If the node is its own parent, return it
        if (parent[node] == node) {
            return node;
        }
        // Otherwise, recursively find the parent of the node and update it
        parent[node] = find(parent, parent[node]);
        return parent[node];
    }

    // A helper method to union two nodes in a disjoint set
    void union_nodes(std::map<T, T>& parent, std::map<T, int>& rank, T x, T y) {
        // Find the representatives of x and y
        T xroot = find(parent, x);
        T yroot = find(parent, y);

        // If they are already in the same set, do nothing
        if (xroot == yroot) {
            return;
        }

        // Otherwise, compare their ranks and attach the smaller one to the larger one
        if (rank[xroot] < rank[yroot]) {
            parent[xroot] = yroot;
        }
        else if (rank[xroot] > rank[yroot]) {
            parent[yroot] = xroot;
        }
        else {
            // If they have the same rank, arbitrarily choose one and increment its rank
            parent[yroot] = xroot;
            rank[xroot]++;
        }
    }

    // A method to find the minimum spanning tree of the graph using Kruskal's algorithm
    void kruskalMST() {
        // Create a vector to store all the edges in the graph
        std::vector<Edge> edges;

        // Iterate over the adjacency list and add each edge to the vector
        for (auto it = adj.begin(); it != adj.end(); it++) {
            for (auto e : it->second) {
                edges.push_back(e);
            }
        }

        // Sort the edges by their weights in ascending order
        std::sort(edges.begin(), edges.end(), [](Edge a, Edge b) {
            return a.weight < b.weight;
        });

        // Create a map to store the parent of each node in a disjoint set
        std::map<T, T> parent;

        // Create a map to store the rank of each node in a disjoint set
        std::map<T, int> rank;

        // Initialize the parent map with each node as its own parent
        for (auto it = adj.begin(); it != adj.end(); it++) {
            parent[it->first] = it->first;
        }

        // Initialize the rank map with 0 for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            rank[it->first] = 0;
        }

        // Create a vector to store the edges of the MST
        std::vector<Edge> mst;

        // Create a counter for the number of edges in the MST
        int count = 0;

        // Iterate over the sorted edges and add them to the MST if they don't form a cycle
        for (auto e : edges) {
            // Find the representatives of the source and destination nodes
            T x = find(parent, e.src);
            T y = find(parent, e.dest);

            // If they are different, they don't form a cycle and can be added to the MST
            if (x != y) {
                mst.push_back(e);
                count++;

                // Union the two nodes in the disjoint set
                union_nodes(parent, rank, x, y);
            }

            // If the number of edges in the MST is equal to the number of nodes minus one, break the loop
            if (count == adj.size() - 1) {
                break;
            }
        }

        // Print the edges and weights of the MST
        int total_weight = 0;
        for (auto e : mst) {
            std::cout << "Edge (" << e.src << ", " << e.dest << ") has weight " << e.weight << "\n";
            total_weight += e.weight;
        }
        std::cout << "The total weight of the MST is " << total_weight << "\n";
    }

    // A method to perform iterative depth first search (DFS) from a given node
    void iterativeDFS(T start) {
        // Create a stack to store the nodes to be visited
        std::stack<T> s;

        // Create a map to store the visited status of each node
        std::map<T, bool> visited;

        // Initialize the visited map with false for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            visited[it->first] = false;
        }

        // Push the start node to the stack and mark it as visited
        s.push(start);
        visited[start] = true;

        // While the stack is not empty, pop a node and visit its adjacent nodes
        while (!s.empty()) {
            T u = s.top();
            s.pop();

            // Print the node
            std::cout << u << " ";

            // Iterate over the adjacency list of the node in reverse order
            for (auto it = adj[u].rbegin(); it != adj[u].rend(); it++) {
                // If the destination node is not visited, push it to the stack and mark it as visited
                if (!visited[it->dest]) {
                    s.push(it->dest);
                    visited[it->dest] = true;
                }
            }
        }

        // Print a new line
        std::cout << "\n";
    }

    // A method to perform uniform cost search (UCS) from a given source node to a given destination node
    void ucs(T src, T dest) {
        // Create a priority queue to store the nodes to be visited along with their costs
        std::priority_queue<std::pair<int, T>, std::vector<std::pair<int, T>>, std::greater<std::pair<int, T>>> pq;

        // Create a map to store the cost of each node from the source node
        std::map<T, int> cost;

        // Create a map to store the parent of each node in the shortest path tree
        std::map<T, T> parent;

        // Initialize the cost map with a large value for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            cost[it->first] = INT_MAX;
        }

        // Push the source node to the priority queue with cost 0 and update its cost and parent
        pq.push({0, src});
        cost[src] = 0;
        parent[src] = -1;

        // While the priority queue is not empty, pop a node and visit its adjacent nodes
        while (!pq.empty()) {
            int u_cost = pq.top().first;
            T u = pq.top().second;
            pq.pop();

            // If the node is the destination node, break the loop
            if (u == dest) {
                break;
            }

            // Iterate over the adjacency list of the node
            for (auto e : adj[u]) {
                // If the cost of the destination node can be reduced by going through this edge, update its cost and parent and push it to the priority queue
                if (cost[e.dest] > u_cost + e.weight) {
                    cost[e.dest] = u_cost + e.weight;
                    parent[e.dest] = u;
                    pq.push({cost[e.dest], e.dest});
                }
            }
        }

            // If the destination node is not reachable from the source node, print a message
        if (cost[dest] == INT_MAX) {
            std::cout << "There is no path from " << src << " to " << dest << "\n";
            return;
        }

        // Create a vector to store the shortest path from the source node to the destination node
        std::vector<T> path;

        // Trace the path from the destination node to the source node using the parent map
        T node = dest;
        while (node != -1) {
            path.push_back(node);
            node = parent[node];
        }

        // Reverse the path vector
        std::reverse(path.begin(), path.end());

        // Print the path and its cost
        std::cout << "The shortest path from " << src << " to " << dest << " is:\n";
        for (int i = 0; i < path.size(); i++) {
            std::cout << path[i];
            if (i < path.size() - 1) {
                std::cout << " -> ";
            }
        }
        std::cout << "\n";
        std::cout << "The cost of the path is " << cost[dest] << "\n";
    }

    // A method to perform A* search from a given source node to a given destination node
    void aStar(T src, T dest) {
        // Create a priority queue to store the nodes to be visited along with their costs and heuristics
        std::priority_queue<std::tuple<int, int, T>, std::vector<std::tuple<int, int, T>>, std::greater<std::tuple<int, int, T>>> pq;

        // Create a map to store the cost of each node from the source node
        std::map<T, int> cost;

        // Create a map to store the heuristic of each node to the destination node
        std::map<T, int> heuristic;

        // Create a map to store the parent of each node in the shortest path tree
        std::map<T, T> parent;

        // Initialize the cost map with a large value for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            cost[it->first] = INT_MAX;
        }

        // Initialize the heuristic map with some estimate for each node
        for (auto it = adj.begin(); it != adj.end(); it++) {
            //heuristic[it->first] = someEstimate(it->first, dest);
            heuristic[it->first] = (it->first + dest);
        }

        // Push the source node to the priority queue with cost 0 and heuristic 0 and update its cost and parent
        pq.push({0, 0, src});
        cost[src] = 0;
        parent[src] = -1;

        // While the priority queue is not empty, pop a node and visit its adjacent nodes
        while (!pq.empty()) {
            int u_cost = std::get<0>(pq.top());
            int u_heuristic = std::get<1>(pq.top());
            T u = std::get<2>(pq.top());
            pq.pop();

            // If the node is the destination node, break the loop
            if (u == dest) {
                break;
            }

            // Iterate over the adjacency list of the node
            for (auto e : adj[u]) {
                // If the cost of the destination node can be reduced by going through this edge, update its cost and parent and push it to the priority queue
                if (cost[e.dest] > u_cost + e.weight) {
                    cost[e.dest] = u_cost + e.weight;
                    parent[e.dest] = u;
                    pq.push({cost[e.dest], heuristic[e.dest], e.dest});
                }
            }
        }

            // If the destination node is not reachable from the source node, print a message
        if (cost[dest] == INT_MAX) {
            std::cout << "There is no path from " << src << " to " << dest << "\n";
            return;
        }

        // Create a vector to store the shortest path from the source node to the destination node
        std::vector<T> path;

        // Trace the path from the destination node to the source node using the parent map
        T node = dest;
        while (node != -1) {
            path.push_back(node);
            node = parent[node];
        }

        // Reverse the path vector
        std::reverse(path.begin(), path.end());

        // Print the path and its cost
        std::cout << "The shortest path from " << src << " to " << dest << " is:\n";
        for (int i = 0; i < path.size(); i++) {
            std::cout << path[i];
            if (i < path.size() - 1) {
                std::cout << " -> ";
            }
        }
        std::cout << "\n";
        std::cout << "The cost of the path is " << cost[dest] << "\n";
    }


};
//Adjacenct Nodes in a directed graph in this case is a symmetric relation ie, nodes are adjacent only if the edges are both ways

//Prim's Algorithm gives us the MST starting from the root node and then does it's calculations whereas
//Kruskal's Algorithm gives us the MST starting from the least weighted edge and then does it's calculations

int main() {
    // Create a directed graph object
    Graph<int> g(true);

    // Add some edges to the graph
    g.addEdge(1, 2, 10);
    g.addEdge(1, 3, 5);
    g.addEdge(3, 1, 6);
    g.addEdge(2, 4, 7);
    g.addEdge(3, 4, 3);
    g.addEdge(4, 5, 2);

    // Print the graph
    g.print();

    g.hasCycle();

    g.edgeColouring();
    g.nodeColouring();

    std::cout<<"Number of connected components: "<<g.connectedComponents()<<std::endl;

    g.katzCentrality();

    std::cout<<"Weights of the Edges in the MST (using Prim's Algorithm) is as below:"<<std::endl;
    g.primMST();

    std::cout<<"Weights of the Edges in the MST (using Kruskal's Algorithm) is as below:"<<std::endl;
    g.kruskalMST();

    std::cout << "Iterative DFS from node 1 (below is the list of the nodes in the order in which the DFS goes through ):\n";
    g.iterativeDFS(1);

    // Perform UCS from node 1 to node 5
    std::cout << "UCS from node 1 to node 5:\n";
    g.ucs(1, 5);


    std::cout << "A* Search from node 2 to node 5:\n";
    g.aStar(2, 5);

    std::cout<<"After detecting and completing edges, we get the following: "<<std::endl;
    g.detectAndCompleteEdges();

    return 0;
}
