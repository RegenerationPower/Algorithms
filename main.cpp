#include <iostream>
#include <list>
#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <algorithm>
#include <fstream>
#include <sstream>
using namespace std;

// !!! IMPORTANT !!!
// There are some comments inside method blocks. You can use them for testing
// Last updated: 17.07.2023

// Class that represents graph and include some methods that do various operations on vertices, edges and nodes
class Graph 
{
    private:
        string name;
        int numVertices;
        vector<list<int>> adjList;
        vector<double> nodeValues;
        vector<vector<double>> edgeValues;

    public:
         Graph(string graphName, int vertices): name(graphName), numVertices(vertices), adjList(vertices), nodeValues(vertices), edgeValues(vertices, vector<double>(vertices)) {}
        // Add custom edge to the graph
        void add_edge(int x, int y)
        {
            //cout << name << ": Edge " << y <<" to Vertex " << x << " added" <<endl;
            adjList[x].push_front(y);
            adjList[y].push_front(x); 
        }
        // Delete some edge
        void delete_edge(int x, int y)
        {
            //cout << name << ": Edge " << y <<" from Vertex " << x << " deleted" <<endl;
            adjList[x].remove(y);
            adjList[y].remove(x);
        }
        // Check if x and y edges are adjacent ***Not used, only for testing
        void adjacent_edges(int x, int y)
        {
            bool found = false;
            for (auto node : adjList[x]) {
                if (node == y) 
                {
                    found = true;
                    break;
                }
            }
            
            if (found)
                cout << name  << ": Edge " << x << " is adjacent to edge " << y << endl;
            else
                cout << name  << ": Edge " << x << " is not adjacent to edge " << y << endl;
        }
        // Return the list of all edges that are neighbours to the x
        list<int> neighbour_edges(int x)
        {
            // cout << name << ": Neighbors of vertex " << x << ": ";
            // for (auto neighbor : adjList[x]) 
            //     cout << neighbor << " ";
            // cout << endl;
            return adjList[x];
        }
        // Set some custom node value a ***Not used, only for testing
        void set_node_value(int x, double a)
        {
            //cout << name << ": Value " << a <<" to node " << x << " added" <<endl;
            nodeValues[x] = a;
        }
        // Set some custom edge value v
        void set_edge_value(int x, int y, double v)
        {
            //cout << name << ": Value " << v << " to edge " << y <<" from vertex " << x << " added" <<endl;
            edgeValues[x][y] = v;
            edgeValues[y][x] = v;
        }
        // Return some node value ***Not used, only for testing
        double get_node_value(int x)
        {
            //cout << name << ": Value of node " << x << ": " << nodeValues[x] << endl;
            return nodeValues[x];
        }
        // Return some edge value
        double get_edge_value(int x, int y)
        {
            //cout << name << ": Value of edge[" << x << "][" << y << "]: " << edgeValues[x][y] << endl;
            return edgeValues[x][y];
        }
        // Return total number of Vertices
        int get_vertices()
        {
            //cout << name << ": Number of vertices: " << numVertices << endl;
            return numVertices;
        }
        // Return total number of edges ***Not used, only for testing
        int get_edges()
        {
            int count = 0;
            for (auto list : adjList) 
                count += list.size();
            //cout << name << ": Number of edges: " << count << endl;
            // Maybe it is better to use count without divide sometimes
            return count / 2;
        }
        // Used to print graphs with Vertex and edges connect to them
        void print_graph() 
        {
            for (int i = 0; i < numVertices; i++) 
            {
                cout << name << ": Vertex " << i << ": ";
                for (auto x : adjList[i]) 
                    cout << x << " ";
                cout << endl;
                }
        }
        // Prim algorithm
        pair<double, vector<pair<int, int>>> prim_algorithm() 
        {
            vector<bool> visited(numVertices, false); 
            vector<pair<int, int>> tree; 
            double weight = 0;
            visited[0] = true;

            while (tree.size() < numVertices - 1) 
            {
                double minWeight = numeric_limits<double>::infinity();
                int minU, minV;

                for (int u = 0; u < numVertices; u++) 
                {
                    if (visited[u]) 
                    {
                        for (int v : adjList[u]) 
                        {
                            if (!visited[v] && edgeValues[u][v] < minWeight) 
                            {
                                minWeight = edgeValues[u][v];
                                minU = u;
                                minV = v;
                            }
                        }
                    }
                }
                visited[minV] = true;
                tree.push_back({ minU, minV });
                weight += minWeight;
            }
            return { weight, tree };
        }
        // Kruskal algorithm
        pair<double, vector<pair<int, int>>> kruskal_algorithm()
        {
            vector<pair<double, pair<int, int>>> edges;

            for (int u = 0; u < numVertices; u++)
            {
                for (int v : adjList[u])
                {
                    if (u < v)
                    {
                        edges.push_back({ edgeValues[u][v], {u, v} });
                    }
                }
            }

            sort(edges.begin(), edges.end());
            vector<int> parent(numVertices);

            for (int i = 0; i < numVertices; i++)
                parent[i] = i;

            vector<pair<int, int>> tree;
            double weight = 0;

            for (auto edge : edges)
            {
                double edgeWeight = edge.first;
                int u = edge.second.first;
                int v = edge.second.second;
                int parentU = find_parent(parent, u);
                int parentV = find_parent(parent, v);

                if (parentU != parentV)
                {
                    tree.push_back({ u, v });
                    weight += edgeWeight;
                    parent[parentU] = parentV;
                }
            }
            return { weight, tree };
        }
        // Additional method for Kruskal's algorithm to find parent
        int find_parent(vector<int>& parent, int vertex)
        {
            if (parent[vertex] == vertex)
                return vertex;
            return find_parent(parent, parent[vertex]);
        }

        // operator << overloading for clearer output
        friend ostream &operator<<(ostream& os, Graph &graph)
        {
            os << graph.name;
            return os;
        }

};
// Class that represents priority queue and include some methods that do various operations with it
class PriorityQueue
{
    private:
        priority_queue<int, vector<int>, greater<int>> pq;
    public:
        // Delete top element of the priority queue
        void delete_top()
        {
            if (!pq.empty())
            {
                int topElement = pq.top();
                pq.pop();
                //cout << "Removed element with priority " << topElement << " from the queue" << endl;
            }
            // else
            //     cout << "Priority queue is empty" << endl;
        }
        // Check if there is an element in queue ***Not used, only for testing
        void contains_element(int element)
        {
            priority_queue<int, vector<int>, greater<int>>  tempPQ = pq;
            bool found = false;

            while (!tempPQ.empty())
            {
                int currentElement = tempPQ.top();
                tempPQ.pop();
                if (currentElement == element)
                {
                    found = true;
                    break;
                }
            }

            if (found)
                cout << "Element " << element << " is in the priority queue" << endl;
            else
                cout << "Element " << element << " is not in the priority queue" << endl;
        }
        // Add element to the priority queue
        void insert_element(int element)
        {
            pq.push(element);
            //cout << "Inserted element with priority " << element << " into the queue" << endl;
        }
        // Return the top element of the priority queue
        int get_top()
        {
            if (!pq.empty())
            {
                int topElement = pq.top();
                //cout << "Top element of the queue: " << topElement << endl;
                return topElement;
            }
            else
            {
                //cout << "Priority queue is empty" << endl;
                return 0;
            }
        }
        // Return size of the queue
        int get_size()
        {
            int queueSize = pq.size();
            //cout << "Size of the queue: " << queueSize << endl;
            return queueSize;
        }

};
// Class that represents shortest path and dijkstra algorithm
class ShortestPath 
{
    private:
        Graph &graph;
        PriorityQueue &priorityQueue;
        vector<double> distance;
        vector<int> prev;

    public:
         ShortestPath(Graph& g, PriorityQueue& pq) : graph(g), priorityQueue(pq), distance(g.get_vertices(), numeric_limits<double>::infinity()), prev(g.get_vertices(), -1) {}
        // Return the shortest path from u to w
        list<int> path(int u, int w) 
        {
            dijkstra(u); 
            return buildPath(w);
        }
        // Return the shortest path size from u to w
        int path_size(int u, int w) 
        {
            dijkstra(u);
            return distance[w];
        }
        // Used for printing the shortest path from u to w
        void print_shortest_path(int u, int w) 
        {
            list<int> shortestPath = path(u, w);
            cout << "Shortest path size: " << path_size(u, w) << endl;
            cout << "Shortest path: ";
            for (int vertex : shortestPath)
                cout << vertex << " ";
            cout << endl;
        }
    // Main logic of the dijkstra algorithm
    private:
        void dijkstra(int topVertex) 
        {
            distance[topVertex] = 0;
            priorityQueue.insert_element(topVertex);

            while (priorityQueue.get_size() > 0) 
            {
                int u = priorityQueue.get_top();
                priorityQueue.delete_top();

                list<int> neighbors = graph.neighbour_edges(u);

                for (int v : neighbors) 
                {
                    double edgeWeight = graph.get_edge_value(u, v);
                    double newdistance = distance[u] + edgeWeight;

                    if (newdistance < distance[v]) 
                    {
                        distance[v] = newdistance;
                        prev[v] = u;
                        priorityQueue.insert_element(v);
                    }
                }
            }
        }

        list<int> buildPath(int w) 
        {
            list<int> shortestPath;
            int currentVertex = w;

            while (currentVertex != -1) 
            {
                shortestPath.push_front(currentVertex);
                currentVertex = prev[currentVertex];
            }
            return shortestPath;
        }
};
// Generate random graph using #include <random> and choose edges that are < density in further calculations
Graph generate_random_graph(int numVertices, double density, double minDistance, double maxDistance)
{
    Graph graph("RandomGraph", numVertices);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist(minDistance, maxDistance);

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = i + 1; j < numVertices; j++)
        {
            double randomProbability = static_cast<double>(rand()) / RAND_MAX;
            if (randomProbability < density)
            {
                double randomDistance = dist(gen);
                graph.add_edge(i, j);
                graph.set_edge_value(i, j, randomDistance);
            }
        }
    }

    return graph;
}
// Calculate an average shortest path from i to source
double calculate_average_shortest_path(Graph& graph, int source, PriorityQueue& priorityQueue)
{
    int numVertices = graph.get_vertices();
    double totalShortestPath = 0;
    int validPaths = 0;

    for (int i = 0; i < numVertices; i++)
    {
        if (i == source)
            continue;

        ShortestPath shortestPath(graph, priorityQueue);
        int pathSize = shortestPath.path_size(source, i);

        if (pathSize != numeric_limits<int>::max())
        {
            totalShortestPath += pathSize;
            validPaths++;
        }
    }

    if (validPaths > 0)
        return totalShortestPath / validPaths;
    else
        return 0;
}
// Read graph from file
Graph read_graph_from_file(const string &filename)
{
    ifstream inputFile(filename);
    if (!inputFile)
    {
        cerr << "Failed to open file: " << filename << endl;
        exit(1);
    }

    string line;
    if (!getline(inputFile, line))
    {
        cerr << "Failed to read graph size from file: " << filename << endl;
        exit(1);
    }

    int numVertices = stoi(line);
    Graph graph("graph", numVertices);

    while (getline(inputFile, line))
    {
        if (line.empty())
            continue;

        int x, y, weight;
        stringstream iss(line);
        if (!(iss >> x >> y >> weight))
        {
            cerr << "Invalid line format in file: " << filename << endl;
            exit(1);
        }

        graph.add_edge(x, y);
        graph.set_edge_value(x, y, weight);
    }

    inputFile.close();
    return graph;
}


int main(int argc, char const *argv[])
{
    cout << "During this assignment, I learned some basic aspects of the C++ programming language. Such as:" << endl
     << "Classes and Objects: I created Graph, PriorityQueue and ShortestPath classes and learned how to work with objects of these classes, create them, use constructors and class methods." << endl
     << "Reference Passing: I used object reference passing in constructors and class methods to avoid copying and working with the same object." << endl
     << "Using standard libraries: You used various classes and functions from C++ standard libraries, such as list, vector, queue, numeric_limits, random_device, mt19937, etc." << endl
     << "Memory management: I worked with dynamic memory when I created graphs and a queue using new and delete operators." << endl
     << "Working with arrays and vectors: I created and manipulated arrays and vectors to store data such as list of adjacent vertices, edge weights, etc." << endl
     << "Using Loops and Conditional Statements: I used for and while loops to iterate through the elements, as well as conditional if statements to implement various conditional checks." << endl
     << "Data output: I used streaming output (iostream) to output data to the console." << endl
     << "Random Number Generation: I used functions from the C++ standard library to generate random numbers and use them to generate random graphs." << endl
     << "I also began to understand graphs, their structure and algorithms for finding the fastest paths." << endl << endl;
    // Create custom graph and find its shortest path
    Graph graph("graph", 5);
    graph.add_edge(0, 1);
    graph.add_edge(0, 2);
    graph.add_edge(1, 2);
    graph.add_edge(1, 3);
    graph.add_edge(2, 3);
    graph.add_edge(2, 4);
    graph.add_edge(3, 4);
    graph.print_graph();
    graph.set_edge_value(0, 1, 2.0);
    graph.set_edge_value(0, 2, 1.0);
    graph.set_edge_value(1, 2, 3.0);
    graph.set_edge_value(1, 3, 2.0);
    graph.set_edge_value(2, 3, 1.0);
    graph.set_edge_value(2, 4, 4.0);
    graph.set_edge_value(3, 4, 3.0);

    PriorityQueue priorityQueue;

    ShortestPath shortestPath(graph, priorityQueue);
    shortestPath.print_shortest_path(0, 4);

    // Create random graph and find its average shortest path for some density
    int numVertices = 50;
    double density1 = 0.2;
    double density2 = 0.4;
    double minDistance = 1.0;
    double maxDistance = 10.0;

    Graph graph1 = generate_random_graph(numVertices, density1, minDistance, maxDistance);
    Graph graph2 = generate_random_graph(numVertices, density2, minDistance, maxDistance);

    cout << "Average shortest path for density 0.2: " << calculate_average_shortest_path(graph1, 0, priorityQueue) << endl;
    cout << "Average shortest path for density 0.4: " << calculate_average_shortest_path(graph2, 0, priorityQueue) << endl;

    // Read from file and using prim's and kruskal's algorithms
    Graph graphFile = read_graph_from_file("graph.txt");
    graphFile.print_graph();

    auto [primWeight, primTree] = graphFile.prim_algorithm();
    auto [kruskalWeight, kruskalTree] = graphFile.kruskal_algorithm();

    cout << "Prim's Algorithm - Minimal Spanning Tree:" << endl;
    cout << "Weight: " << primWeight << endl;
    cout << "Edges: ";
    for (auto [u, v] : primTree)
    {
        cout << "(" << u << ", " << v << ") ";
    }
    cout << endl;

    cout << "Kruskal's Algorithm - Minimal Spanning Tree:" << endl;
    cout << "Weight: " << kruskalWeight << endl;
    cout << "Edges: ";
    for (auto [u, v] : kruskalTree)
    {
        cout << "(" << u << ", " << v << ") ";
    }
    cout << endl;

    return 0;
}