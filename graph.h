/***************************************************************************************\
* Class: directed and undirected graph using adjacency list representation 				*
* based on https://www.geeksforgeeks.org/graph-and-its-representations/					*
\***************************************************************************************/

#include<iostream> 
#include <list> 
using namespace std; 

class Graph 
{ 
	private:
		int V; 									// Number of vertices 	
		list<int> *adj; 						// Points to adjacency list
		void DFSUtil(int v, bool visited[], int *connected_v, int comp_index); 
	public: 
		Graph(int V); 
		~Graph(); 
		void addUndirectedEdge(int v, int w); 
		void addDirectedEdge(int v, int w); 
		void connectedComponents(int *connected_v); 
		int isThereEdge(int v, int w);
		void printGraph(void); 
}; 

/******************************************************************************\
*								Constructor									   *
\******************************************************************************/
Graph::Graph(int V) 
{ 
	this->V = V; 
	adj = new list<int>[V]; 
} 

/******************************************************************************\
*								 Destructor									   *
\******************************************************************************/
Graph::~Graph() 
{ 
	delete [] adj; 
} 

/******************************************************************************\
*					Add an undirected edge									   *
\******************************************************************************/
void Graph::addUndirectedEdge(int v, int w) 
{ 
	adj[v].push_back(w); 
	adj[w].push_back(v); 
} 

/******************************************************************************\
*					Add directed edge									   *
\******************************************************************************/
void Graph::addDirectedEdge(int v, int w) 
{ 
	adj[v].push_back(w); 
} 

/******************************************************************************\
*	Finding the connected components in an undirected graph					   *
\******************************************************************************/
void Graph::connectedComponents(int *connected_v) 
{ 
	int v, comp_index=0;		// index of the component
	
	// Mark all the vertices as not visited 
	bool *visited = new bool[V]; 
	for(v = 0; v < V; v++) 
		visited[v] = false; 

	for (v=0; v<V; v++) 
	{ 
		if (visited[v] == false) 
		{ 
			// find all reachable vertices from v 
			DFSUtil(v, visited, connected_v, comp_index); 			
			comp_index++;
			//cout <<endl; 
		} 
	} 
	delete [] visited; 
} 

/******************************************************************************\
*					Visiting vertices in a connected component				   *
\******************************************************************************/
void Graph::DFSUtil(int v, bool visited[], int *connected_v, int comp_index) 
{ 
	list<int>::iterator i; 

	// Mark the current node as visited 
	visited[v] = true; 
	//cout << v << " "; 
	connected_v[v]=comp_index;

	// Recurrence for all the vertices adjacent to vertex v
	for(i = adj[v].begin(); i != adj[v].end(); i++) 
		if(!visited[*i]) 
			DFSUtil(*i, visited, connected_v, comp_index); 
} 

/******************************************************************************\
*					Check if there is a given edge in the Graph				   *
\******************************************************************************/
int Graph::isThereEdge(int v, int w) 
{ 
	list<int>::iterator i;
	
	i = adj[v].begin();
	while ( i != adj[v].end() && w != *i ) 
		i++;

	if (i != adj[v].end()) 
		return 1;
	else
		return 0;
} 

/******************************************************************************\
*					Print a graph											   *
\******************************************************************************/
void Graph::printGraph(void) 
{ 
	int x, v;
	list<int>::iterator i;
	
	cout <<"Graph:"<<endl;  
    for (v = 0; v < V; ++v) 
    { 
        cout << v << "->";    	 
		for(i = adj[v].begin(); i != adj[v].end(); ++i) {
			x=*i;
		    cout << x << " ";
		}
        cout<<endl; 
    } 
} 

