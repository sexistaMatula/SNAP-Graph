#include <deque>
#include <chrono>
#include <functional>
#include <algorithm>
#include <queue>
#include <vector>
#include <iostream>
#include <utility>
#include <limits>
#include <iomanip>
#include "stdafx.h"
#undef min
#undef max
#define inf std::numeric_limits<float>::infinity()

typedef TNodeEDatNet<TInt,TFlt>::TNodeI SnapNode;
typedef TNodeEDatNet<TInt,TFlt>::TEdgeI SnapEdge;
typedef TPt<TNodeEDatNet<TInt, TFlt> > Graph;
typedef std::chrono::microseconds micros;
bool vertexExists(TPt<TNodeEDatNet<TInt, TFlt> >  & G, int id)
{
    for (SnapNode NI = G->BegNI(); NI < G->EndNI(); NI++)
    {
        if(NI.GetDat() == id)
        {
            return true;
        }
    }
    return false;
}


bool edgeExists(TPt<TNodeEDatNet<TInt, TFlt> >  & G, int id1, int id2)
{
    for (SnapEdge EI = G->BegEI(); EI < G->EndEI(); EI++)
    {
        if(EI.GetDstNDat() == id2 && EI.GetSrcNDat() == id1)
            return true;
    }
    return false;
}

bool edgeExists(TPt<TNodeEDatNet<TInt, TFlt> >  & G, int id1, int id2, TFlt weight)
{
    for (SnapEdge EI = G->BegEI(); EI < G->EndEI(); EI++)
    {
        if(EI.GetDstNDat() == id2 && EI.GetSrcNDat() == id1 && EI.GetDat() == weight)
            return true;
    }
    return false;
}

class myVis {
public:
    myVis() { }
    myVis(const int& Nodes){ }
    void DiscoverNode(int NId) { std::cout << NId << std::endl; }
    void FinishNode(const int& NId) { }
    void ExamineEdge(const int& NId1, const int& NId2) { }
    void TreeEdge(const int& NId1, const int& NId2) { }
    void BackEdge(const int& NId1, const int& NId2) { }
    void FwdEdge(const int& NId1, const int& NId2) { }
};

class myVisBench {
public:
    myVisBench() { }
    myVisBench(const int& Nodes){ }
    void DiscoverNode(int NId) {}
    void FinishNode(const int& NId) { }
    void ExamineEdge(const int& NId1, const int& NId2) { }
    void TreeEdge(const int& NId1, const int& NId2) { }
    void BackEdge(const int& NId1, const int& NId2) { }
    void FwdEdge(const int& NId1, const int& NId2) { }
};

void printTable(std::vector<std::vector< float > > graph)
{
    for(int i = 0; i<graph.size(); i++)
    {
        for(int j = 0; j<graph.size(); j++)
        {
            std::cout << "|"<< j << "|";
        }
    }
    std::cout << std::endl;
}

void floydWarshall( std::vector< std::vector<float > > & graph)
{
    for(int k = 0; k<graph.size(); k++)
    {
        for(int i = 0; i<graph.size(); i++)
        {
            for(int j = 0; j<graph.size(); j++)
            {
                if(graph[i][k] + graph[k][j] < graph[i][j])
                {
                    graph[i][j] = graph[i][k] + graph[k][j];
                }
            }
        }
    }
}

void printFloydWarshall(std::vector< std::vector<float > > graph)
{
    for(int i = 0; i<graph.size(); i++)
    {
        for(int j = 0; j<graph.size(); j++)
        {
            std::cout << "From vertex " << i+1 << " to " << j+1 << " : ";
            if(graph[i][j] == inf)
                std::cout << "inf" << std::endl;
            else
                std::cout << graph[i][j] << std::endl;
        }
        std::cout << std::endl;
    }
}

void printDetail(Graph & G)
{
    for (SnapNode NI = G->BegNI(); NI != G->EndNI(); NI++)
    {
        std::cout << "Id: " << NI.GetId() << " Out Degree: " << NI.GetOutDeg() << " In Degree: " << NI.GetInDeg() << std::endl;
    }

    for (SnapEdge NI = G->BegEI(); NI < G->EndEI(); NI++)
    {
        std::cout << "Edge from " << (int)NI.GetSrcNDat() << " --> " << (int)NI.GetDstNDat() << " with weight " << (int)NI.GetDstNDat() << std::endl;
    }
}

struct QueueItem{
    int VertexID;
    float distance;
};


bool sortQueue(QueueItem a, QueueItem b)
{
    return a.distance < b.distance;
}

bool sortQueueInt(int a, int b, std::vector<int>distances)
{
    return distances[a-1] < distances[b-1];
}

bool sortClusters(SnapEdge a, SnapEdge b)
{
    return (float) a.GetDat() < (float) b.GetDat();
}

class Compare
{
public:
    bool operator() (SnapEdge a, SnapEdge b)
    {
        return (float) a.GetDat() > (float) b.GetDat();
    }
};



class disjointTreeNode{
public:
    int id;
    disjointTreeNode* parent;
    int rank;
    disjointTreeNode():id(0),parent(nullptr),rank(0){}
    disjointTreeNode(int val){
        id = val;
        parent = nullptr;
        rank = 0;
    }
};


void makeSet(disjointTreeNode* x)
{
    x->parent = x;
    x->rank = 0;
}

disjointTreeNode* Find(disjointTreeNode* x)
{
    if(x->parent != x)
        x->parent = Find(x->parent);
    return x->parent;
}

bool sortFunction(int a, int b, std::vector<float> distances)
{
    return distances[a-1] < distances[b-1];
}

class sorter {
    std::vector<float> type_;
public:
    sorter(std::vector<float> type) : type_(type) {}
    bool operator()(int o1, int o2) const {
        return sortFunction(o1, o2, type_ );
    }
};

void Union(disjointTreeNode* x, disjointTreeNode* y)
{
    disjointTreeNode* xRoot = Find(x);
    disjointTreeNode* yRoot = Find(y);
    if(xRoot == yRoot)
        return;
    if(xRoot->rank < yRoot->rank)
        xRoot->parent = yRoot;
    else if(xRoot->rank > yRoot->rank)
        yRoot->parent = xRoot;
    else
    {
        yRoot->parent = xRoot;
        xRoot->rank = xRoot->rank+1;
    }
}

int getLargestNode(Graph & G)
{
    int largest = G->BegNI().GetDat();
    for (SnapNode NI = G->BegNI(); NI != G->EndNI(); NI++)
    {
        if(NI.GetDat() > largest)
            largest = NI.GetDat();
    }
    return largest;
}

void kruskal(Graph & G, bool bench)
{
    int largestNode = getLargestNode(G);
    std::vector<disjointTreeNode*> disjointNodes(largestNode);
    std::vector<SnapEdge> result;
    std::deque<SnapEdge> Q;

    for(SnapNode NI = G->BegNI(); NI < G->EndNI(); NI++)
    {
        disjointTreeNode* x =  new disjointTreeNode();
        x->id = NI.GetDat();
        makeSet(x);
        disjointNodes[x->id-1] = x;
    }


    for (SnapEdge NI = G->BegEI(); NI < G->EndEI(); NI++)
    {
        Q.push_back(NI);
    }

    while (! Q.empty())
    {
        std::sort(Q.begin(), Q.end(), sortClusters);
        SnapEdge e = Q.front();
        Q.pop_front();
        disjointTreeNode* u = Find(disjointNodes[e.GetSrcNDat()-1]);
        disjointTreeNode* v = Find(disjointNodes[e.GetDstNDat()()-1]);
        if ( u != v )
        {
            result.push_back(e);
            Union(disjointNodes[e.GetSrcNDat()-1], disjointNodes[e.GetDstNDat()-1]);
        }

    }
    if(!bench)
    {
        for(int i = 0; i<result.size(); i++)
        {
            std::cout << result[i].GetSrcNDat() << " --> " << result[i].GetDstNDat() << " peso: " << result[i].GetDat() << std::endl;
        }
    }
}

struct Cluster{
    bool empty;
    std::vector<int> clusterElements;
};

bool compareCluster(int vert1, int vert2, std::vector<Cluster> & clusters)
{
    int v1pos = 0;
    int v2pos = 0;
    for(int i=0; i<clusters.size(); i++)
    {
        for(int j = 0; j < clusters[i].clusterElements.size() ; j++)
        {
            if(clusters[i].clusterElements[j] == vert1)
            {
                v1pos = i;
            }
            if(clusters[i].clusterElements[j] == vert2)
            {
                v2pos = i;
            }
        }
    }
    if(v1pos == v2pos)
    {
        return true;
    }
    else
    {
        Cluster newCluster;
        clusters[v1pos].empty = true;
        clusters[v2pos].empty = true;
        for(int i = 0; i < clusters[v1pos].clusterElements.size(); i++)
            newCluster.clusterElements.push_back(clusters[v1pos].clusterElements[i]);
        for(int i = 0; i < clusters[v2pos].clusterElements.size(); i++)
            newCluster.clusterElements.push_back(clusters[v2pos].clusterElements[i]);
        for(int i = 0; i<clusters.size(); i++)
        {
            if(clusters[i].empty)
            {
                clusters.erase(clusters.begin()+i);
                i--;
            }
        }
        newCluster.empty = false;
        clusters.push_back(newCluster);
        return false;
    }
}

void kruskalCluster(Graph & G)
{
    std::vector<Cluster> clusters;
    for (SnapNode NI = G->BegNI(); NI < G->EndNI(); NI++)
    {
        Cluster c;
        c.empty = false;
        c.clusterElements.push_back(NI.GetDat());
        clusters.push_back(c);
    }

    std::vector <SnapEdge> result;
    std::deque <SnapEdge> Q;

    for (SnapEdge NI = G->BegEI(); NI < G->EndEI(); NI++)
    {
        Q.push_back(NI);
    }

    while(!Q.empty())
    {
        std::sort(Q.begin(), Q.end(), sortClusters);
        SnapEdge edge = Q.front();
        Q.pop_front();
        if(!compareCluster((int)edge.GetSrcNDat(), (int)edge.GetDstNDat(), clusters))
        {
            result.push_back(edge);
        }
    }

    for(int i = 0; i<result.size(); i++)
    {
        std::cout << "(" << result[i].GetSrcNDat() << ", " << result[i].GetDstNDat() << ", " << result[i].GetDat() << ")" << std::endl;
    }
}

void floydWarshallPaths(Graph & G, int size, bool bench)
{
    std::vector< std::vector<float > > VectorGraph(size, std::vector<float>(size, inf));
    for(int i=0; i<VectorGraph.size(); i++)
    {
        for(int j = 0; j<VectorGraph.size(); j++)
        {
            if(i == j)
                VectorGraph[i][j] = 0;
        }
    }
    for (TNodeEDatNet<TInt, TFlt>::TEdgeI NI = G->BegEI(); NI < G->EndEI(); NI++)
    {
        int vert1 = NI.GetSrcNDat();
        int vert2 = NI.GetDstNDat();
        VectorGraph[vert1-1][vert2-1] = (float) NI.GetDat();
    }
    floydWarshall(VectorGraph);
    if(!bench)
        printFloydWarshall(VectorGraph);
}

void dijkstra(Graph & G, int v, bool bench)
{
    int size = G->GetNodes();
    std::vector<float> distances (size,inf);
    std::vector<int> parents (size,-1);
    std::priority_queue<std::pair<float, Graph::TObj::TNodeI> > queue;

    distances[v-1] = 0;
	parents[v-1] = -1;

    for(SnapNode NI = G->BegNI(); NI != G->EndNI(); NI++)
    {
		queue.push(std::make_pair(distances[NI.GetDat()-1], NI));
    }

    while(!queue.empty())
    {
        SnapNode u = queue.top().second;
        float dist = queue.top().first;
        queue.pop();

        int sourceNode = u.GetDat();
        for (int ed = 0; ed < u.GetOutDeg(); ed++)//For each adjacent vertex to node NI
        {
            int destNode = u.GetOutNDat(ed);
            SnapEdge e =  G->GetEI(sourceNode,destNode);
			float alt = dist + e.GetDat();
            if(alt < distances[destNode-1])//relax
            {
                distances[destNode-1] = alt;
                parents[destNode-1] = sourceNode;
                queue.push(std::make_pair(alt, G->GetNI(destNode)));
            }
        }
    }
    if(!bench)
    {
        int node = 1;
        for (int i = 0 ; i< parents.size(); i++)
        {
            if(parents[node-1] == -1)
                std::cout << "start -> " << node << " distance:  " << distances[node-1] << std::endl;
            else
                std::cout << parents[node-1] << " -> " << node << " distance:  " << distances[node-1] << std::endl;
            node++;
        }
    }
}

std::vector<int> primPrincipal(Graph & G, float & length, int v)
{
    int size = G->GetNodes();
    std::vector<float> distances (size,inf);
    std::vector<bool> visitedNodes (size,false);
    std::vector<int> parents (size,-1);
    std::deque<QueueItem> queue;
    for(int i=0; i<size; i++)
    {
        QueueItem q;
        q.VertexID = i+1;
        q.distance = inf;
        queue.push_back(q);
    }
    distances[v-1] = 0;
    visitedNodes[v-1] = true;
    while(!queue.empty())
    {
        std::sort(queue.begin(), queue.end(), sortQueue);
        QueueItem u = queue.front();
        queue.pop_front();
        SnapEdge EI;
        if(parents[u.VertexID-1] != -1)
        {
            TNodeEDatNet<TInt, TFlt>::TEdgeI EI = G->GetEI(parents[u.VertexID-1],u.VertexID);
            length+=(float)EI.GetDat();
        }
        visitedNodes[u.VertexID-1] = true;
        SnapNode NI = G->GetNI(u.VertexID);
        for (int e = 0; e < NI.GetOutDeg(); e++)
        {
            SnapEdge EI = G->GetEI(u.VertexID,NI.GetOutNDat(e));
            int outNode = NI.GetOutNDat(e);
            float edgeWeight = (float)EI.GetDat();
            if(!visitedNodes[outNode-1] && distances[outNode-1] > edgeWeight)
            {
                parents[outNode-1] = u.VertexID;
                distances[outNode-1] = edgeWeight;
                for(int i=0; i< queue.size(); i++)
                {
                    if(queue[i].VertexID == outNode)
                        queue[i].distance = edgeWeight;
                }
            }
        }
    }
    return parents;
}

void prim(Graph & G, int bench)
{
    if(!bench)
    {
        int v;
        float length = 0;
        std::cout << "Vértice de origen: ";
        while(!(std::cin >> v)){
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Entrada inválida, intenta otra vez: ";
        }
        if(vertexExists(G, v))
        {
            std::vector<int> parents = primPrincipal(G,length,v);
            int node = 1;
            for (int i = 0 ; i< parents.size(); i++)
            {
                if(parents[node-1] == -1)
                    std::cout << "null -> " << node << std::endl;
                else
                    std::cout << parents[node-1] << " -> " << node << std::endl;
                node++;
            }
            std::cout << "Length: " << length << std::endl;
        }
        else
        {
            std::cout << "El vértice " << v << " no existe" << std::endl;
        }
    }
    else
    {
        float t;
        primPrincipal(G,t,1);
    }
}

void printAlgorithm(std::string alg, micros dur)
{
   std::cout << alg << ": " << dur.count()<< " micros" << std::endl;
   std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    TPt<TNodeEDatNet<TInt, TFlt> >  G = TNodeEDatNet<TInt, TFlt>::New();
    bool negativeEdge = false;
   for(int i=1; i<15; i++)
   {
       G->AddNode(i,i);
   }

   G->AddEdge(1,4,8);
   G->AddEdge(1,3,8);
   G->AddEdge(2,5,7);
   G->AddEdge(3,10,4);
   G->AddEdge(3,2,7);
   G->AddEdge(3,5,8);
   G->AddEdge(4,7,3);
   G->AddEdge(4,5,1);
   G->AddEdge(4,8,2);
   G->AddEdge(5,6,9);
   G->AddEdge(6,13,4);
   G->AddEdge(7,4,6);
   G->AddEdge(8,9,3);
   G->AddEdge(8,7,3);
   G->AddEdge(9,10,2);
   G->AddEdge(9,12,4);
   G->AddEdge(10,3,10);
   G->AddEdge(10,6,6);
   G->AddEdge(11,12,6);
   G->AddEdge(12,11,8);
   G->AddEdge(12,9,2);
   G->AddEdge(12,14,9);
   G->AddEdge(13,14,6);
   G->AddEdge(14,13,2);

    int choice;
    bool done = false;
    while(!done)
    {
        std::cout << "===== SNAP Library Implementation ===== " << std::endl;
        std::cout << "1) Insertar vértice en el grafo" << std::endl;
        std::cout << "2) Insertar arista en el grafo" << std::endl;
        std::cout << "3) Eliminar vértice del grafo" << std::endl;
        std::cout << "4) Eliminar arista del grafo" << std::endl;
        std::cout << "5) Depth First Traversal" << std::endl;
        std::cout << "6) Breadth First Traversal" << std::endl;
        std::cout << "7) Algoritmo de Prim" << std::endl;
        std::cout << "8) Algoritmo de Kruskal" << std::endl;
        std::cout << "9) Ruta mínima con Dijkstra" << std::endl;
        std::cout << "10) Ruta mínima con Floyd-Warshall" << std::endl;
        std::cout << "11) Imprimir grafo" << std::endl;
        std::cout << "12) Tiempos" << std::endl;
        std::cout << "13) Salir" << std::endl;
        std::cout << "Elige: ";

        while(!(std::cin >> choice)){
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Entrada inválida, intenta otra vez: ";
        }

        switch(choice)
        {
            case 1:{
                int v;
                std::cout << "Id del nodo a agregar: ";
                while(!(std::cin >> v)){
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "Entrada inválida, intenta otra vez: ";
                }

                if(vertexExists(G, v))
                {
                    std::cout << "El vértice con id: " << v << " ya existe" << std::endl;
                }
                else
                {
                    G->AddNode(v, v);
                    std::cout << "Vértice con id: " << v << " insertado correctamente" << std::endl;
                }
                break;}
            case 2:{

                int id1, id2;
                std::cout << "Vértice 1: ";
                while(!(std::cin >> id1)){
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "Entrada inválida, intenta otra vez: ";
                }

                if(!vertexExists(G, id1))
                    std::cout << "El vértice con id: " << id1 << " no existe" << std::endl;
                else
                {
                    std::cout << "Vértice 2: ";
                    while(!(std::cin >> id2)){
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::cout << "Entrada inválida, intenta otra vez: ";
                    }

                    if(!vertexExists(G, id2))
                        std::cout << "El vértice con id: " << id2 << " no existe" << std::endl;
                    else
                    {
                        float weight;
                        std::cout << "Peso de la arista: ";
                        while(!(std::cin >> weight)){
                            std::cin.clear();
                            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                            std::cout << "Entrada inválida, intenta otra vez: ";
                        }
                        if(edgeExists(G,id1,id2,weight))
                        {
                            std::cout << "La arista de " << id1 << " a " << id2 << " con peso " << weight << " ya existe" << std::endl;
                        }
                        else
                        {
                            if(weight < 0)
                            {
                                negativeEdge = true;
                                std::cout << "No se pueden agregar aristas negativas" << std::endl;
                                break;
                            }
                            G->AddEdge(id1, id2, weight);
                            std::cout << "Se añadió la arista de " << id1 << " a " << id2 << " con peso " << weight << std::endl;
                        }
                    }
                }
                break;}
            case 3:{
                int v;
                std::cout << "Id del nodo a borrar: ";
                while(!(std::cin >> v)){
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "Entrada inválida, intenta otra vez: ";
                }

                SnapNode Vertex;
                if (vertexExists(G,v))
                    G->DelNode(v);
                else
                    std::cout << "El vértice con id " << v << " no existe" << std::endl;
                break;}
            case 4:{
                int id1, id2;
                std::cout << "Id del vértice 1 conectado por la arista: ";
                while(!(std::cin >> id1)){
                    std::cin.clear();
                    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                    std::cout << "Entrada inválida, intenta otra vez: ";
                }

                if(vertexExists(G,id1))
                {
                    std::cout << "Id del vértice 2 conectado por la arista: ";
                    while(!(std::cin >> id2)){
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::cout << "Entrada inválida, intenta otra vez: ";
                    }

                    if(vertexExists(G,id2))
                    {
                        if(edgeExists(G,id1,id2))
                        {
                            G->DelEdge(id1, id2);
                            std::cout << "La arista de " << id1 << " a " << id2 << " fue borrada exitosamente" << std::endl;
                        }
                        else
                        {
                            std::cout << "La arista de " << id1 << " a " << id2 << " no existe" << std::endl;
                        }
                    }
                    else
                    {
                        std::cout << "El vértice con id " << id2 << " no existe" << std::endl;
                    }
                }
                else
                {
                    std::cout << "El vértice con id " << id1 << " no existe" << std::endl;
                }
                break;}
            case 5:{
                if(G->GetNodes() == 0)
                {
                    std::cout << "No hay nodos en el grafo" << std::endl;
                }
                else
                {
                    myVis vis(G->GetNodes());
                    TCnCom::GetDfsVisitor(G, vis);
                }
                break;}
            case 6:{
                if(G->GetNodes() == 0)
                {
                    std::cout << "No hay nodos en el grafo" << std::endl;
                }
                else
                {
                    PNGraph GBFS = TSnap::GetBfsTree(G, 1, true, true);
                    for (TNGraph::TNodeI NI = GBFS->BegNI(); NI < GBFS->EndNI(); NI++)
                    {
                        std::cout << NI.GetId() << std::endl;
                    }
                }
                break;}
            case 7:{
                prim(G, false);
                break;}
            case 8:{
                if(G->GetNodes() != 0)
                    kruskal(G, false);
                else
                    std::cout << "El grafo está vacío" << std::endl;
                break;}
            case 9:{
                if(G->GetNodes() != 0)
                {
                    int v;
                    std::cout << "Vértice de incicio: ";
                    while(!(std::cin >> v)){
                        std::cin.clear();
                        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                        std::cout << "Entrada inválida, intenta otra vez: ";
                    }
                    if(!negativeEdge)
                    {
                        if(vertexExists(G, v))
                            dijkstra(G, v, false);
                        else
                            std::cout << "El vértice con id " << v << " no existe" << std::endl;
                    }
                    else
                        std::cout << "El grafo tiene aristas con peso negativo, no se puede ejecutar el algoritmo de Dijkstra correctamente" << std::endl;
                }
                else
                    std::cout << "El grafo está vacío" << std::endl;
                break;}
            case 10:{
                int size = G->GetNodes();
                if(size != 0)
                    floydWarshallPaths(G, size, false);
                else
                    std::cout << "El grafo está vacío" << std::endl;
                break;}
            case 11:{
                printDetail(G);
                break;}
            case 12:{
                TPt<TNodeEDatNet<TInt, TFlt> >  g = TNodeEDatNet<TInt, TFlt>::New();
                for(int i=1; i<15; i++)
                {
                    g->AddNode(i,i);
                }

                g->AddEdge(1,4,8);
                g->AddEdge(1,3,8);
                g->AddEdge(2,5,7);
                g->AddEdge(3,10,4);
                g->AddEdge(3,2,7);
                g->AddEdge(3,5,8);
                g->AddEdge(4,7,3);
                g->AddEdge(4,5,1);
                g->AddEdge(4,8,2);
                g->AddEdge(5,6,9);
                g->AddEdge(6,13,4);
                g->AddEdge(7,4,6);
                g->AddEdge(8,9,3);
                g->AddEdge(8,7,3);
                g->AddEdge(9,10,2);
                g->AddEdge(9,12,4);
                g->AddEdge(10,3,10);
                g->AddEdge(10,6,6);
                g->AddEdge(11,12,6);
                g->AddEdge(12,11,8);
                g->AddEdge(12,9,2);
                g->AddEdge(12,14,9);
                g->AddEdge(13,14,6);
                g->AddEdge(14,13,2);


                 //--------------------- DFS ---------------------//
                 std::cout << "==============DFS===============" << std::endl;
                 auto begin = std::chrono::high_resolution_clock::now();

                 myVisBench vis(g->GetNodes());
                 TCnCom::GetDfsVisitor(G, vis);

                 auto end = std::chrono::high_resolution_clock::now();
                 auto alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("DFS", alg);

                 //--------------------- BFS ---------------------//
                 std::cout << "==============BFS===============" << std::endl;
                 begin = std::chrono::high_resolution_clock::now();

                 TSnap::GetBfsTree(g, 1, true, true);

                 end = std::chrono::high_resolution_clock::now();
                 alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("BFS", alg);

                 //--------------------- Prim ---------------------//
                 std::cout << "==============PRIM===============" << std::endl;
                 begin = std::chrono::high_resolution_clock::now();

                 prim(g, true);

                 end = std::chrono::high_resolution_clock::now();
                 alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("Prim", alg);


                 //--------------------- Kruskal ---------------------//
                 std::cout << "==============Kruskal===============" << std::endl;
                 begin = std::chrono::high_resolution_clock::now();

                 kruskal(g,true);

                 end = std::chrono::high_resolution_clock::now();
                 alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("Kruskal", alg);


                 //--------------------- Dijkstra ---------------------//

                 std::cout << "==============Dijkstra===============" << std::endl;
                 begin = std::chrono::high_resolution_clock::now();

                 dijkstra(g, 1, true);

                 end = std::chrono::high_resolution_clock::now();
                 alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("Dijkstra", alg);


                 //--------------------- Floyd Warshall ---------------------//

                 std::cout << "==============Floyd Warshall===============" << std::endl;
                 begin = std::chrono::high_resolution_clock::now();

                 int size = g->GetNodes();
                 floydWarshallPaths(g, size, true);

                 end = std::chrono::high_resolution_clock::now();
                 alg = std::chrono::duration_cast<micros>(end-begin);
                 printAlgorithm("Floyd Warshall", alg);


                break;}
            case 13:{
                printf("Adios!");
                done = true;
                break;}
            default:
                printf("Opción inválida, elige otra");
                break;
        }
    }
}
