#ifndef GRAPHALGORITHMS_H
#define GRAPHALGORITHMS_H

#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <limits>
#include <cmath>
#include <chrono>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>

class Graph {
private:
    std::vector<std::vector<int>> adjacencyMatrix;
    std::vector<std::vector<std::pair<int, int>>> adjacencyList;
    bool useMatrix;

public:
    Graph(bool matrixForm = true) : useMatrix(matrixForm) {}

    void loadFromMatrix(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        adjacencyMatrix.clear();
        std::string line;
        while (std::getline(file, line)) {
            std::vector<int> row;
            std::istringstream iss(line);
            int value;
            while (iss >> value) {
                row.push_back(value);
            }
            adjacencyMatrix.push_back(row);
        }
        file.close();
        useMatrix = true;
    }

    void loadFromEdges(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // First pass: find maximum vertex
        int maxVertex = 0;
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int u, v, w;
            if (iss >> u >> v >> w) {
                maxVertex = std::max({ maxVertex, u, v });
            }
        }

        // Initialize adjacency list
        adjacencyList.resize(maxVertex + 1);

        // Second pass: build graph
        file.clear();
        file.seekg(0);
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int u, v, w;
            if (iss >> u >> v >> w) {
                adjacencyList[u].push_back({ v, w });
                adjacencyList[v].push_back({ u, w }); // Undirected
            }
        }
        file.close();
        useMatrix = false;
    }

    int getVerticesCount() const {
        return useMatrix ? adjacencyMatrix.size() : adjacencyList.size();
    }

    std::vector<int> getNeighbors(int vertex) const {
        std::vector<int> neighbors;
        if (useMatrix) {
            for (int i = 0; i < adjacencyMatrix[vertex].size(); i++) {
                if (adjacencyMatrix[vertex][i] > 0) {
                    neighbors.push_back(i);
                }
            }
        }
        else {
            for (const auto& neighbor : adjacencyList[vertex]) {
                neighbors.push_back(neighbor.first);
            }
        }
        return neighbors;
    }

    int getWeight(int u, int v) const {
        if (useMatrix) {
            return adjacencyMatrix[u][v];
        }
        else {
            for (const auto& neighbor : adjacencyList[u]) {
                if (neighbor.first == v) {
                    return neighbor.second;
                }
            }
            return 0;
        }
    }

    const std::vector<std::vector<int>>& getMatrix() const { return adjacencyMatrix; }
    const std::vector<std::vector<std::pair<int, int>>>& getList() const { return adjacencyList; }
    bool isMatrixForm() const { return useMatrix; }
};

class FordFulkerson {
private:
    std::vector<std::vector<int>> capacity;
    std::vector<std::vector<int>> flow;
    int source, sink;

    bool bfs(std::vector<int>& parent) {
        std::vector<bool> visited(capacity.size(), false);
        std::queue<int> q;
        q.push(source);
        visited[source] = true;
        parent[source] = -1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v = 0; v < capacity.size(); v++) {
                if (!visited[v] && capacity[u][v] - flow[u][v] > 0) {
                    parent[v] = u;
                    visited[v] = true;
                    q.push(v);
                    if (v == sink) return true;
                }
            }
        }
        return false;
    }

public:
    FordFulkerson(const std::vector<std::vector<int>>& cap, int s, int t)
        : capacity(cap), source(s), sink(t) {
        flow.resize(capacity.size(), std::vector<int>(capacity.size(), 0));
    }

    int computeMaxFlow() {
        int maxFlow = 0;
        std::vector<int> parent(capacity.size());

        while (bfs(parent)) {
            int pathFlow = std::numeric_limits<int>::max();

            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = std::min(pathFlow, capacity[u][v] - flow[u][v]);
            }

            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                flow[u][v] += pathFlow;
                flow[v][u] -= pathFlow;
            }

            maxFlow += pathFlow;
        }
        return maxFlow;
    }

    const std::vector<std::vector<int>>& getFlow() const { return flow; }
    const std::vector<std::vector<int>>& getCapacity() const { return capacity; }
    int getSource() const { return source; }
    int getSink() const { return sink; }
};

class PathFinding {
private:
    std::vector<std::vector<int>> grid;
    int rows, cols;

    struct Node {
        int x, y;
        double f, g, h;
        std::shared_ptr<Node> parent;

        Node(int x, int y) : x(x), y(y), f(0), g(0), h(0) {}

        bool operator>(const Node& other) const {
            return f > other.f;
        }
    };

    double chebyshevHeuristic(int x1, int y1, int x2, int y2) {
        return std::max(std::abs(x1 - x2), std::abs(y1 - y2));
    }

    double euclideanHeuristic(int x1, int y1, int x2, int y2) {
        return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
    }

    bool isValid(int x, int y) {
        return x >= 0 && x < rows && y >= 0 && y < cols && grid[x][y] > 0;
    }

public:
    PathFinding(const std::vector<std::vector<int>>& g) : grid(g) {
        rows = grid.size();
        cols = grid[0].size();
    }

    struct PathResult {
        double length;
        std::vector<std::pair<int, int>> path;
        double exploredPercent;
        double timeMs;
    };

    PathResult dijkstra(int startX, int startY, int endX, int endY) {
        auto startTime = std::chrono::high_resolution_clock::now();

        std::vector<std::vector<double>> dist(rows, std::vector<double>(cols, std::numeric_limits<double>::max()));
        std::vector<std::vector<bool>> visited(rows, std::vector<bool>(cols, false));
        std::vector<std::vector<std::pair<int, int>>> parent(rows, std::vector<std::pair<int, int>>(cols, { -1, -1 }));

        auto cmp = [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return dist[a.first][a.second] > dist[b.first][b.second];
            };
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, decltype(cmp)> pq(cmp);

        dist[startX][startY] = 0;
        pq.push({ startX, startY });

        int explored = 0;
        int totalValid = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (grid[i][j] > 0) totalValid++;
            }
        }

        while (!pq.empty()) {
            auto [x, y] = pq.top();
            pq.pop();

            if (visited[x][y]) continue;
            visited[x][y] = true;
            explored++;

            if (x == endX && y == endY) break;

            // 4-directional movement
            int dx[] = { -1, 1, 0, 0 };
            int dy[] = { 0, 0, -1, 1 };

            for (int i = 0; i < 4; i++) {
                int nx = x + dx[i];
                int ny = y + dy[i];

                if (isValid(nx, ny)) {
                    double newDist = dist[x][y] + std::abs(grid[nx][ny] - grid[x][y]) + 1;
                    if (newDist < dist[nx][ny]) {
                        dist[nx][ny] = newDist;
                        parent[nx][ny] = { x, y };
                        pq.push({ nx, ny });
                    }
                }
            }
        }

        // Reconstruct path
        std::vector<std::pair<int, int>> path;
        if (dist[endX][endY] != std::numeric_limits<double>::max()) {
            int x = endX, y = endY;
            while (x != -1 && y != -1) {
                path.push_back({ x, y });
                auto [px, py] = parent[x][y];
                x = px; y = py;
            }
            std::reverse(path.begin(), path.end());
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        double timeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();

        return { dist[endX][endY], path, (explored * 100.0) / totalValid, timeMs };
    }

    PathResult aStar(int startX, int startY, int endX, int endY) {
        auto startTime = std::chrono::high_resolution_clock::now();

        auto cmp = [](const std::shared_ptr<Node>& a, const std::shared_ptr<Node>& b) {
            return a->f > b->f;
            };
        std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, decltype(cmp)> openSet(cmp);
        std::vector<std::vector<bool>> closedSet(rows, std::vector<bool>(cols, false));

        std::vector<std::vector<std::shared_ptr<Node>>> nodes(rows, std::vector<std::shared_ptr<Node>>(cols));
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                nodes[i][j] = std::make_shared<Node>(i, j);
            }
        }

        nodes[startX][startY]->g = 0;
        nodes[startX][startY]->h = chebyshevHeuristic(startX, startY, endX, endY);
        nodes[startX][startY]->f = nodes[startX][startY]->g + nodes[startX][startY]->h;
        openSet.push(nodes[startX][startY]);

        int explored = 0;
        int totalValid = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (grid[i][j] > 0) totalValid++;
            }
        }

        while (!openSet.empty()) {
            auto current = openSet.top();
            openSet.pop();

            int x = current->x;
            int y = current->y;

            if (closedSet[x][y]) continue;
            closedSet[x][y] = true;
            explored++;

            if (x == endX && y == endY) break;

            // 4-directional movement
            int dx[] = { -1, 1, 0, 0 };
            int dy[] = { 0, 0, -1, 1 };

            for (int i = 0; i < 4; i++) {
                int nx = x + dx[i];
                int ny = y + dy[i];

                if (isValid(nx, ny) && !closedSet[nx][ny]) {
                    double tentative_g = current->g + std::abs(grid[nx][ny] - grid[x][y]) + 1;

                    if (tentative_g < nodes[nx][ny]->g) {
                        nodes[nx][ny]->parent = current;
                        nodes[nx][ny]->g = tentative_g;
                        nodes[nx][ny]->h = chebyshevHeuristic(nx, ny, endX, endY);
                        nodes[nx][ny]->f = nodes[nx][ny]->g + nodes[nx][ny]->h;
                        openSet.push(nodes[nx][ny]);
                    }
                }
            }
        }

        // Reconstruct path
        std::vector<std::pair<int, int>> path;
        if (nodes[endX][endY]->parent != nullptr) {
            auto current = nodes[endX][endY];
            while (current != nullptr) {
                path.push_back({ current->x, current->y });
                current = current->parent;
            }
            std::reverse(path.begin(), path.end());
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        double timeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
        double pathLength = nodes[endX][endY]->g;

        return { pathLength, path, (explored * 100.0) / totalValid, timeMs };
    }
};

class BipartiteMatching {
private:
    std::vector<std::vector<int>> graph;
    std::vector<int> match;
    std::vector<bool> used;
    std::vector<int> color;

    bool dfs(int u) {
        if (used[u]) return false;
        used[u] = true;

        for (int v : graph[u]) {
            if (match[v] == -1 || dfs(match[v])) {
                match[v] = u;
                match[u] = v;
                return true;
            }
        }
        return false;
    }

    bool isBipartiteUtil(int u, int c) {
        color[u] = c;
        for (int v : graph[u]) {
            if (color[v] == -1) {
                if (!isBipartiteUtil(v, 1 - c)) {
                    return false;
                }
            }
            else if (color[v] == c) {
                return false;
            }
        }
        return true;
    }

public:
    BipartiteMatching(const std::vector<std::vector<int>>& g) : graph(g) {
        int n = g.size();
        color.assign(n, -1);
        match.assign(n, -1);
    }

    bool isBipartite() {
        for (int i = 0; i < graph.size(); i++) {
            if (color[i] == -1) {
                if (!isBipartiteUtil(i, 0)) {
                    return false;
                }
            }
        }
        return true;
    }

    int findMaxMatching() {
        if (!isBipartite()) return 0;

        int n = graph.size();
        match.assign(n, -1);

        // Kuhn's algorithm
        for (int u = 0; u < n; u++) {
            if (color[u] == 0) { // Left side
                used.assign(n, false);
                dfs(u);
            }
        }

        int matching = 0;
        for (int u = 0; u < n; u++) {
            if (color[u] == 1 && match[u] != -1) { // Right side with match
                matching++;
            }
        }
        return matching;
    }

    const std::vector<int>& getMatching() const { return match; }
    const std::vector<int>& getColor() const { return color; }
};

#endif