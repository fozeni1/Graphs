#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <limits>
#include <cmath>
#include <chrono>
#include <memory>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>

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
                nodes[i][j]->g = std::numeric_limits<double>::max();
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

        for (int u = 0; u < n; u++) {
            if (color[u] == 0) {
                used.assign(n, false);
                dfs(u);
            }
        }

        int matching = 0;
        for (int u = 0; u < n; u++) {
            if (color[u] == 1 && match[u] != -1) {
                matching++;
            }
        }
        return matching;
    }

    const std::vector<int>& getMatching() const { return match; }
    const std::vector<int>& getColor() const { return color; }
};

std::vector<std::vector<int>> loadMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::vector<std::vector<int>> matrix;
    std::string line;
    while (std::getline(file, line)) {
        std::vector<int> row;
        std::istringstream iss(line);
        int value;
        while (iss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
    }
    file.close();
    return matrix;
}

std::vector<std::vector<int>> loadEdges(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    int maxVertex = 0;
    std::vector<std::tuple<int, int, int>> edges;
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int u, v, w;
        if (iss >> u >> v >> w) {
            edges.push_back({ u, v, w });
            maxVertex = std::max({ maxVertex, u, v });
        }
    }
    file.close();

    std::vector<std::vector<int>> graph(maxVertex + 1);
    for (const auto& edge : edges) {
        int u = std::get<0>(edge);
        int v = std::get<1>(edge);
        graph[u].push_back(v);
        graph[v].push_back(u);
    }

    return graph;
}

void printUsage() {
    std::cout << "Usage:\n";
    std::cout << "1. Max flow: program -f input.txt [-o output.txt]\n";
    std::cout << "2. Path finding: program -p input.txt -a startX,startY -b endX,endY [-o output.txt]\n";
    std::cout << "3. Bipartite matching: program -m input.txt [-o output.txt]\n";
}

void taskMaxFlow(const std::string& inputFile, const std::string& outputFile) {
    try {
        auto matrix = loadMatrix(inputFile);
        int n = matrix.size();

        FordFulkerson ff(matrix, 0, n - 1);
        int maxFlow = ff.computeMaxFlow();

        std::ofstream out(outputFile);
        out << "Source: 0\n";
        out << "Sink: " << n - 1 << "\n";
        out << "Max Flow: " << maxFlow << "\n\n";
        out << "Edge Flows:\n";

        auto flow = ff.getFlow();
        auto capacity = ff.getCapacity();

        for (int u = 0; u < n; u++) {
            for (int v = 0; v < n; v++) {
                if (capacity[u][v] > 0) {
                    out << u << " -> " << v << ": Flow = " << flow[u][v]
                        << ", Capacity = " << capacity[u][v] << "\n";
                }
            }
        }

        std::cout << "Max flow computed successfully. Results saved to " << outputFile << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void taskPathFinding(const std::string& inputFile, const std::string& startStr,
    const std::string& endStr, const std::string& outputFile) {
    try {
        size_t comma1 = startStr.find(',');
        size_t comma2 = endStr.find(',');

        if (comma1 == std::string::npos || comma2 == std::string::npos) {
            throw std::invalid_argument("Invalid coordinate format");
        }

        int startX = std::stoi(startStr.substr(0, comma1));
        int startY = std::stoi(startStr.substr(comma1 + 1));
        int endX = std::stoi(endStr.substr(0, comma2));
        int endY = std::stoi(endStr.substr(comma2 + 1));

        auto grid = loadMatrix(inputFile);
        PathFinding pf(grid);

        auto dijkstraResult = pf.dijkstra(startX, startY, endX, endY);
        auto aStarResult = pf.aStar(startX, startY, endX, endY);

        std::ofstream out(outputFile);

        out << "DIJKSTRA ALGORITHM:\n";
        out << "Path Length: " << dijkstraResult.length << "\n";
        out << "Path: ";
        for (const auto& point : dijkstraResult.path) {
            out << "(" << point.first << "," << point.second << ") ";
        }
        out << "\nExplored: " << dijkstraResult.exploredPercent << "%\n";
        out << "Time: " << dijkstraResult.timeMs << "ms\n\n";

        out << "A* ALGORITHM:\n";
        out << "Path Length: " << aStarResult.length << "\n";
        out << "Path: ";
        for (const auto& point : aStarResult.path) {
            out << "(" << point.first << "," << point.second << ") ";
        }
        out << "\nExplored: " << aStarResult.exploredPercent << "%\n";
        out << "Time: " << aStarResult.timeMs << "ms\n";

        std::cout << "Path finding completed. Results saved to " << outputFile << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void taskBipartiteMatching(const std::string& inputFile, const std::string& outputFile) {
    try {
        auto graph = loadEdges(inputFile);
        BipartiteMatching bm(graph);

        std::ofstream out(outputFile);

        if (bm.isBipartite()) {
            out << "Graph is bipartite\n";
            int maxMatching = bm.findMaxMatching();
            out << "Max Matching Size: " << maxMatching << "\n";

            auto matching = bm.getMatching();
            auto color = bm.getColor();

            out << "Matching Edges:\n";
            for (int u = 0; u < graph.size(); u++) {
                if (color[u] == 0 && matching[u] != -1) {
                    out << u << " - " << matching[u] << "\n";
                }
            }
        }
        else {
            out << "Graph is not bipartite\n";
        }

        std::cout << "Bipartite matching completed. Results saved to " << outputFile << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage();
        return 1;
    }

    std::map<std::string, std::string> args;
    std::string inputFile, outputFile = "ans.txt";
    std::string taskType;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-f" || arg == "-p" || arg == "-m") {
            taskType = arg;
        }
        else if (arg == "-o" && i + 1 < argc) {
            outputFile = argv[++i];
        }
        else if (arg == "-a" && i + 1 < argc) {
            args["-a"] = argv[++i];
        }
        else if (arg == "-b" && i + 1 < argc) {
            args["-b"] = argv[++i];
        }
        else if (arg[0] != '-') {
            inputFile = arg;
        }
    }

    if (inputFile.empty()) {
        std::cerr << "Error: Input file is required\n";
        printUsage();
        return 1;
    }

    try {
        if (taskType == "-f") {
            taskMaxFlow(inputFile, outputFile);
        }
        else if (taskType == "-p") {
            if (args["-a"].empty() || args["-b"].empty()) {
                std::cerr << "Error: Start and end points are required for path finding\n";
                printUsage();
                return 1;
            }
            taskPathFinding(inputFile, args["-a"], args["-b"], outputFile);
        }
        else if (taskType == "-m") {
            taskBipartiteMatching(inputFile, outputFile);
        }
        else {
            std::cerr << "Error: Unknown task type\n";
            printUsage();
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}