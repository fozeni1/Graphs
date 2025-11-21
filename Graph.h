#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <queue>
#include <algorithm>
#include <limits>
#include <set>

class Graph {
private:
    bool directed;
    std::map<int, std::vector<std::pair<int, double>>> adjacencyList;

public:
    // Конструктор для создания пустого графа
    Graph(bool isDirected = false) : directed(isDirected) {}

    // Конструктор из файла со списком рёбер
    Graph(const std::string& filename, bool isDirected) : directed(isDirected) {
        loadFromEdgeList(filename);
    }

    // Загрузка графа из файла со списком рёбер
    void loadFromEdgeList(const std::string& filename) {
        adjacencyList.clear();
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << filename << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            int from, to;
            double weight;

            if (iss >> from >> to >> weight) {
                add_edge(from, to, weight);
            }
        }
        file.close();
    }

    // Количество вершин
    int size() const {
        return adjacencyList.size();
    }

    // Вес ребра между двумя вершинами
    double weight(int vertex1, int vertex2) {
        if (!is_edge(vertex1, vertex2)) {
            return std::numeric_limits<double>::infinity();
        }

        for (const auto& neighbor : adjacencyList[vertex1]) {
            if (neighbor.first == vertex2) {
                return neighbor.second;
            }
        }
        return std::numeric_limits<double>::infinity();
    }

    // Проверка наличия ребра
    bool is_edge(int vertex1, int vertex2) {
        if (adjacencyList.find(vertex1) == adjacencyList.end()) {
            return false;
        }

        for (const auto& neighbor : adjacencyList[vertex1]) {
            if (neighbor.first == vertex2) {
                return true;
            }
        }

        // Для неориентированного графа проверяем в обе стороны
        if (!directed) {
            if (adjacencyList.find(vertex2) == adjacencyList.end()) {
                return false;
            }

            for (const auto& neighbor : adjacencyList[vertex2]) {
                if (neighbor.first == vertex1) {
                    return true;
                }
            }
        }

        return false;
    }

    // Добавление вершины
    void add_vertex(int vertex) {
        if (adjacencyList.find(vertex) == adjacencyList.end()) {
            adjacencyList[vertex] = std::vector<std::pair<int, double>>();
        }
    }

    // Добавление ребра
    void add_edge(int vertex1, int vertex2, double weight = 1.0) {
        add_vertex(vertex1);
        add_vertex(vertex2);

        // Проверяем, существует ли уже такое ребро
        for (auto& neighbor : adjacencyList[vertex1]) {
            if (neighbor.first == vertex2) {
                neighbor.second = weight;
                return;
            }
        }

        adjacencyList[vertex1].push_back({ vertex2, weight });

        // Для неориентированного графа добавляем обратное ребро
        if (!directed) {
            for (auto& neighbor : adjacencyList[vertex2]) {
                if (neighbor.first == vertex1) {
                    neighbor.second = weight;
                    return;
                }
            }
            adjacencyList[vertex2].push_back({ vertex1, weight });
        }
    }

    // Удаление вершины
    void remove_vertex(int vertex) {
        // Удаляем все рёбра, инцидентные вершине
        if (adjacencyList.find(vertex) != adjacencyList.end()) {
            // Удаляем все исходящие рёбра
            adjacencyList.erase(vertex);

            // Удаляем все входящие рёбра
            for (auto& pair : adjacencyList) {
                auto& neighbors = pair.second;
                neighbors.erase(
                    std::remove_if(neighbors.begin(), neighbors.end(),
                        [vertex](const std::pair<int, double>& neighbor) {
                            return neighbor.first == vertex;
                        }),
                    neighbors.end());
            }
        }
    }

    // Удаление ребра
    void remove_edge(int vertex1, int vertex2) {
        if (adjacencyList.find(vertex1) != adjacencyList.end()) {
            auto& neighbors = adjacencyList[vertex1];
            neighbors.erase(
                std::remove_if(neighbors.begin(), neighbors.end(),
                    [vertex2](const std::pair<int, double>& neighbor) {
                        return neighbor.first == vertex2;
                    }),
                neighbors.end());
        }

        if (!directed && adjacencyList.find(vertex2) != adjacencyList.end()) {
            auto& neighbors = adjacencyList[vertex2];
            neighbors.erase(
                std::remove_if(neighbors.begin(), neighbors.end(),
                    [vertex1](const std::pair<int, double>& neighbor) {
                        return neighbor.first == vertex1;
                    }),
                neighbors.end());
        }
    }

    // Список всех вершин
    std::vector<int> list_of_vertices() const {
        std::vector<int> vertices;
        for (const auto& pair : adjacencyList) {
            vertices.push_back(pair.first);
        }
        return vertices;
    }

    // Список всех рёбер
    std::vector<std::tuple<int, int, double>> list_of_edges() {
        std::vector<std::tuple<int, int, double>> edges;
        std::set<std::pair<int, int>> added_edges;

        for (const auto& pair : adjacencyList) {
            int from = pair.first;
            for (const auto& neighbor : pair.second) {
                int to = neighbor.first;
                double weight = neighbor.second;

                // Для неориентированного графа добавляем каждое ребро только один раз
                if (!directed) {
                    int min_vertex = std::min(from, to);
                    int max_vertex = std::max(from, to);
                    if (added_edges.find({ min_vertex, max_vertex }) == added_edges.end()) {
                        edges.push_back({ from, to, weight });
                        added_edges.insert({ min_vertex, max_vertex });
                    }
                }
                else {
                    edges.push_back({ from, to, weight });
                }
            }
        }
        return edges;
    }

    // Список рёбер, инцидентных вершине / дуг, исходящих из вершины
    std::vector<std::pair<int, double>> list_of_edges(int vertex) {
        if (adjacencyList.find(vertex) != adjacencyList.end()) {
            return adjacencyList[vertex];
        }
        return {};
    }

    // Проверка связности и поиск компонент связности
    std::vector<std::vector<int>> connected_components() {
        std::vector<std::vector<int>> components;
        std::set<int> visited;

        for (const auto& pair : adjacencyList) {
            int vertex = pair.first;
            if (visited.find(vertex) == visited.end()) {
                std::vector<int> component;
                dfs_undirected(vertex, visited, component);
                components.push_back(component);
            }
        }

        return components;
    }

    // Алгоритм Беллмана-Форда для поиска кратчайших путей
    std::pair<std::vector<double>, std::vector<int>> bellman_ford(int start) {
        std::vector<int> vertices = list_of_vertices();
        int n = vertices.size();

        // Создаем mapping вершин к индексам
        std::map<int, int> vertex_to_index;
        for (int i = 0; i < n; i++) {
            vertex_to_index[vertices[i]] = i;
        }

        std::vector<double> dist(n, std::numeric_limits<double>::infinity());
        std::vector<int> prev(n, -1);

        int start_index = vertex_to_index[start];
        dist[start_index] = 0;

        // Релаксация рёбер n-1 раз
        for (int i = 0; i < n - 1; i++) {
            for (const auto& edge : list_of_edges()) {
                int u = std::get<0>(edge);
                int v = std::get<1>(edge);
                double w = std::get<2>(edge);

                int u_idx = vertex_to_index[u];
                int v_idx = vertex_to_index[v];

                if (dist[u_idx] != std::numeric_limits<double>::infinity() &&
                    dist[u_idx] + w < dist[v_idx]) {
                    dist[v_idx] = dist[u_idx] + w;
                    prev[v_idx] = u_idx;
                }
            }
        }

        // Проверка на отрицательные циклы
        for (const auto& edge : list_of_edges()) {
            int u = std::get<0>(edge);
            int v = std::get<1>(edge);
            double w = std::get<2>(edge);

            int u_idx = vertex_to_index[u];
            int v_idx = vertex_to_index[v];

            if (dist[u_idx] != std::numeric_limits<double>::infinity() &&
                dist[u_idx] + w < dist[v_idx]) {
                std::cerr << "Граф содержит цикл отрицательного веса!" << std::endl;
                return { {}, {} };
            }
        }

        return { dist, prev };
    }

    // Кратчайший путь между двумя вершинами
    std::pair<std::vector<int>, double> shortest_path(int start, int end) {
        auto [dist, prev] = bellman_ford(start);

        std::vector<int> vertices = list_of_vertices();
        std::map<int, int> vertex_to_index;
        for (int i = 0; i < vertices.size(); i++) {
            vertex_to_index[vertices[i]] = i;
        }

        int start_idx = vertex_to_index[start];
        int end_idx = vertex_to_index[end];

        if (dist[end_idx] == std::numeric_limits<double>::infinity()) {
            return { {}, std::numeric_limits<double>::infinity() };
        }

        // Восстановление пути
        std::vector<int> path;
        int current = end_idx;
        while (current != -1) {
            path.push_back(vertices[current]);
            current = prev[current];
        }
        std::reverse(path.begin(), path.end());

        return { path, dist[end_idx] };
    }

    // Алгоритм Краскала для минимального остовного дерева
    std::vector<std::tuple<int, int, double>> kruskal_mst() {
        // Создаем неориентированную версию графа для MST
        std::vector<std::tuple<int, int, double>> edges = list_of_edges();

        // Сортируем рёбра по весу
        std::sort(edges.begin(), edges.end(),
            [](const std::tuple<int, int, double>& a, const std::tuple<int, int, double>& b) {
                return std::get<2>(a) < std::get<2>(b);
            });

        std::vector<int> vertices = list_of_vertices();
        int n = vertices.size();

        // Создаем mapping вершин к индексам
        std::map<int, int> vertex_to_index;
        for (int i = 0; i < n; i++) {
            vertex_to_index[vertices[i]] = i;
        }

        // Система непересекающихся множеств
        std::vector<int> parent(n);
        std::vector<int> rank(n, 0);
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }

        std::vector<std::tuple<int, int, double>> mst;

        for (const auto& edge : edges) {
            int u = std::get<0>(edge);
            int v = std::get<1>(edge);
            double weight = std::get<2>(edge);

            int u_idx = vertex_to_index[u];
            int v_idx = vertex_to_index[v];

            int root_u = find_set(u_idx, parent);
            int root_v = find_set(v_idx, parent);

            if (root_u != root_v) {
                mst.push_back(edge);
                union_sets(root_u, root_v, parent, rank);
            }
        }

        return mst;
    }

    // Сохранение графа в файл
    void save_to_file(const std::string& filename) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла для записи: " << filename << std::endl;
            return;
        }

        auto edges = list_of_edges();
        for (const auto& edge : edges) {
            file << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << std::endl;
        }

        file.close();
    }

private:
    // DFS для неориентированного графа
    void dfs_undirected(int vertex, std::set<int>& visited, std::vector<int>& component) {
        visited.insert(vertex);
        component.push_back(vertex);

        for (const auto& neighbor : adjacencyList[vertex]) {
            if (visited.find(neighbor.first) == visited.end()) {
                dfs_undirected(neighbor.first, visited, component);
            }
        }
    }

    // Поиск корня множества для DSU
    int find_set(int v, std::vector<int>& parent) {
        if (v == parent[v]) {
            return v;
        }
        return parent[v] = find_set(parent[v], parent);
    }

    // Объединение множеств для DSU
    void union_sets(int a, int b, std::vector<int>& parent, std::vector<int>& rank) {
        a = find_set(a, parent);
        b = find_set(b, parent);

        if (a != b) {
            if (rank[a] < rank[b]) {
                std::swap(a, b);
            }
            parent[b] = a;
            if (rank[a] == rank[b]) {
                rank[a]++;
            }
        }
    }
};

#endif