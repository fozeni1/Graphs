#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include "Graph.h"

void print_menu() {
    std::cout << "\n=== МЕНЮ УПРАВЛЕНИЯ ГРАФОМ ===\n";
    std::cout << "1. Создать пустой граф\n";
    std::cout << "2. Считать граф из файла\n";
    std::cout << "3. Добавить вершину\n";
    std::cout << "4. Добавить ребро\n";
    std::cout << "5. Удалить вершину\n";
    std::cout << "6. Удалить ребро\n";
    std::cout << "7. Вывести список всех вершин\n";
    std::cout << "8. Вывести список всех рёбер\n";
    std::cout << "9. Проверить связность и вывести компоненты\n";
    std::cout << "10. Найти кратчайший путь между вершинами\n";
    std::cout << "11. Найти расстояния от вершины до всех остальных\n";
    std::cout << "12. Вывести минимальное остовное дерево\n";
    std::cout << "13. Сохранить граф в файл\n";
    std::cout << "14. Информация о графе\n";
    std::cout << "0. Выход\n";
    std::cout << "Выберите действие: ";
}

int main() {
    setlocale(LC_ALL, "Russian");
    std::unique_ptr<Graph> graph;
    bool graph_created = false;

    while (true) {
        print_menu();
        int choice;
        std::cin >> choice;

        if (choice == 0) {
            std::cout << "Выход из программы.\n";
            break;
        }

        switch (choice) {
        case 1: {
            std::cout << "Создать ориентированный граф? (1 - да, 0 - нет): ";
            bool directed;
            std::cin >> directed;
            graph = std::make_unique<Graph>(directed);
            graph_created = true;
            std::cout << "Пустой граф создан.\n";
            break;
        }

        case 2: {
            std::string filename;
            std::cout << "Введите путь к файлу: ";
            std::cin >> filename;
            std::cout << "Ориентированный граф? (1 - да, 0 - нет): ";
            bool directed;
            std::cin >> directed;
            graph = std::make_unique<Graph>(filename, directed);
            graph_created = true;
            std::cout << "Граф загружен из файла.\n";
            break;
        }

        case 3: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int vertex;
            std::cout << "Введите номер вершины: ";
            std::cin >> vertex;
            graph->add_vertex(vertex);
            std::cout << "Вершина " << vertex << " добавлена.\n";
            break;
        }

        case 4: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int v1, v2;
            double weight;
            std::cout << "Введите первую вершину: ";
            std::cin >> v1;
            std::cout << "Введите вторую вершину: ";
            std::cin >> v2;
            std::cout << "Введите вес ребра: ";
            std::cin >> weight;
            graph->add_edge(v1, v2, weight);
            std::cout << "Ребро (" << v1 << ", " << v2 << ") добавлено.\n";
            break;
        }

        case 5: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int vertex;
            std::cout << "Введите номер вершины для удаления: ";
            std::cin >> vertex;
            graph->remove_vertex(vertex);
            std::cout << "Вершина " << vertex << " удалена.\n";
            break;
        }

        case 6: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int v1, v2;
            std::cout << "Введите первую вершину: ";
            std::cin >> v1;
            std::cout << "Введите вторую вершину: ";
            std::cin >> v2;
            graph->remove_edge(v1, v2);
            std::cout << "Ребро (" << v1 << ", " << v2 << ") удалено.\n";
            break;
        }

        case 7: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            auto vertices = graph->list_of_vertices();
            std::cout << "Вершины графа: ";
            for (int v : vertices) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
            std::cout << "Всего вершин: " << vertices.size() << std::endl;
            break;
        }

        case 8: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            auto edges = graph->list_of_edges();
            std::cout << "Рёбра графа:\n";
            for (const auto& edge : edges) {
                std::cout << "(" << std::get<0>(edge) << ", " << std::get<1>(edge)
                    << ") вес: " << std::get<2>(edge) << std::endl;
            }
            std::cout << "Всего рёбер: " << edges.size() << std::endl;
            break;
        }

        case 9: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            auto components = graph->connected_components();
            std::cout << "Компоненты связности:\n";
            for (size_t i = 0; i < components.size(); i++) {
                std::cout << "Компонента " << i + 1 << ": ";
                for (int v : components[i]) {
                    std::cout << v << " ";
                }
                std::cout << std::endl;
            }
            std::cout << "Всего компонент: " << components.size() << std::endl;
            break;
        }

        case 10: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int start, end;
            std::cout << "Введите начальную вершину: ";
            std::cin >> start;
            std::cout << "Введите конечную вершину: ";
            std::cin >> end;

            auto [path, distance] = graph->shortest_path(start, end);

            if (path.empty()) {
                std::cout << "Путь не существует!\n";
            }
            else {
                std::cout << "Кратчайший путь: ";
                for (size_t i = 0; i < path.size(); i++) {
                    std::cout << path[i];
                    if (i < path.size() - 1) {
                        std::cout << " -> ";
                    }
                }
                std::cout << "\nДлина пути: " << distance << std::endl;
            }
            break;
        }

        case 11: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            int start;
            std::cout << "Введите начальную вершину: ";
            std::cin >> start;

            auto [distances, prev] = graph->bellman_ford(start);
            auto vertices = graph->list_of_vertices();

            std::cout << "Расстояния от вершины " << start << ":\n";
            for (size_t i = 0; i < vertices.size(); i++) {
                std::cout << "До " << vertices[i] << ": ";
                if (distances[i] == std::numeric_limits<double>::infinity()) {
                    std::cout << "нет пути";
                }
                else {
                    std::cout << distances[i];
                }
                std::cout << std::endl;
            }
            break;
        }

        case 12: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            auto mst = graph->kruskal_mst();
            std::cout << "Минимальное остовное дерево:\n";
            double total_weight = 0;
            for (const auto& edge : mst) {
                std::cout << "(" << std::get<0>(edge) << ", " << std::get<1>(edge)
                    << ") вес: " << std::get<2>(edge) << std::endl;
                total_weight += std::get<2>(edge);
            }
            std::cout << "Общий вес MST: " << total_weight << std::endl;
            break;
        }

        case 13: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            std::string filename;
            std::cout << "Введите имя файла для сохранения: ";
            std::cin >> filename;
            graph->save_to_file(filename);
            std::cout << "Граф сохранен в файл " << filename << std::endl;
            break;
        }

        case 14: {
            if (!graph_created) {
                std::cout << "Сначала создайте граф!\n";
                break;
            }
            std::cout << "Информация о графе:\n";
            std::cout << "Количество вершин: " << graph->size() << std::endl;
            auto edges = graph->list_of_edges();
            std::cout << "Количество рёбер: " << edges.size() << std::endl;
            break;
        }

        default:
            std::cout << "Неверный выбор. Попробуйте снова.\n";
        }
    }

    return 0;
}