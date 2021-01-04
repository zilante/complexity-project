#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <memory>
#include <ostream>
#include <unordered_set>
#include <vector>

typedef std::shared_ptr<std::unordered_set<size_t>> size_t_unordered_set_ptr;

class Graph {
 public:
    explicit Graph(size_t vertex_count);

    void AddEdge(size_t vertex1, size_t vertex2);

    void ReadEdges(std::istream& istream);

    void ReadWeightFunction(std::istream& istream);

    // функция, возвращающая 2-аппроксимацию вершинного покрытия минимального веса.
    // может менять в процессе выполнения
    // adjacency_matrix_, weight_function_, degrees_
    size_t_unordered_set_ptr GetApproximatedVertexCover();

    // функция, находящая вес вершинного покрытия минимального веса
    double GetMinVertexCoverWeight() const;

    bool IsVertexCover(const std::unordered_set<size_t>& vertices) const;

    // функция, находящая суммарный вес всех вершин графа
    double GetWeight() const;

    // функция, находящая суммарный вес всех вершин заданного множества
    double GetWeight(const std::unordered_set<size_t>& vertices) const;

 private:
    // возвращает коэффициент C в наибольшей
    // пропорционально-степенной функции
    double GetMaxProportionalDegreeFuncCoefficient() const;

    // вершина считается удаленной <=> ее степень равна 0
    bool IsRemovedVertex(size_t vertex) const;

    // обновление степеней вершин
    void RenewDegrees();

    void RemoveAllEdgesFrom(size_t vertex);

    size_t vertex_count_;

    // матрица смежности
    std::vector<std::vector<bool>> adjacency_matrix_;

    // функция весов вершин
    std::vector<double> weight_function_;

    // степени вершин
    std::vector<size_t> degrees_;
};

Graph::Graph(size_t vertex_count) {
    vertex_count_ = vertex_count;
    weight_function_.resize(vertex_count, 1);
    degrees_.resize(vertex_count, 0);

    adjacency_matrix_.resize(vertex_count, std::vector<bool>());
    for (auto&& row : adjacency_matrix_) {
        row.resize(vertex_count, false);
    }
}

void Graph::AddEdge(size_t vertex1, size_t vertex2) {
    adjacency_matrix_[vertex1][vertex2] = true;
    adjacency_matrix_[vertex2][vertex1] = true;

    if (vertex1 != vertex2) {
        ++degrees_[vertex1];
        ++degrees_[vertex2];
    }
}

void Graph::ReadEdges(std::istream& istream) {
    size_t vertex1 = 0;
    size_t vertex2 = 0;

    while (istream >> vertex1 >> vertex2) {
        if (vertex1 == -1) { // to be deleted
            break;
        }

        AddEdge(vertex1, vertex2);
    }
}

void Graph::ReadWeightFunction(std::istream& istream) {
    for (auto&& vertex_weight : weight_function_) {
        istream >> vertex_weight;
    }
}

bool Graph::IsRemovedVertex(size_t vertex) const {
    return degrees_[vertex] == 0;
}

double Graph::GetMaxProportionalDegreeFuncCoefficient() const {
    double coefficient = 0;
    for (size_t i = 0; i < vertex_count_; ++i) {
        if (!IsRemovedVertex(i)) {
            if (coefficient == 0 ||
                    coefficient > weight_function_[i] / degrees_[i]) {
                coefficient = weight_function_[i] / degrees_[i];
            }
        }
    }

    return coefficient;
}

void Graph::RenewDegrees() {
    for (size_t i = 0; i < vertex_count_; ++i) {
        degrees_[i] = 0;

        for (size_t j = 0; j < vertex_count_; ++j) {
            if (adjacency_matrix_[i][j]) {
                ++degrees_[i];
            }
        }
    }
}

void Graph::RemoveAllEdgesFrom(size_t vertex) {
    for (size_t j = 0; j < vertex_count_; ++j) {
        if (adjacency_matrix_[vertex][j]) {
            adjacency_matrix_[vertex][j] = false;
            adjacency_matrix_[j][vertex] = false;
        }
    }
}

double Graph::GetWeight() const {
    double sum_weight = 0;
    for (auto weight : weight_function_) {
        sum_weight += weight;
    }

    return sum_weight;
}

double Graph::GetWeight(const std::unordered_set<size_t>& vertices) const {
    double sum_weight = 0;
    for (auto vertex : vertices) {
        sum_weight += weight_function_[vertex];
    }

    return sum_weight;
}

size_t_unordered_set_ptr Graph::GetApproximatedVertexCover() {
    size_t_unordered_set_ptr vertex_cover =
            std::make_shared<std::unordered_set<size_t>>();

    size_t new_solution_vertex_count = 0;
    do {
        new_solution_vertex_count = 0;

        double coefficient = GetMaxProportionalDegreeFuncCoefficient();

        for (size_t i = 0; i < vertex_count_; ++i) {
            if (!IsRemovedVertex(i)) {
                weight_function_[i] -=
                        coefficient * static_cast<double>(degrees_[i]);

                const double eps = 0.0001;
                if (fabs(weight_function_[i]) < eps) {
                    vertex_cover->insert(i);
                    ++new_solution_vertex_count;

                    RemoveAllEdgesFrom(i);
                }
            }
        }

        RenewDegrees();
    } while (new_solution_vertex_count != 0);

    return vertex_cover;
}

// если количество подмножеств множества вершин не вмещается
// в size_t, то поведение функции не определено
double Graph::GetMinVertexCoverWeight() const {
    size_t subset_count = static_cast<size_t>(1 << vertex_count_);

    double min_weight = GetWeight();
    for (size_t i = 0; i < subset_count; ++i) {
        size_t mask = i;
        std::unordered_set<size_t> subset = std::unordered_set<size_t>();

        size_t vertex = 0;
        while (mask > 0) {
            if (mask % 2 == 1) {
                subset.insert(vertex);
            }

            vertex++;
            mask /= 2;
        }

        if (IsVertexCover(subset)) {
            double weight = GetWeight(subset);
            if (weight < min_weight) {
                min_weight = weight;
            }
        }
    }

    return min_weight;
}

bool Graph::IsVertexCover(const std::unordered_set<size_t>& vertices) const {
    for (size_t i = 0; i < vertex_count_; ++i) {
        if (vertices.count(i) == 0) {
            for (size_t j = 0; j < vertex_count_; ++j) {
                if (adjacency_matrix_[i][j]) {
                    if (vertices.count(j) == 0) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

template <typename T>
void PrintContainer(const T& container, std::ostream& ostream) {
    for (auto&& element : container) {
        ostream << element << std::endl;
    }
}

template <typename T>
void PrintVector(const std::vector<T>& to_print_vector, std::ostream& ostream) {
    for (auto&& element : to_print_vector) {
        ostream << element << std::endl;
    }
}

void run_test(const std::string& infile_name, const std::string& outfile_name) {
    std::ifstream infile;
    infile.open(infile_name);

    size_t vertex_count = 0;
    infile >> vertex_count;
    Graph graph = Graph(vertex_count);
    graph.ReadWeightFunction(infile);
    graph.ReadEdges(infile);

    infile.close();


    std::ofstream outfile;
    outfile.open(outfile_name);

    Graph copy_graph = Graph(graph);

    auto start = std::chrono::high_resolution_clock::now();
    auto approximated_vertex_cover = copy_graph.GetApproximatedVertexCover();
    auto end = std::chrono::high_resolution_clock::now();
    auto approximation_duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end - start).count();

    bool is_vertex_cover = graph.IsVertexCover(*approximated_vertex_cover);
    double approximated_vertex_cover_weight =
            graph.GetWeight(*approximated_vertex_cover);

    start = std::chrono::high_resolution_clock::now();
    double min_vertex_cover_weight = graph.GetMinVertexCoverWeight();
    end = std::chrono::high_resolution_clock::now();
    auto exhaustive_search_duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end - start).count();

    outfile << "Is correct vertex cover found: " << is_vertex_cover << std::endl;
    outfile << "Found vertex cover weight: " << approximated_vertex_cover_weight
                                             << std::endl;
    outfile << "Min vertex cover weight: " << min_vertex_cover_weight << std::endl;
    outfile << "Approximation searching duration: " << approximation_duration
                                                    << std::endl;
    outfile << "Exhaustive searching duration: " << exhaustive_search_duration
                                                 << std::endl;

    outfile.close();
}

int main() {
    std::string input_dir = "./tests_input/";
    std::string output_dir = "./tests_output/";

    run_test(input_dir + "in.txt", output_dir + "out.txt");

    return 0;
}
