#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <limits>
#include <iomanip>

using namespace std;

class Roads {
public:
    vector<vector<double>> distances;
    vector<vector<double>> pheromones;
    vector<vector<bool>> precedence;

    void PrintPrecedence()
    {
        for (const auto& direction : precedence)
        {
            for (const auto& edge : direction)
            {
                std::cout << edge << ' ';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    void PrintPheromones()
    {
        int counter = 0;
        std::cout << std::fixed << std::setprecision(2);
        for (const auto& direction : pheromones)
        {
            std::cout << counter++ << ": " << std::flush;
            for (const auto& edge : direction)
            {
                std::cout << edge << ' ';
            }
            std::cout << '\n';
        }
        std::cout << std::endl;

    }

    void PrintRoads()
    {
        for (const auto& direction : distances)
        {
            for (const auto& edge : direction)
            {
                std::cout << edge << ' ';
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    Roads(int n) {
        distances.resize(n, vector<double>(n, 0.0));
        pheromones.resize(n, vector<double>(n, 1.0));
        precedence.resize(n, vector<bool>(n, false));
    }

    int size() const { return distances.size(); }

    void generateRandomDistances(unsigned seed, double minDist = 1.0, double maxDist = 10.0) {
        mt19937 rng(seed);
        uniform_int_distribution<int> dist(minDist, maxDist);
        int n = size();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j)
                    distances[i][j] = dist(rng);
                else
                    distances[i][j] = 0;
            }
        }
    }
};

class Ant {
public:
    vector<int> path;
    double cost;

    Ant() : cost(0.0) {}

    void reset() {
        path.clear();
        cost = 0.0;
    }
};

class Solver {
private:
    Roads& roads;
    vector<Ant> ants;
    int numAnts;
    int numIterations;
    int topW;
    double alpha, beta, rho;
    mt19937 rng;
    int beg, end;

public:
    Solver(Roads& r, int nAnts, int nIter, int w, double a, double b, double evaporation, unsigned seed, int start, int finish)
        : roads(r), numAnts(nAnts), numIterations(nIter), topW(w),
        alpha(a), beta(b), rho(evaporation), rng(seed), beg(start), end(finish)
    {
        ants.resize(numAnts);
    }

    bool isFeasible(int node, const vector<int>& visited) {
        for (int pred = 0; pred < roads.size(); ++pred) {
            if (roads.precedence[pred][node] &&
                find(visited.begin(), visited.end(), pred) == visited.end())
                return false;
        }
        return true;
    }

    vector<int> getFeasibleNodes(int current, const vector<int>& visited) {
        vector<int> feasible;
        int n = roads.size();

        for (int j = 0; j < n; ++j) {
            if (find(visited.begin(), visited.end(), j) != visited.end())
                continue;

            if (j == end) {
                if (visited.size() == n - 1) {
                    feasible.push_back(j);
                }
                continue;
            }
            if (isFeasible(j, visited))
                feasible.push_back(j);
        }
        return feasible;
    }

    int selectNextNode(int current, const vector<int>& visited) {
        auto candidates = getFeasibleNodes(current, visited);
        if (candidates.empty()) return -1;

        vector<double> probs;
        double sum = 0.0;
        for (int j : candidates) {
            double tau = pow(roads.pheromones[current][j], alpha);
            double eta = pow(1.0 / roads.distances[current][j], beta);
            probs.push_back(tau * eta);
            sum += tau * eta;
        }
        for (double& p : probs) p /= sum;

        uniform_real_distribution<double> dist(0.0, 1.0);
        double r = dist(rng);
        double cum = 0.0;
        for (size_t i = 0; i < candidates.size(); ++i) {
            cum += probs[i];
            if (r <= cum) return candidates[i];
        }
        return candidates.back();
    }

    double calculateCost(const vector<int>& path) {
        double cost = 0.0;
        for (size_t i = 0; i < path.size() - 1; ++i)
            cost += roads.distances[path[i]][path[i + 1]];
        return cost;
    }

    void constructSolutions() {
        for (Ant& ant : ants) {
            ant.reset();
            int start = beg;
            ant.path.push_back(start);

            while (ant.path.size() < roads.size()) {
                int nextNode = selectNextNode(ant.path.back(), ant.path);
                if (nextNode == -1) break;
                ant.path.push_back(nextNode);
            }

            ant.cost = calculateCost(ant.path);
        }
    }

    void updatePheromones() {
        for (int i = 0; i < roads.size(); ++i)
            for (int j = 0; j < roads.size(); ++j)
                roads.pheromones[i][j] *= (1.0 - rho);

        sort(ants.begin(), ants.end(), [](const Ant& a, const Ant& b) { return a.cost < b.cost; });

        for (int r = 0; r < min(topW, static_cast<int>(ants.size())); ++r) {
            double delta = (topW - r) / ants[r].cost;
            for (size_t i = 0; i < ants[r].path.size() - 1; ++i) {
                int u = ants[r].path[i];
                int v = ants[r].path[i + 1];
                roads.pheromones[u][v] += delta;
            }
        }
    }

    pair<vector<int>, double> solve() {
        vector<int> bestPath;
        double bestCost = numeric_limits<double>::max();

        for (int iter = 0; iter < numIterations; ++iter) {
            constructSolutions();
            updatePheromones();
            roads.PrintPheromones();


            for (const Ant& ant : ants) {
                if (ant.cost < bestCost) {
                    bestCost = ant.cost;
                    bestPath = ant.path;
                }
            }

            std::cout << iter << ": ";
            for (auto const& el : bestPath)
            {
                std::cout << el << ' ';
            }
            std::cout << "Price: " << bestCost << std::endl;
            std::cout << std::endl;
        }

        return { bestPath, bestCost };
    }
};

int main() {
    int nNodes = 20;
    Roads roads(nNodes);
    roads.generateRandomDistances(42, 1, 100);
    roads.precedence[1][3] = true;
    roads.precedence[2][4] = true;
    roads.precedence[3][5] = true;
    roads.PrintRoads();
    roads.PrintPrecedence();

    Solver solver(roads, 10, 50, 3, 3.0, 1.0, 0.1, 123, 0, 4);
    auto [path, cost] = solver.solve();

    cout << "Best path: ";
    for (int node : path) cout << node << " ";
    cout << "\nCost: " << cost << endl;

    return 0;
}
