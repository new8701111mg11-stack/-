# include <iostream>
#include <vector>
#include <chrono>
#include "data.h"
#include "data.cpp"
#include "ga.cpp"

using namespace std;

bool isBetter(const Individual& a, const Individual& b) {
    const auto& fa = a.fitness;
    const auto& fb = b.fitness;

    // 先比 fitness[1]（租用成本），越小越好
    if (fa[1] != fb[1]) {
        return fa[1] < fb[1];
    }
    // 再比 fitness[0]（體積差距），越小越好
    return fa[0] < fb[0];
}

int main(){

    srand(time(NULL));
    int noImproveCount = 0;
    const int patience = 1;

    // 讀檔
    Data parameters;
    readParameters("customerInfo.csv", "goods.csv", "serviceArea.csv", "routes.csv", parameters);
    
    // 編碼 & 初始母體生成
    vector<Individual> population = initializePopulation(populationSize, parameters);
    // printChromosomeInfo(population[0]);
    Individual globalBest;
    bool hasGlobalBest = false;

    for (int generation = 0; generation < maxGenerations; ++generation) {
        vector<Individual> undecodedPopulation = population;

        // 建立貨物對照表，方便進行解碼和貨物裝載對應
        auto cargoLookUp = createCargoLookup(parameters);

        // 解碼
        decodePopulation(population,parameters,cargoLookUp);

        for (size_t i = 0; i < population.size(); ++i) {
            evaluateFitness(population[i], parameters);
        }

        // ====== 這一代的最佳解 ======
        size_t bestIdx = 0;
        for (size_t i = 1; i < population.size(); ++i) {
            if (isBetter(population[i], population[bestIdx])) {
                bestIdx = i;
            }
        }
        const Individual& genBest = population[bestIdx];

        // 印出「這一代」的最佳 fitness（只印值，不印路線）
        // cout << "Generation " << generation
        //      << " best fitness -> "
        //      << "f[1] (rented cost) = " << genBest.fitness[1]
        //      << ", f[0] (volume diff) = " << genBest.fitness[0] << '\n';

        // 更新「全程最佳解」
        bool improvedThisGen = false;
        if (!hasGlobalBest || isBetter(genBest, globalBest)) {
            globalBest = genBest;
            hasGlobalBest = true;
            noImproveCount = 0;
            improvedThisGen = true;
        } else {
            ++noImproveCount;
        }

        cout << "Generation " << generation
        << " global best fitness so far -> "
        << "f[1] (rented cost) = " << globalBest.fitness[1]
        << ", f[0] (volume diff) = " << globalBest.fitness[0] << '\n';

        if (noImproveCount >= patience) {
            cout << "No improvement in " << patience
                << " generations. Early stopping at generation "
                << generation << ".\n";
            break;
        }

        vector<Individual> selectedPopulation = selection(undecodedPopulation, population);
        vector<Individual> crossoveredPopulation = crossoverPopulation(selectedPopulation, crossoverRate);
        for (int i = 0; i < populationSize; ++i) {
            mutateServiceArea(crossoveredPopulation[i], parameters, mutationRate);   
            mutateRotation(crossoveredPopulation[i], mutationRate);
        }
        population = crossoveredPopulation;
        // printChromosomeInfo(population[0]);
    }

    // ====== 所有世代跑完後，印「全程最佳染色體」的完整解 ======
    if (hasGlobalBest) {
        cout << "\n===== Global Best Solution =====\n";
        cout << "Best fitness -> "
             << "f[1] (rented cost) = " << globalBest.fitness[1]
             << ", f[0] (volume diff) = " << globalBest.fitness[0] << '\n';

        // ===== 自有車：印路線 =====
        for (int i = 1; i <= regionNum; ++i) {
            const Truck& truck = globalBest.selfOwnedTrucks[i];
            cout << "\nSelf-owned Truck (Area " << i << ") route: ";

            // 根據 assignedCargo 推出「拜訪客戶順序」（同一客戶只顯示一次）
            vector<int> route;
            unordered_set<int> seen;
            for (const auto& g : truck.assignedCargo) {
                if (!seen.count(g.customerId)) {
                    seen.insert(g.customerId);
                    route.push_back(g.customerId);
                }
            }

            if (route.empty()) {
                cout << "(no customers)\n";
            } else {
                for (size_t j = 0; j < route.size(); ++j) {
                    cout << route[j];
                    if (j + 1 < route.size()) cout << " -> ";
                }
                cout << '\n';
            }
        }

        // 印出自有車的裝載／路線（用你之前的寫法）
        for (int i = 1; i <= regionNum; ++i) {
            const Truck& truck = globalBest.selfOwnedTrucks[i];

            cout << "Truck " << i << " cargos:\n";
            for (const auto& g : truck.assignedCargo) {
                cout << "  Customer: " << g.customerId
                     << " CargoID: " << g.cargoId
                     << " Position: (" << g.position[0] << ", "
                                       << g.position[1] << ", "
                                       << g.position[2] << ")\n";
            }
        }
    
        // 如果也想看租用車，可以再加：
        
        for (size_t k = 0; k < globalBest.rentedTrucks.size(); ++k) {
            const Truck& truck = globalBest.rentedTrucks[k];
            cout << "Rented Truck " << k << " cargos:\n";
            for (const auto& g : truck.assignedCargo) {
                cout << "  Customer: " << g.customerId
                     << " CargoID: " << g.cargoId
                     << " Position: (" << g.position[0] << ", "
                                       << g.position[1] << ", "
                                       << g.position[2] << ")\n";
            }
        }
        
    }
    return 0;
}
