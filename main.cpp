#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <omp.h>
#include "data.h"
#include "data.cpp"
#include "ga.cpp"
#include <iomanip>

using namespace std;

struct Box {
    int x, y, z;
    int l, w, h;
    int customerId, cargoId;
};




// 全程保留最好的 10 個 better-but-failed 解


static bool inContainer(const Box& b, int L, int W, int H) {
    return b.x >= 0 && b.y >= 0 && b.z >= 0 &&
           b.x + b.l <= L && b.y + b.w <= W && b.z + b.h <= H;
}

static bool collide(const Box& a, const Box& p) {
    return !(a.x + a.l <= p.x || p.x + p.l <= a.x ||
             a.y + a.w <= p.y || p.y + p.w <= a.y ||
             a.z + a.h <= p.z || p.z + p.h <= a.z);
}


template <class CargoGeneT>
static Box geneToBox(const CargoGeneT& g,
                     const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup) {
    const Cargo& c = cargoLookup.at(g.customerId).at(g.cargoId);

    int l = c.lwh[0], w = c.lwh[1], h = c.lwh[2];

    // 如果你的欄位不是 decodedRotation，請把下面這行改成你的欄位名
    int rot = g.decodedRotation;

    switch (rot) {
        case 1: l = c.lwh[0]; w = c.lwh[1]; h = c.lwh[2]; break;
        case 2: l = c.lwh[0]; w = c.lwh[2]; h = c.lwh[1]; break;
        case 3: l = c.lwh[1]; w = c.lwh[0]; h = c.lwh[2]; break;
        case 4: l = c.lwh[1]; w = c.lwh[2]; h = c.lwh[0]; break;
        case 5: l = c.lwh[2]; w = c.lwh[0]; h = c.lwh[1]; break;
        case 6: l = c.lwh[2]; w = c.lwh[1]; h = c.lwh[0]; break;
        default: break;
    }

    Box b;
    b.customerId = g.customerId;
    b.cargoId = g.cargoId;
    b.x = g.position[0];
    b.y = g.position[1];
    b.z = g.position[2];
    b.l = l; b.w = w; b.h = h;
    return b;
}

static int checkTruckBoxes(const vector<Box>& boxes, int L, int W, int H, int truckTag, bool verbose=true) {
    int viol = 0;

    // bounds + collision
    for (int i = 0; i < (int)boxes.size(); ++i) {
        if (!inContainer(boxes[i], L, W, H)) {
            ++viol;
            if (verbose) {
                cerr << "[VIOL] truck " << truckTag << " out of bounds: cust "
                     << boxes[i].customerId << " cargo " << boxes[i].cargoId << "\n";
            }
        }
        for (int j = i + 1; j < (int)boxes.size(); ++j) {
            if (collide(boxes[i], boxes[j])) {
                ++viol;
                if (verbose) {
                    cerr << "[VIOL] truck " << truckTag << " collision: ("
                         << boxes[i].customerId << "," << boxes[i].cargoId << ") vs ("
                         << boxes[j].customerId << "," << boxes[j].cargoId << ")\n";
                }
            }
        }
    }
    return viol;
}
static int checkIndividualPlacement(const Individual& ind,
                                    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
                                    bool verbose=true) {
    const int L = 300, W = 170, H = 165;
    int viol = 0;

    // self-owned trucks: 1..regionNum
    for (int i = 1; i <= regionNum; ++i) {
        vector<Box> boxes;
        const Truck& t = ind.selfOwnedTrucks[i];
        boxes.reserve(t.assignedCargo.size());
        for (const auto& g : t.assignedCargo) {
            boxes.push_back(geneToBox(g, cargoLookup));
        }
        viol += checkTruckBoxes(boxes, L, W, H, i, verbose);
    }

    // rented trucks: 0..size-1 (tag 用 1000+idx 避免跟自有車混淆)
    for (int k = 0; k < (int)ind.rentedTrucks.size(); ++k) {
        vector<Box> boxes;
        const Truck& t = ind.rentedTrucks[k];
        boxes.reserve(t.assignedCargo.size());
        for (const auto& g : t.assignedCargo) {
            boxes.push_back(geneToBox(g, cargoLookup));
        }
        viol += checkTruckBoxes(boxes, L, W, H, 1000 + k, verbose);
    }

    return viol;
}


bool isBetterCostOnly(const Individual& a, const Individual& b) {
    return a.fitness[0] < b.fitness[0];
}

static SingleTruckGurobiCase buildSingleTruckGurobiCaseFromRecord(
    const BetterButFailedRecord& rec,
    int caseId
) {
    SingleTruckGurobiCase out;

    out.caseId = caseId;
    out.fitnessValue = rec.fitnessValue;

    out.area = rec.area;
    out.truckId = rec.truckId;

    out.tryingCustomer = rec.tryingCustomer;
    out.failedCustomer = rec.failedCustomer;
    out.failedCargo = rec.failedCargo;

    out.containerL = rec.containerL;
    out.containerW = rec.containerW;
    out.containerH = rec.containerH;

    out.candidateCount = rec.candidateCount;
    out.boundaryFailCount = rec.boundaryFailCount;
    out.collisionFailCount = rec.collisionFailCount;
    out.bestBoundaryViolation = rec.bestBoundaryViolation;
    out.bestTotalOverlap = rec.bestTotalOverlap;
    out.collisionDominant = rec.collisionDominant ? 1 : 0;

    out.loadedCustomerList = rec.loadedCustomerList;
    out.truckCaseCustomerList = rec.truckCaseCustomerList;
    out.fullPlannedCustomerList = rec.fullPlannedCustomerList;
    out.remainingCustomerListAfterFail = rec.remainingCustomerListAfterFail;
    out.failIndexInPlannedRoute = rec.failIndexInPlannedRoute;

    // 1 = loaded before fail
    for (const auto& c : rec.baseLoadedCargoes) {
        SingleTruckGurobiCase::Item item;
        item.customerId = c.customerId;
        item.cargoId = c.cargoId;
        item.l = c.l;
        item.w = c.w;
        item.h = c.h;
        item.isPreloaded = 1;
        item.itemRole = 1;
        for (int r = 0; r < 6; ++r) item.orientation[r] = c.orientation[r];
        out.items.push_back(item);
    }

    // 2 = trying customer
    for (const auto& c : rec.tryingCustomerCargoes) {
        SingleTruckGurobiCase::Item item;
        item.customerId = c.customerId;
        item.cargoId = c.cargoId;
        item.l = c.l;
        item.w = c.w;
        item.h = c.h;
        item.isPreloaded = 0;
        item.itemRole = 2;
        for (int r = 0; r < 6; ++r) item.orientation[r] = c.orientation[r];
        out.items.push_back(item);
    }

    // 3 = remaining customers after fail
    for (const auto& c : rec.remainingCustomerCargoes) {
        SingleTruckGurobiCase::Item item;
        item.customerId = c.customerId;
        item.cargoId = c.cargoId;
        item.l = c.l;
        item.w = c.w;
        item.h = c.h;
        item.isPreloaded = 0;
        item.itemRole = 3;
        for (int r = 0; r < 6; ++r) item.orientation[r] = c.orientation[r];
        out.items.push_back(item);
    }

    return out;
}



static void exportSingleTruckGurobiCaseToTxt(
    const SingleTruckGurobiCase& c,
    const std::string& filePath
) {
    std::ofstream fout(filePath);
    if (!fout) {
        std::cerr << "[WARN] cannot open file: " << filePath << "\n";
        return;
    }

    fout << "CASE_ID " << c.caseId << "\n";
    fout << "FITNESS " << std::fixed << std::setprecision(6) << c.fitnessValue << "\n";

    fout << "AREA " << c.area << "\n";
    fout << "TRUCK_ID " << c.truckId << "\n";

    fout << "TRYING_CUSTOMER " << c.tryingCustomer << "\n";
    fout << "FAILED_CUSTOMER " << c.failedCustomer << "\n";
    fout << "FAILED_CARGO " << c.failedCargo << "\n";

    fout << "CONTAINER " << c.containerL << " " << c.containerW << " " << c.containerH << "\n";

    fout << "CANDIDATE_COUNT " << c.candidateCount << "\n";
    fout << "BOUNDARY_FAIL_COUNT " << c.boundaryFailCount << "\n";
    fout << "COLLISION_FAIL_COUNT " << c.collisionFailCount << "\n";
    fout << "BEST_BOUNDARY_VIOLATION " << std::fixed << std::setprecision(6) << c.bestBoundaryViolation << "\n";
    fout << "BEST_TOTAL_OVERLAP " << std::fixed << std::setprecision(6) << c.bestTotalOverlap << "\n";
    fout << "COLLISION_DOMINANT " << c.collisionDominant << "\n";

    fout << "LOADED_CUSTOMERS_BEFORE_FAIL";
    for (int cid : c.loadedCustomerList) {
        fout << " " << cid;
    }
    fout << "\n";

    fout << "TRUCK_CASE_CUSTOMERS";
    for (int cid : c.truckCaseCustomerList) {
        fout << " " << cid;
    }
    fout << "\n";

    fout << "FULL_PLANNED_CUSTOMERS";
    for (int cid : c.fullPlannedCustomerList) {
        fout << " " << cid;
    }
    fout << "\n";

    fout << "REMAINING_CUSTOMERS_AFTER_FAIL";
    for (int cid : c.remainingCustomerListAfterFail) {
        fout << " " << cid;
    }
    fout << "\n";

    fout << "FAIL_INDEX_IN_PLANNED_ROUTE " << c.failIndexInPlannedRoute << "\n";

    fout << "ITEM_COUNT " << c.items.size() << "\n";
    fout << "ITEMS_BEGIN\n";
    fout << "# customerId cargoId l w h isPreloaded itemRole ori1 ori2 ori3 ori4 ori5 ori6\n";

    for (const auto& item : c.items) {
        fout << item.customerId << " "
             << item.cargoId << " "
             << item.l << " "
             << item.w << " "
             << item.h << " "
             << item.isPreloaded << " "
             << item.itemRole << " "
             << item.orientation[0] << " "
             << item.orientation[1] << " "
             << item.orientation[2] << " "
             << item.orientation[3] << " "
             << item.orientation[4] << " "
             << item.orientation[5] << "\n";
    }

    fout << "ITEMS_END\n";
}

static void exportSavedViolationSolutionsToTxt(
    const std::vector<SavedViolationSolution>& sols,
    const std::string& prefix = "saved_sol"
) {
    int globalCaseId = 0;

    for (size_t si = 0; si < sols.size(); ++si) {
        const auto& sol = sols[si];

        for (size_t ri = 0; ri < sol.failedRecords.size(); ++ri) {
            SingleTruckGurobiCase c =
                buildSingleTruckGurobiCaseFromRecord(sol.failedRecords[ri], globalCaseId);

            std::ostringstream oss;
            oss << prefix
                << "_" << si
                << "_case_" << ri
                << "_area" << c.area
                << "_truck" << c.truckId
                << "_tryCust" << c.tryingCustomer
                << ".txt";

            exportSingleTruckGurobiCaseToTxt(c, oss.str());
            ++globalCaseId;
        }
    }
}


static void printFinalBestSummary(const Individual& indiv, const Data& parameters) {
    std::cout << "\n==============================\n";
    std::cout << "===== Global Best Feasible Solution =====\n";
    std::cout << "Best feasible cost = " << indiv.fitness[0] << "\n\n";

    // Self-owned trucks
    for (int area = 1; area <= regionNum; ++area) {
        const auto& truck = indiv.selfOwnedTrucks[area];

        std::cout << "Self-owned Truck (Area " << area << ") route: ";
        if (truck.route.empty()) {
            std::cout << "(empty)";
        } else {
            for (int i = 0; i < (int)truck.route.size(); ++i) {
                std::cout << truck.route[i];
                if (i + 1 < (int)truck.route.size()) std::cout << " -> ";
            }
        }
        std::cout << "\n";

        if (!truck.assignedCargo.empty()) {
            std::cout << "Truck " << area << " cargos:\n";
            for (const auto& g : truck.assignedCargo) {
                std::cout << "  Customer: " << g.customerId
                          << " CargoID: " << g.cargoId
                          << " Position: ("
                          << g.position[0] << ", "
                          << g.position[1] << ", "
                          << g.position[2] << ")"
                          << " Rotation: " << g.decodedRotation
                          << "\n";
            }
        } else {
            std::cout << "Truck " << area << " cargos: (empty)\n";
        }

        std::cout << "\n";
    }

    // Rented trucks
    if (indiv.rentedTrucks.empty()) {
        std::cout << "Rented trucks: none\n";
    } else {
        for (int r = 0; r < (int)indiv.rentedTrucks.size(); ++r) {
            const auto& truck = indiv.rentedTrucks[r];
            std::cout << "Rented Truck " << r << " cargos:\n";
            for (const auto& g : truck.assignedCargo) {
                std::cout << "  Customer: " << g.customerId
                          << " CargoID: " << g.cargoId
                          << " Position: ("
                          << g.position[0] << ", "
                          << g.position[1] << ", "
                          << g.position[2] << ")"
                          << " Rotation: " << g.decodedRotation
                          << "\n";
            }
        }
    }

    // Customer summary
    std::cout << "\n=== [DEBUG] Customer Summary (globalBest) ===\n";
    for (int cid = 1; cid <= Customer; ++cid) {
        bool inSelf = false;
        bool inRented = false;
        int selfArea = -1;

        for (int area = 1; area <= regionNum; ++area) {
            for (const auto& g : indiv.selfOwnedTrucks[area].assignedCargo) {
                if (g.customerId == cid) {
                    inSelf = true;
                    selfArea = area;
                    break;
                }
            }
            if (inSelf) break;
        }

        for (const auto& rt : indiv.rentedTrucks) {
            for (const auto& g : rt.assignedCargo) {
                if (g.customerId == cid) {
                    inRented = true;
                    break;
                }
            }
            if (inRented) break;
        }

        std::cout << "cust " << cid
                  << " | decoded=" << indiv.chromosome[cid - 1].decodedServiceArea
                  << " | selfOwned=" << (inSelf ? "Y" : "N")
                  << " | rented=" << (inRented ? "Y" : "N")
                  << " | selfArea=" << selfArea
                  << " | totalVol=" << parameters.totalVolume[cid - 1]
                  << "\n";
    }

    // Outsourcing compare table
    std::cout << "\n=== [DEBUG] Outsourcing Compare Table (GA side) ===\n";
    for (int cid = 1; cid <= Customer; ++cid) {
        bool inSelf = false;
        bool inRented = false;

        for (int area = 1; area <= regionNum; ++area) {
            for (const auto& g : indiv.selfOwnedTrucks[area].assignedCargo) {
                if (g.customerId == cid) {
                    inSelf = true;
                    break;
                }
            }
            if (inSelf) break;
        }

        for (const auto& rt : indiv.rentedTrucks) {
            for (const auto& g : rt.assignedCargo) {
                if (g.customerId == cid) {
                    inRented = true;
                    break;
                }
            }
            if (inRented) break;
        }

        double vol = parameters.totalVolume[cid - 1];
        int ceilVol = (int)std::ceil(vol);
        double fee = 100 + 30 * std::max(0, ceilVol - 3);

        std::cout << "cust " << cid
                  << " | totalVol=" << vol
                  << " | ceil=" << ceilVol
                  << " | fee=" << fee
                  << " | self=" << (inSelf ? "Y" : "N")
                  << " | rented=" << (inRented ? "Y" : "N")
                  << "\n";
    }

    std::cout << "==============================\n";
}

int main(){
    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();
    srand(time(NULL));
    int noImproveCount = 0;
    const int patience =500;
    double C0 = 6;
    int noImproveObj1Count = 0;
    

    // 讀檔
    string folder = "datasets/N12_A4_S20250102";
    Data parameters;
    readParameters(
        folder + "/customerInfo.csv",
        folder + "/goods.csv",
        folder + "/serviceArea.csv",
        folder + "/routes_star.csv",
        parameters
    );

    if (!parameters.loadGeneratedXY("datasets/N12_A4_S20250102/generatedCoordinates.csv")) {
        return 1;
    }

    // 編碼 & 初始母體生成
    vector<Individual> population = initializePopulation(populationSize, parameters);

std::vector<std::pair<int, double>> obj1_improve_points;
std::vector<std::pair<int, double>> obj1_best_so_far_series;
std::vector<SavedViolationSolution> savedViolationSolutions;

double bestObj1SoFar = std::numeric_limits<double>::infinity();
LoadingCache loadingCache;   // ⭐ 放這裡
Individual globalBest;
bool hasGlobalBest = false;
for (int generation = 0; generation < maxGenerations; ++generation) {
    vector<Individual> undecodedPopulation = population;

    // 建立貨物對照表
    auto cargoLookUp = createCargoLookup(parameters);

    // 解碼
    decodePopulation(population, parameters, cargoLookUp);

//#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)population.size(); ++i) {

        population[i].fitness.clear();
        evaluateFitness(population[i], parameters, loadingCache);
        if (!population[i].fitness.empty() &&
            population[i].fitness[0] < 1e12 &&!population[i].failedTruckRecords.empty() &&population[i].failedTruckRecords.size() <= 2) {

            SavedViolationSolution cand = buildSavedViolationSolution(population[i], 2);
            tryAddSavedViolationSolution(savedViolationSolutions, cand, 5);
        }
        if (!population[i].rentedTrucks.empty()) {
            localSearchImproveByRealObjective(population[i], parameters, cargoLookUp, 5);

            //population[i].fitness.clear();
            //evaluateFitness(population[i], parameters, loadingCache);
            if (!population[i].fitness.empty() &&
            population[i].fitness[0] < 1e12 &&!population[i].failedTruckRecords.empty() &&population[i].failedTruckRecords.size() <= 2) {

            SavedViolationSolution cand = buildSavedViolationSolution(population[i], 2);
            tryAddSavedViolationSolution(savedViolationSolutions, cand, 5);
        
        }
        }
    }

    // ==============================
    // 找本代最佳解（只比 cost）
    // ==============================
    size_t bestFeasibleIdx = population.size();
    double genBestFeasibleCost = std::numeric_limits<double>::infinity();

    for (size_t i = 0; i < population.size(); ++i) {
        if (!population[i].fitness.empty() &&
            population[i].fitness[0] < 1e12 &&
            population[i].fitness[0] < genBestFeasibleCost) {
            genBestFeasibleCost = population[i].fitness[0];
            bestFeasibleIdx = i;
        }
    }

    if (bestFeasibleIdx == population.size()) {
        ++noImproveCount;
        cout << "Generation " << generation
             << " no feasible solution yet"
             << " | noImproveCount = " << noImproveCount
             << '\n';
    } else {
        const Individual& genBest = population[bestFeasibleIdx];

        if (!hasGlobalBest || genBest.fitness[0] < globalBest.fitness[0]) {
            globalBest = genBest;
            hasGlobalBest = true;
            noImproveCount = 0;

            auto cargoLookupExploit = createCargoLookup(parameters);

            Individual baseline = globalBest;
            Individual candidate = globalBest;

            localSearchImproveByRealObjective(candidate, parameters, cargoLookupExploit, 30);

            //candidate.fitness.clear();
            //evaluateFitness(candidate, parameters, loadingCache);
            if (!candidate.fitness.empty() &&
            candidate.fitness[0] < 1e12 &&!candidate.failedTruckRecords.empty() &&candidate.failedTruckRecords.size() <= 2) {

            SavedViolationSolution cand = buildSavedViolationSolution(candidate, 2);
            tryAddSavedViolationSolution(savedViolationSolutions, cand, 5);
            
        }

            int viol2 = checkIndividualPlacement(candidate, cargoLookupExploit, false);

            if (viol2 == 0 && candidate.fitness[0] < baseline.fitness[0]) {
                globalBest = candidate;
            }
        } else {
            ++noImproveCount;
        }

        cout << "[DEBUG] generation " << generation
             << " | gen best feasible = " << genBestFeasibleCost
             << " | global best feasible = " << globalBest.fitness[0]
             << " | noImproveCount = " << noImproveCount
             << '\n';
    }

    if (noImproveCount >= patience) {
        cout << "No improvement in " << patience
             << " generations. Early stopping at generation "
             << generation << ".\n";
        break;
    }

    // ==============================
    // GA reproduction
    // ==============================
    vector<Individual> selectedPopulation = selection(undecodedPopulation, population);
    vector<Individual> crossoveredPopulation = crossoverPopulation(selectedPopulation, crossoverRate);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)crossoveredPopulation.size(); ++i) {
        mutateServiceArea(crossoveredPopulation[i], parameters, mutationRate);
        mutateRotation(crossoveredPopulation[i], mutationRate);
    }

    population = crossoveredPopulation;
}

       // ==============================
    // 最後印 better-but-infeasible top10
    // ==============================
    

  
    // ==============================
    // 印 best feasible solution
    // ==============================
    if (hasGlobalBest) {

// ⭐ FULL packing 多跑幾次，挑最好的穩定結果
{
    Individual rebuilt = globalBest;

    bool ok = reconstructFullPackingFromBestAssignment(
        rebuilt,
        parameters,
        300,   // selfOwnedTrialsPerArea
        100    // rentedTrialsPerCustomer
    );

    if (ok && !rebuilt.fitness.empty() &&
        rebuilt.fitness[0] < 1e12) {
        globalBest = rebuilt;
    }
}
printFinalBestSummary(globalBest, parameters);
cout << "\n===== GA BEST FIXED COST =====\n";

double distSum = 0.0;
for (int area = 1; area <= regionNum; ++area) {
    const auto& seq = globalBest.selfOwnedTrucks[area].route;
    if (seq.empty()) continue;

    distSum += parameters.getDistance(0, seq.front());
    for (int k = 0; k + 1 < (int)seq.size(); ++k) {
        distSum += parameters.getDistance(seq[k], seq[k + 1]);
    }
    distSum += parameters.getDistance(seq.back(), 0);
}
        std::cerr << "[Obj1-B] self-owned fuel (distance) = " << distSum << "\n";

        // ===== 自有車：印路線 =====
        for (int i = 1; i <= regionNum; ++i) {
            const Truck& truck = globalBest.selfOwnedTrucks[i];
            cout << "\nSelf-owned Truck (Area " << i << ") route: ";

            const auto& route = truck.route;
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

        // ===== 印出自有車裝載內容 =====
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

        // ===== 印出租用車裝載內容 =====
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

        auto cargoLookUp2 = createCargoLookup(parameters);
        
        vector<Individual> one{globalBest};
        decodePopulation(one, parameters, cargoLookUp2);
        one[0].fitness.clear();
        evaluateFitness(one[0], parameters, loadingCache);
        
        int viol = checkIndividualPlacement(one[0], cargoLookUp2, true);
        cout << "\n[CHECK] placement violations = " << viol << "\n";

        unordered_set<int> omegaSet;
        for (const auto& t : globalBest.rentedTrucks) {
            for (const auto& g : t.assignedCargo) {
                omegaSet.insert(g.customerId);
            }
        }
        
        cout << "\n=== [DEBUG] Customer Summary (globalBest) ===\n";
        cout << "\n[FINAL] savedViolationSolutions count = "
        << savedViolationSolutions.size() << "\n";

        exportSavedViolationSolutionsToTxt(
            savedViolationSolutions,
            "saved_sol"
        );
        unordered_map<int,int> custArea;
        unordered_set<int> selfOwnedCust, rentedCust;

        for (int a = 1; a <= regionNum; ++a) {
            for (const auto& g : globalBest.selfOwnedTrucks[a].assignedCargo) {
                selfOwnedCust.insert(g.customerId);
                custArea[g.customerId] = a;
            }
        }
        for (size_t k = 0; k < globalBest.rentedTrucks.size(); ++k) {
            for (const auto& g : globalBest.rentedTrucks[k].assignedCargo) {
                rentedCust.insert(g.customerId);
            }
        }

        for (int cid = 1; cid <= Customer; ++cid) {
            cout << "cust " << cid << " | allowed=[";
            for (int b = 1; b <= regionNum; ++b) {
                if (parameters.serviceRegion[cid - 1][b - 1] == 1) cout << b << " ";
            }
            cout << "]";

            int decodedArea = -1;
            for (const auto& gene : globalBest.chromosome) {
                if (gene.customerId == cid) {
                    decodedArea = gene.decodedServiceArea;
                    break;
                }
            }
            cout << " | decoded=" << decodedArea;
            cout << " | selfOwned=" << (selfOwnedCust.count(cid) ? "Y" : "N");
            cout << " | rented=" << (rentedCust.count(cid) ? "Y" : "N");

            if (custArea.count(cid)) cout << " | selfArea=" << custArea[cid];
            else cout << " | selfArea=-";

            cout << " | totalVol=" << parameters.totalVolume[cid - 1];
            cout << "\n";
        }

        cout << "\n=== [DEBUG] Outsourcing Compare Table (GA side) ===\n";
        for (int cid = 1; cid <= Customer; ++cid) {
            double volume = parameters.totalVolume[cid - 1];
            int v = (int)ceil(volume);
            double fee = 100 + 30 * max(0, v - 3);

            cout << "cust " << cid
                 << " | totalVol=" << std::fixed << std::setprecision(6) << volume
                 << " | ceil=" << v
                 << " | fee=" << std::fixed << std::setprecision(6) << fee
                 << " | self=" << (selfOwnedCust.count(cid) ? "Y" : "N")
                 << " | rented=" << (rentedCust.count(cid) ? "Y" : "N")
                 << "\n";
        }
    }

    auto t_end = Clock::now();
    double total_sec = std::chrono::duration<double>(t_end - t_start).count();
    cout << "\nTotal runtime: " << total_sec << " seconds\n";
    
    return 0;
}