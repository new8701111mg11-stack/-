#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#include "data.h"
#include "data.cpp"
#include "ga.cpp"

using namespace std;

struct Box {
    int x, y, z;
    int l, w, h;
    int customerId, cargoId;
};

static bool inContainer(const Box& b, int L, int W, int H) {
    return b.x >= 0 && b.y >= 0 && b.z >= 0 &&
           b.x + b.l <= L && b.y + b.w <= W && b.z + b.h <= H;
}

static bool collide(const Box& a, const Box& p) {
    return !(a.x + a.l <= p.x || p.x + p.l <= a.x ||
             a.y + a.w <= p.y || p.y + p.w <= a.y ||
             a.z + a.h <= p.z || p.z + p.h <= a.z);
}

// 80% 支撐（你的模型/程式用的規則）
static bool supported80(const Box& b, const vector<Box>& boxes) {
    if (b.z == 0) return true;
    long long base = 1LL * b.l * b.w;
    long long sup  = 0;

    for (const auto& p : boxes) {
        if (p.z + p.h == b.z) {
            int xo = max(0, min(b.x + b.l, p.x + p.l) - max(b.x, p.x));
            int yo = max(0, min(b.y + b.w, p.y + p.w) - max(b.y, p.y));
            sup += 1LL * xo * yo;
        }
    }
    return sup >= (long long)floor(0.8 * base);
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

    // support 80%
    for (const auto& b : boxes) {
        if (!supported80(b, boxes)) {
            ++viol;
            if (verbose) {
                cerr << "[VIOL] truck " << truckTag << " not supported: cust "
                     << b.customerId << " cargo " << b.cargoId << "\n";
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
        viol += checkTruckBoxes(boxes, L, W, H, /*truckTag=*/i, verbose);
    }

    // rented trucks: 0..size-1 (tag 用 1000+idx 避免跟自有車混淆)
    for (int k = 0; k < (int)ind.rentedTrucks.size(); ++k) {
        vector<Box> boxes;
        const Truck& t = ind.rentedTrucks[k];
        boxes.reserve(t.assignedCargo.size());
        for (const auto& g : t.assignedCargo) {
            boxes.push_back(geneToBox(g, cargoLookup));
        }
        viol += checkTruckBoxes(boxes, L, W, H, /*truckTag=*/1000 + k, verbose);
    }

    return viol;
}
// GA 比較（你原本那個，保留）
static bool isBetter(const Individual& a, const Individual& b) {
    const auto& fa = a.fitness;
    const auto& fb = b.fitness;
    if (fa[1] != fb[1]) return fa[1] < fb[1];
    return fa[0] < fb[0];
}
/*
// seed 分區 area(1..regionNum) -> undecodedServiceArea code（讓 decode 後等於 targetArea）
static int serviceAreaToUndecodedCode(const Data& p, int customerId, int targetArea) {
    int idx = customerId - 1;
    vector<int> feasible;
    for (int r = 0; r < regionNum; ++r) {
        if (p.serviceRegion[idx][r] == 1) feasible.push_back(r + 1);
    }
    sort(feasible.begin(), feasible.end());

    if (feasible.size() == 1) return 1;

    for (int k = 0; k < (int)feasible.size(); ++k) {
        if (feasible[k] == targetArea) return k + 1;
    }
    cerr << "[ERROR] customer " << customerId
         << " cannot be assigned to area " << targetArea << "\n";
    return 1;
}

// 如果你之後要用 seed 的 rot 去推 undecodedRotation 才需要；
// 若你走「直接用 seed 的 placedL/W/H」驗算，這個其實可以不必用。
static int rotationToUndecodedRotation(const Cargo& c, int desiredDecodedRot) {
    vector<int> feasible;
    for (int ori = 0; ori < 6; ++ori) if (c.orientation[ori] == 1) feasible.push_back(ori + 1);
    if (feasible.empty()) return 1;

    int k = (int)feasible.size();
    int pos = -1;
    for (int i = 0; i < k; ++i) if (feasible[i] == desiredDecodedRot) { pos = i; break; }

    if (pos < 0) {
        cerr << "[WARN] desired rotation " << desiredDecodedRot << " not feasible\n";
        return 1;
    }
    // 要讓 undecodedRotation % k == pos
    return (pos == 0) ? k : pos;
}
*/
int main(){
    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();
    srand(time(NULL));
    int noImproveCount = 0;
    const int patience = 1000;

    // 讀檔
    string folder = "datasets/N11_A4_S20250102";
    Data parameters;
    readParameters(
    folder + "/customerInfo.csv",
    folder + "/goods.csv",
    folder + "/serviceArea.csv",
    folder + "/routes.csv",
    parameters
);
    
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
    
    auto cargoLookUp2 = createCargoLookup(parameters);

// 重要：globalBest 可能是 undecoded，先 decode 成「可檢查」版本
    vector<Individual> one{globalBest};
    decodePopulation(one, parameters, cargoLookUp2);
    one[0].fitness.clear();
    evaluateFitness(one[0], parameters);

    int viol = checkIndividualPlacement(one[0], cargoLookUp2, /*verbose=*/true);
    cout << "\n[CHECK] placement violations = " << viol << "\n";
   /*
   auto cargoLookup = createCargoLookup(parameters);
 
   // ===== 2) read seed =====
      auto seedCustomer = readSeedCustomerCSV(folder + "/gurobi_seed_customer.csv");
    auto seedCargo    = readSeedCargoCSV(folder + "/gurobi_seed_cargo.csv");

    cerr << "[seed] customer rows=" << seedCustomer.size()
         << ", cargo rows=" << seedCargo.size() << "\n";

    if (seedCustomer.empty() || seedCargo.empty()) {
        cerr << "[ERROR] seed csv not loaded. Please check file paths.\n";
        return 1;
    }

    // （可留可不留：你這支驗算程式其實不用 omegaCustomers）
    unordered_set<int> omegaCustomers;
    for (auto& [cid, sc] : seedCustomer) {
        if (sc.omega == 1) omegaCustomers.insert(cid);
    }

    // ===== 3) build boxes by area（用 seed 的 placedL/W/H）=====
    vector<Box> boxesByArea[regionNum + 1];
    unordered_set<int> customerSeen[regionNum + 1];

    for (const auto& [key, sg] : seedCargo) {
        int cid = (int)(key / 1000LL);
        int tid = (int)(key % 1000LL);

        if (sg.omega == 1) continue; // 外包貨物不檢查裝載

        int area = sg.area;
        if (area < 1 || area > regionNum) continue;

        customerSeen[area].insert(cid);

        Box b;
        b.customerId = cid;
        b.cargoId = tid;
        b.x = sg.x; b.y = sg.y; b.z = sg.z;
        b.l = sg.placedL; b.w = sg.placedW; b.h = sg.placedH;
        boxesByArea[area].push_back(b);
    }

    // ===== 4) fitness[0] volume gap =====
    long long volByArea[regionNum + 1] = {0};
    for (int area = 1; area <= regionNum; ++area) {
        for (int cid : customerSeen[area]) {
            volByArea[area] += parameters.totalVolume[cid - 1];
        }
    }
    long long mx = 0, mn = (1LL << 60);
    for (int area = 1; area <= regionNum; ++area) {
        mx = max(mx, volByArea[area]);
        mn = min(mn, volByArea[area]);
    }
    long long gap = mx - mn;

    // ===== 5) fitness[1] rented cost（照你 GA evaluateFitness 的算法）=====
    long long rentedCost = 0;
    for (auto& [cid, sc] : seedCustomer) {
        if (sc.omega != 1) continue;

        int idx = cid - 1;
        int cargoCnt = parameters.cargoNumber[idx];
        for (int t = 1; t <= cargoCnt; ++t) {
            const Cargo& c = cargoLookup.at(cid).at(t);
            long long chargeUnits = c.volume / 27000;
            rentedCost += chargeUnits * 6;
        }
    }

    long long rentedCost_GA = 0;     // 原本 GA: sum floor(v_it/27000)*6
double obj1_GurobiStyle = 0.0;   // 新增: sum C0*vi*omega

double C0 = 6;

for (const auto& [cid, sc] : seedCustomer) {
    if (sc.omega != 1) continue;

    // (A) GA 原本算法：按貨物離散計價
    int idx = cid - 1;
    int cargoCnt = parameters.cargoNumber[idx];
    for (int t = 1; t <= cargoCnt; ++t) {
        const Cargo& c = cargoLookup.at(cid).at(t);
        long long chargeUnits = c.volume / 27000;
        rentedCost_GA += chargeUnits * 6;
    }

    // (B) Gurobi-style：按客戶總體積線性計價
    obj1_GurobiStyle += C0 * (double)parameters.totalVolume[idx];
}

    // ===== 6) feasibility checks =====
    int L = 300, W = 170, H = 165;
    int viol = 0;

    for (int area = 1; area <= regionNum; ++area) {
        auto& boxes = boxesByArea[area];

        // bounds + collision
        for (int i = 0; i < (int)boxes.size(); ++i) {
            if (!inContainer(boxes[i], L, W, H)) {
                ++viol;
                cerr << "[VIOL] area " << area << " out of bounds: cust "
                     << boxes[i].customerId << " cargo " << boxes[i].cargoId << "\n";
            }
            for (int j = i + 1; j < (int)boxes.size(); ++j) {
                if (collide(boxes[i], boxes[j])) {
                    ++viol;
                    cerr << "[VIOL] area " << area << " collision: ("
                         << boxes[i].customerId << "," << boxes[i].cargoId << ") vs ("
                         << boxes[j].customerId << "," << boxes[j].cargoId << ")\n";
                }
            }
        }

        // support 80%
        for (const auto& b : boxes) {
            if (!supported80(b, boxes)) {
                ++viol;
                cerr << "[VIOL] area " << area << " not supported: cust "
                     << b.customerId << " cargo " << b.cargoId << "\n";
            }
        }
    }

    cout << "=== Gurobi-seed fixed-placement evaluation ===\n";
    cout << "fitness[0] (volume gap) = " << gap << "\n";
    cout << "fitness[1] (rented cost)= " << rentedCost << "\n";
    cout << "placement violations     = " << viol << "\n";
    cout << "fitness[1] (rented cost, GA style)= " << rentedCost_GA << "\n";
cout << "obj1 (Gurobi style: Σ C0*vi*omega)= " << obj1_GurobiStyle << "\n";
*/
double C0 = 6;

unordered_set<int> omegaSet; // GA 端視為 omega=1 的客戶：出現在 rentedTrucks 的客戶
for (const auto& t : globalBest.rentedTrucks) {
    for (const auto& g : t.assignedCargo) {
        omegaSet.insert(g.customerId);
    }
}

// 依客戶總體積計價（線性）
double obj1_gurobi_style = 0.0;
for (int cid : omegaSet) {
    obj1_gurobi_style += C0 * (double)parameters.totalVolume[cid - 1];
}

cout << "\n[Gurobi-style obj1] Σ C0 * v_i * ω_i = " << obj1_gurobi_style << "\n";
cout << "[omega customers] count = " << omegaSet.size() << " : ";
for (int cid : omegaSet) cout << cid << " ";
cout << "\n";
    auto t_end = Clock::now();
    double total_sec = std::chrono::duration<double>(t_end - t_start).count();
    cout << "\nTotal runtime: " << total_sec << " seconds\n";
    return 0;

}
