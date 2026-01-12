// 處理GA相關程式，包括編碼、解碼、計算適應度等等
#include "data.h"
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <unordered_map>
#include <chrono>
#include <utility>
#include <cassert>

using namespace std;

// 產生初始染色體，包含編碼
vector<Individual> initializePopulation(const int population_size, const Data& parameters) {

    // 1. 先生成初始種群，Individual分別代表的是一條染色體，種群內會有 populationSize 個染色體
    vector<Individual> population;
    for (int i = 0; i < population_size; ++i) {
        Individual ind;
        vector <int> remembered;
        unordered_map<int, vector<Gene>> twoRegionCustomerBank;
        // 依照每個區域的路線順序建立基因序列
        for (int region = 0; region < regionNum; ++region) {
            const auto& route = parameters.route[region];
            for (int customerId : route) {
                if (customerId <= 0) continue;
                int idx = customerId - 1;

                int cargoCount = parameters.cargoNumber[idx];
                if (cargoCount <= 0) cargoCount = 1;

                // 計算該顧客可服務的區域數
                int regionCount = 0;
                for (int r = 0; r < regionNum; ++r)
                    if (parameters.serviceRegion[idx][r] == 1) regionCount++;

                if (regionCount == 2) {
                    // 兩區顧客：模板確保一致，但插入前要覆寫 routeArea
                    if (!twoRegionCustomerBank.count(customerId)) {
                        int assignedServiceArea = rand() % 2 + 1; // 1 or 2
                        vector<Gene> templ;
                        templ.reserve(cargoCount);
                        for (int c = 0; c < cargoCount; ++c) {
                            Gene g;
                            g.cargoId     = c + 1;
                            g.customerId  = customerId;
                            g.undecodedServiceArea = assignedServiceArea;
                            g.undecodedRotation    = rand() % 6 + 1;
                            g.routeArea   = -1;                // 先留空，插入時再依 region 設定
                            templ.push_back(g);
                        }
                        twoRegionCustomerBank[customerId] = std::move(templ);
                    }
                    // 插入時把 routeArea 設為當前 region+1
                    const auto& templ = twoRegionCustomerBank[customerId];
                    for (auto g : templ) {
                        g.routeArea = region + 1;             // ★ 關鍵：來源路線區域
                        ind.chromosome.push_back(g);
                    }
                } else {
                    // 非兩區顧客：現場生成並直接設 routeArea
                    int assignedServiceArea = 1; // 你的規則
                    for (int c = 0; c < cargoCount; ++c) {
                        Gene gene;
                        gene.cargoId     = c + 1;
                        gene.customerId  = customerId;
                        gene.undecodedServiceArea = assignedServiceArea;
                        gene.undecodedRotation    = rand() % 6 + 1;
                        gene.routeArea   = region + 1;        // ★ 關鍵：來源路線區域
                        ind.chromosome.push_back(gene);
                    }
                }
            }
        }
        population.push_back(ind);
    }
    return population;
}

// 進行服務區域的解碼
void decodeServiceArea(Individual &indiv, const Data &parameters) {
    for (auto &gene : indiv.chromosome) {
        int customerIdx = gene.customerId - 1; //這是為了對應parameters內的矩陣格式

        vector<int> feasible_regions;
        // 取得客戶的可行服務區域
        for (int area = 0; area < regionNum; area++) {
            if (parameters.serviceRegion[customerIdx][area] == 1) {
                feasible_regions.push_back(area + 1); // 轉換為 {1, 2, 3}
            }
        }

        // 排序確保由小到大
        sort(feasible_regions.begin(), feasible_regions.end());

        if (gene.undecodedServiceArea == 1) {
            // 可行區域必定只有一個，直接使用
            gene.decodedServiceArea = feasible_regions[0];
        } else {
            int idx = gene.undecodedServiceArea - 1; // 編碼值從1開始
            if (idx < feasible_regions.size()) {
                gene.decodedServiceArea = feasible_regions[idx];
            } else {
                cerr << "customer " << gene.customerId << " area wrong" << endl;
            }
        }    
    }
    // 解碼完區域後，可以知道客戶真正的服務區域，接下來對各區域不屬於該路線的客戶進行刪除，留下只屬於該區域的
    indiv.chromosome.erase(
        remove_if(
            indiv.chromosome.begin(),
            indiv.chromosome.end(),
            [](const Gene& g) {
                return g.routeArea != g.decodedServiceArea;
            }),
        indiv.chromosome.end()
    );
}

void decodeCargoRotation(Individual &indiv, const Data &parameters,const unordered_map<int, unordered_map<int, Cargo>> &cargoLookup) {
    for (auto &gene : indiv.chromosome) {
        // 正確透過 customerId 和 cargoId 找到貨物資訊
        const Cargo &cargo = cargoLookup.at(gene.customerId).at(gene.cargoId);
        vector<int> feasibleOrientations;
        for (int ori = 0; ori < 6; ori++) {
            if (cargo.orientation[ori] == 1) {
                feasibleOrientations.push_back(ori + 1);  // 方向1~6
            }
        }

        int orientationCount = feasibleOrientations.size();
        if (orientationCount == 0) {
            cerr << "error " << gene.customerId << "cargo " << gene.cargoId 
            << " no feasible orientation" << endl;
            continue;
        }
        // 依照你原本的解碼規則
        int decodedIndex = (gene.undecodedRotation  % orientationCount);
        gene.decodedRotation = feasibleOrientations[decodedIndex];
    }
}

// 建立貨物對應表
unordered_map<int, unordered_map<int, Cargo>> createCargoLookup(const Data &parameters) {
    unordered_map<int, unordered_map<int, Cargo>> cargoLookup;
    for (const auto &cargo : parameters.cargoInformation) {
        cargoLookup[cargo.customerId][cargo.cargoId] = cargo;
    }
    return cargoLookup;
}

void decodePopulation(vector<Individual>& decodedPopulation, const Data &parameters, const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup){
   for (int i = 0; i < decodedPopulation.size(); ++i) {
        decodeServiceArea(decodedPopulation[i], parameters);
        decodeCargoRotation(decodedPopulation[i], parameters, cargoLookup);
    }
    // printChromosomeInfo(decodedPopulation[0]);
}

void evaluateFitness(Individual &indiv, const Data &parameters) {

    long long rentedVehicleCargoCost = 0;
    long long vechicleLoadedMaxGap = 0;

    Truck selfOwnedTrucks[regionNum + 1]; // 每一區域皆有自己的自有車輛
    vector<Truck> rentedTrucks;
    unordered_map<int, vector<Gene>> regionMap;
    unordered_map<int, bool> isLoadedGlobal; // 紀錄所有裝過的客戶
    vector<int> notLoadedCustomer;           // 紀錄裝不了的客戶
    
    for (int i = 1; i <= regionNum; ++i) { 
        selfOwnedTrucks[i].truckId = i;
    }

    // 先把染色體中出現過的所有客戶初始化為 false（預設都沒裝到）
    for (const auto &gene : indiv.chromosome) {
        isLoadedGlobal[gene.customerId] = false;
        regionMap[gene.decodedServiceArea].push_back(gene);  // 順便建 regionMap
    }

    for (int area = 1; area <= regionNum; ++area) {
        if (regionMap.find(area) == regionMap.end()) continue;
        vector<Gene> &cargoList = regionMap[area];

        // 按染色體順序處理（從後往前）
        unordered_map<int, vector<Gene>> customerGrouped;
        for (int i = (int)cargoList.size() - 1; i >= 0; --i) {
            const Gene& g = cargoList[i];
            customerGrouped[g.customerId].push_back(g);
        }

        // === 自有車輛處理===
        Truck& truck = selfOwnedTrucks[area];
        unordered_set<int> seen;
        
        BLPlacement3D loader(truck.length, truck.width, truck.height);
        loader.setCargoLookup(createCargoLookup(parameters));

        for (const auto& gene : cargoList) {
            if (seen.count(gene.customerId)) continue;
            seen.insert(gene.customerId);

            auto& cargoGroup = customerGrouped[gene.customerId];
            if (loader.tryInsert(cargoGroup,20)) {
                truck.loadedVolume += parameters.totalVolume[gene.customerId - 1];
                truck.assignedCargo.insert(truck.assignedCargo.end(), cargoGroup.begin(), cargoGroup.end());
                isLoadedGlobal[gene.customerId] = true;

                for (const auto& g : cargoGroup) {
                    for (auto& indivGene : indiv.chromosome) {
                        if (indivGene.customerId == g.customerId && indivGene.cargoId == g.cargoId) {
                            indivGene.position[0] = g.position[0];
                            indivGene.position[1] = g.position[1];
                            indivGene.position[2] = g.position[2];
                        }
                    }
                }
            } else {
                // 這個顧客裝不下 → 保持 false，後面交給租用車
                break;
            }
        }
    }

    long long maxVol = 0;
    long long minVol = INT_MAX;
    for (int area = 1; area <= regionNum; ++area) {
        long long v = selfOwnedTrucks[area].loadedVolume;
        if (v > maxVol) maxVol = v;
        if (v < minVol) minVol = v;
    }
    vechicleLoadedMaxGap = maxVol - minVol;
    indiv.fitness.push_back(vechicleLoadedMaxGap);

    // ✅ 這裡一定會包含所有「出現在 chromosome 中」但沒被自有車載到的客戶
    notLoadedCustomer.clear();
    for (const auto& [customerId, loaded] : isLoadedGlobal) {
        if (!loaded) {
            notLoadedCustomer.push_back(customerId);
        }
    }
    sort(notLoadedCustomer.begin(), notLoadedCustomer.end());
    
    size_t cursor = 0;
    int rentedTruckId = 0;
    unordered_set<int> rentedSeen;

    unordered_set<int> notLoadedSet(notLoadedCustomer.begin(), notLoadedCustomer.end());
    
    // 以「顧客」為單位，把沒裝到的客戶分批塞進租用車
    while (cursor < notLoadedCustomer.size()) {

        Truck rentedTruck;
        rentedTruck.truckId = ++rentedTruckId;

        BLPlacement3D loader(rentedTruck.length, rentedTruck.width, rentedTruck.height);
        loader.setCargoLookup(createCargoLookup(parameters));

        bool anyLoadedThisTruck = false;

        for (; cursor < notLoadedCustomer.size(); ++cursor) {
            int custId = notLoadedCustomer[cursor];

            if (rentedSeen.count(custId)) continue;

            // 建 cargoGroup
            vector<Gene> cargoGroup;
            for (const auto& g : indiv.chromosome) {
                if (g.customerId == custId) {
                    cargoGroup.push_back(g);
                }
            }
            if (cargoGroup.empty()) continue;

            bool canLoad = loader.tryInsert(cargoGroup,20);

            if (canLoad) {
                anyLoadedThisTruck = true;
                rentedSeen.insert(custId);

                rentedTruck.assignedCargo.insert(rentedTruck.assignedCargo.end(),
                                                cargoGroup.begin(), cargoGroup.end());
                rentedTruck.route.push_back(custId);

                for (const auto& g : cargoGroup) {
                    const Cargo& c = loader.cargoLookup[g.customerId][g.cargoId];
                    long long chargeUnits = c.volume / 27000;
                    rentedVehicleCargoCost += chargeUnits * 6;
                }
            } else {
                // 分兩種狀況：
                if (!anyLoadedThisTruck) {    
                    rentedSeen.insert(custId);
                    cerr << "Customer " << custId << " cannot be loaded into any rented truck.\n";
                    continue; 
                } else {
                    break;
                }
            }
        }

        if (anyLoadedThisTruck) {
            rentedTrucks.push_back(rentedTruck);
        } else {
            break;
        }
    }

    for (int i = 1; i <= regionNum; ++i) {
        indiv.selfOwnedTrucks[i] = selfOwnedTrucks[i];
    }
    indiv.rentedTrucks = rentedTrucks;
    indiv.fitness.push_back(rentedVehicleCargoCost);
}

vector<Individual> selection(const vector<Individual>& population, const vector<Individual>& decodedPopulation, double eliteRatio = 0.05, int tournamentSize = 2) {
    
    int eliteCount = static_cast<int>(population.size() * eliteRatio);
    vector<int> indices(population.size());
    iota(indices.begin(), indices.end(), 0); // [0, 1, 2, ..., N-1]

    // 依照 fitness 排序 index
    // 這個 fitness 就是最後目標函數值，因為long long太大如果換成小數可能出現錯誤，所以要把fitness越小的擺在最前面
    sort(indices.begin(), indices.end(), [&](int a, int b) {
        const auto& fa = decodedPopulation[a].fitness;
        const auto& fb = decodedPopulation[b].fitness;

        // 先比第二個目標：fitness[1](外用租用成本)
        if (fa[1] != fb[1]) {
            return fa[1] < fb[1];   // 越小排越前面
        }
        // 若 fitness[1] 一樣，再比第一個目標：fitness[0](各車裝載材積差距)
        return fa[0] < fb[0];       // 也是越小排越前面
    });

    // 1. 先選出 top N% elite
    vector<Individual> newPopulation;
    for (int i = 0; i < eliteCount; ++i) {
        newPopulation.push_back(population[indices[i]]);
    }

    // 2. Tournament selection 直到補滿
    while (newPopulation.size() < population.size()) {
        // 隨機抽出 tournamentSize 個 index
        vector<int> tournament;
        for (int i = 0; i < tournamentSize; ++i) {
            int r = rand() % population.size();
            tournament.push_back(r);
        }

        // 找出 tournament 中 fitness 最好的
        int bestIdx = *min_element(tournament.begin(), tournament.end(), [&](int a, int b) {
            const auto& fa = decodedPopulation[a].fitness;
            const auto& fb = decodedPopulation[b].fitness;
            // 先看 fitness[1]，小的視為「比較小」
            if (fa[1] != fb[1]) {
                return fa[1] < fb[1];  // fa[1] 比 fb[1] 小 >> a 比 b 小
            }
            // 若 fitness[1] 相同，再比 fitness[0]
            return fa[0] < fb[0];
        });

        newPopulation.push_back(population[bestIdx]);
    }

    return newPopulation;
}

void crossoverServiceArea(Individual& child1, Individual& child2, double swapProb = 0.5) {
    assert(child1.chromosome.size() == child2.chromosome.size());
    size_t n = child1.chromosome.size();
    if (n == 0) return;

    std::unordered_set<int> visited;   // 確保每個 customer 只處理一次

    for (size_t i = 0; i < n; ++i) {
        int cid = child1.chromosome[i].customerId;

        // 同一個 customer 只決定一次要不要交換
        if (visited.count(cid)) continue;
        visited.insert(cid);

        // 如果設計保證同一位置 customerId 一樣，可以檢查一下：
        if (child1.chromosome[i].customerId != child2.chromosome[i].customerId) {
            // 如果有可能不一樣，也可以改成 continue;
            continue;
        }

        int code1 = child1.chromosome[i].undecodedServiceArea;
        int code2 = child2.chromosome[i].undecodedServiceArea;

        // 若兩條染色體中此客戶的編碼一樣，就完全不交換
        if (code1 == code2) continue;

        // 可選：用機率決定要不要交換這個客戶
        double r = static_cast<double>(rand()) / RAND_MAX;
        if (r > swapProb) continue;

        // 對這個顧客的「所有位置」做交換
        for (size_t j = 0; j < n; ++j) {
            if (child1.chromosome[j].customerId == cid) {
                std::swap(child1.chromosome[j].undecodedServiceArea,
                          child2.chromosome[j].undecodedServiceArea);
            }
        }
    }
}

void crossoverLoadingRotation(Individual& child1, Individual& child2) {
    int N = child1.chromosome.size();

    // 隨機選一個切斷點
    int cutIdx = rand() % (N - 1) + 1; 

    for (int i = cutIdx; i < N; ++i) {
        swap(child1.chromosome[i].undecodedRotation, child2.chromosome[i].undecodedRotation);
    }
}

pair<Individual, Individual> crossover(const Individual& parent1, const Individual& parent2) {
    Individual child1 = parent1;
    Individual child2 = parent2;
    crossoverServiceArea(child1, child2);
    crossoverLoadingRotation(child1, child2);
    return {child1, child2};
}

vector<Individual> crossoverPopulation(const vector<Individual>& selectedPopulation, double crossoverRate) {
    vector<Individual> newPopulation;
    int popSize = selectedPopulation.size();

    vector<int> indices(popSize);
    iota(indices.begin(), indices.end(),0);
    random_device rd;
    mt19937 g(rd());
    shuffle(indices.begin(), indices.end(),g);
    for (int i = 0; i < popSize; i += 2) {
        if (i + 1 >= popSize) break;

        const Individual& parent1 = selectedPopulation[indices[i]];
        const Individual& parent2 = selectedPopulation[indices[i+1]];

        // 交配機率
        double r = (double)rand() / RAND_MAX;
        if (r < crossoverRate) {
            Individual child1, child2;
            tie(child1, child2) = crossover(parent1, parent2);
            newPopulation.push_back(child1);
            newPopulation.push_back(child2);
        } else {
            // 不交配，直接複製父母
            newPopulation.push_back(parent1);
            newPopulation.push_back(parent2);
        }
    }

    // 如果新族群數量多於原來，切掉多餘的
    if (newPopulation.size() > popSize) {
        newPopulation.resize(popSize);
    }

    return newPopulation;
}

void mutateServiceArea(Individual& indiv, const Data& parameters, double mutationRate) {
    // 1. 先找出「可行服務區域數量 > 1」的顧客
    unordered_set<int> multiRegionCustomers; // 只記 "可在多個區域服務" 的客戶ID

    for (int i = 0; i < Customer; ++i) { // i: 0-based index for parameters.serviceRegion
        int cid = i + 1;                 // 客戶ID從 1 開始
        int cnt = 0;
        for (int area = 0; area < regionNum; ++area) {
            if (parameters.serviceRegion[i][area] == 1) {
                ++cnt;
                if (cnt > 1) {
                    multiRegionCustomers.insert(cid);
                    break; // 已經確定>1就可以不用再數了
                }
            }
        }
    }

    // 2. 以顧客為單位處理：同一顧客只決定一次要不要突變
    unordered_set<int> visitedCustomers;

    for (auto& gene : indiv.chromosome) {
        int cid = gene.customerId;
        if (!multiRegionCustomers.count(cid)) continue;

        // 同一個顧客只決定一次
        if (visitedCustomers.count(cid)) continue;
        visitedCustomers.insert(cid);

        double r = static_cast<double>(rand()) / RAND_MAX;
        if (r >= mutationRate) continue;

        int currentEncoding = gene.undecodedServiceArea;
        int newEncoding = (currentEncoding == 1 ? 2 : 1);

        // 將這個顧客在整條 chromosome 中的所有基因一起改
        for (auto& g2 : indiv.chromosome) {
            if (g2.customerId == cid) {
                g2.undecodedServiceArea = newEncoding;
            }
        }
    }
}

void mutateRotation(Individual& indiv, double mutationRate) {
    for (auto& gene : indiv.chromosome) {
        if ((double)rand() / RAND_MAX < mutationRate) {
            int originalRotation = gene.undecodedRotation;
            int newRotation = rand() % 6 + 1;  // 產生 1~6
            while (newRotation == originalRotation) {
                newRotation = rand() % 6 + 1;  // 避免跟原本一樣
            }
            gene.undecodedRotation = newRotation;
        }
    }
}