// 處理GA相關程式，包括編碼、解碼、計算適應度等等
#include "data.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <utility>
#include <cassert>
#include <omp.h>
#include <limits>
#include <cmath>
#include <array>
#include <map>
#include <utility>

using namespace std;
const std::unordered_set<int> TEST_FORCED_OUTSOURCE = {};
const std::unordered_set<int> TEST_FORCED_SELFOWNED = {};

unordered_set<string> successCache;


enum ReduceOp {
    opReinsertion = 0,   // rented -> self
    opConsolidation = 1, // merge rented trucks
    opExchange = 2       // swap rented<->self customer packs
};

// ===== Helper: 從 rentedTrucks 移除某客戶 cid 的所有 Gene =====
static void removeCustomerFromRented(Individual& ind, int cid)
{
    for (auto& rt : ind.rentedTrucks) {
        auto& v = rt.assignedCargo;
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&](const Gene& g){ return g.customerId == cid; }),
                v.end());
        // route 也順便移掉（如果你 route 只是顯示用）
        rt.route.erase(std::remove(rt.route.begin(), rt.route.end(), cid), rt.route.end());
    }

    // 刪掉空的 rented truck（可選，但通常要做）
    ind.rentedTrucks.erase(
        std::remove_if(ind.rentedTrucks.begin(), ind.rentedTrucks.end(),
                       [&](const Truck& t){ return t.assignedCargo.empty(); }),
        ind.rentedTrucks.end()
    );
}
static std::vector<Gene> collectCustomerGenesFromRented(Individual& ind, int cid) {
    std::vector<Gene> pack;

    for (auto itTruck = ind.rentedTrucks.begin(); itTruck != ind.rentedTrucks.end(); ) {
        Truck& rt = *itTruck;

        std::vector<Gene> kept;
        kept.reserve(rt.assignedCargo.size());

        bool foundInThisTruck = false;

        for (const auto& g : rt.assignedCargo) {
            if (g.customerId == cid) {
                pack.push_back(g);
                foundInThisTruck = true;
            } else {
                kept.push_back(g);
            }
        }

        if (foundInThisTruck) {
            rt.assignedCargo.swap(kept);

            // 重建 route：只保留還有貨的 customer，且不重複
            std::vector<int> newRoute;
            std::unordered_set<int> seen;
            for (const auto& g : rt.assignedCargo) {
                if (!seen.count(g.customerId)) {
                    newRoute.push_back(g.customerId);
                    seen.insert(g.customerId);
                }
            }
            rt.route = newRoute;

            // 重新計算 loadedVolume
            rt.loadedVolume = 0;
           

            // 空車就刪掉
            if (rt.assignedCargo.empty()) {
                itTruck = ind.rentedTrucks.erase(itTruck);
            } else {
                ++itTruck;
            }

            // 一個 customer 理論上只會在一台 rented truck
            break;
        } else {
            ++itTruck;
        }
    }

    return pack;
}
// ==============================
// Helper 2: 把某台車的 assignedCargo 全部重新裝箱
// 成功就更新每個 Gene 的 position / rotation
// ==============================
static bool repackTruck(Truck& t,
                        const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
                        int maxTriesInsert = 500)
{
    BLPlacement3D loader(t.length, t.width, t.height);
    loader.setCargoLookup(cargoLookup);
    loader.placedBoxes.clear();

    // 把原本貨物抓出來（整包一起 repack）
    vector<Gene> group = t.assignedCargo;

    // 清掉舊位置（避免用到舊座標）
    for (auto& g : group) {
        g.position[0] = g.position[1] = g.position[2] = -1;
    }

    // tryInsert 會回寫 position / rotation（你已經做了 rotation 探索版）
    bool ok = loader.tryInsert(group, maxTriesInsert);
    if (!ok) return false;

    t.assignedCargo = std::move(group);
    return true;
}
static int encodeServiceAreaChoice(int customerId, int targetArea, const Data& parameters)
{
    int customerIdx = customerId - 1;

    vector<int> feasible_regions;
    for (int area = 0; area < regionNum; ++area) {
        if (parameters.serviceRegion[customerIdx][area] == 1) {
            feasible_regions.push_back(area + 1);
        }
    }

    sort(feasible_regions.begin(), feasible_regions.end());

    for (int i = 0; i < (int)feasible_regions.size(); ++i) {
        if (feasible_regions[i] == targetArea) {
            return i + 1;  // undecodedServiceArea 是從 1 開始
        }
    }

    return -1; // 找不到，代表 targetArea 不可行
}
// ===== 核心 Move：嘗試把 cid 從 rented -> 自有車 area =====
static bool tryMoveRentedCustomerToSelfTruck(Individual& ind,
                                             int cid, int area,
                                             const Data& parameters,
                                             const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    // area 必須可行（服務區限制）
    if (parameters.serviceRegion[cid - 1][area - 1] != 1) return false;

    int encodedAreaChoice = encodeServiceAreaChoice(cid, area, parameters);
    if (encodedAreaChoice == -1) return false;

    // 取出整包貨
    vector<Gene> pack = collectCustomerGenesFromRented(ind, cid);
    if (pack.empty()) return false;

    Truck& selfT = ind.selfOwnedTrucks[area];

    // 暫存，失敗要回復
    auto backupSelf = selfT.assignedCargo;
    auto backupRented = ind.rentedTrucks;
    auto backupChromosome = ind.chromosome;

    const auto originalCargo = selfT.assignedCargo;
    const int n = (int)originalCargo.size();

    for (int pos = 0; pos <= n; ++pos) {
        selfT.assignedCargo = originalCargo;
        ind.rentedTrucks = backupRented;
        ind.chromosome = backupChromosome;

        selfT.assignedCargo.insert(
            selfT.assignedCargo.begin() + pos,
            pack.begin(),
            pack.end()
        );

        if (repackTruck(selfT, cargoLookup, 80)) {
            removeCustomerFromRented(ind, cid);

            // 真正重要的是改 raw encoding，而不是只改 decoded 結果
            for (auto& g : ind.chromosome) {
                if (g.customerId == cid) {
                    g.undecodedServiceArea = encodedAreaChoice;
                    g.routeArea = area;
                    g.decodedServiceArea = area; // 先同步，之後 decode 也會再算一次
                }
            }

            return true;
        }
    }

    selfT.assignedCargo = std::move(backupSelf);
    ind.rentedTrucks = std::move(backupRented);
    ind.chromosome = std::move(backupChromosome);
    return false;
}
static vector<Gene> collectCustomerGenesFromSelf(const Truck& t, int cid) {
    vector<Gene> out;
    for (const auto& g : t.assignedCargo)
        if (g.customerId == cid) out.push_back(g);
    return out;
}

static void removeCustomerFromSelf(Truck& t, int cid) {
    auto& v = t.assignedCargo;
    v.erase(remove_if(v.begin(), v.end(),
                      [&](const Gene& g){ return g.customerId == cid; }),
            v.end());
    t.route.erase(remove(t.route.begin(), t.route.end(), cid), t.route.end());
}



static int rentedTruckCount(const Individual& ind) {
    return (int)ind.rentedTrucks.size();
}


struct SeedCustomer { int area; int omega; };

struct SeedCargo {
    int area;   // 0 if outsourced
    int omega;  // 0/1
    int rot;    // 1..6
    int x, y, z;
    int placedL, placedW, placedH; // ← 新增：用 CSV 的旋轉後尺寸
};

static inline void stripBOM(std::string& s) {
    // UTF-8 BOM: EF BB BF
    if (s.size() >= 3 &&
        (unsigned char)s[0] == 0xEF &&
        (unsigned char)s[1] == 0xBB &&
        (unsigned char)s[2] == 0xBF) {
        s.erase(0, 3);
    }
}

std::unordered_map<int, SeedCustomer>
readSeedCustomerCSV(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open " + path);

    std::unordered_map<int, SeedCustomer> mp;

    std::string line;
    if (!std::getline(fin, line)) return mp;
    stripBOM(line); // header

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string a,b,c;

        std::getline(ss, a, ',');
        std::getline(ss, b, ',');
        std::getline(ss, c, ',');

        int customer = std::stoi(a);
        int area     = std::stoi(b);
        int omega    = std::stoi(c);

        mp[customer] = SeedCustomer{area, omega};
    }
    return mp;
}

std::unordered_map<long long, SeedCargo>
readSeedCargoCSV(const std::string& path) {
    std::ifstream fin(path);
    if (!fin) throw std::runtime_error("Cannot open " + path);

    // key = customer*1000 + cargo  (簡單做法，cargoId 不大)
    std::unordered_map<long long, SeedCargo> mp;

    std::string line;
    if (!std::getline(fin, line)) return mp;
    stripBOM(line); // header

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string sCust,sCargo,sArea,sOmega,sRot,sx,sy,sz;

        std::getline(ss, sCust, ',');
        std::getline(ss, sCargo, ',');
        std::getline(ss, sArea, ',');
        std::getline(ss, sOmega, ',');
        std::getline(ss, sRot, ',');
        std::getline(ss, sx, ',');
        std::getline(ss, sy, ',');
        std::getline(ss, sz, ',');

        int customer = std::stoi(sCust);
        int cargo    = std::stoi(sCargo);

        long long key = (long long)customer * 1000LL + cargo;
        mp[key] = SeedCargo{
            std::stoi(sArea),
            std::stoi(sOmega),
            std::stoi(sRot),
            std::stoi(sx),
            std::stoi(sy),
            std::stoi(sz)
        };
    }
    return mp;
}
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
#pragma omp parallel for
    for (int i = 0; i < decodedPopulation.size(); ++i) {
        decodeServiceArea(decodedPopulation[i], parameters);
        decodeCargoRotation(decodedPopulation[i], parameters, cargoLookup);
    }
    // printChromosomeInfo(decodedPopulation[0]);
}

//local search improve
static long long usedVolumeFromPlaced(const BLPlacement3D& loader) {
    long long used = 0;
    for (const auto& b : loader.placedBoxes) used += 1LL * b.l * b.w * b.h;
    return used;
}

static long long customerGroupVolume(const vector<Gene>& group,
                                     const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup) {
    long long v = 0;
    for (const auto& g : group) {
        v += cargoLookup.at(g.customerId).at(g.cargoId).volume;
    }
    return v;
}

// 回傳：是否成功塞入；同時輸出 bestOverflow (越小越好)
// 若成功：group 會被替換成「成功那次」的擺放結果（position/rotation 也會被 loader 更新）
// 若失敗：group 會被替換成「overflow 最小那次」的 rotation（但仍然沒塞進去）
static bool localSearchRotateAllAndTryInsert(
    BLPlacement3D& loader,
    vector<Gene>& group,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int maxOuterTries,
    int maxInsertTries,
    long long& bestOverflow
) {
    // 計算「剩餘可用體積」(用已放入的 box 估)
    const long long containerVol = 1LL * loader.containerL * loader.containerW * loader.containerH;
    const long long remainVol = max(0LL, containerVol - usedVolumeFromPlaced(loader));

    const long long groupVol = customerGroupVolume(group, cargoLookup);

    // overflow = max(0, groupVol - remainVol)
    bestOverflow = max(0LL, groupVol - remainVol);

    // 留一份最佳 rotation 版本
    vector<Gene> bestGroup = group;

    // 多次嘗試：每次把整包貨物的 rotation 全部改掉
    for (int it = 0; it < maxOuterTries; ++it) {
        vector<Gene> cand = group;

        // ★ 整包換方向：每件貨都換一個隨機 undecodedRotation，然後重解碼成 decodedRotation
        for (auto& g : cand) {
            g.undecodedRotation = rand() % 6 + 1;
        }
        // cand 的 decodedRotation 需要更新（你原本 decodeCargoRotation 就能做）
        // 這裡借用你現成的邏輯：用一個暫時 Individual 包起來解碼
        {
            Individual tmp;
            tmp.chromosome = cand;
            decodeCargoRotation(tmp, Data{}, cargoLookup); // 這行不能用 Data{}，所以下面我給你「正確寫法」
        }

        // ↑ 上面那段不能直接 Data{}，所以更穩的做法是「直接照你的 decodeCargoRotation 寫一個小版」：
        for (auto& g : cand) {
            const Cargo& c = cargoLookup.at(g.customerId).at(g.cargoId);
            vector<int> feasible;
            for (int ori = 0; ori < 6; ++ori) if (c.orientation[ori] == 1) feasible.push_back(ori + 1);
            if (feasible.empty()) { g.decodedRotation = 1; continue; }
            int k = (int)feasible.size();
            int idx = (g.undecodedRotation % k);
            g.decodedRotation = feasible[idx];
        }

        bool ok = loader.tryInsert(cand, maxInsertTries);
        if (ok) {
            group = std::move(cand);
            return true;
        }

        // 仍失敗：overflow 仍用「體積超出」衡量（這裡不因 rotation 改變，值會一樣）
        // 但你也可以改成「若 remainVol 變動（不同車、不同已塞貨）就會不同」，所以仍有意義
        long long overflow = max(0LL, groupVol - remainVol);
        if (overflow < bestOverflow) {
            bestOverflow = overflow;
            bestGroup = cand;
        }
    }

    group = std::move(bestGroup);
    return false;
}
void evaluateFitness(Individual& indiv,
                     const Data& parameters,
                     LoadingCache& loadingCache) {

    const double C1 = 5.41;

    double outsourcedCost = 0.0;
    double selfOwnedFuelCost = 0.0;

    Truck selfOwnedTrucks[regionNum + 1];
    std::unordered_map<int, std::vector<Gene>> regionMap;
    std::unordered_map<int, bool> isLoadedGlobal;

    auto cargoLookup = createCargoLookup(parameters);

    indiv.failedTruckRecords.clear();
    indiv.rentedTrucks.clear();

    for (int area = 1; area <= regionNum; ++area) {
        selfOwnedTrucks[area] = Truck();
        selfOwnedTrucks[area].truckId = area;
        selfOwnedTrucks[area].assignedCargo.clear();
        selfOwnedTrucks[area].route.clear();
        selfOwnedTrucks[area].loadedVolume = 0.0;
    }

    // build region map
    for (const auto& g : indiv.chromosome) {
        isLoadedGlobal[g.customerId] = false;
        regionMap[g.decodedServiceArea].push_back(g);
    }

    // self-owned loading
    for (int area = 1; area <= regionNum; ++area) {

        auto it = regionMap.find(area);
        if (it == regionMap.end()) continue;

        auto& cargoList = it->second;

        std::unordered_map<int, std::vector<Gene>> customerGrouped;
        for (auto& g : cargoList) {
            customerGrouped[g.customerId].push_back(g);
        }

        std::vector<int> customerSet;
        for (auto& kv : customerGrouped) {
            customerSet.push_back(kv.first);
        }

        std::sort(customerSet.begin(), customerSet.end());
        customerSet.erase(
            std::unique(customerSet.begin(), customerSet.end()),
            customerSet.end()
        );

        std::string setKey = buildSetKey(customerSet);

        std::vector<int> routeOrder = buildRouteFromFixedOrder(
            customerSet,
            parameters,
            area
        );

        std::string orderKey = buildOrderKey(routeOrder);

        loadingCache.saveSeen(area, setKey, orderKey);

        Truck& truck = selfOwnedTrucks[area];

        // cache hit: only reuse feasibility, route must follow routeOrder
        if (loadingCache.hasSuccess(area, setKey, orderKey)) {
            for (int cid : routeOrder) {
                isLoadedGlobal[cid] = true;
                truck.loadedVolume += parameters.totalVolume[cid - 1];
                truck.route.push_back(cid);

                const auto& grp = customerGrouped[cid];
                truck.assignedCargo.insert(
                    truck.assignedCargo.end(),
                    grp.begin(),
                    grp.end()
                );
            }
            continue;
        }

        BLPlacement3D loader(truck.length, truck.width, truck.height);
        loader.setCargoLookup(cargoLookup);

        bool allSuccess = true;
        std::vector<int> loadedCustomers;
        std::vector<Gene> loadedCargoSnapshot;

        for (int idx = 0; idx < (int)routeOrder.size(); ++idx) {
            int custId = routeOrder[idx];

            auto cargoGroup = customerGrouped[custId];
            bool success = false;

            bool hasViolationRecord = false;
            PlacementViolationInfo bestTrialInfo;
            std::vector<Gene> successGroup;

            for (int trial = 0; trial < 5 && !success; ++trial) {
                auto temp = cargoGroup;

                if (trial > 0) {
                    static thread_local std::mt19937 rng(std::random_device{}());
                    std::shuffle(temp.begin(), temp.end(), rng);
                }

                if (loader.tryInsert(temp, 50)) {
                    success = true;
                    successGroup = temp;
                } else {
                    if (loader.lastViolationInfo.hasFailedPlacement) {
                        const auto& cur = loader.lastViolationInfo;

                        if (!hasViolationRecord) {
                            bestTrialInfo = cur;
                            hasViolationRecord = true;
                        } else {
                            bool curBetter = false;

                            if (cur.hasBoundaryZeroButCollision &&
                                !bestTrialInfo.hasBoundaryZeroButCollision) {
                                curBetter = true;
                            }
                            else if (cur.hasBoundaryZeroButCollision ==
                                     bestTrialInfo.hasBoundaryZeroButCollision) {

                                if (cur.bestBoundaryViolation <
                                    bestTrialInfo.bestBoundaryViolation) {
                                    curBetter = true;
                                }
                                else if (cur.bestBoundaryViolation ==
                                         bestTrialInfo.bestBoundaryViolation) {

                                    if (cur.bestTotalOverlap <
                                        bestTrialInfo.bestTotalOverlap) {
                                        curBetter = true;
                                    }
                                }
                            }

                            if (curBetter) {
                                bestTrialInfo = cur;
                            }
                        }
                    }
                }
            }

            if (success) {
                loadedCustomers.push_back(custId);
                isLoadedGlobal[custId] = true;
                truck.loadedVolume += parameters.totalVolume[custId - 1];

                truck.assignedCargo.insert(
                    truck.assignedCargo.end(),
                    successGroup.begin(),
                    successGroup.end()
                );

                loadedCargoSnapshot.insert(
                    loadedCargoSnapshot.end(),
                    successGroup.begin(),
                    successGroup.end()
                );
            } else {
                allSuccess = false;

                if (hasViolationRecord) {
                    double currentApproxFitness =
                        outsourcedCost + C1 * selfOwnedFuelCost;

                    std::vector<Gene> remainingCustomerCargoGroup;
                    for (int j = idx + 1; j < (int)routeOrder.size(); ++j) {
                        int futureCid = routeOrder[j];
                        const auto& futureGroup = customerGrouped[futureCid];

                        remainingCustomerCargoGroup.insert(
                            remainingCustomerCargoGroup.end(),
                            futureGroup.begin(),
                            futureGroup.end()
                        );
                    }

                    BetterButFailedRecord rec = buildBetterButFailedRecord(
                        area,
                        truck.truckId,
                        truck.length,
                        truck.width,
                        truck.height,
                        custId,
                        currentApproxFitness,
                        bestTrialInfo,
                        loadedCargoSnapshot,
                        cargoGroup,
                        remainingCustomerCargoGroup,
                        routeOrder,
                        idx,
                        cargoLookup
                    );

                    tryAddFailedRecord(indiv.failedTruckRecords, rec, 10);
                }

                // 保留原本邏輯：該 customer 失敗後，後面的 customer 仍可繼續嘗試
                // 如果你想要「一個 customer 失敗後整台車停止裝載」，這裡才改成 break;
            }
        }

        truck.route = loadedCustomers;

        if (allSuccess) {
            loadingCache.saveSuccess(area, setKey, orderKey);
        }
    }

    // outsourced customers
    std::vector<int> notLoaded;
    for (auto& kv : isLoadedGlobal) {
        if (!kv.second) {
            notLoaded.push_back(kv.first);
        }
    }

    std::sort(notLoaded.begin(), notLoaded.end());
    notLoaded.erase(
        std::unique(notLoaded.begin(), notLoaded.end()),
        notLoaded.end()
    );

    int rentedTruckId = 1001;

    for (int cid : notLoaded) {
        double vol = parameters.totalVolume[cid - 1];
        int v = (int)std::ceil(vol);

        outsourcedCost += 100 + 30 * std::max(0, v - 3);

        Truck rented;
        rented.truckId = rentedTruckId++;
        rented.route.clear();
        rented.route.push_back(cid);
        rented.loadedVolume = vol;

        for (const auto& g : indiv.chromosome) {
            if (g.customerId == cid) {
                rented.assignedCargo.push_back(g);
            }
        }

        indiv.rentedTrucks.push_back(rented);
    }

    // self-owned distance cost
    for (int area = 1; area <= regionNum; ++area) {
        auto& route = selfOwnedTrucks[area].route;
        if (route.empty()) continue;

        double dist = 0.0;

        dist += parameters.getDistance(0, route.front());

        for (int i = 0; i + 1 < (int)route.size(); ++i) {
            dist += parameters.getDistance(route[i], route[i + 1]);
        }

        dist += parameters.getDistance(route.back(), 0);

        selfOwnedFuelCost += dist;
    }

    double obj = outsourcedCost + C1 * selfOwnedFuelCost;

    for (auto& rec : indiv.failedTruckRecords) {
        rec.fitnessValue = obj;
    }

    std::sort(
        indiv.failedTruckRecords.begin(),
        indiv.failedTruckRecords.end(),
        isBetterFailedRecord
    );

    for (int area = 1; area <= regionNum; ++area) {
        indiv.selfOwnedTrucks[area] = selfOwnedTrucks[area];
    }

    indiv.fitness.clear();
    indiv.fitness.push_back(obj);
}

vector<Individual> selection(const vector<Individual>& population,
                             const vector<Individual>& decodedPopulation,
                             double eliteRatio = 0.05,
                             int tournamentSize = 2) {

    int eliteCount = static_cast<int>(population.size() * eliteRatio);
    vector<int> indices(population.size());
    iota(indices.begin(), indices.end(), 0);

    // 依照單一目標 fitness[0]（總成本）排序：越小越好
    sort(indices.begin(), indices.end(), [&](int a, int b) {
        const auto& fa = decodedPopulation[a].fitness;
        const auto& fb = decodedPopulation[b].fitness;
        return fa[0] < fb[0];
    });

    // 1) Elite
    vector<Individual> newPopulation;
    newPopulation.reserve(population.size());

    for (int i = 0; i < eliteCount; ++i) {
        newPopulation.push_back(population[indices[i]]);
    }

    // 2) Tournament selection 補滿
    while (newPopulation.size() < population.size()) {
        vector<int> tournament;
        tournament.reserve(tournamentSize);

        for (int i = 0; i < tournamentSize; ++i) {
            int r = rand() % population.size();
            tournament.push_back(r);
        }

        int bestIdx = *min_element(tournament.begin(), tournament.end(),
                                   [&](int a, int b) {
            const auto& fa = decodedPopulation[a].fitness;
            const auto& fb = decodedPopulation[b].fitness;
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
static bool isBetter(const Individual& a, const Individual& b) {
    return a.fitness[0] < b.fitness[0];
}
static vector<int> getCustomerOrderFromChromosome(const Individual& ind)
{
    vector<int> order;
    unordered_set<int> seen;
    order.reserve(ind.chromosome.size());

    for (const auto& g : ind.chromosome) {
        if (!seen.count(g.customerId)) {
            seen.insert(g.customerId);
            order.push_back(g.customerId);
        }
    }
    return order;
}
static void swapCustomerBlocksInChromosome(Individual& ind, int c1, int c2)
{
    if (c1 == c2) return;

    vector<Gene> block1, block2;
    vector<Gene> newChrom;
    newChrom.reserve(ind.chromosome.size());

    for (const auto& g : ind.chromosome) {
        if (g.customerId == c1) block1.push_back(g);
        else if (g.customerId == c2) block2.push_back(g);
    }

    bool inserted = false;
    for (size_t i = 0; i < ind.chromosome.size(); ++i) {
        int cid = ind.chromosome[i].customerId;

        if (cid == c1 || cid == c2) {
            if (!inserted) {
                newChrom.insert(newChrom.end(), block2.begin(), block2.end());
                newChrom.insert(newChrom.end(), block1.begin(), block1.end());
                inserted = true;
            }

            while (i + 1 < ind.chromosome.size() &&
                   (ind.chromosome[i + 1].customerId == c1 ||
                    ind.chromosome[i + 1].customerId == c2)) {
                ++i;
            }
        } else {
            newChrom.push_back(ind.chromosome[i]);
        }
    }

    if (!newChrom.empty() && newChrom.size() == ind.chromosome.size()) {
        ind.chromosome = std::move(newChrom);
    }
}
static void moveCustomerBlockInChromosome(Individual& ind, int movingCid, int beforeCid)
{
    if (movingCid == beforeCid) return;

    vector<Gene> movingBlock;
    vector<Gene> rest;
    rest.reserve(ind.chromosome.size());

    for (const auto& g : ind.chromosome) {
        if (g.customerId == movingCid) movingBlock.push_back(g);
        else rest.push_back(g);
    }

    if (movingBlock.empty()) return;

    vector<Gene> newChrom;
    newChrom.reserve(ind.chromosome.size());

    bool inserted = false;
    size_t i = 0;
    while (i < rest.size()) {
        int cid = rest[i].customerId;

        if (!inserted && cid == beforeCid) {
            newChrom.insert(newChrom.end(), movingBlock.begin(), movingBlock.end());
            inserted = true;
        }

        newChrom.push_back(rest[i]);
        ++i;
    }

    if (!inserted) {
        newChrom.insert(newChrom.end(), movingBlock.begin(), movingBlock.end());
    }

    if (newChrom.size() == ind.chromosome.size()) {
        ind.chromosome = std::move(newChrom);
    }
}
// ==============================
// Helper 3: 重新解碼並重算 fitness
// ==============================
static void reDecodeAndEvaluate(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    std::vector<Individual> one{ind};
    LoadingCache loadingCache;   // ⭐ 放這裡
    decodePopulation(one, parameters, cargoLookup);

    one[0].fitness.clear();
    evaluateFitness(one[0], parameters,loadingCache);

    ind = one[0];
}
static bool betterObj(const Individual& a, const Individual& b, double eps = 1e-9)
{
    return a.fitness[0] < b.fitness[0] - eps;
}

static bool doReinsertionOnceStable(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    unordered_set<int> omegaSet;
    for (const auto& rt : ind.rentedTrucks)
        for (const auto& g : rt.assignedCargo)
            omegaSet.insert(g.customerId);

    vector<int> omega(omegaSet.begin(), omegaSet.end());
    sort(omega.begin(), omega.end());

    for (int cid : omega) {
        for (int area = 1; area <= regionNum; ++area) {
            if (parameters.serviceRegion[cid - 1][area - 1] != 1) continue;

            Individual backup = ind;
            if (tryMoveRentedCustomerToSelfTruck(ind, cid, area, parameters, cargoLookup)) {
                return true;
            }
            ind = std::move(backup);
        }
    }
    return false;
}

static bool doAreaRebuildForRentedOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    const double EPS = 1e-9;

    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);
    if (base.fitness.empty()) return false;

    double baseObj = base.fitness[0];

    // 1) 找 rented customers，優先挑外包成本高的
    unordered_set<int> rentedSet;
    for (const auto& rt : base.rentedTrucks) {
        for (const auto& g : rt.assignedCargo) {
            rentedSet.insert(g.customerId);
        }
    }
    if (rentedSet.empty()) return false;

    auto outsourceFee = [&](int cid) {
        double vol = parameters.totalVolume[cid - 1];
        int v = (int)std::ceil(vol);
        return 100 + 30 * std::max(0, v - 3);
    };

    vector<int> rentedCustomers(rentedSet.begin(), rentedSet.end());
    sort(rentedCustomers.begin(), rentedCustomers.end(),
         [&](int a, int b) {
             return outsourceFee(a) > outsourceFee(b);
         });

    // helper: 把 cid 移到另一個 feasible area；若做不到就回傳 false
    auto moveCustomerToAnotherArea =
        [&](Individual& cand, int cid, int forbiddenArea) -> bool
    {
        for (int altArea = 1; altArea <= regionNum; ++altArea) {
            if (altArea == forbiddenArea) continue;
            if (parameters.serviceRegion[cid - 1][altArea - 1] != 1) continue;

            int encOut = encodeServiceAreaChoice(cid, altArea, parameters);
            if (encOut == -1) continue;

            for (auto& g : cand.chromosome) {
                if (g.customerId == cid) {
                    g.undecodedServiceArea = encOut;
                    g.routeArea = altArea;
                    g.decodedServiceArea = altArea;
                }
            }
            return true;
        }
        return false;
    };

    // 2) 逐個 rented customer 試
    for (int cidR : rentedCustomers) {

        // cidR 的可行區域
        vector<int> feasibleAreasR;
        for (int area = 1; area <= regionNum; ++area) {
            if (parameters.serviceRegion[cidR - 1][area - 1] == 1) {
                feasibleAreasR.push_back(area);
            }
        }

        // 3) 對每個可行區域做區域重組
        for (int targetArea : feasibleAreasR) {

            unordered_set<int> areaSet;
            for (const auto& g : base.selfOwnedTrucks[targetArea].assignedCargo) {
                areaSet.insert(g.customerId);
            }

            vector<int> areaCustomers(areaSet.begin(), areaSet.end());
            sort(areaCustomers.begin(), areaCustomers.end(),
                 [&](int a, int b) {
                     return parameters.totalVolume[a - 1] < parameters.totalVolume[b - 1];
                 });

            // -------- Case A: 原區域 + cidR --------
            {
                Individual cand = base;

                int encR = encodeServiceAreaChoice(cidR, targetArea, parameters);
                if (encR != -1) {
                    for (auto& g : cand.chromosome) {
                        if (g.customerId == cidR) {
                            g.undecodedServiceArea = encR;
                            g.routeArea = targetArea;
                            g.decodedServiceArea = targetArea;
                        }
                    }

                    reDecodeAndEvaluate(cand, parameters, cargoLookup);
                    if (!cand.fitness.empty() && cand.fitness[0] < baseObj - EPS) {
                        ind = std::move(cand);
                        return true;
                    }
                }
            }

            // -------- Case B: 移掉 1 個 customer，再把 cidR 放進來 --------
            for (int cidOut1 : areaCustomers) {

                Individual cand = base;

                // 搬不走就直接放棄這個 candidate，不要讓客戶消失
                if (!moveCustomerToAnotherArea(cand, cidOut1, targetArea)) {
                    continue;
                }

                int encR = encodeServiceAreaChoice(cidR, targetArea, parameters);
                if (encR == -1) continue;

                for (auto& g : cand.chromosome) {
                    if (g.customerId == cidR) {
                        g.undecodedServiceArea = encR;
                        g.routeArea = targetArea;
                        g.decodedServiceArea = targetArea;
                    }
                }

                reDecodeAndEvaluate(cand, parameters, cargoLookup);
                if (!cand.fitness.empty() && cand.fitness[0] < baseObj - EPS) {
                    ind = std::move(cand);
                    return true;
                }
            }

            // -------- Case C: 移掉 2 個 customer，再把 cidR 放進來 --------
            for (int i = 0; i < (int)areaCustomers.size(); ++i) {
                for (int j = i + 1; j < (int)areaCustomers.size(); ++j) {
                    int cidOut1 = areaCustomers[i];
                    int cidOut2 = areaCustomers[j];

                    Individual cand = base;

                    if (!moveCustomerToAnotherArea(cand, cidOut1, targetArea)) {
                        continue;
                    }
                    if (!moveCustomerToAnotherArea(cand, cidOut2, targetArea)) {
                        continue;
                    }

                    int encR = encodeServiceAreaChoice(cidR, targetArea, parameters);
                    if (encR == -1) continue;

                    for (auto& g : cand.chromosome) {
                        if (g.customerId == cidR) {
                            g.undecodedServiceArea = encR;
                            g.routeArea = targetArea;
                            g.decodedServiceArea = targetArea;
                        }
                    }

                    reDecodeAndEvaluate(cand, parameters, cargoLookup);
                    if (!cand.fitness.empty() && cand.fitness[0] < baseObj - EPS) {
                        ind = std::move(cand);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}
static bool doExchangeOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    const double EPS = 1e-9;

    // 先把 base 解整理成最新狀態
    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);
    if (base.fitness.empty()) return false;

    double baseObj = base.fitness[0];

    // 找 rented customers
    unordered_set<int> rentedSet;
    for (const auto& rt : base.rentedTrucks) {
        for (const auto& g : rt.assignedCargo) {
            rentedSet.insert(g.customerId);
        }
    }
    if (rentedSet.empty()) return false;

    vector<int> rentedCustomers(rentedSet.begin(), rentedSet.end());

    // 可選：優先試外包成本高的
    auto outsourceFee = [&](int cid) {
        double vol = parameters.totalVolume[cid - 1];
        int v = (int)std::ceil(vol);
        return 100 + 30 * std::max(0, v - 3);
    };

    sort(rentedCustomers.begin(), rentedCustomers.end(),
         [&](int a, int b) {
             return outsourceFee(a) > outsourceFee(b);
         });

    for (int cidR : rentedCustomers) {

        // cidR 的所有可行區域
        vector<int> feasibleAreasR;
        for (int area = 1; area <= regionNum; ++area) {
            if (parameters.serviceRegion[cidR - 1][area - 1] == 1) {
                feasibleAreasR.push_back(area);
            }
        }

        for (int targetArea : feasibleAreasR) {

            // 目標區域目前有哪些 self-owned customers
            unordered_set<int> selfSet;
            for (const auto& g : base.selfOwnedTrucks[targetArea].assignedCargo) {
                selfSet.insert(g.customerId);
            }
            if (selfSet.empty()) continue;

            vector<int> selfCustomers(selfSet.begin(), selfSet.end());

            // 可選：先搬比較小 / 比較便宜的
            sort(selfCustomers.begin(), selfCustomers.end(),
                 [&](int a, int b) {
                     return parameters.totalVolume[a - 1] < parameters.totalVolume[b - 1];
                 });

            for (int cidS : selfCustomers) {
                if (cidS == cidR) continue;

                // cidS 必須還有其他 feasible area 可搬
                vector<int> altAreasS;
                for (int area = 1; area <= regionNum; ++area) {
                    if (area == targetArea) continue;
                    if (parameters.serviceRegion[cidS - 1][area - 1] == 1) {
                        altAreasS.push_back(area);
                    }
                }
                if (altAreasS.empty()) continue;

                for (int altArea : altAreasS) {
                    int encodedChoiceS = encodeServiceAreaChoice(cidS, altArea, parameters);
                    if (encodedChoiceS == -1) continue;

                    Individual cand = base;

                    // 第一步：先把 cidS 從 targetArea 搬到 altArea
                    for (auto& g : cand.chromosome) {
                        if (g.customerId == cidS) {
                            g.undecodedServiceArea = encodedChoiceS;
                            g.routeArea = altArea;
                            g.decodedServiceArea = altArea;
                        }
                    }

                    reDecodeAndEvaluate(cand, parameters, cargoLookup);
                    if (cand.fitness.empty()) continue;

                    // 第二步：再把 rented 的 cidR 塞進 targetArea
                    if (!tryMoveRentedCustomerToSelfTruck(cand, cidR, targetArea, parameters, cargoLookup)) {
                        continue;
                    }

                    reDecodeAndEvaluate(cand, parameters, cargoLookup);
                    if (cand.fitness.empty()) continue;

                    if (cand.fitness[0] < baseObj - EPS) {
                        ind = std::move(cand);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}
bool doKickOutAndInsertOnce(
    Individual& ind,
    const Data& parameters,
    const std::unordered_map<int, std::unordered_map<int, Cargo>>& cargoLookup
) {
    if (ind.rentedTrucks.empty()) return false;

    const double EPS = 1e-9;

    Individual base = ind;
    Individual best = ind;

    double baseCost = ind.fitness.empty() ? 1e18 : ind.fitness[0];
    double bestCost = baseCost;

    bool improved = false;

    // 1. collect rented customers
    std::vector<int> rentedCustomers;
    for (const auto& rt : ind.rentedTrucks) {
        if (!rt.assignedCargo.empty()) {
            int cid = rt.assignedCargo.front().customerId;
            if (std::find(rentedCustomers.begin(), rentedCustomers.end(), cid)
                == rentedCustomers.end()) {
                rentedCustomers.push_back(cid);
            }
        }
    }

    if (rentedCustomers.empty()) return false;

    // 2. collect self-owned customers by area
    std::vector<int> selfCustomersByArea[regionNum + 1];

    for (int area = 1; area <= regionNum; ++area) {
        for (int cid : ind.selfOwnedTrucks[area].route) {
            selfCustomersByArea[area].push_back(cid);
        }
    }

    // 3. Try insert + kick
    for (int insertCid : rentedCustomers) {

        for (int area = 1; area <= regionNum; ++area) {

            if (!parameters.serviceRegion[insertCid - 1][area - 1]) {
                continue;
            }

            for (int kickCid : selfCustomersByArea[area]) {

                if (kickCid == insertCid) continue;

                Individual cand = base;

                for (auto& g : cand.chromosome) {

                    if (g.customerId == insertCid) {
                        g.decodedServiceArea = area;
                    }
                    else if (g.customerId == kickCid) {
                        g.decodedServiceArea = -1;
                    }
                }

                LoadingCache tempCache;

                evaluateFitness(cand, parameters, tempCache);

                if (cand.fitness.empty()) continue;

                double candCost = cand.fitness[0];

                if (candCost + EPS < bestCost) {
                    bestCost = candCost;
                    best = cand;
                    improved = true;
                }
            }
        }
    }

    if (improved) {
        ind = best;
        return true;
    }

    return false;
}
static bool doCrossAreaRelocateOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    const double EPS = 1e-9;

    // baseline
    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);
    if (base.fitness.empty()) return false;

    double baseObj = base.fitness[0];

    // 先找目前在 self-owned 的 customer 及其所在 area
    vector<pair<int, int>> selfCustomers; // (cid, currentArea)
    unordered_set<int> seenCustomer;

    for (int area = 1; area <= regionNum; ++area) {
        for (const auto& g : ind.selfOwnedTrucks[area].assignedCargo) {
            if (!seenCustomer.count(g.customerId)) {
                selfCustomers.push_back({g.customerId, area});
                seenCustomer.insert(g.customerId);
            }
        }
    }

    if (selfCustomers.empty()) return false;

    // 可加排序，讓結果較穩定
    sort(selfCustomers.begin(), selfCustomers.end());

    for (const auto& [cid, currentArea] : selfCustomers) {

        // 找此 customer 的所有 feasible areas
        vector<int> feasibleAreas;
        int customerIdx = cid - 1;
        for (int area = 1; area <= regionNum; ++area) {
            if (parameters.serviceRegion[customerIdx][area - 1] == 1) {
                feasibleAreas.push_back(area);
            }
        }

        // 只有一個 feasible area 就沒得搬
        if ((int)feasibleAreas.size() <= 1) continue;

        for (int targetArea : feasibleAreas) {
            if (targetArea == currentArea) continue;

            int encodedChoice = encodeServiceAreaChoice(cid, targetArea, parameters);
            if (encodedChoice == -1) continue;

            Individual cand = ind;

            // 關鍵：改 raw encoding，讓 reDecodeAndEvaluate 真的吃到
            for (auto& g : cand.chromosome) {
                if (g.customerId == cid) {
                    g.undecodedServiceArea = encodedChoice;
                    g.routeArea = targetArea;
                    g.decodedServiceArea = targetArea; // 先同步，之後 decode 還會重算
                }
            }

            reDecodeAndEvaluate(cand, parameters, cargoLookup);
            if (cand.fitness.empty()) continue;

            if (cand.fitness[0] < baseObj - EPS) {
                ind = std::move(cand);
                return true;
            }
        }
    }

    return false;
}
static std::vector<OperatorStats> gLocalSearchOpStats(4);

static void printLocalSearchOperatorStats()
{
    cout << "\n===== Local Search Operator Stats "
         << "(OP0 + CrossAreaRelocate + AreaRebuild + KickOutInsert) =====\n";

    for (int i = 0; i < 4; ++i) {
        double returnRate = 0.0;
        double acceptRate = 0.0;
        double avgImprovement = 0.0;

        if (gLocalSearchOpStats[i].selectedCount > 0) {
            returnRate =
                static_cast<double>(gLocalSearchOpStats[i].returnedTrueCount)
                / gLocalSearchOpStats[i].selectedCount;

            acceptRate =
                static_cast<double>(gLocalSearchOpStats[i].acceptedCount)
                / gLocalSearchOpStats[i].selectedCount;
        }

        if (gLocalSearchOpStats[i].acceptedCount > 0) {
            avgImprovement =
                gLocalSearchOpStats[i].totalImprovement
                / gLocalSearchOpStats[i].acceptedCount;
        }

        cout << "Operator " << i
             << " | selected = " << gLocalSearchOpStats[i].selectedCount
             << " | returnedTrue = " << gLocalSearchOpStats[i].returnedTrueCount
             << " | accepted = " << gLocalSearchOpStats[i].acceptedCount
             << " | returnRate = " << returnRate
             << " | acceptRate = " << acceptRate
             << " | totalImprovement = " << gLocalSearchOpStats[i].totalImprovement
             << " | avgImprovement = " << avgImprovement
             << " | bestImprovement = " << gLocalSearchOpStats[i].bestImprovement
             << "\n";
    }
}
static void localSearchImproveByRealObjective(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int budgetMoves = 100)
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pickOp(0, 3);

    const double EPS = 1e-9;

    reDecodeAndEvaluate(ind, parameters, cargoLookup);
    if (ind.fitness.empty()) return;

    double bestObj = ind.fitness[0];

    for (int iter = 0; iter < budgetMoves; ++iter) {
        int op = pickOp(rng);
        gLocalSearchOpStats[op].selectedCount++;

        Individual backup = ind;
        bool moved = false;
        double oldObj = bestObj;

        switch (op) {
        case 0:
            moved = doReinsertionOnceStable(ind, parameters, cargoLookup);
            break;

        case 1:
            moved = doCrossAreaRelocateOnce(ind, parameters, cargoLookup);
            break;

        case 2:
            moved = doAreaRebuildForRentedOnce(ind, parameters, cargoLookup);
            break;

        case 3:
            moved = doKickOutAndInsertOnce(ind, parameters, cargoLookup);
            break;

        default:
            break;
        }

        if (!moved) {
            ind = std::move(backup);
            continue;
        }

        gLocalSearchOpStats[op].returnedTrueCount++;

        if (op != 3) {
            reDecodeAndEvaluate(ind, parameters, cargoLookup);
        }

        if (!ind.fitness.empty() && ind.fitness[0] < bestObj - EPS) {
            double improvement = oldObj - ind.fitness[0];

            bestObj = ind.fitness[0];

            gLocalSearchOpStats[op].acceptedCount++;
            gLocalSearchOpStats[op].totalImprovement += improvement;

            if (improvement > gLocalSearchOpStats[op].bestImprovement) {
                gLocalSearchOpStats[op].bestImprovement = improvement;
            }
        } else {
            ind = std::move(backup);
        }
    }
}
void evaluateFitnessFullPacking(Individual& indiv, const Data& parameters) {
    const double C1 = 5.41;
    const double BIG = 1e12;

    const auto& forcedOutsource = TEST_FORCED_OUTSOURCE;
    const auto& forcedSelfOwned = TEST_FORCED_SELFOWNED;

    double outsourcedCost = 0.0;
    double selfOwnedFuelCost = 0.0;

    Truck selfOwnedTrucks[regionNum + 1];
    vector<Truck> rentedTrucks;

    unordered_map<int, vector<Gene>> regionMap;
    unordered_map<int, bool> isLoadedGlobal;
    vector<int> notLoadedCustomer;

    auto cargoLookup = createCargoLookup(parameters);

    for (int area = 1; area <= regionNum; ++area) {
        selfOwnedTrucks[area].truckId = area;
        selfOwnedTrucks[area].loadedVolume = 0;
        selfOwnedTrucks[area].assignedCargo.clear();
        selfOwnedTrucks[area].route.clear();
    }

    regionMap.clear();
    isLoadedGlobal.clear();

    for (const auto& gene : indiv.chromosome) {
        isLoadedGlobal[gene.customerId] = false;
        regionMap[gene.decodedServiceArea].push_back(gene);
    }

    // ==============================
    // Self-owned loading
    // ==============================
    for (int area = 1; area <= regionNum; ++area) {
        auto itRM = regionMap.find(area);
        if (itRM == regionMap.end()) continue;

        vector<Gene>& cargoList = itRM->second;

        unordered_map<int, vector<Gene>> customerGrouped;
        customerGrouped.reserve(cargoList.size() * 2 + 1);

        for (int idx = (int)cargoList.size() - 1; idx >= 0; --idx) {
            const Gene& g = cargoList[idx];
            customerGrouped[g.customerId].push_back(g);
        }

        vector<int> customerOrder;
        customerOrder.reserve(customerGrouped.size());

        unordered_map<int, int> routePos;
        const auto& baseRoute = parameters.route[area - 1];
        for (int pos = 0; pos < (int)baseRoute.size(); ++pos) {
            int node = baseRoute[pos];
            if (node > 0) {
                routePos[node] = pos;
            }
        }

        for (const auto& kv : customerGrouped) {
            int custId = kv.first;
            if (forcedOutsource.count(custId)) continue;
            customerOrder.push_back(custId);
        }

        sort(customerOrder.begin(), customerOrder.end(),
            [&](int a, int b) {
                bool inA = routePos.count(a);
                bool inB = routePos.count(b);

                if (inA && inB) return routePos[a] < routePos[b];
                if (inA != inB) return inA > inB;
                return a < b;
            });

        Truck& truck = selfOwnedTrucks[area];
        BLPlacement3D loader(truck.length, truck.width, truck.height);
        loader.setCargoLookup(cargoLookup);

        for (int custId : customerOrder) {
            auto& cargoGroup = customerGrouped[custId];

            bool success = false;
            vector<Gene> bestGroup = cargoGroup;

            for (int trial = 0; trial < 5 && !success; ++trial) {
                auto tempGroup = cargoGroup;

                if (trial > 0) {
                    static thread_local std::mt19937 rng(std::random_device{}());
                    std::shuffle(tempGroup.begin(), tempGroup.end(), rng);
                }

                if (loader.tryInsert(tempGroup, 50)) {
                    bestGroup = tempGroup;
                    success = true;
                } else {
                    long long bestOverflow = 0;
                    bool repaired = localSearchRotateAllAndTryInsert(
                        loader,
                        tempGroup,
                        loader.cargoLookup,
                        20,
                        30,
                        bestOverflow
                    );

                    if (repaired) {
                        bestGroup = tempGroup;
                        success = true;
                    }
                }
            }

            if (success) {
                cargoGroup = bestGroup;

                truck.loadedVolume += parameters.totalVolume[custId - 1];
                truck.assignedCargo.insert(
                    truck.assignedCargo.end(),
                    cargoGroup.begin(),
                    cargoGroup.end()
                );

                isLoadedGlobal[custId] = true;

                for (const auto& g : cargoGroup) {
                    for (auto& indivGene : indiv.chromosome) {
                        if (indivGene.customerId == g.customerId &&
                            indivGene.cargoId == g.cargoId) {
                            indivGene.position[0] = g.position[0];
                            indivGene.position[1] = g.position[1];
                            indivGene.position[2] = g.position[2];
                            indivGene.undecodedRotation = g.undecodedRotation;
                            indivGene.decodedRotation = g.decodedRotation;
                        }
                    }
                }
            }
            // self-owned 裝不下，後面交給 rented 補救
        }
    }

    // ==============================
    // Collect not-loaded customers
    // ==============================
    notLoadedCustomer.clear();

    for (const auto& kv : isLoadedGlobal) {
        if (!kv.second) {
            notLoadedCustomer.push_back(kv.first);
        }
    }

    for (int custId : forcedOutsource) {
        if (!isLoadedGlobal[custId]) {
            notLoadedCustomer.push_back(custId);
        }
    }

    sort(notLoadedCustomer.begin(), notLoadedCustomer.end());
    notLoadedCustomer.erase(
        unique(notLoadedCustomer.begin(), notLoadedCustomer.end()),
        notLoadedCustomer.end()
    );

    // ==============================
    // Self-owned fuel cost
    // ==============================
    selfOwnedFuelCost = 0.0;

    for (int area = 1; area <= regionNum; ++area) {
        const auto& baseRoute = parameters.route[area - 1];

        unordered_set<int> servedSet;
        servedSet.reserve(selfOwnedTrucks[area].assignedCargo.size() * 2 + 1);

        for (const auto& g : selfOwnedTrucks[area].assignedCargo) {
            servedSet.insert(g.customerId);
        }

        if (servedSet.empty()) {
            selfOwnedTrucks[area].route.clear();
            continue;
        }

        vector<int> servedSeq;
        servedSeq.reserve(servedSet.size());

        for (int node : baseRoute) {
            if (node > 0 && servedSet.count(node)) {
                servedSeq.push_back(node);
            }
        }

        for (int cid : servedSet) {
            bool found = false;
            for (int node : servedSeq) {
                if (node == cid) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                servedSeq.push_back(cid);
            }
        }

        if (servedSeq.empty()) {
            selfOwnedTrucks[area].route.clear();
            continue;
        }

        selfOwnedTrucks[area].route = servedSeq;

        double distSum = 0.0;
        distSum += parameters.getDistance(0, servedSeq.front());

        for (int k = 0; k + 1 < (int)servedSeq.size(); ++k) {
            distSum += parameters.getDistance(servedSeq[k], servedSeq[k + 1]);
        }

        distSum += parameters.getDistance(servedSeq.back(), 0);

        selfOwnedFuelCost += distSum;
    }

    // ==============================
    // Rented loading
    // ==============================
    size_t cursor = 0;
    int rentedTruckId = 0;
    unordered_set<int> rentedSeen;

    while (cursor < notLoadedCustomer.size()) {
        Truck rentedTruck;
        rentedTruck.truckId = ++rentedTruckId;
        rentedTruck.loadedVolume = 0;
        rentedTruck.assignedCargo.clear();
        rentedTruck.route.clear();

        BLPlacement3D loader(rentedTruck.length, rentedTruck.width, rentedTruck.height);
        loader.setCargoLookup(cargoLookup);

        bool anyLoadedThisTruck = false;

        for (; cursor < notLoadedCustomer.size(); ++cursor) {
            int custId = notLoadedCustomer[cursor];

            if (rentedSeen.count(custId)) continue;
            if (forcedSelfOwned.count(custId)) continue;

            vector<Gene> cargoGroup;
            cargoGroup.reserve(16);

            for (const auto& g : indiv.chromosome) {
                if (g.customerId == custId) {
                    cargoGroup.push_back(g);
                }
            }

            if (cargoGroup.empty()) continue;

            bool canLoad = loader.tryInsert(cargoGroup, 50);

            if (!canLoad) {
                long long bestOverflow = 0;
                canLoad = localSearchRotateAllAndTryInsert(
                    loader,
                    cargoGroup,
                    loader.cargoLookup,
                    20,
                    30,
                    bestOverflow
                );
            }

            if (canLoad) {
                anyLoadedThisTruck = true;
                rentedSeen.insert(custId);

                rentedTruck.loadedVolume += parameters.totalVolume[custId - 1];
                rentedTruck.assignedCargo.insert(
                    rentedTruck.assignedCargo.end(),
                    cargoGroup.begin(),
                    cargoGroup.end()
                );
                rentedTruck.route.push_back(custId);

                double volume = parameters.totalVolume[custId - 1];
                int v = (int)ceil(volume);
                double cost = 100 + 30 * max(0, v - 3);

                outsourcedCost += cost;
            } else {
                indiv.fitness.clear();
                indiv.fitness.push_back(BIG);
                return;
            }
        }

        if (anyLoadedThisTruck) {
            rentedTrucks.push_back(rentedTruck);
        } else {
            break;
        }
    }

    for (int area = 1; area <= regionNum; ++area) {
        indiv.selfOwnedTrucks[area] = selfOwnedTrucks[area];
    }
    indiv.rentedTrucks = rentedTrucks;

    for (int cid : forcedSelfOwned) {
        bool inSelfOwned = false;

        for (int area = 1; area <= regionNum; ++area) {
            for (const auto& g : indiv.selfOwnedTrucks[area].assignedCargo) {
                if (g.customerId == cid) {
                    inSelfOwned = true;
                    break;
                }
            }
            if (inSelfOwned) break;
        }

        if (!inSelfOwned) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return;
        }
    }

    double obj = outsourcedCost + C1 * selfOwnedFuelCost;

    indiv.fitness.clear();
    indiv.fitness.push_back(obj);
}

bool reconstructFullPackingFromBestAssignment(
    Individual& indiv,
    const Data& parameters,
    int selfOwnedTrialsPerArea = 200,
    int rentedTrialsPerCustomer = 100
) {
    const double C1 = 5.41;
    const double BIG = 1e12;

    const auto& forcedOutsource = TEST_FORCED_OUTSOURCE;
    const auto& forcedSelfOwned = TEST_FORCED_SELFOWNED;

    auto cargoLookup = createCargoLookup(parameters);

    double outsourcedCost = 0.0;
    double selfOwnedFuelCost = 0.0;

    Truck newSelfOwnedTrucks[regionNum + 1];
    vector<Truck> newRentedTrucks;

    unordered_map<int, bool> isLoadedGlobal;

    // ========= init =========
    for (int area = 1; area <= regionNum; ++area) {
        newSelfOwnedTrucks[area] = Truck();
        newSelfOwnedTrucks[area].truckId = area;
        newSelfOwnedTrucks[area].loadedVolume = 0.0;
        newSelfOwnedTrucks[area].assignedCargo.clear();
        newSelfOwnedTrucks[area].route.clear();
    }

    for (const auto& g : indiv.chromosome) {
        isLoadedGlobal[g.customerId] = false;
    }

    // 儲存原本 FAST 找到的「固定自有車 assignment」
    vector<vector<int>> targetRoutes(regionNum + 1);
    for (int area = 1; area <= regionNum; ++area) {
        targetRoutes[area] = indiv.selfOwnedTrucks[area].route;

        // 清掉 0 / 去重
        vector<int>& r = targetRoutes[area];
        r.erase(remove(r.begin(), r.end(), 0), r.end());
        sort(r.begin(), r.end());
        r.erase(unique(r.begin(), r.end()), r.end());

        // 用固定路線順序重建一次，避免順序漂移
        r = buildRouteFromFixedOrder(r, parameters, area);
    }

    // 先建立 gene lookup：customer -> cargo genes
    unordered_map<int, vector<Gene>> genesByCustomer;
    genesByCustomer.reserve(indiv.chromosome.size() * 2 + 1);
    for (const auto& g : indiv.chromosome) {
        genesByCustomer[g.customerId].push_back(g);
    }

    // ========= Self-owned reconstruction =========
    for (int area = 1; area <= regionNum; ++area) {
        vector<int> customerOrder = targetRoutes[area];

        // forcedSelfOwned 若被分到這個 area，但 route 裡沒有，補進去
        for (int cid : forcedSelfOwned) {
            bool belongsToArea = false;
            for (const auto& g : indiv.chromosome) {
                if (g.customerId == cid && g.decodedServiceArea == area) {
                    belongsToArea = true;
                    break;
                }
            }
            if (belongsToArea &&
                find(customerOrder.begin(), customerOrder.end(), cid) == customerOrder.end()) {
                customerOrder.push_back(cid);
            }
        }

        customerOrder.erase(
            remove_if(customerOrder.begin(), customerOrder.end(),
                      [&](int cid) { return forcedOutsource.count(cid) > 0; }),
            customerOrder.end()
        );

        if (customerOrder.empty()) {
            newSelfOwnedTrucks[area].route.clear();
            continue;
        }

        // 用 base route 順序整理
        customerOrder = buildRouteFromFixedOrder(customerOrder, parameters, area);

        bool areaSuccess = false;
        Truck bestTruck;
        unordered_map<int, vector<Gene>> bestPackedGenesByCustomer;

        for (int trial = 0; trial < selfOwnedTrialsPerArea && !areaSuccess; ++trial) {
            Truck truck;
            truck = Truck();
            truck.truckId = area;
            truck.loadedVolume = 0.0;
            truck.assignedCargo.clear();
            truck.route.clear();

            BLPlacement3D loader(truck.length, truck.width, truck.height);
            loader.setCargoLookup(cargoLookup);

            bool allLoaded = true;
            unordered_map<int, vector<Gene>> packedGenesByCustomer;

            for (int custId : customerOrder) {
                auto it = genesByCustomer.find(custId);
                if (it == genesByCustomer.end() || it->second.empty()) {
                    allLoaded = false;
                    break;
                }

                vector<Gene> baseGroup = it->second;
                bool success = false;
                vector<Gene> bestGroup = baseGroup;

                for (int inner = 0; inner < 5 && !success; ++inner) {
                    vector<Gene> tempGroup = baseGroup;

                    if (inner > 0) {
                        static thread_local std::mt19937 rng(std::random_device{}());
                        std::shuffle(tempGroup.begin(), tempGroup.end(), rng);
                    }

                    if (loader.tryInsert(tempGroup, 50)) {
                        bestGroup = tempGroup;
                        success = true;
                    } else {
                        long long bestOverflow = 0;
                        bool repaired = localSearchRotateAllAndTryInsert(
                            loader,
                            tempGroup,
                            loader.cargoLookup,
                            20,
                            30,
                            bestOverflow
                        );

                        if (repaired) {
                            bestGroup = tempGroup;
                            success = true;
                        }
                    }
                }

                if (!success) {
                    allLoaded = false;
                    break;
                }

                packedGenesByCustomer[custId] = bestGroup;
                truck.loadedVolume += parameters.totalVolume[custId - 1];
                truck.assignedCargo.insert(
                    truck.assignedCargo.end(),
                    bestGroup.begin(),
                    bestGroup.end()
                );
            }

            if (allLoaded) {
                truck.route = customerOrder;
                bestTruck = truck;
                bestPackedGenesByCustomer = packedGenesByCustomer;
                areaSuccess = true;
            }
        }

        if (!areaSuccess) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return false;
        }

        newSelfOwnedTrucks[area] = bestTruck;

        // 寫回 chromosome positions / rotations
        for (const auto& kv : bestPackedGenesByCustomer) {
            int custId = kv.first;
            isLoadedGlobal[custId] = true;

            for (const auto& packedGene : kv.second) {
                for (auto& indivGene : indiv.chromosome) {
                    if (indivGene.customerId == packedGene.customerId &&
                        indivGene.cargoId == packedGene.cargoId) {
                        indivGene.position[0] = packedGene.position[0];
                        indivGene.position[1] = packedGene.position[1];
                        indivGene.position[2] = packedGene.position[2];
                        indivGene.undecodedRotation = packedGene.undecodedRotation;
                        indivGene.decodedRotation = packedGene.decodedRotation;
                    }
                }
            }
        }
    }

    // ========= collect outsourced customers =========
    vector<int> notLoadedCustomers;
    for (const auto& kv : isLoadedGlobal) {
        if (!kv.second) {
            notLoadedCustomers.push_back(kv.first);
        }
    }

    for (int cid : forcedOutsource) {
        if (!isLoadedGlobal[cid]) {
            notLoadedCustomers.push_back(cid);
        }
    }

    sort(notLoadedCustomers.begin(), notLoadedCustomers.end());
    notLoadedCustomers.erase(
        unique(notLoadedCustomers.begin(), notLoadedCustomers.end()),
        notLoadedCustomers.end()
    );

    // ========= Rented reconstruction =========
    // 簡化成「每個 outsourced customer 一台 rented truck」
    // 因為你現在的 outsourced cost 是按 customer 算，和合併到同台 rented truck 無關
    int rentedTruckId = 0;

    for (int custId : notLoadedCustomers) {
        if (forcedSelfOwned.count(custId)) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return false;
        }

        auto it = genesByCustomer.find(custId);
        if (it == genesByCustomer.end() || it->second.empty()) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return false;
        }

        bool success = false;
        Truck rentedTruck;
        unordered_map<int, vector<Gene>> packedGenesByCustomer;

        for (int trial = 0; trial < rentedTrialsPerCustomer && !success; ++trial) {
            Truck trialTruck;
            trialTruck = Truck();
            trialTruck.truckId = ++rentedTruckId;
            trialTruck.loadedVolume = 0.0;
            trialTruck.assignedCargo.clear();
            trialTruck.route.clear();

            BLPlacement3D loader(trialTruck.length, trialTruck.width, trialTruck.height);
            loader.setCargoLookup(cargoLookup);

            vector<Gene> baseGroup = it->second;
            vector<Gene> bestGroup = baseGroup;
            bool packed = false;

            for (int inner = 0; inner < 5 && !packed; ++inner) {
                vector<Gene> tempGroup = baseGroup;

                if (inner > 0) {
                    static thread_local std::mt19937 rng(std::random_device{}());
                    std::shuffle(tempGroup.begin(), tempGroup.end(), rng);
                }

                if (loader.tryInsert(tempGroup, 50)) {
                    bestGroup = tempGroup;
                    packed = true;
                } else {
                    long long bestOverflow = 0;
                    bool repaired = localSearchRotateAllAndTryInsert(
                        loader,
                        tempGroup,
                        loader.cargoLookup,
                        20,
                        30,
                        bestOverflow
                    );

                    if (repaired) {
                        bestGroup = tempGroup;
                        packed = true;
                    }
                }
            }

            if (packed) {
                trialTruck.loadedVolume += parameters.totalVolume[custId - 1];
                trialTruck.assignedCargo = bestGroup;
                trialTruck.route.push_back(custId);

                rentedTruck = trialTruck;
                packedGenesByCustomer[custId] = bestGroup;
                success = true;
            }
        }

        if (!success) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return false;
        }

        newRentedTrucks.push_back(rentedTruck);

        for (const auto& packedGene : packedGenesByCustomer[custId]) {
            for (auto& indivGene : indiv.chromosome) {
                if (indivGene.customerId == packedGene.customerId &&
                    indivGene.cargoId == packedGene.cargoId) {
                    indivGene.position[0] = packedGene.position[0];
                    indivGene.position[1] = packedGene.position[1];
                    indivGene.position[2] = packedGene.position[2];
                    indivGene.undecodedRotation = packedGene.undecodedRotation;
                    indivGene.decodedRotation = packedGene.decodedRotation;
                }
            }
        }

        double volume = parameters.totalVolume[custId - 1];
        int v = (int)ceil(volume);
        outsourcedCost += 100 + 30 * max(0, v - 3);
    }

    // ========= self-owned fuel cost =========
    selfOwnedFuelCost = 0.0;
    for (int area = 1; area <= regionNum; ++area) {
        const auto& route = newSelfOwnedTrucks[area].route;
        if (route.empty()) continue;

        double dist = 0.0;
        dist += parameters.getDistance(0, route.front());

        for (int i = 0; i + 1 < (int)route.size(); ++i) {
            dist += parameters.getDistance(route[i], route[i + 1]);
        }

        dist += parameters.getDistance(route.back(), 0);
        selfOwnedFuelCost += dist;
    }

    // ========= forced self-owned check =========
    for (int cid : forcedSelfOwned) {
        bool inSelfOwned = false;
        for (int area = 1; area <= regionNum; ++area) {
            for (const auto& g : newSelfOwnedTrucks[area].assignedCargo) {
                if (g.customerId == cid) {
                    inSelfOwned = true;
                    break;
                }
            }
            if (inSelfOwned) break;
        }

        if (!inSelfOwned) {
            indiv.fitness.clear();
            indiv.fitness.push_back(BIG);
            return false;
        }
    }

    // ========= final write back =========
    for (int area = 1; area <= regionNum; ++area) {
        indiv.selfOwnedTrucks[area] = newSelfOwnedTrucks[area];
    }
    indiv.rentedTrucks = newRentedTrucks;

    double obj = outsourcedCost + C1 * selfOwnedFuelCost;
    indiv.fitness.clear();
    indiv.fitness.push_back(obj);

    return true;
}