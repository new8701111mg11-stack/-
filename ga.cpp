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

// ===== 核心 Move：嘗試把 cid 從 rented -> 自有車 area =====
static bool tryMoveRentedCustomerToSelfTruck(Individual& ind,
                                            int cid, int area,
                                            const Data& parameters,
                                            const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    // area 必須可行（服務區限制）
    if (parameters.serviceRegion[cid - 1][area - 1] != 1) return false;

    // 取出整包貨
    vector<Gene> pack = collectCustomerGenesFromRented(ind, cid);
    if (pack.empty()) return false;

    Truck& selfT = ind.selfOwnedTrucks[area];

    // 暫存，失敗要回復
    auto backupSelf = selfT.assignedCargo;
    auto backupRented = ind.rentedTrucks;

    const auto originalCargo = selfT.assignedCargo;
    const int n = (int)originalCargo.size();

    // 嘗試不同插入位置，而不是只 append 到尾端
    for (int pos = 0; pos <= n; ++pos) {
        selfT.assignedCargo = originalCargo;

        selfT.assignedCargo.insert(
            selfT.assignedCargo.begin() + pos,
            pack.begin(),
            pack.end()
        );

        if (repackTruck(selfT, cargoLookup, /*maxTriesInsert=*/80)) {
            removeCustomerFromRented(ind, cid);
            return true;
        }
    }

    // 全部插入位置都失敗才 rollback
    selfT.assignedCargo = std::move(backupSelf);
    ind.rentedTrucks = std::move(backupRented);
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

static bool doReinsertionOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup
){
    unordered_set<int> omega;
    for (const auto& rt : ind.rentedTrucks)
        for (const auto& g : rt.assignedCargo)
            omega.insert(g.customerId);

    for (int cid : omega) {
        for (int area = 1; area <= regionNum; ++area) {
            if (tryMoveRentedCustomerToSelfTruck(ind, cid, area, parameters, cargoLookup)) {
                return true;
            }
        }
    }
    return false;
}


static bool doExchangeOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup
){
    // 1) 找 rented 客戶集合
    unordered_set<int> omegaSet;
    for (const auto& rt : ind.rentedTrucks)
        for (const auto& g : rt.assignedCargo)
            omegaSet.insert(g.customerId);

    if (omegaSet.empty()) return false;

    vector<int> omega(omegaSet.begin(), omegaSet.end());

    // 👉 關鍵改動：全部試，不只第一個
    for (int cidR : omega) {

        vector<Gene> packR = collectCustomerGenesFromRented(ind, cidR);
        if (packR.empty()) continue;

        // 2) 找 self truck + cidS
        for (int area = 1; area <= regionNum; ++area) {
            Truck& selfT = ind.selfOwnedTrucks[area];

            unordered_set<int> selfCustomers;
            for (auto& g : selfT.assignedCargo)
                selfCustomers.insert(g.customerId);

            if (selfCustomers.empty()) continue;

            for (int cidS : selfCustomers) {

                if (parameters.serviceRegion[cidR - 1][area - 1] != 1) continue;

                Individual backup = ind;

                vector<Gene> packS = collectCustomerGenesFromSelf(selfT, cidS);
                if (packS.empty()) { ind = std::move(backup); continue; }

                // self: 移除 S，加入 R
                removeCustomerFromSelf(selfT, cidS);
                selfT.assignedCargo.insert(
                    selfT.assignedCargo.end(),
                    packR.begin(),
                    packR.end()
                );

                // rented: 移除 R
                removeCustomerFromRented(ind, cidR);

                // 👉 嘗試放到不同 rented truck
                bool success = false;

                for (int rtIdx = 0; rtIdx < (int)ind.rentedTrucks.size(); ++rtIdx) {
                    ind.rentedTrucks[rtIdx].assignedCargo.insert(
                        ind.rentedTrucks[rtIdx].assignedCargo.end(),
                        packS.begin(),
                        packS.end()
                    );

                    bool ok1 = repackTruck(ind.selfOwnedTrucks[area], cargoLookup, 80);
                    bool ok2 = repackTruck(ind.rentedTrucks[rtIdx], cargoLookup, 100);

                    if (ok1 && ok2) {
                        success = true;
                        break;
                    }

                    // rollback 該台 rented truck
                    ind = backup;
                }

                if (success) {
                    return true;
                } else {
                    ind = std::move(backup);
                }
            }
        }
    }
    return false;
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
    unordered_map<int, vector<Gene>> regionMap;
    unordered_map<int, bool> isLoadedGlobal;

    auto cargoLookup = createCargoLookup(parameters);

    long long dbgHit = 0;
    long long dbgReal = 0;

    // ==============================
    // init
    // ==============================
    for (int area = 1; area <= regionNum; ++area) {
        selfOwnedTrucks[area] = Truck();
        selfOwnedTrucks[area].truckId = area;
    }

    // ==============================
    // build regionMap
    // ==============================
    for (const auto& g : indiv.chromosome) {
        regionMap[g.decodedServiceArea].push_back(g);
        isLoadedGlobal[g.customerId] = false;
    }

    // ==============================
    // self-owned loading
    // ==============================
    for (int area = 1; area <= regionNum; ++area) {

        auto it = regionMap.find(area);
        if (it == regionMap.end()) continue;

        auto& cargoList = it->second;

        unordered_map<int, vector<Gene>> customerGrouped;
        for (auto& g : cargoList) {
            customerGrouped[g.customerId].push_back(g);
        }

        vector<int> customerSet;
        customerSet.reserve(customerGrouped.size());
        for (auto& kv : customerGrouped) {
            customerSet.push_back(kv.first);
        }

        sort(customerSet.begin(), customerSet.end());
        customerSet.erase(unique(customerSet.begin(), customerSet.end()), customerSet.end());

        string key = buildSortedCustomerKey(customerSet);

        // exact cache hit: 整組都成功過
        if (loadingCache.selfOwnedSuccessCache[area].count(key)) {
            dbgHit++;

            for (int cid : customerSet) {
                isLoadedGlobal[cid] = true;
                selfOwnedTrucks[area].loadedVolume += parameters.totalVolume[cid - 1];
            }

            selfOwnedTrucks[area].route =
                buildRouteFromFixedOrder(customerSet, parameters, area);

            continue;
        }

        dbgReal++;

        Truck& truck = selfOwnedTrucks[area];
        BLPlacement3D loader(truck.length, truck.width, truck.height);
        loader.setCargoLookup(cargoLookup);

        bool allSuccess = true;
        vector<int> loadedCustomers;

        // 逐 customer 試，fail 不 break
        for (int custId : customerSet) {

            auto cargoGroup = customerGrouped[custId];
            bool success = false;

            for (int trial = 0; trial < 5 && !success; ++trial) {
                auto temp = cargoGroup;

                if (trial > 0) {
                    static thread_local std::mt19937 rng(std::random_device{}());
                    shuffle(temp.begin(), temp.end(), rng);
                }

                if (loader.tryInsert(temp, 50)) {
                    success = true;
                }
            }

            if (success) {
                loadedCustomers.push_back(custId);
                isLoadedGlobal[custId] = true;
                truck.loadedVolume += parameters.totalVolume[custId - 1];
            } else {
                allSuccess = false;
            }
        }

        sort(loadedCustomers.begin(), loadedCustomers.end());
        loadedCustomers.erase(unique(loadedCustomers.begin(), loadedCustomers.end()), loadedCustomers.end());

        if (!loadedCustomers.empty()) {
            truck.route = buildRouteFromFixedOrder(loadedCustomers, parameters, area);
        } else {
            truck.route.clear();
        }

        // 只有整組都成功才 cache
        if (allSuccess) {
            loadingCache.selfOwnedSuccessCache[area].insert(key);
        }
    }

    // ==============================
    // collect not loaded
    // ==============================
    vector<int> notLoaded;
    for (auto& kv : isLoadedGlobal) {
        if (!kv.second) notLoaded.push_back(kv.first);
    }

    // ==============================
    // self-owned fuel cost
    // ==============================
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

    // ==============================
    // rented cost
    // ==============================
    for (int cid : notLoaded) {
        double vol = parameters.totalVolume[cid - 1];
        int v = (int)ceil(vol);
        outsourcedCost += 100 + 30 * max(0, v - 3);
    }

    // ==============================
    // final
    // ==============================
    for (int area = 1; area <= regionNum; ++area) {
        indiv.selfOwnedTrucks[area] = selfOwnedTrucks[area];
    }

    double obj = outsourcedCost + C1 * selfOwnedFuelCost;
    indiv.fitness.clear();
    indiv.fitness.push_back(obj);

    cout << "[FAST] hit=" << dbgHit << " real=" << dbgReal << endl;
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
static bool doBlockSwapSampled(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int samplePairs = 8)
{
    const double EPS = 1e-9;

    // base 解先確定是最新的
    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);

    vector<int> order = getCustomerOrderFromChromosome(base);
    int n = (int)order.size();
    if (n < 2) return false;

    // 先產生要測的 pair，避免在 parallel 區塊裡共用 rng
    vector<pair<int, int>> sampledPairs;
    sampledPairs.reserve(samplePairs);

    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pick(0, n - 1);

    for (int s = 0; s < samplePairs; ++s) {
        int i = pick(rng);
        int j = pick(rng);
        if (i == j) {
            --s;
            continue;
        }
        sampledPairs.push_back({i, j});
    }

    // 每個 sample 一個結果槽
    vector<Individual> candidates(samplePairs);
    vector<double> objs(samplePairs, 1e18);
    vector<char> valid(samplePairs, 0);

#pragma omp parallel for schedule(dynamic)
    for (int s = 0; s < samplePairs; ++s) {
        int i = sampledPairs[s].first;
        int j = sampledPairs[s].second;

        Individual cand = ind;  // thread-local copy
        swapCustomerBlocksInChromosome(cand, order[i], order[j]);
        reDecodeAndEvaluate(cand, parameters, cargoLookup);

        if (!cand.fitness.empty()) {
            candidates[s] = std::move(cand);
            objs[s] = candidates[s].fitness[0];
            valid[s] = 1;
        }
    }

    // 單執行緒選 best
    int bestIdx = -1;
    double bestObj = base.fitness[0];

    for (int s = 0; s < samplePairs; ++s) {
        if (valid[s] && objs[s] < bestObj - EPS) {
            bestObj = objs[s];
            bestIdx = s;
        }
    }

    if (bestIdx != -1) {
        ind = std::move(candidates[bestIdx]);
        return true;
    }

    return false;
}
static bool doBlockInsertSampled(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int sampleMoves = 10)
{
    const double EPS = 1e-9;

    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);

    vector<int> order = getCustomerOrderFromChromosome(base);
    int n = (int)order.size();
    if (n < 2) return false;

    vector<pair<int, int>> sampledMoves;
    sampledMoves.reserve(sampleMoves);

    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pick(0, n - 1);

    for (int s = 0; s < sampleMoves; ++s) {
        int i = pick(rng);
        int j = pick(rng);
        if (i == j) {
            --s;
            continue;
        }
        sampledMoves.push_back({i, j});
    }

    vector<Individual> candidates(sampleMoves);
    vector<double> objs(sampleMoves, 1e18);
    vector<char> valid(sampleMoves, 0);

#pragma omp parallel for schedule(dynamic)
    for (int s = 0; s < sampleMoves; ++s) {
        int i = sampledMoves[s].first;
        int j = sampledMoves[s].second;

        Individual cand = ind;  // thread-local copy
        moveCustomerBlockInChromosome(cand, order[i], order[j]);
        reDecodeAndEvaluate(cand, parameters, cargoLookup);

        if (!cand.fitness.empty()) {
            candidates[s] = std::move(cand);
            objs[s] = candidates[s].fitness[0];
            valid[s] = 1;
        }
    }

    int bestIdx = -1;
    double bestObj = base.fitness[0];

    for (int s = 0; s < sampleMoves; ++s) {
        if (valid[s] && objs[s] < bestObj - EPS) {
            bestObj = objs[s];
            bestIdx = s;
        }
    }

    if (bestIdx != -1) {
        ind = std::move(candidates[bestIdx]);
        return true;
    }

    return false;
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
    sort(omega.begin(), omega.end());  // 穩定順序，至少先排除 unordered_set 的亂序影響

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

static bool doForceSelfAttempt(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup)
{
    // 找目前所有 rented customers
    unordered_set<int> rentedSet;
    for (const auto& rt : ind.rentedTrucks) {
        for (const auto& g : rt.assignedCargo) {
            rentedSet.insert(g.customerId);
        }
    }

    if (rentedSet.empty()) return false;

    // baseline objective
    Individual base = ind;
    reDecodeAndEvaluate(base, parameters, cargoLookup);
    if (base.fitness.empty()) return false;

    const double baseObj = base.fitness[0];
    const double EPS = 1e-9;

    vector<int> rentedCustomers(rentedSet.begin(), rentedSet.end());
    sort(rentedCustomers.begin(), rentedCustomers.end());

    for (int cid : rentedCustomers) {

        // 從 rented 取出整包貨
        vector<Gene> pack = collectCustomerGenesFromRented(ind, cid);
        if (pack.empty()) continue;

        for (int area = 1; area <= regionNum; ++area) {

            // 服務區限制
            if (parameters.serviceRegion[cid - 1][area - 1] != 1) continue;

            Individual backup = ind;

            // 先從 rented 移掉這個 customer
            removeCustomerFromRented(ind, cid);

            // 加到對應 self-owned truck
            Truck& selfT = ind.selfOwnedTrucks[area];
            const auto originalCargo = selfT.assignedCargo;
            const int n = (int)originalCargo.size();

            bool accepted = false;

            // 嘗試不同插入位置，不只放尾端
            for (int pos = 0; pos <= n; ++pos) {
                selfT.assignedCargo = originalCargo;

                selfT.assignedCargo.insert(
                    selfT.assignedCargo.begin() + pos,
                    pack.begin(),
                    pack.end()
                );

                // 先看這台車能不能重排成功
                if (!repackTruck(selfT, cargoLookup, 80)) {
                    continue;
                }

                // 整體重新解碼與評估
                reDecodeAndEvaluate(ind, parameters, cargoLookup);

                if (!ind.fitness.empty() && ind.fitness[0] < baseObj - EPS) {
                    accepted = true;
                    return true;
                }
            }

            // 全部位置都不行，rollback
            if (!accepted) {
                ind = std::move(backup);
            }
        }
    }

    return false;
}
static bool doConsolidationOnce(
    Individual& ind,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup
){
    for (int i = 0; i < (int)ind.rentedTrucks.size(); ++i) {
        for (int j = i + 1; j < (int)ind.rentedTrucks.size(); ++j) {
            auto backup = ind.rentedTrucks;

            auto& A = ind.rentedTrucks[i];
            auto& B = ind.rentedTrucks[j];
            A.assignedCargo.insert(A.assignedCargo.end(), B.assignedCargo.begin(), B.assignedCargo.end());

            if (repackTruck(A, cargoLookup, 100)) {
                ind.rentedTrucks.erase(ind.rentedTrucks.begin() + j);
                return true;
            } else {
                ind.rentedTrucks = std::move(backup);
            }
        }
    }
    return false;
}

static void localSearchImproveByRealObjective(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int budgetMoves = 100)
{
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pickOp(0, 4); // 0..4 共五種操作
    const double EPS = 1e-9;

    reDecodeAndEvaluate(ind, parameters, cargoLookup);
    double bestObj = ind.fitness[0];

    for (int iter = 0; iter < budgetMoves; ++iter) {
        int op = pickOp(rng);

        Individual backup = ind;
        bool moved = false;

        switch (op) {
        case 0:
            moved = doReinsertionOnceStable(ind, parameters, cargoLookup);
            break;
        case 1:
            moved = doForceSelfAttempt(ind, parameters, cargoLookup);
            break;
        case 2:
            moved = doExchangeOnce(ind, parameters, cargoLookup);
            break;
        case 3:
            moved = doBlockSwapSampled(ind, parameters, cargoLookup);
            break;
        case 4:
            moved = doBlockInsertSampled(ind, parameters, cargoLookup);
            break;
        default:
            break;
        }

        if (!moved) {
            ind = std::move(backup);
            continue;
        }

        reDecodeAndEvaluate(ind, parameters, cargoLookup);

        if (ind.fitness[0] < bestObj - EPS) {
            bestObj = ind.fitness[0];
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