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

using namespace std;

enum ReduceOp {
    opReinsertion = 0,   // rented -> self
    opConsolidation = 1, // merge rented trucks
    opExchange = 2       // swap rented<->self customer packs
};

// ===== Helper: 把一台 Truck 內的 assignedCargo 全部重排（repack） =====
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

// ===== Helper: 從 ind.rentedTrucks 裡收集某個客戶 cid 的所有 Gene（整包） =====
static vector<Gene> collectCustomerGenesFromRented(const Individual& ind, int cid)
{
    vector<Gene> out;
    for (const auto& rt : ind.rentedTrucks) {
        for (const auto& g : rt.assignedCargo) {
            if (g.customerId == cid) out.push_back(g);
        }
    }
    return out;
}

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

    // 加到自有車（先加再 repack）
    selfT.assignedCargo.insert(selfT.assignedCargo.end(), pack.begin(), pack.end());

    // repack 自有車
    if (!repackTruck(selfT, cargoLookup, /*maxTriesInsert=*/80)) {
        // rollback
        selfT.assignedCargo = std::move(backupSelf);
        ind.rentedTrucks    = std::move(backupRented);
        return false;
    }

    // 成功：從 rented 移除 cid
    removeCustomerFromRented(ind, cid);

    return true;
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

static bool doExchangeOnce(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup
){
    // 1) 找 rented 客戶集合
    unordered_set<int> omega;
    for (const auto& rt : ind.rentedTrucks)
        for (const auto& g : rt.assignedCargo)
            omega.insert(g.customerId);
    if (omega.empty()) return false;

    // 2) 找一個 rented 客戶（你也可以改成隨機）
    int cidR = *omega.begin();
    vector<Gene> packR = collectCustomerGenesFromRented(ind, cidR);
    if (packR.empty()) return false;

    // 3) 找一台自有車 + 一個 self 客戶 cidS
    for (int area = 1; area <= regionNum; ++area) {
        Truck& selfT = ind.selfOwnedTrucks[area];

        // selfT 裡有哪些客戶
        unordered_set<int> selfCustomers;
        for (auto& g : selfT.assignedCargo) selfCustomers.insert(g.customerId);
        if (selfCustomers.empty()) continue;

        for (int cidS : selfCustomers) {
            // cidR 必須允許放到 area（服務限制）
            if (parameters.serviceRegion[cidR - 1][area - 1] != 1) continue;

            // ===== 交換開始：先備份 =====
            Individual backup = ind;

            // packS
            vector<Gene> packS = collectCustomerGenesFromSelf(selfT, cidS);
            if (packS.empty()) { ind = std::move(backup); continue; }

            // (a) self 車移除 cidS，加入 cidR
            removeCustomerFromSelf(selfT, cidS);
            selfT.assignedCargo.insert(selfT.assignedCargo.end(), packR.begin(), packR.end());

            // (b) rented 移除 cidR，加入 cidS（丟回 rented：先丟到「第一台有 cidR 的 rented truck」）
            // 你也可以更聰明：找原本裝 cidR 的那台 rented truck
            removeCustomerFromRented(ind, cidR);
            if (ind.rentedTrucks.empty()) {
                // 沒有 rented truck 就代表 cidR 被移掉了，但我們還要把 cidS 丟去哪？=> rollback
                ind = std::move(backup);
                continue;
            }
            ind.rentedTrucks[0].assignedCargo.insert(ind.rentedTrucks[0].assignedCargo.end(), packS.begin(), packS.end());

            // (c) repack 兩台車（self + rented[0]）
            bool ok1 = repackTruck(ind.selfOwnedTrucks[area], cargoLookup, 80);
            bool ok2 = repackTruck(ind.rentedTrucks[0], cargoLookup, 100);

            if (ok1 && ok2) {
                // 成功一次 exchange
                return true;
            } else {
                ind = std::move(backup);
            }
        }
    }
    return false;
}

static int rentedTruckCount(const Individual& ind) {
    return (int)ind.rentedTrucks.size();
}

static void localSearchReduceRentedBudget(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int budgetMoves = 1800   // 老師說的總次數（例如 90）
){
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> pickOp(0, 2); // 0..2 => 1/3 each

    // 你可以用「租用車數」當主要改善判斷（老師要看 Obj1）
    int bestRentedCnt = rentedTruckCount(ind);

    for (int iter = 0; iter < budgetMoves; ++iter) {
        int op = pickOp(rng);

        // 備份：失敗就 rollback
        Individual backup = ind;

        bool moved = false;

        if (op == opReinsertion) {
            // 你原本 LS-A 的內容：tryMoveRentedCustomerToSelfTruck(...)
            moved = doReinsertionOnce(ind, parameters, cargoLookup);   // 你等下照我第3段做
        } 
        else if (op == opConsolidation) {
            // 你原本 LS-B：merge rented trucks
            moved = doConsolidationOnce(ind, cargoLookup);            // 你等下照我第3段做
        } 
        else if (op == opExchange) {
            // 新增：交換 rented 客戶 與 self 客戶
            moved = doExchangeOnce(ind, parameters, cargoLookup);     // 你等下照我第4段做
        }

        // 接受準則：租用車數下降才算改善（最直觀、也符合老師主軸）
        int nowCnt = rentedTruckCount(ind);
        if (moved && nowCnt <= bestRentedCnt) {
            bestRentedCnt = nowCnt;
            break;
        } else {
            ind = std::move(backup); // rollback
        }

        // 已經 0 台租用車就提前停
        if (bestRentedCnt == 0) break;
    }
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

void evaluateFitness(Individual &indiv, const Data &parameters) {
    indiv.fitness.clear();
    long long vechicleLoadedMaxGap = 0;
    long long rentedVehicleCargoCost = 0;

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
            if (loader.tryInsert(cargoGroup,50)) {
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
                continue;
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

            bool canLoad = loader.tryInsert(cargoGroup,50);
            //試錯
            if (!canLoad) {
                long long bestOverflow = 0;
                // ★ 老師要的：cannot be loaded -> 整包換方向做 local search
                canLoad = localSearchRotateAllAndTryInsert(
                loader,
                cargoGroup,
                loader.cargoLookup,
                /*maxOuterTries=*/20,
                /*maxInsertTries=*/30,
                bestOverflow
            );

    // 你可以順便印 overflow 讓你知道是不是「體積真的爆了」
    /*if (!canLoad) {
        cerr << "[LS FAIL] cust " << custId
             << " bestOverflowVolume=" << bestOverflow << "\n";
    }*/
}
            if (canLoad) {
                anyLoadedThisTruck = true;
                rentedSeen.insert(custId);

                rentedTruck.assignedCargo.insert(rentedTruck.assignedCargo.end(),
                                                cargoGroup.begin(), cargoGroup.end());
                rentedTruck.route.push_back(custId);

                for (const auto& g : cargoGroup) {
                    const Cargo& c = loader.cargoLookup[g.customerId][g.cargoId];
                    int chargeUnits = c.volume;
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
    indiv.rentedVehicleLoadingCost = (double)rentedVehicleCargoCost;

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
// ===================== 你已有/你之前用過的 helper（若已存在可略過） =====================

static bool evalCandidate_NoWorseObj1(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    long long baselineF1
){
    vector<Individual> tmp{ ind };
    decodePopulation(tmp, parameters, cargoLookup);
    ind = tmp[0];

    ind.fitness.clear();
    evaluateFitness(ind, parameters);

    if (ind.fitness.size() < 2) return false;

    // Obj1 不變差（允許變更好）
    if (ind.fitness[1] > baselineF1) return false;

    return true;
}


static inline bool canServe(const Data& parameters, int customerId, int area) {
    return parameters.serviceRegion[customerId - 1][area - 1] == 1;
}

static unordered_set<int> collectOmegaRentedCustomers(const Individual& ind) {
    unordered_set<int> omega;
    for (const auto& rt : ind.rentedTrucks)
        for (const auto& g : rt.assignedCargo)
            omega.insert(g.customerId);
    return omega;
}

// 每客戶一次，排除Ω
static void buildCustomerAreaMap_SelfOnly(
    const Individual& ind,
    const unordered_set<int>& omega,
    unordered_map<int,int>& custArea,
    vector<vector<int>>& customersInArea
){
    custArea.clear();
    customersInArea.assign(regionNum + 1, {});
    for (const auto& g : ind.chromosome) {
        int cid = g.customerId;
        if (omega.count(cid)) continue;
        if (custArea.find(cid) != custArea.end()) continue;
        int a = g.routeArea;
        if (a < 1 || a > regionNum) continue;
        custArea[cid] = a;
        customersInArea[a].push_back(cid);
    }
}

static void applyMoveCustomerArea(Individual& ind, int cid, int newArea) {
    for (auto& g : ind.chromosome)
        if (g.customerId == cid)
            g.routeArea = newArea;
}

static void applySwapCustomersArea(Individual& ind, int cidA, int areaA, int cidB, int areaB) {
    for (auto& g : ind.chromosome) {
        if (g.customerId == cidA) g.routeArea = areaB;
        else if (g.customerId == cidB) g.routeArea = areaA;
    }
}

// ===================== Teacher-style Obj2 比較（字典序） =====================
struct Obj2State {
    long long maxLoad;
    long long minLoad;
    long long gap() const { return maxLoad - minLoad; }
};

static Obj2State getObj2State(const Individual& ind) {
    long long mx = 0;
    long long mn = LLONG_MAX;
    for (int a = 1; a <= regionNum; ++a) {
        mx = max(mx, ind.selfOwnedTrucks[a].loadedVolume);
        mn = min(mn, ind.selfOwnedTrucks[a].loadedVolume);
    }
    return {mx, mn};
}

static bool betterObj2Lex(const Obj2State& cur, const Obj2State& base) {
    if (cur.maxLoad != base.maxLoad) return cur.maxLoad < base.maxLoad;
    if (cur.minLoad != base.minLoad) return cur.minLoad > base.minLoad;
    return cur.gap() < base.gap();
}


static bool equalObj2(const Obj2State& a, const Obj2State& b) {
    return a.maxLoad == b.maxLoad && a.minLoad == b.minLoad && a.gap() == b.gap();
}

// cur 不比 base 差（允許平手 / 或在 lex 序上不後退）
static bool notWorseObj2Lex(const Obj2State& cur, const Obj2State& base) {
    // 只要 base 不比 cur 好，就代表 cur >= base（lex）
    return !betterObj2Lex(base, cur);
}
// maxLoad 越小越好；同 maxLoad 時 minLoad 越大越好；再同看 gap




// ===================== OSR 候選：serviceRegion 交集 + 相鄰區域 =====================
// 從 areaFrom 移到 areaTo 的 OSR 客戶：目前在 areaFrom，且同時可由 areaTo 服務（也可由 areaFrom 服務通常成立）
static vector<int> getOSRCustomers_FromTo(
    const vector<vector<int>>& customersInArea,
    const Data& parameters,
    const unordered_set<int>& omega,
    int areaFrom, int areaTo
){
    vector<int> res;
    for (int cid : customersInArea[areaFrom]) {
        if (omega.count(cid)) continue;
        if (canServe(parameters, cid, areaFrom) && canServe(parameters, cid, areaTo)) {
            res.push_back(cid);
        }
    }
    return res;
}

// 用 totalVolume 當作啟發式排序（大件優先）
static void sortByVolumeDesc(vector<int>& cands, const Data& parameters) {
    sort(cands.begin(), cands.end(), [&](int x, int y){
        return parameters.totalVolume[x - 1] > parameters.totalVolume[y - 1];
    });
}

// ===================== 核心：對指定 area 嘗試 Reduce-Max（先 insertion 再 exchange） =====================
static bool decodeEvalKeepObj1Baseline_BestOfN(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    long long obj1Baseline,
    int N
){
    Individual best = ind;
    bool hasBest = false;

    long long bestF1 = LLONG_MAX;
    Obj2State bestObj2{LLONG_MAX, 0};

    for (int t = 0; t < N; ++t) {
        Individual cand = ind;

        vector<Individual> tmp{ cand };
        decodePopulation(tmp, parameters, cargoLookup);
        cand = tmp[0];

        cand.fitness.clear();
        evaluateFitness(cand, parameters);

        if (cand.fitness.size() < 2) continue;
        if (cand.fitness[1] > obj1Baseline) continue;  // Obj1 不變差

        Obj2State st = getObj2State(cand);

        // 先比 Obj1，再比 Obj2Lex（你也可以只挑 Obj2 最好）
        if (!hasBest ||
            cand.fitness[1] < bestF1 ||
            (cand.fitness[1] == bestF1 && betterObj2Lex(st, bestObj2)))
        {
            hasBest = true;
            best = cand;
            bestF1 = cand.fitness[1];
            bestObj2 = st;
        }
    }

    if (!hasBest) return false;

    ind = best;
    return true;
}


static bool tryReduceMaxOnArea_KeepObj1Baseline(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    long long obj1Baseline,
    int areaA,
    int repackRestarts
){
    Obj2State base = getObj2State(ind);

    vector<int> neighbors;
    if (areaA - 1 >= 1) neighbors.push_back(areaA - 1);
    if (areaA + 1 <= regionNum) neighbors.push_back(areaA + 1);

    // 這裡 omegaBaseline 已移除：改成 self-only map（你原本的 build 需要 omega，我用空集合表示「沒有人禁止動」，
    // 但真正禁止動的是「不能讓 Obj1 變大」，由 obj1Baseline 保護。
    unordered_set<int> omegaEmpty;
    unordered_map<int,int> custArea;
    vector<vector<int>> customersInArea;
    buildCustomerAreaMap_SelfOnly(ind, omegaEmpty, custArea, customersInArea);

    // A1) Insertion
    for (int areaB : neighbors) {
        auto Sa = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaA, areaB);
        if (Sa.empty()) continue;
        sortByVolumeDesc(Sa, parameters);

        for (int cid : Sa) {
            Individual backup = ind;

            for (int r = 0; r < repackRestarts; ++r) {
                ind = backup;
                applyMoveCustomerArea(ind, cid, areaB);

                if (!decodeEvalKeepObj1Baseline_BestOfN(ind, parameters, cargoLookup, obj1Baseline, /*N=*/10)
) {
                    continue;
                }

                Obj2State cur = getObj2State(ind);
                if (betterObj2Lex(cur, base)) return true;
            }
            ind = backup;
        }
    }

    // A2) Pair-exchange
    for (int areaB : neighbors) {
        auto Sa = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaA, areaB);
        auto Sb = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaB, areaA);
        if (Sa.empty() || Sb.empty()) continue;

        struct Cand { int i, j; long long improveA; };
        vector<Cand> cands;
        long long loadA = ind.selfOwnedTrucks[areaA].loadedVolume;

        for (int i : Sa) {
            long long vi = parameters.totalVolume[i - 1];
            for (int j : Sb) {
                long long vj = parameters.totalVolume[j - 1];
                long long newA = loadA - vi + vj;
                long long improve = loadA - newA;
                if (improve > 0) cands.push_back({i, j, improve});
            }
        }
        if (cands.empty()) continue;

        sort(cands.begin(), cands.end(), [&](const Cand& x, const Cand& y){
            return x.improveA > y.improveA;
        });

        for (const auto& pc : cands) {
            Individual backup = ind;

            for (int r = 0; r < repackRestarts; ++r) {
                ind = backup;
                applySwapCustomersArea(ind, pc.i, areaA, pc.j, areaB);

                if (!decodeEvalKeepObj1Baseline_BestOfN(ind, parameters, cargoLookup, obj1Baseline, /*N=*/10)
) {
                    continue;
                }

                Obj2State cur = getObj2State(ind);
                if (betterObj2Lex(cur, base)) return true;
            }
            ind = backup;
        }
    }

    return false;
}

static bool tryIncreaseMinOnArea_KeepObj1Baseline(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    long long obj1Baseline,
    int areaA,
    int repackRestarts
){
    Obj2State base = getObj2State(ind);

    vector<int> neighbors;
    if (areaA - 1 >= 1) neighbors.push_back(areaA - 1);
    if (areaA + 1 <= regionNum) neighbors.push_back(areaA + 1);

    unordered_set<int> omegaEmpty;
    unordered_map<int,int> custArea;
    vector<vector<int>> customersInArea;
    buildCustomerAreaMap_SelfOnly(ind, omegaEmpty, custArea, customersInArea);

    // B1) Insertion
    for (int areaB : neighbors) {
        auto Sb_to_A = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaB, areaA);
        if (Sb_to_A.empty()) continue;
        sortByVolumeDesc(Sb_to_A, parameters);

        for (int cid : Sb_to_A) {
            Individual backup = ind;

            for (int r = 0; r < repackRestarts; ++r) {
                ind = backup;
                applyMoveCustomerArea(ind, cid, areaA);

                if (!decodeEvalKeepObj1Baseline_BestOfN(ind, parameters, cargoLookup, obj1Baseline, /*N=*/10)
) {
                    continue;
                }

                Obj2State cur = getObj2State(ind);
                if (betterObj2Lex(cur, base)) return true;
            }
            ind = backup;
        }
    }

    // B2) Pair-exchange
    for (int areaB : neighbors) {
        auto A_to_B = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaA, areaB);
        auto B_to_A = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaB, areaA);
        if (A_to_B.empty() || B_to_A.empty()) continue;

        struct Cand { int iFromB, jFromA; long long improveMin; };
        vector<Cand> cands;
        long long loadA = ind.selfOwnedTrucks[areaA].loadedVolume;

        for (int i : B_to_A) {
            long long vi = parameters.totalVolume[i - 1];
            for (int j : A_to_B) {
                long long vj = parameters.totalVolume[j - 1];
                long long newA = loadA + vi - vj;
                long long improve = newA - loadA;
                if (improve > 0) cands.push_back({i, j, improve});
            }
        }
        if (cands.empty()) continue;

        sort(cands.begin(), cands.end(), [&](const Cand& x, const Cand& y){
            return x.improveMin > y.improveMin;
        });

        for (const auto& pc : cands) {
            Individual backup = ind;

            for (int r = 0; r < repackRestarts; ++r) {
                ind = backup;
                applySwapCustomersArea(ind, pc.iFromB, areaB, pc.jFromA, areaA);

                if (!decodeEvalKeepObj1Baseline_BestOfN(ind, parameters, cargoLookup, obj1Baseline, /*N=*/10)
) {
                    continue;
                }

                Obj2State cur = getObj2State(ind);
                if (betterObj2Lex(cur, base)) return true;
            }
            ind = backup;
        }
    }

    return false;
}

static bool tryTwoStepImproveOnArea_Obj2Lex(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    long long baselineF1,
    int areaA,
    int repackRestarts,
    bool isReduceMaxPhase
){
    Obj2State base2 = getObj2State(ind);

    // neighbors
    vector<int> neighbors;
    if (areaA - 1 >= 1) neighbors.push_back(areaA - 1);
    if (areaA + 1 <= regionNum) neighbors.push_back(areaA + 1);

    // 用你原本的 OSR 產生法（這裡不排除 omega，靠 Obj1 baseline 擋）
    unordered_set<int> omegaEmpty;
    unordered_map<int,int> custArea;
    vector<vector<int>> customersInArea;
    buildCustomerAreaMap_SelfOnly(ind, omegaEmpty, custArea, customersInArea);

    struct Move {
        int type; // 0=move, 1=swap
        int c1, from1, to1;
        int c2, from2, to2;
        long long pri;
    };

    vector<Move> step1;
    step1.reserve(2000);

    // --- build step1 moves (move + swap) ---
    for (int areaB : neighbors) {
        if (isReduceMaxPhase) {
            // move: areaA -> areaB
            auto Sa = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaA, areaB);
            sortByVolumeDesc(Sa, parameters);
            for (int cid : Sa) {
                step1.push_back({0, cid, areaA, areaB, -1, -1, -1, parameters.totalVolume[cid - 1]});
            }

            // swap: areaA <-> areaB (prioritize reducing load of areaA)
            auto Sb = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaB, areaA);
            long long loadA = ind.selfOwnedTrucks[areaA].loadedVolume;
            for (int i : Sa) {
                long long vi = parameters.totalVolume[i - 1];
                for (int j : Sb) {
                    long long vj = parameters.totalVolume[j - 1];
                    long long improveA = vi - vj; // >0 means reduce A
                    if (improveA > 0) step1.push_back({1, i, areaA, areaB, j, areaB, areaA, improveA});
                }
            }
        } else {
            // increase min: move areaB -> areaA
            auto Sb_to_A = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaB, areaA);
            sortByVolumeDesc(Sb_to_A, parameters);
            for (int cid : Sb_to_A) {
                step1.push_back({0, cid, areaB, areaA, -1, -1, -1, parameters.totalVolume[cid - 1]});
            }

            // swap: bring bigger into A
            auto A_to_B = getOSRCustomers_FromTo(customersInArea, parameters, omegaEmpty, areaA, areaB);
            long long loadA = ind.selfOwnedTrucks[areaA].loadedVolume;
            for (int i : Sb_to_A) {
                long long vi = parameters.totalVolume[i - 1];
                for (int j : A_to_B) {
                    long long vj = parameters.totalVolume[j - 1];
                    long long improve = vi - vj;
                    if (improve > 0) step1.push_back({1, i, areaB, areaA, j, areaA, areaB, improve});
                }
            }
        }
    }

    sort(step1.begin(), step1.end(), [](const Move& a, const Move& b){
        return a.pri > b.pri;
    });

    // --- try step1 then step2 ---
    for (const auto& m1 : step1) {
        Individual s0 = ind;

        for (int r1 = 0; r1 < repackRestarts; ++r1) {
            ind = s0;

            if (m1.type == 0) applyMoveCustomerArea(ind, m1.c1, m1.to1);
            else applySwapCustomersArea(ind, m1.c1, m1.from1, m1.c2, m1.from2);

            if (!evalCandidate_NoWorseObj1(ind, parameters, cargoLookup, baselineF1)) continue;

            Obj2State st1 = getObj2State(ind);

            // ✅ step1 允許不變差（plateau）
            if (!notWorseObj2Lex(st1, base2)) continue;

            Individual after1 = ind;

            // step2: 再做一個 move（這裡先用 insertion 就很有效；你要也可再加 swap）
            unordered_map<int,int> custArea2;
            vector<vector<int>> customersInArea2;
            buildCustomerAreaMap_SelfOnly(after1, omegaEmpty, custArea2, customersInArea2);

            for (int areaB : neighbors) {
                vector<int> cand2;
                if (isReduceMaxPhase)
                    cand2 = getOSRCustomers_FromTo(customersInArea2, parameters, omegaEmpty, areaA, areaB);
                else
                    cand2 = getOSRCustomers_FromTo(customersInArea2, parameters, omegaEmpty, areaB, areaA);

                if (cand2.empty()) continue;
                sortByVolumeDesc(cand2, parameters);

                for (int cid2 : cand2) {
                    Individual s1 = after1;

                    for (int r2 = 0; r2 < repackRestarts; ++r2) {
                        ind = s1;

                        if (isReduceMaxPhase) applyMoveCustomerArea(ind, cid2, areaB);
                        else applyMoveCustomerArea(ind, cid2, areaA);

                        if (!evalCandidate_NoWorseObj1(ind, parameters, cargoLookup, baselineF1)) continue;

                        Obj2State st2 = getObj2State(ind);

                        // ✅ step2 要求嚴格改善（真正進步）
                        if (betterObj2Lex(st2, base2)) return true;
                    }
                }
            }
        }

        ind = s0;
    }

    return false;
}

void localSearch_Obj2_TeacherStyle(
    Individual& ind,
    const Data& parameters,
    const unordered_map<int, unordered_map<int, Cargo>>& cargoLookup,
    int repackRestarts
){
    // 入口假設：ind 已經有有效的 fitness & trucks 狀態
    if (ind.fitness.size() < 2) return;
    const long long obj1Baseline = ind.fitness[1]; // ✅ Obj1 不變差（或你要完全不變就改 evalCandidate）

    bool improved = true;
    while (improved) {
        improved = false;

        // ===================== Phase A: Reduce Max =====================
        {
            vector<int> areas(regionNum);
            iota(areas.begin(), areas.end(), 1);
            sort(areas.begin(), areas.end(), [&](int x, int y){
                return ind.selfOwnedTrucks[x].loadedVolume > ind.selfOwnedTrucks[y].loadedVolume;
            });

            for (int a : areas) {
                // ✅ 先用兩步 look-ahead（允許 step1 plateau，step2 必須嚴格改善 Obj2Lex）
                if (tryTwoStepImproveOnArea_Obj2Lex(
                        ind, parameters, cargoLookup,
                        obj1Baseline, a, repackRestarts,
                        /*isReduceMaxPhase=*/true
                    )) {
                    ind.fitness.clear();
                    evaluateFitness(ind, parameters);
                    improved = true;
                    break;
                }

                // （可選）fallback：原本一步 operator（如果你想保留）
                if (tryReduceMaxOnArea_KeepObj1Baseline(
                        ind, parameters, cargoLookup,
                        obj1Baseline, a, repackRestarts
                    )) {
                    ind.fitness.clear();
                    evaluateFitness(ind, parameters);
                    improved = true;
                    break;
                }
            }
        }

        // ===================== Phase B: Increase Min =====================
        {
            vector<int> areas(regionNum);
            iota(areas.begin(), areas.end(), 1);
            sort(areas.begin(), areas.end(), [&](int x, int y){
                return ind.selfOwnedTrucks[x].loadedVolume < ind.selfOwnedTrucks[y].loadedVolume;
            });

            for (int a : areas) {
                if (tryTwoStepImproveOnArea_Obj2Lex(
                        ind, parameters, cargoLookup,
                        obj1Baseline, a, repackRestarts,
                        /*isReduceMaxPhase=*/false
                    )) {
                    ind.fitness.clear();
                    evaluateFitness(ind, parameters);
                    improved = true;
                    break;
                }

                // （可選）fallback：原本一步 operator
                if (tryIncreaseMinOnArea_KeepObj1Baseline(
                        ind, parameters, cargoLookup,
                        obj1Baseline, a, repackRestarts
                    )) {
                    ind.fitness.clear();
                    evaluateFitness(ind, parameters);
                    improved = true;
                    break;
                }
            }
        }
    }
}
