#include "data.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <cmath>
#include <cctype>

using namespace std;

// 處理讀參數相關資料
const double SCALE = 4.5;          // 改成 1.0, 1.25, 1.5, 1.75, 2.0
const double LEN_SCALE = cbrt(SCALE);

static bool isCollisionDominantViolation(const PlacementViolationInfo& v) {
    return
        v.hasFailedPlacement &&
        v.bestBoundaryViolation <= 1e-9 &&
        v.hasBoundaryZeroCandidate &&
        v.hasBoundaryZeroButCollision &&
        v.collisionFailCount > 0;
}

static RecordCargoInfo makeRecordCargoInfo(
    const Gene& g,
    const std::unordered_map<int, std::unordered_map<int, Cargo>>& cargoLookup
) {
    RecordCargoInfo info;
    info.customerId = g.customerId;
    info.cargoId = g.cargoId;

    auto it1 = cargoLookup.find(g.customerId);
    if (it1 != cargoLookup.end()) {
        auto it2 = it1->second.find(g.cargoId);
        if (it2 != it1->second.end()) {
            const Cargo& c = it2->second;
            info.l = c.lwh[0];
            info.w = c.lwh[1];
            info.h = c.lwh[2];
            for (int r = 0; r < 6; ++r) {
                info.orientation[r] = c.orientation[r];
            }
        }
    }

    return info;
}

static std::vector<RecordCargoInfo> makeRecordCargoInfoList(
    const std::vector<Gene>& genes,
    const std::unordered_map<int, std::unordered_map<int, Cargo>>& cargoLookup
) {
    std::vector<RecordCargoInfo> out;
    out.reserve(genes.size());
    for (const auto& g : genes) {
        out.push_back(makeRecordCargoInfo(g, cargoLookup));
    }
    return out;
}

static std::vector<int> extractUniqueCustomerList(const std::vector<Gene>& genes) {
    std::vector<int> out;
    std::unordered_set<int> seen;
    for (const auto& g : genes) {
        if (!seen.count(g.customerId)) {
            seen.insert(g.customerId);
            out.push_back(g.customerId);
        }
    }
    return out;
}

static std::vector<int> mergeCustomerListWithTryingCustomer(
    const std::vector<int>& loadedCustomerList,
    int tryingCustomer
) {
    std::vector<int> out = loadedCustomerList;
    bool exists = false;
    for (int x : out) {
        if (x == tryingCustomer) {
            exists = true;
            break;
        }
    }
    if (!exists) out.push_back(tryingCustomer);
    return out;
}


static bool isSameTruckCase(
    const BetterButFailedRecord& a,
    const BetterButFailedRecord& b
) {
    if (a.area != b.area) return false;
    if (a.truckId != b.truckId) return false;
    if (a.tryingCustomer != b.tryingCustomer) return false;
    if (a.truckCaseCustomerList.size() != b.truckCaseCustomerList.size()) return false;

    for (size_t i = 0; i < a.truckCaseCustomerList.size(); ++i) {
        if (a.truckCaseCustomerList[i] != b.truckCaseCustomerList[i]) return false;
    }
    return true;
}

static bool isBetterFailedRecord(
    const BetterButFailedRecord& a,
    const BetterButFailedRecord& b
) {
    // 1. 先比 fitness，越小越好
    if (a.fitnessValue != b.fitnessValue) {
        return a.fitnessValue < b.fitnessValue;
    }

    // 2. collision-dominant 優先
    if (a.collisionDominant != b.collisionDominant) {
        return a.collisionDominant && !b.collisionDominant;
    }

    // 3. boundary violation 越小越好
    if (a.bestBoundaryViolation != b.bestBoundaryViolation) {
        return a.bestBoundaryViolation < b.bestBoundaryViolation;
    }

    // 4. overlap 越小越好
    if (a.bestTotalOverlap != b.bestTotalOverlap) {
        return a.bestTotalOverlap < b.bestTotalOverlap;
    }

    // 5. collisionFailCount 越小越好
    if (a.collisionFailCount != b.collisionFailCount) {
        return a.collisionFailCount < b.collisionFailCount;
    }

    // 6. candidateCount 越小越好
    return a.candidateCount < b.candidateCount;
}

void tryAddFailedRecord(
    std::vector<BetterButFailedRecord>& pool,
    const BetterButFailedRecord& rec,
    size_t maxKeep
) {
    // 先做同類 case 去重：若同 case 已存在，只保留更好的那筆
    for (auto& oldRec : pool) {
        if (isSameTruckCase(oldRec, rec)) {
            if (isBetterFailedRecord(rec, oldRec)) {
                oldRec = rec;
            }
            return;
        }
    }

    pool.push_back(rec);

    std::sort(pool.begin(), pool.end(), isBetterFailedRecord);

    if (pool.size() > maxKeep) {
        pool.resize(maxKeep);
    }
}

double BLPlacement3D::computeBoundaryViolation(
    const Box& box,
    double& overflowX,
    double& overflowY,
    double& overflowZ
) {
    overflowX = std::max(0.0, double(box.x + box.l - containerL));
    overflowY = std::max(0.0, double(box.y + box.w - containerW));
    overflowZ = std::max(0.0, double(box.z + box.h - containerH));

    return overflowX + overflowY + overflowZ;
}
double BLPlacement3D::computeTotalOverlapVolume(
    const Box& b,
    const vector<Box>& boxes
) {
    double totalOverlap = 0.0;

    for (const auto& p : boxes) {
        int overlapX = max(0, min(b.x + b.l, p.x + p.l) - max(b.x, p.x));
        int overlapY = max(0, min(b.y + b.w, p.y + p.w) - max(b.y, p.y));
        int overlapZ = max(0, min(b.z + b.h, p.z + p.h) - max(b.z, p.z));

        totalOverlap += 1.0 * overlapX * overlapY * overlapZ;
    }

    return totalOverlap;
}

void readGoodsCSV(const string& filename, Data& data) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open goods file: " << filename << endl;
        return;
    }
    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        Cargo c;

        getline(ss, value, ','); c.customerId = stoi(value);
        getline(ss, value, ','); c.cargoId = stoi(value);

        getline(ss, value, ',');
        c.volume = stod(value) * SCALE;

        getline(ss, value, ',');
        c.lwh[0] = static_cast<int>(round(stod(value) * LEN_SCALE));

        getline(ss, value, ',');
        c.lwh[1] = static_cast<int>(round(stod(value) * LEN_SCALE));

        getline(ss, value, ',');
        c.lwh[2] = static_cast<int>(round(stod(value) * LEN_SCALE));

        for (int i = 0; i < 6; ++i) {
            getline(ss, value, ',');
            c.orientation[i] = stoi(value);
        }

        getline(ss, value, ',');
        c.fragility = stoi(value);

        data.cargoInformation.push_back(c);
    }
    file.close();
}

// 讀取服務區域
void readServiceAreaCSV(const string& filename, Data& data) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open service area file: " << filename << endl;
        return;
    }
    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        int customer_id;

        getline(ss, value, ',');
        customer_id = stoi(value);

        for (int i = 0; i < regionNum; ++i) {
            getline(ss, value, ',');
            data.serviceRegion[customer_id - 1][i] = stoi(value);
        }
    }
    file.close();
}

// 讀取貨物數量
void readCustomerInfoCSV(const string& filename, Data& data) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open cargo count file: " << filename << endl;
        return;
    }
    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        int customer_id, count;
        double totalVolume;

        getline(ss, value, ',');
        customer_id = stoi(value);

        getline(ss, value, ',');
        count = stoi(value);

        getline(ss, value, ',');
        totalVolume = stod(value);

        data.cargoNumber[customer_id - 1] = count;
        data.totalVolume[customer_id - 1] = totalVolume * SCALE;
    }
    file.close();
}
void readRouteToCSV(const string& filename, Data& data) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open cargo count file: " << filename << endl;
        return;
    }
    string line;
    getline(file, line); // Skip header
    while (getline(file, line)) {
        stringstream ss(line);
        string value;
        int region;
        vector <int> route;
        getline(ss, value, ',');
        region = stoi(value);
        while (getline(ss, value, ',')) {
            if (value.empty()) continue;  // 空欄略過
            int node = stoi(value);
            if (node > 0) route.push_back(node); // 忽略 depot
        }
        if (region >= 0 && region < regionNum)
            data.route[region] = route;
        else
            cerr << "Invalid region index: " << region << endl;
    }
    file.close();
}
void readParameters(const string& customerInfo, const string& goods, const string& serviceArea, const string& routes, Data& parameter) {
    readCustomerInfoCSV(customerInfo, parameter);
    readGoodsCSV(goods, parameter);
    readServiceAreaCSV(serviceArea, parameter);
    readRouteToCSV(routes, parameter);
}

inline std::string trim2(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace((unsigned char)s[start])) {
        ++start;
    }

    size_t end = s.size();
    while (end > start && std::isspace((unsigned char)s[end - 1])) {
        --end;
    }

    return s.substr(start, end - start);
}

void printData(const Data& data) {
    cout << "===== Cargo Info =====" << endl;
    for (size_t i = 0; i < data.cargoInformation.size(); ++i) {
        cout << "CustomerId=" << data.cargoInformation[i].customerId
             << ", CargoId= " << data.cargoInformation[i].cargoId
             << ", volume=" << data.cargoInformation[i].volume << endl;  
    }

    cout << "\n===== Customer Cargo Number & Volume =====" << endl;
    for (int i = 0; i < Customer; ++i) {
        cout << "Customer " << i
            << " | CargoNum=" << data.cargoNumber[i]
            << " | TotalVol=" << data.totalVolume[i]
            << endl;
    }

    cout << "\n===== Service Region Assignment =====" << endl;
    for (int i = 0; i < Customer; ++i) {
        cout << "Customer " << i << ": ";
        for (int r = 0; r < regionNum; ++r) {
            cout << data.serviceRegion[i][r] << " ";
        }
        cout << endl;
    }

    cout << "\n===== Routes =====" << endl;
    for (int r = 0; r < regionNum; ++r) {
        if (data.route[r].empty()) continue;
        cout << "Region " << r + 1 << ": ";
        for (int node : data.route[r]) {
            cout << node << " ";
        }
        cout << endl;
    }
    cout << "=================================================" << endl;
}

void printChromosomeInfo(const Individual& indiv){
    cout << "======== Individual ========\n";
    for (size_t j = 0; j < indiv.chromosome.size(); ++j) {
        const Gene& g = indiv.chromosome[j];
        cout << "Gene " << j + 1
             << " | Customer: " << g.customerId
             << " | CargoID: "  << g.cargoId
             << " | RouteArea: "  << g.routeArea
             << " | undecodedServiceArea: " << g.undecodedServiceArea
             << " | decodedServiceArea: "   << g.decodedServiceArea
             << " | undecodedRotation: " << g.undecodedRotation 
             << " | decodedRotation: " << g.decodedRotation << "\n";
    }
    cout << endl;
}

void BLPlacement3D::setCargoLookup(const unordered_map<int, unordered_map<int, Cargo>>& lookup) {
    cargoLookup = lookup;
}
bool BLPlacement3D::tryInsert(std::vector<Gene>& group, int maxTries) {
    static std::mt19937 rng(std::random_device{}());

    resetViolationInfo();

    auto tryOneOrder = [&](std::vector<Gene> cand) -> bool {
        std::vector<Box> tempPlaced = placedBoxes;

        for (auto& g : cand) {
            int originalRot = g.undecodedRotation;
            if (originalRot < 1 || originalRot > 6) originalRot = 1;

            int rotList[6] = {1, 2, 3, 4, 5, 6};
            std::swap(rotList[0], rotList[originalRot - 1]);

            bool placed = false;
            Box placedBox;

            PlacementViolationInfo bestFailInfo;
            bool hasBestFailInfo = false;

            for (int k = 0; k < 6; ++k) {
                g.undecodedRotation = rotList[k];
                g.decodedRotation = rotList[k];

                Box b = getBoxFromGene(g);

                if (placeBox(b, tempPlaced)) {
                    placed = true;
                    placedBox = b;
                    break;
                } else {
                    PlacementViolationInfo cur;
                    cur.hasFailedPlacement = true;
                    cur.failedCustomerId = g.customerId;
                    cur.failedCargoId = g.cargoId;

                    cur.bestBoundaryViolation = lastPlaceBoxSummary.bestBoundaryViolation;
                    cur.overflowX = lastPlaceBoxSummary.bestOverflowX;
                    cur.overflowY = lastPlaceBoxSummary.bestOverflowY;
                    cur.overflowZ = lastPlaceBoxSummary.bestOverflowZ;

                    cur.candidateCount = lastPlaceBoxSummary.candidateCount;
                    cur.boundaryFailCount = lastPlaceBoxSummary.boundaryFailCount;
                    cur.collisionFailCount = lastPlaceBoxSummary.collisionFailCount;
                    cur.hasBoundaryZeroCandidate = lastPlaceBoxSummary.hasBoundaryZeroCandidate;
                    cur.hasBoundaryZeroButCollision = lastPlaceBoxSummary.hasBoundaryZeroButCollision;
                    cur.bestTotalOverlap = lastPlaceBoxSummary.bestTotalOverlap;

                    if (!hasBestFailInfo) {
                        bestFailInfo = cur;
                        hasBestFailInfo = true;
                    } else {
                        bool curBetter = false;

                        // 規則 1：優先記錄「沒超界但撞到」的 case
                        if (cur.hasBoundaryZeroButCollision && !bestFailInfo.hasBoundaryZeroButCollision) {
                            curBetter = true;
                        }
                        else if (cur.hasBoundaryZeroButCollision == bestFailInfo.hasBoundaryZeroButCollision) {
                            // 規則 2：boundary violation 較小者優先
                            if (cur.bestBoundaryViolation < bestFailInfo.bestBoundaryViolation) {
                                curBetter = true;
                            }
                            else if (cur.bestBoundaryViolation == bestFailInfo.bestBoundaryViolation) {
                                // 規則 3：overlap 較小者優先
                                if (cur.bestTotalOverlap < bestFailInfo.bestTotalOverlap) {
                                    curBetter = true;
                                }
                            }
                        }

                        if (curBetter) {
                            bestFailInfo = cur;
                        }
                    }
                }
            }

            if (!placed) {
                if (hasBestFailInfo) {
                    lastViolationInfo = bestFailInfo;
                } else {
                    lastViolationInfo.hasFailedPlacement = true;
                    lastViolationInfo.failedCustomerId = g.customerId;
                    lastViolationInfo.failedCargoId = g.cargoId;
                }
                return false;
            }

            g.position[0] = placedBox.x;
            g.position[1] = placedBox.y;
            g.position[2] = placedBox.z;
            tempPlaced.push_back(placedBox);
        }

        group = std::move(cand);
        placedBoxes = std::move(tempPlaced);
        return true;
    };

    if ((int)group.size() <= 4) {
        std::vector<Gene> base = group;

        if (tryOneOrder(base)) return true;

        std::sort(base.begin(), base.end(),
                  [](const Gene& a, const Gene& b) {
                      if (a.customerId != b.customerId) return a.customerId < b.customerId;
                      return a.cargoId < b.cargoId;
                  });

        do {
            std::vector<Gene> cand = base;
            if (tryOneOrder(cand)) return true;
        } while (std::next_permutation(
            base.begin(), base.end(),
            [](const Gene& a, const Gene& b) {
                if (a.customerId != b.customerId) return a.customerId < b.customerId;
                return a.cargoId < b.cargoId;
            }
        ));

        return false;
    }

    for (int t = 0; t < maxTries; ++t) {
        std::vector<Gene> cand = group;
        if (t > 0) std::shuffle(cand.begin(), cand.end(), rng);
        if (tryOneOrder(cand)) return true;
    }

    return false;
}
bool BLPlacement3D::placeBox(Box& box, const vector<Box>& currentBoxes) {
    resetPlaceBoxSummary();

    vector<tuple<int, int, int>> anchorPoints;
    anchorPoints.reserve(1 + 3 * currentBoxes.size());
    anchorPoints.push_back({0, 0, 0});

    for (const auto& b : currentBoxes) {
        anchorPoints.push_back({b.x + b.l, b.y, b.z});
        anchorPoints.push_back({b.x, b.y + b.w, b.z});
        anchorPoints.push_back({b.x, b.y, b.z + b.h});
    }

    sort(anchorPoints.begin(), anchorPoints.end(),
         [](const auto& a, const auto& b) {
             if (get<2>(a) != get<2>(b)) return get<2>(a) < get<2>(b); // z
             if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b); // y
             return get<0>(a) < get<0>(b);                             // x
         });

    anchorPoints.erase(unique(anchorPoints.begin(), anchorPoints.end()), anchorPoints.end());

    for (const auto& [ax, ay, az] : anchorPoints) {
        box.x = ax;
        box.y = ay;
        box.z = az;

        lastPlaceBoxSummary.candidateCount++;

        double ox = 0.0, oy = 0.0, oz = 0.0;
        double boundaryV = computeBoundaryViolation(box, ox, oy, oz);

        bool within = isWithinContainer(box);
        bool collision = false;
        double overlapV = 0.0;

        if (within) {
            lastPlaceBoxSummary.hasBoundaryZeroCandidate = true;
            collision = hasCollision(box, currentBoxes);

            if (collision) {
                lastPlaceBoxSummary.hasBoundaryZeroButCollision = true;
                lastPlaceBoxSummary.collisionFailCount++;

                overlapV = computeTotalOverlapVolume(box, currentBoxes);
                if (overlapV < lastPlaceBoxSummary.bestTotalOverlap) {
                    lastPlaceBoxSummary.bestTotalOverlap = overlapV;
                }
            }
        } else {
            lastPlaceBoxSummary.boundaryFailCount++;
        }

        if (boundaryV < lastPlaceBoxSummary.bestBoundaryViolation) {
            lastPlaceBoxSummary.bestBoundaryViolation = boundaryV;
            lastPlaceBoxSummary.bestOverflowX = ox;
            lastPlaceBoxSummary.bestOverflowY = oy;
            lastPlaceBoxSummary.bestOverflowZ = oz;
        }

        if (within && !collision) {
            return true;
        }
    }

    if (lastPlaceBoxSummary.bestBoundaryViolation == 1e18) {
        lastPlaceBoxSummary.bestBoundaryViolation = 0.0;
        lastPlaceBoxSummary.bestOverflowX = 0.0;
        lastPlaceBoxSummary.bestOverflowY = 0.0;
        lastPlaceBoxSummary.bestOverflowZ = 0.0;
    }

    if (lastPlaceBoxSummary.bestTotalOverlap == 1e18) {
        lastPlaceBoxSummary.bestTotalOverlap = 0.0;
    }

    return false;
}
bool BLPlacement3D::isWithinContainer(const Box& b) {
    return b.x >= 0 && b.y >= 0 && b.z >= 0 &&
           b.x + b.l <= containerL &&
           b.y + b.w <= containerW &&
           b.z + b.h <= containerH;
}

bool BLPlacement3D::hasCollision(const Box& b, const vector<Box>& boxes) {
    for (const auto& p : boxes) {
        bool separated =
            (b.x + b.l <= p.x) || (p.x + p.l <= b.x) ||
            (b.y + b.w <= p.y) || (p.y + p.w <= b.y) ||
            (b.z + b.h <= p.z) || (p.z + p.h <= b.z);

        if (!separated) {
            return true;
        }
    }
    return false;
}



bool Data::loadGeneratedXY(const string& filepath) {
    solomonXY.clear();

    ifstream fin(filepath);
    if (!fin.is_open()) {
        cerr << "Cannot open generated coordinate file: " << filepath << endl;
        return false;
    }

    string line;

    // 跳過表頭
    if (!getline(fin, line)) {
        cerr << "Empty coordinate file: " << filepath << endl;
        return false;
    }

    // depot 手動補上
    solomonXY[0] = {50.0, 25.0};

    while (getline(fin, line)) {
        string t = trim2(line);
        if (t.empty()) continue;

        stringstream ss(t);
        string token;

        vector<string> cols;
        while (getline(ss, token, ',')) {
            cols.push_back(trim2(token));
        }

        if (cols.size() < 3) continue;

        try {
            int id = stoi(cols[0]);
            double x = stod(cols[1]);
            double y = stod(cols[2]);
            solomonXY[id] = {x, y};
        }
        catch (...) {
            continue;
        }
    }

    if (solomonXY.empty()) {
        cerr << "Generated XY empty, check file format: " << filepath << endl;
        return false;
    }

    return true;
}

double Data::getDistance(int i, int j) const {
    auto it1 = solomonXY.find(i);
    auto it2 = solomonXY.find(j);
    if (it1 == solomonXY.end() || it2 == solomonXY.end()) {
        // 如果你的 route/customer id 沒有對到 Solomon，會走到這裡
        return 1e100;
    }
    double dx = it1->second.first  - it2->second.first;
    double dy = it1->second.second - it2->second.second;
    return sqrt(dx*dx + dy*dy);
}

inline std::string buildSetKey(std::vector<int> customers) {
    std::sort(customers.begin(), customers.end());

    std::string key;
    for (int i = 0; i < (int)customers.size(); ++i) {
        key += std::to_string(customers[i]);
        if (i + 1 < (int)customers.size()) key += '-';
    }
    return key;
}


inline std::string buildOrderKey(const std::vector<int>& route) {
    std::string key;
    for (size_t i = 0; i < route.size(); ++i) {
        if (i) key += "-";
        key += std::to_string(route[i]);
    }
    return key;
}

static double getTotalOverlapScore(const SavedViolationSolution& s) {
    double sum = 0.0;
    for (const auto& r : s.failedRecords) {
        sum += r.bestTotalOverlap;
    }
    return sum;
}

static bool isBetterSavedViolationSolution(
    const SavedViolationSolution& a,
    const SavedViolationSolution& b
) {
    if (a.fitnessValue != b.fitnessValue) {
        return a.fitnessValue < b.fitnessValue;
    }

    if (a.failTruckCount != b.failTruckCount) {
        return a.failTruckCount < b.failTruckCount;
    }

    return getTotalOverlapScore(a) < getTotalOverlapScore(b);
}

static SavedViolationSolution buildSavedViolationSolution(
    const Individual& indiv,
    size_t maxCasesPerSolution = 2
) {
    SavedViolationSolution s;
    s.fitnessValue = (!indiv.fitness.empty() ? indiv.fitness[0] : 1e18);

    size_t takeN = std::min(maxCasesPerSolution, indiv.failedTruckRecords.size());
    s.failedRecords.assign(
        indiv.failedTruckRecords.begin(),
        indiv.failedTruckRecords.begin() + takeN
    );

    s.failTruckCount = (int)s.failedRecords.size();
    return s;
}




static bool isSameSavedViolationSolution(
    const SavedViolationSolution& a,
    const SavedViolationSolution& b
) {
    if (a.failTruckCount != b.failTruckCount) return false;
    if (a.failedRecords.size() != b.failedRecords.size()) return false;

    for (size_t i = 0; i < a.failedRecords.size(); ++i) {
        const auto& ra = a.failedRecords[i];
        const auto& rb = b.failedRecords[i];

        if (ra.area != rb.area) return false;
        if (ra.truckId != rb.truckId) return false;
        if (ra.tryingCustomer != rb.tryingCustomer) return false;
        if (ra.truckCaseCustomerList != rb.truckCaseCustomerList) return false;
    }
    return true;
}

static void tryAddSavedViolationSolution(
    std::vector<SavedViolationSolution>& pool,
    const SavedViolationSolution& cand,
    size_t maxKeep = 5
) {
    if (cand.failedRecords.empty()) return;

    for (auto& oldSol : pool) {
        if (isSameSavedViolationSolution(oldSol, cand)) {
            if (isBetterSavedViolationSolution(cand, oldSol)) {
                oldSol = cand;
            }
            return;
        }
    }

    pool.push_back(cand);

    std::sort(pool.begin(), pool.end(), isBetterSavedViolationSolution);

    if (pool.size() > maxKeep) {
        pool.resize(maxKeep);
    }
}

static string buildSortedCustomerKey(vector<int> customerSet) {
    sort(customerSet.begin(), customerSet.end());
    customerSet.erase(unique(customerSet.begin(), customerSet.end()), customerSet.end());

    string key;
    for (size_t i = 0; i < customerSet.size(); ++i) {
        if (i) key += "-";
        key += to_string(customerSet[i]);
    }
    return key;
}

static vector<int> buildRouteFromFixedOrder(
    const vector<int>& customerSet,
    const Data& parameters,
    int area
) {
    vector<int> route;
    unordered_set<int> inSet(customerSet.begin(), customerSet.end());

    const auto& baseRoute = parameters.route[area - 1];

    for (int cid : baseRoute) {
        if (cid != 0 && inSet.count(cid)) {
            route.push_back(cid);
        }
    }

    if (route.size() != customerSet.size()) {
        route = customerSet;
    }

    return route;
}



static BetterButFailedRecord buildBetterButFailedRecord(
    int area,
    int truckId,
    int containerL,
    int containerW,
    int containerH,
    int tryingCustomer,
    double fitnessValue,
    const PlacementViolationInfo& v,
    const std::vector<Gene>& loadedCargoSnapshot,
    const std::vector<Gene>& tryingCustomerCargoGroup,
    const std::vector<Gene>& remainingCustomerCargoGroup,
    const std::vector<int>& routeOrder,
    int failIndex,
    const std::unordered_map<int, std::unordered_map<int, Cargo>>& cargoLookup
) {
    BetterButFailedRecord rec;

    rec.fitnessValue = fitnessValue;

    rec.area = area;
    rec.truckId = truckId;

    rec.tryingCustomer = tryingCustomer;
    rec.failedCustomer = v.failedCustomerId;
    rec.failedCargo = v.failedCargoId;

    rec.candidateCount = v.candidateCount;
    rec.boundaryFailCount = v.boundaryFailCount;
    rec.collisionFailCount = v.collisionFailCount;

    rec.hasBoundaryZeroCandidate = v.hasBoundaryZeroCandidate;
    rec.hasBoundaryZeroButCollision = v.hasBoundaryZeroButCollision;
    rec.collisionDominant = isCollisionDominantViolation(v);

    rec.bestBoundaryViolation = v.bestBoundaryViolation;
    rec.overflowX = v.overflowX;
    rec.overflowY = v.overflowY;
    rec.overflowZ = v.overflowZ;
    rec.bestTotalOverlap = v.bestTotalOverlap;

    rec.containerL = containerL;
    rec.containerW = containerW;
    rec.containerH = containerH;

    // local rescue case
    rec.loadedCustomerList = extractUniqueCustomerList(loadedCargoSnapshot);
    rec.truckCaseCustomerList = mergeCustomerListWithTryingCustomer(
        rec.loadedCustomerList,
        tryingCustomer
    );

    // full planned route info
    rec.fullPlannedCustomerList = routeOrder;
    rec.failIndexInPlannedRoute = failIndex;

    rec.remainingCustomerListAfterFail.clear();
    for (int i = failIndex + 1; i < (int)routeOrder.size(); ++i) {
        rec.remainingCustomerListAfterFail.push_back(routeOrder[i]);
    }

    // cargo info
    rec.baseLoadedCargoes = makeRecordCargoInfoList(loadedCargoSnapshot, cargoLookup);
    rec.tryingCustomerCargoes = makeRecordCargoInfoList(tryingCustomerCargoGroup, cargoLookup);
    rec.remainingCustomerCargoes = makeRecordCargoInfoList(remainingCustomerCargoGroup, cargoLookup);

    // local case = loaded + trying
    rec.allCargoesForTruckCase = rec.baseLoadedCargoes;
    rec.allCargoesForTruckCase.insert(
        rec.allCargoesForTruckCase.end(),
        rec.tryingCustomerCargoes.begin(),
        rec.tryingCustomerCargoes.end()
    );

    // full case = loaded + trying + remaining
    rec.fullPlannedCargoes = rec.allCargoesForTruckCase;
    rec.fullPlannedCargoes.insert(
        rec.fullPlannedCargoes.end(),
        rec.remainingCustomerCargoes.begin(),
        rec.remainingCustomerCargoes.end()
    );

    return rec;
}
