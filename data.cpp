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

    auto tryOneOrder = [&](std::vector<Gene> cand) -> bool {
        std::vector<Box> tempPlaced = placedBoxes;

        for (auto& g : cand) {
            int originalRot = g.undecodedRotation;
            if (originalRot < 1 || originalRot > 6) originalRot = 1;

            int rotList[6] = {1, 2, 3, 4, 5, 6};
            std::swap(rotList[0], rotList[originalRot - 1]);

            bool placed = false;
            Box placedBox;

            for (int k = 0; k < 6; ++k) {
                g.undecodedRotation = rotList[k];
                g.decodedRotation = rotList[k];

                Box b = getBoxFromGene(g);
                if (placeBox(b, tempPlaced)) {
                    placed = true;
                    placedBox = b;
                    break;
                }
            }

            if (!placed) {
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

        if (isWithinContainer(box) && !hasCollision(box, currentBoxes)) {
            return true;
        }
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

    // 假設 parameters.route 是 0-based: route[0], route[1], ...
    const auto& baseRoute = parameters.route[area - 1];

    for (int cid : baseRoute) {
        if (cid != 0 && inSet.count(cid)) {
            route.push_back(cid);
        }
    }

    // fallback：如果 base route 沒有完整涵蓋 customerSet
    if (route.size() != customerSet.size()) {
        route = customerSet;
    }

    return route;
}