#include "data.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <unordered_map>

using namespace std;

// 處理讀參數相關資料

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
        getline(ss, value, ','); c.volume = stoi(value);
        getline(ss, value, ','); c.lwh[0] = stoi(value);
        getline(ss, value, ','); c.lwh[1] = stoi(value);
        getline(ss, value, ','); c.lwh[2] = stoi(value);
        for (int i = 0; i < 6; ++i) {
            getline(ss, value, ','); c.orientation[i] = stoi(value);
        }
        getline(ss, value, ','); c.fragility = stoi(value);

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
        int customer_id, count, totalVolume;

        getline(ss, value, ',');
        customer_id = stoi(value);
        getline(ss, value, ',');
        count = stoi(value);
        getline(ss, value, ',');
        totalVolume = stoi(value);

        data.cargoNumber[customer_id - 1] = count;
        data.totalVolume[customer_id - 1] = totalVolume;
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
bool BLPlacement3D::tryInsert(vector<Gene>& group) {
    vector<Box> tempPlaced = placedBoxes;
    for (auto& g : group) {
        Box b = getBoxFromGene(g);
        if (!placeBox(b, tempPlaced)) {
            return false;
        }
        g.position[0] = b.x;
        g.position[1] = b.y;
        g.position[2] = b.z;
        tempPlaced.push_back(b);
    }

    placedBoxes = tempPlaced; 
    return true;
} 

bool BLPlacement3D::placeBox(Box& box, const vector<Box>& currentBoxes) {
    vector<tuple<int, int, int>> anchorPoints = {{0, 0, 0}};
    for (const auto& b : currentBoxes) {
        anchorPoints.push_back({b.x + b.l, b.y, b.z});
        anchorPoints.push_back({b.x, b.y + b.w, b.z});
        anchorPoints.push_back({b.x, b.y, b.z + b.h});
    }

    for (const auto& [ax, ay, az] : anchorPoints) {
        box.x = ax;
        box.y = ay;
        box.z = az;
        if (isWithinContainer(box) && !hasCollision(box, currentBoxes) ) {
            return true;
        }
    }
    return false;
}

bool BLPlacement3D::isWithinContainer(const Box& b) {
    return b.x + b.l <= containerL && b.y + b.w <= containerW && b.z + b.h <= containerH;
}

bool BLPlacement3D::hasCollision(const Box& b, const vector<Box>& boxes) {
    for (const auto& p : boxes) {
        if (!(b.x + b.l <= p.x || p.x + p.l <= b.x ||
                b.y + b.w <= p.y || p.y + p.w <= b.y ||
                b.z + b.h <= p.z || p.z + p.h <= b.z)) {
            return true;
        }
    }
    return false;
}

bool BLPlacement3D::isSupported(const Box& b, const vector<Box>& boxes) {
    if (b.z == 0) return true;
    int supportArea = 0, baseArea = b.l * b.w;
    for (const auto& p : boxes) {
        if (p.z + p.h == b.z) {
            int xOverlap = max(0, min(b.x + b.l, p.x + p.l) - max(b.x, p.x));
            int yOverlap = max(0, min(b.y + b.w, p.y + p.w) - max(b.y, p.y));
            supportArea += xOverlap * yOverlap;
        }
    }
    return supportArea >= 0.8 * baseArea; 
}