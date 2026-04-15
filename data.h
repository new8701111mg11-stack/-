// data.h
#ifndef DATA_H
#define DATA_H
#include <vector>
#include <string>
#include <string>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

using namespace std;
const int N = 61;
const int Customer = 12;
const int regionNum = 4;
const int selfOwnedTruck = regionNum;
const int rentedTruck = regionNum;
const int populationSize = 2;
const int maxGenerations = 10000;
const double crossoverRate = 0.8;
const double mutationRate = 0.1;
const double eliteRatio = 0.4;


struct Cargo{
    int customerId;
    int cargoId;
    int volume;
    int lwh[3];
    int orientation[6];
    int fragility;
};

struct Data{
    int cargoNumber[Customer];
    int totalVolume[Customer];
    int serviceRegion[Customer][regionNum];
    vector <Cargo> cargoInformation;
    vector <int> route[regionNum];
    unordered_map<int, pair<double,double>> solomonXY;
    bool loadGeneratedXY(const std::string& filepath);
    double getDistance(int i, int j) const;
};

struct Gene {
    int customerId;
    int cargoId;
    int routeArea;   
    int undecodedServiceArea;
    int decodedServiceArea = 0;       
    int undecodedRotation; 
    int decodedRotation = 0;       
    int position[3] = { -1, -1, -1 };
};
struct LoadingCache {
    // area -> setKey -> orderKey 成功
    std::vector<std::unordered_map<std::string, std::unordered_set<std::string>>> selfOwnedSuccessCache;

    LoadingCache() {
        selfOwnedSuccessCache.resize(regionNum + 1);
    }

    void clear() {
        for (auto& mp : selfOwnedSuccessCache) {
            mp.clear();
        }
    }

    bool hasSuccess(int area, const std::string& setKey, const std::string& orderKey) const {
        if (area < 1 || area > regionNum) return false;

        const auto& mp = selfOwnedSuccessCache[area];
        auto it = mp.find(setKey);
        if (it == mp.end()) return false;

        return it->second.count(orderKey) > 0;
    }

    void saveSuccess(int area, const std::string& setKey, const std::string& orderKey) {
        if (area < 1 || area > regionNum) return;
        selfOwnedSuccessCache[area][setKey].insert(orderKey);
    }
};
struct Truck{
    int truckId;
    int length = 300;
    int width = 170;
    int height = 165;
    vector<int> route;
    long long loadedVolume = 0;
    vector<Gene> assignedCargo; 
};
struct RecordCargoInfo {
    int customerId = -1;
    int cargoId = -1;

    int l = 0;
    int w = 0;
    int h = 0;

    int orientation[6] = {0, 0, 0, 0, 0, 0};
};
struct BetterButFailedRecord {
    // solution level info
    double fitnessValue = 0.0;

    // truck / violation meta
    int area = -1;
    int truckId = -1;

    int tryingCustomer = -1;
    int failedCustomer = -1;
    int failedCargo = -1;

    int candidateCount = 0;
    int boundaryFailCount = 0;
    int collisionFailCount = 0;

    bool hasBoundaryZeroCandidate = false;
    bool hasBoundaryZeroButCollision = false;
    bool collisionDominant = false;

    double bestBoundaryViolation = 0.0;
    double overflowX = 0.0;
    double overflowY = 0.0;
    double overflowZ = 0.0;
    double bestTotalOverlap = 0.0;

    // container
    int containerL = 0;
    int containerW = 0;
    int containerH = 0;

    // customer sets
    std::vector<int> loadedCustomerList;              // fail 前已成功放入的 customer
    std::vector<int> truckCaseCustomerList;           // 目前 local case = loaded + tryingCustomer

    std::vector<int> fullPlannedCustomerList;         // fixed route 下，這台車原本完整要處理的 customer
    std::vector<int> remainingCustomerListAfterFail;  // fail 之後，routeOrder 中尚未處理的 customer

    int failIndexInPlannedRoute = -1;                 // tryingCustomer 在 fullPlannedCustomerList 中的位置

    // cargo sets
    std::vector<RecordCargoInfo> baseLoadedCargoes;      // fail 前已在車上的貨
    std::vector<RecordCargoInfo> tryingCustomerCargoes;  // 當前 trying customer 的貨
    std::vector<RecordCargoInfo> allCargoesForTruckCase; // baseLoaded + tryingCustomer
    std::vector<RecordCargoInfo> remainingCustomerCargoes; // fail 後 remaining customers 的貨
    std::vector<RecordCargoInfo> fullPlannedCargoes;   
};


struct Individual {
    vector<Gene> chromosome;
    Truck selfOwnedTrucks[regionNum + 1];
    vector<Truck> rentedTrucks;
    vector<double> fitness; //第一個放才積差距，第二個放租用的成本
    double rentedVehicleLoadingCost = 0.0;
    double maxVolumeDifferenceOfEachCar = 0.0;
    std::vector<std::vector<int>> lastCustomerSet; 

    std::vector<BetterButFailedRecord> failedTruckRecords;
};

struct SavedViolationSolution {
    double fitnessValue = 0.0;
    int failTruckCount = 0;
    std::vector<BetterButFailedRecord> failedRecords;
};

struct SingleTruckGurobiCase {
    struct Item {
        int customerId = -1;
        int cargoId = -1;
        int l = 0;
        int w = 0;
        int h = 0;

        int orientation[6] = {0, 0, 0, 0, 0, 0};

        int isPreloaded = 0;
        // 1 = fail 前已在車上的貨
        // 0 = 其他（保留舊欄位相容性）

        int itemRole = 0;
        // 1 = loaded before fail
        // 2 = trying customer
        // 3 = remaining customers after fail
    };

    int caseId = -1;
    double fitnessValue = 0.0;

    int area = -1;
    int truckId = -1;

    int tryingCustomer = -1;
    int failedCustomer = -1;
    int failedCargo = -1;

    int containerL = 0;
    int containerW = 0;
    int containerH = 0;

    int candidateCount = 0;
    int boundaryFailCount = 0;
    int collisionFailCount = 0;
    double bestBoundaryViolation = 0.0;
    double bestTotalOverlap = 0.0;
    int collisionDominant = 0;

    // customer-level info
    std::vector<int> loadedCustomerList;
    std::vector<int> truckCaseCustomerList;          // local rescue = loaded + trying
    std::vector<int> fullPlannedCustomerList;        // fixed route 下完整規劃 customer
    std::vector<int> remainingCustomerListAfterFail; // fail 後還沒處理的 customer
    int failIndexInPlannedRoute = -1;

    // item-level info
    std::vector<Item> items;              // full planned cargoes 全部放進來
};


struct PlacementViolationInfo {
    bool hasFailedPlacement = false;

    int failedCustomerId = -1;
    int failedCargoId = -1;

    double bestBoundaryViolation = 0.0;
    double overflowX = 0.0;
    double overflowY = 0.0;
    double overflowZ = 0.0;

    int candidateCount = 0;
    int boundaryFailCount = 0;
    int collisionFailCount = 0;

    bool hasBoundaryZeroCandidate = false;
    bool hasBoundaryZeroButCollision = false;

    double bestTotalOverlap = 0.0;
};
class BLPlacement3D {
public:
    struct Box {
        int x, y, z, l, w, h;
        int customerId, cargoId;
    };

    struct CandidateCheckSummary {
        int candidateCount = 0;
        int boundaryFailCount = 0;
        int collisionFailCount = 0;

        bool hasBoundaryZeroCandidate = false;
        bool hasBoundaryZeroButCollision = false;

        double bestBoundaryViolation = 1e18;
        double bestOverflowX = 0.0;
        double bestOverflowY = 0.0;
        double bestOverflowZ = 0.0;

        double bestTotalOverlap = 1e18;
    };

    vector<Box> placedBoxes;
    unordered_map<int, unordered_map<int, Cargo>> cargoLookup;
    int containerL, containerW, containerH;

    PlacementViolationInfo lastViolationInfo;
    CandidateCheckSummary lastPlaceBoxSummary;

    BLPlacement3D(int L, int W, int H) : containerL(L), containerW(W), containerH(H) {}

    void setCargoLookup(const unordered_map<int, unordered_map<int, Cargo>>& lookup);
    bool tryInsert(vector<Gene>& group, int maxTries);

    void resetViolationInfo() {
        lastViolationInfo = PlacementViolationInfo{};
    }

    void resetPlaceBoxSummary() {
        lastPlaceBoxSummary = CandidateCheckSummary{};
    }

private:
    Box getBoxFromGene(const Gene& g) {
        const Cargo& c = cargoLookup[g.customerId][g.cargoId];
        int l = c.lwh[0], w = c.lwh[1], h = c.lwh[2];
        switch (g.decodedRotation) {
            case 1: l = c.lwh[0]; w = c.lwh[1]; h = c.lwh[2]; break;
            case 2: l = c.lwh[0]; w = c.lwh[2]; h = c.lwh[1]; break;
            case 3: l = c.lwh[1]; w = c.lwh[0]; h = c.lwh[2]; break;
            case 4: l = c.lwh[1]; w = c.lwh[2]; h = c.lwh[0]; break;
            case 5: l = c.lwh[2]; w = c.lwh[0]; h = c.lwh[1]; break;
            case 6: l = c.lwh[2]; w = c.lwh[1]; h = c.lwh[0]; break;
            default: break;
        }
        return Box{0, 0, 0, l, w, h, g.customerId, g.cargoId};
    }

    bool placeBox(Box& box, const vector<Box>& currentBoxes);
    bool isWithinContainer(const Box& b);
    bool hasCollision(const Box& b, const vector<Box>& boxes);

    double computeBoundaryViolation(
        const Box& box,
        double& overflowX,
        double& overflowY,
        double& overflowZ
    );

    double computeTotalOverlapVolume(
        const Box& b,
        const vector<Box>& boxes
    );
};
void readParameters(const string& customerInfo, const string& goods, const string& serviceArea, const string& routes, Data& parameter);
void printData(const Data& data);
void printChromosomeInfo(const Individual& indiv);

#endif