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

struct Truck{
    int truckId;
    int length = 300;
    int width = 170;
    int height = 165;
    vector<int> route;
    long long loadedVolume = 0;
    vector<Gene> assignedCargo; 
};

struct Individual {
    vector<Gene> chromosome;
    Truck selfOwnedTrucks[regionNum + 1];
    vector<Truck> rentedTrucks;
    vector<long long> fitness; //第一個放才積差距，第二個放租用的成本
    double rentedVehicleLoadingCost = 0.0;
    double maxVolumeDifferenceOfEachCar = 0.0;
};

class BLPlacement3D {
    public:
        struct Box {
            int x, y, z, l, w, h;
            int customerId, cargoId;
        };
        vector<Box> placedBoxes;
        unordered_map<int, unordered_map<int, Cargo>> cargoLookup;
        int containerL, containerW, containerH;

        BLPlacement3D(int L, int W, int H) : containerL(L), containerW(W), containerH(H) {}
        void setCargoLookup(const unordered_map<int, unordered_map<int, Cargo>>& lookup);       
        bool tryInsert(vector<Gene>& group, int  maxTries); 

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
        bool isSupported(const Box& b, const vector<Box>& boxes);
};
void readParameters(const string& customerInfo, const string& goods, const string& serviceArea, const string& routes, Data& parameter);
void printData(const Data& data);
void printChromosomeInfo(const Individual& indiv);
#endif