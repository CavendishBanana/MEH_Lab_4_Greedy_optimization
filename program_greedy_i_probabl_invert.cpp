//Krzysztof Gryko - Lab 4 -Greedy optimization methods
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <fstream>
#include <bits/stdc++.h>

using namespace std;

double getNoConnectionValue()
{
    return std::numeric_limits<double>::max();
}

void buildMesh(vector<vector<double> >& graph,const vector<vector< double> > &coordinates)
{
    unsigned int nodesCount = coordinates.size();
    double noConnectionValue = getNoConnectionValue();
    for(unsigned int i=0; i < nodesCount; i++)
    {
        for(unsigned int j = i + 1; j < nodesCount; j++)
        {
            double distance = sqrt(pow(coordinates[j][0] - coordinates[i][0], 2.0) + pow(coordinates[j][1] - coordinates[i][1], 2.0));   
            graph[i][j] = distance;
            graph[j][i] = distance;
        }
        
    }
    for(unsigned int i=0u; i < nodesCount; i++)
    {
        graph[i][i]= noConnectionValue;
    }
    
}

unsigned int getClosestUnvisitedCityIdx(const vector<double> & neighboursDistances, const vector<bool> visitedFlags)
{
    double smallestDistance = getNoConnectionValue();
    unsigned int closestCityIdx = 0u;
    double noConnectVal = getNoConnectionValue();
    for(unsigned int i =0u; i < neighboursDistances.size(); i++)
    {
        if(neighboursDistances[i] != noConnectVal && neighboursDistances[i] < smallestDistance)
        {
            smallestDistance = neighboursDistances[i];
            closestCityIdx = i;
        }
    }
    return closestCityIdx;
}

vector<unsigned int> greedyMethodShortestPathInTSP(const vector<vector<double> >& graph, bool generateRandomStart=true, unsigned int startingCityIdx=0u)
{
    unsigned int startingCity = startingCityIdx;
    if(generateRandomStart)
    {
        default_random_engine engine = default_random_engine();
        uniform_int_distribution<unsigned int> distribution = uniform_int_distribution<unsigned int>(0, graph.size() - 1u);
        startingCity = distribution(engine);
    }
    vector<bool> visitedFlags = vector<bool>(graph.size(), false);
    visitedFlags[startingCity] = true;
    double totalDistance = 0.0;
    vector<unsigned int> citiesVisitingOrder(graph.size(),0u);
    citiesVisitingOrder[0u] = startingCity;
    unsigned int lastVisitedCityIdx = 0u;
    for(unsigned int i=0u; i < graph.size() -1u; i++)
    {
        unsigned int closestNeighbourIdx = getClosestUnvisitedCityIdx(graph[citiesVisitingOrder[lastVisitedCityIdx]], visitedFlags);
        visitedFlags[closestNeighbourIdx] = true;
        totalDistance += graph[citiesVisitingOrder[lastVisitedCityIdx]][closestNeighbourIdx];
        lastVisitedCityIdx++;
        citiesVisitingOrder[lastVisitedCityIdx] = closestNeighbourIdx;
        
    }
    return citiesVisitingOrder;
}

unsigned int selectRandomUnvisitedCityWithProbabsInvertlyProportionalToDistance(const vector<double>& neighbours, vector<bool> &visitedFlags, default_random_engine& engine)
{
    double inversionsSum = 0.0;
    vector<double> inversions(neighbours.size(), 0.0);
    unsigned int lastFoundUnvisitedIdx =0u;
    for(unsigned int i=0; i < neighbours.size(); i++)
    {
        if(!visitedFlags[i] && neighbours[i] != getNoConnectionValue())
        {
            inversions[i] = 1.0/neighbours[i];
            inversionsSum+=inversions[i];
            lastFoundUnvisitedIdx = i;
        }
    }
    double inversionOfInversionsSum = 1.0/inversionsSum;
    for(unsigned int i = 0; i < neighbours.size(); i++)
    {
        inversions[i]*=inversionOfInversionsSum;
    }
    uniform_real_distribution<double> distribution = uniform_real_distribution<double>(0.0,1.0);
    double randomValue = distribution(engine);
    double encounteredProbablsSum = 0.0;
    for(unsigned int i = 0; i < neighbours.size(); i++)
    {
        if(!visitedFlags[i] && neighbours[i] != getNoConnectionValue())
        {
            if(encounteredProbablsSum <= randomValue && randomValue < encounteredProbablsSum + inversions[i])
            {
                return i;
            }
            encounteredProbablsSum += inversions[i];
        }
    }
    return lastFoundUnvisitedIdx;
}

vector<unsigned int> greedyRandomMethodShortestPathInTSP(const vector<vector<double> >& graph, bool generateRandomStart=true, unsigned int startingCityIdx=0u)
{
    unsigned int startingCity = startingCityIdx;
    default_random_engine engine = default_random_engine();
    if(generateRandomStart)
    {
        
        uniform_int_distribution<unsigned int> distribution = uniform_int_distribution<unsigned int>(0, graph.size() - 1u);
        startingCity = distribution(engine);
    }
    vector<bool> visitedFlags = vector<bool>(graph.size(), false);
    visitedFlags[startingCity] = true;
    double totalDistance = 0.0;
    vector<unsigned int> citiesVisitingOrder(graph.size(),0u);
    citiesVisitingOrder[0u] = startingCity;
    unsigned int lastVisitedCityIdx = 0u;
    for(unsigned int i=0u; i < graph.size() -1u; i++)
    {
        unsigned int closestNeighbourIdx = selectRandomUnvisitedCityWithProbabsInvertlyProportionalToDistance(graph[citiesVisitingOrder[lastVisitedCityIdx]], visitedFlags, engine);
        visitedFlags[closestNeighbourIdx] = true;
        totalDistance += graph[citiesVisitingOrder[lastVisitedCityIdx]][closestNeighbourIdx];
        lastVisitedCityIdx++;
        citiesVisitingOrder[lastVisitedCityIdx] = closestNeighbourIdx;
        
    }
    return citiesVisitingOrder;
}

double getTotalPathLength(const vector<vector<double> > &graph, const vector<unsigned int>& citiesVisitingOrder)
{
    double pathLength = 0.0;
    for(unsigned int i =0; i< graph.size() - 1u; i++)
    {
        pathLength += graph[i][i+1];
    }
    pathLength+=graph[citiesVisitingOrder[citiesVisitingOrder.size() - 1u]][citiesVisitingOrder[0u]];
    return pathLength;
}

vector<string> extractCoordinatesStringRepresentations(const string& str)
{
    stringstream ss(str);
    vector<string> v(2,"");
    string s="";
    getline(ss, s, ' ');
    getline(ss, s, ' ');
    v[0] = s;
    getline(ss, s, '\n');
    v[1] = s;
    return v;
}

vector<vector<double> > readCoordinatesFromFile(const string& inputFilePath)
{
    ifstream fileReader(inputFilePath);
    string lineHolder = "";
    getline(fileReader, lineHolder);
    unsigned int citiesCount = static_cast<unsigned int>(stoi(lineHolder));
    vector<vector<double> > allCoordinates(citiesCount, vector<double>());
    for(unsigned int i =0u; i<citiesCount; i++)
    {
        getline(fileReader, lineHolder);
        vector<double> coords(2, 0.0);
        vector<string> coordsString = extractCoordinatesStringRepresentations(lineHolder);
        coords[0u] = stod(coordsString[0u]);
        coords[1u] = stod(coordsString[1u]);
        allCoordinates[i] = coords;
    }
    fileReader.close();
    return allCoordinates;
}


int main()
{
    
    string inputFilePath = "ch130.tsp";
    vector<vector<double> > coordinates = readCoordinatesFromFile(inputFilePath);
    unsigned int citiesCount = coordinates.size();
    vector <vector<double> > graphMesh(citiesCount, vector<double>(citiesCount, 0.0));
    buildMesh(graphMesh, coordinates);
    unsigned int runsCount = 5;
    double averagePathLengthGreedyAlgorithm = 0.0;
    double averagePathLengthGreedyRandomAlgorithm =0.0;
    default_random_engine engine;
    uniform_int_distribution<unsigned int> distribution(0u, citiesCount -1);
    unsigned int startingCityGreedyRandom = distribution(engine);
    for(unsigned int runIdx = 0u; runIdx < runsCount; runIdx++)
    {
        unsigned int startingCityGreedy = distribution(engine);
        startingCityGreedyRandom = startingCityGreedy;
        vector<unsigned int> greedyAlgorithmSolutionPath = greedyMethodShortestPathInTSP(graphMesh, false, startingCityGreedy);
        vector<unsigned int> greedyRandomAlgorithmSolutionPath = greedyRandomMethodShortestPathInTSP(graphMesh, false, startingCityGreedyRandom);
        double greedyPathLength = getTotalPathLength(graphMesh, greedyAlgorithmSolutionPath);
        double greedyRandomPathLength = getTotalPathLength(graphMesh, greedyRandomAlgorithmSolutionPath);
        averagePathLengthGreedyAlgorithm += greedyPathLength;
        averagePathLengthGreedyRandomAlgorithm += greedyRandomPathLength;
    }
    averagePathLengthGreedyAlgorithm /= static_cast<double>(runsCount);
    averagePathLengthGreedyRandomAlgorithm /= static_cast<double>(runsCount);
    cout<<"Average greedy path length: "<<averagePathLengthGreedyAlgorithm<<std::endl;
    cout<<"Average greedy random path length: "<<averagePathLengthGreedyRandomAlgorithm<<std::endl;
    return 0;
}
