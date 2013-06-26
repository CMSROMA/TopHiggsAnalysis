#! /bin/sh


mkdir -p ../results/Data

make LatinosAnalyzer

echo "Adding weights for eee datasets..."

./LatinosAnalyzer ../results2/Data/DoubleElectron-Run2012-datasetEEE.root 1 1 0 
./LatinosAnalyzer ../results2/Data/DoubleMu-Run2012-datasetEEE.root 1 1 0 
./LatinosAnalyzer ../results2/Data/MuEG-Run2012-datasetEEE.root 1 1 0 

echo "Adding weights for mmm datasets..."

./LatinosAnalyzer ../results2/Data/DoubleElectron-Run2012-datasetMMM.root 1 1 1 
./LatinosAnalyzer ../results2/Data/DoubleMu-Run2012-datasetMMM.root 1 1 1
./LatinosAnalyzer ../results2/Data/MuEG-Run2012-datasetMMM.root 1 1 1 

echo "Adding weights for eem datasets..."

./LatinosAnalyzer ../results2/Data/DoubleElectron-Run2012-datasetEEM.root 1 1 2 
./LatinosAnalyzer ../results2/Data/DoubleMu-Run2012-datasetEEM.root 1 1 2
./LatinosAnalyzer ../results2/Data/MuEG-Run2012-datasetEEM.root 1 1 2

echo "Adding weights for mme datasets..."

./LatinosAnalyzer ../results2/Data/DoubleElectron-Run2012-datasetMME.root 1 1 3 
./LatinosAnalyzer ../results2/Data/DoubleMu-Run2012-datasetMME.root 1 1 3
./LatinosAnalyzer ../results2/Data/MuEG-Run2012-datasetMME.root 1 1 3

echo "done weighting."
