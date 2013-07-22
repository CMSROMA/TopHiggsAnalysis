#! /bin/sh

make LatinosAnalyzer

VERSION="V04"
FOLDER="results_"${VERSION}

FINALFOLDER=/cmsrm/pc24_2/jorda/data/${FOLDER}/

echo "Adding weights for eee datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetEEE.root 1 0 0 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetEEE.root 1 0 0 
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetEEE.root 1 0 0 

echo "Adding weights for mmm datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetMMM.root 1 0 1 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetMMM.root 1 0 1
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetMMM.root 1 0 1 

echo "Adding weights for eem datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetEEM.root 1 0 2 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetEEM.root 1 0 2
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetEEM.root 1 0 2

echo "Adding weights for mme datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetMME.root 1 0 3 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetMME.root 1 0 3
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetMME.root 1 0 3

echo "Adding weights for ee datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetEE.root 1 0 4 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetEE.root 1 0 4
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetEE.root 1 0 4

echo "Adding weights for mm datasets..."

./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleElectron-Run2012-datasetMM.root 1 0 5 
./LatinosAnalyzer ${FINALFOLDER}/Data/DoubleMu-Run2012-datasetMM.root 1 0 5
./LatinosAnalyzer ${FINALFOLDER}/Data/MuEG-Run2012-datasetMM.root 1 0 5

echo "done weighting."
