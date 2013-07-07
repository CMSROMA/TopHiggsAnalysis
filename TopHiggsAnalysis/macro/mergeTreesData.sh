#! /bin/sh 
# ./mergeTrees.sh expects results of the selection in subdirectories such as results/H120, results/bkg1, results/bkg2...
# it creates a merged root file in the self-created results/merged/

mkdir -p results_data/merged

echo "Now merging EEE datasets..."
hadd results_data/merged/dataset_DoubleElectron_eee.root results_data/testV1/DoubleElectron/*datasetEEE.root

echo "Now merging MMM datasets..."
hadd results_data/merged/dataset_DoubleMu_mmm.root results_data/testV1/DoubleMu/*datasetMMM.root

echo "Now merging EEM datasets..."
hadd results_data/merged/dataset_MuEG_eem.root results_data/testV1/MuEG/*datasetEEM.root
hadd results_data/merged/dataset_DoubleElectron_eem.root results_data/testV1/DoubleElectron/*datasetEEM.root

echo "Now merging MME datasets..."
hadd results_data/merged/dataset_MuEG_mme.root results_data/testV1/MuEG/*datasetMME.root
hadd results_data/merged/dataset_DoubleMu_mme.root results_data/testV1/DoubleMu/*datasetMME.root

echo "Now merging EE datasets..."
hadd results_data/merged/dataset_DoubleElectron_ee.root results_data/testV1/DoubleElectron/*datasetEE.root

echo "Now merging MM datasets..."
hadd results_data/merged/dataset_DoubleMu_mm.root results_data/testV1/DoubleMu/*datasetMM.root
