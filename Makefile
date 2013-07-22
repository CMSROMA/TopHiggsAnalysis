ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lMinuit2
GLIBS          = $(filter-out -lNew, $(NGLIBS))

# For MT2:
# MT2LIB_INC_DIR := mT2
# MT2LIB_LIB_DIR := mT2/src
# MT2LIB_CPPFLAGS := -I $(MT2LIB_INC_DIR)
# MT2LIB_LDFLAGS  := $(MT2LIB_LIB_DIR)/libMt2.so
# CXXFLAGS += $(MT2LIB_CPPFLAGS)
# LDFLAGS  += $(MT2LIB_LDFLAGS)

INCLUDEDIR        = ./
INCLUDEDIRCOMMON  = ../
INCLUDEDIRHIGGS   = ../HiggsAnalysisTools
CXX	          += -I$(INCLUDEDIR) -I$(INCLUDEDIRCOMMON) -I$(INCLUDEDIRHIGGS) -I.
OUTLIB	          = $(INCLUDEDIR)/lib/
OUTLIBCOMMON      = $(INCLUDEDIRCOMMON)/CommonTools/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

$(OUTLIB)TopHiggsBase.o: $(INCLUDEDIR)/src/TopHiggsBase.C \
	$(INCLUDEDIR)/src/TopHiggs.cc \
	$(INCLUDEDIR)/src/RochCor2012.cc \
	$(INCLUDEDIRHIGGS)/src/HiggsSelection.cc \
	$(INCLUDEDIRHIGGS)/src/HiggsMLSelection.cc \
	$(INCLUDEDIR)/src/tHSelection.cc \
	$(INCLUDEDIRHIGGS)/src/HiggsEleIdOptimToyMC.cc \
	$(INCLUDEDIRHIGGS)/src/RedEleIDOptimTree.cc \
	$(INCLUDEDIRHIGGS)/src/RedLikeOptimTree.cc \
	$(INCLUDEDIRHIGGS)/src/HiggsIsolationOptimToyMC.cc \
	$(INCLUDEDIRHIGGS)/src/RedIsolationOptimTree.cc \
	$(INCLUDEDIRHIGGS)/src/RedEleIDTree.cc \
	$(INCLUDEDIRHIGGS)/src/ZplusJetsSelection.cc \
	$(INCLUDEDIRHIGGS)/src/LeptonPlusFakeSelection.cc \
	$(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection.cc \
	$(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_fullEE.cc \
	$(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_fullME.cc \
	$(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_ME.cc \
	$(INCLUDEDIRHIGGS)/src/HiggsVertexing.cpp \
	$(INCLUDEDIRHIGGS)/src/VertexTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TopHiggsBase.o $<

$(OUTLIBCOMMON)Conditions.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Conditions.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Conditions.o $<

$(OUTLIBCOMMON)PUWeight.o: $(INCLUDEDIRCOMMON)/CommonTools/src/PUWeight.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)PUWeight.o $<

$(OUTLIBCOMMON)Utils.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Utils.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIBCOMMON)Utils.o $<

$(OUTLIBCOMMON)Counters.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Counters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Counters.o $<

$(OUTLIBCOMMON)Selection.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Selection.o $<

$(OUTLIBCOMMON)EfficiencyEvaluator.o: $(INCLUDEDIRCOMMON)/CommonTools/src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)EfficiencyEvaluator.o $<

$(OUTLIBCOMMON)Monitor.o: $(INCLUDEDIRCOMMON)/CommonTools/src/Monitor.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)Monitor.o $<

$(OUTLIBCOMMON)SprDataFiller.o: $(INCLUDEDIRCOMMON)/CommonTools/src/SprDataFiller.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)SprDataFiller.o $<

$(OUTLIBCOMMON)TriggerMask.o: $(INCLUDEDIRCOMMON)/CommonTools/src/TriggerMask.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIBCOMMON)TriggerMask.o $<

$(OUTLIB)EcalCleaner.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/EcalCleaner.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)EcalCleaner.o $<

$(OUTLIB)CutBasedEleIDSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CutBasedEleIDSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)CutBasedEleIDSelector.o $<

$(OUTLIB)CiCBasedEleSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/CiCBasedEleSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)CiCBasedEleSelector.o $<

$(OUTLIB)TopHiggs.o: $(INCLUDEDIR)/src/TopHiggs.cc $(OUTLIB)JetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)TopHiggs.o $<

$(OUTLIB)JetCorrectorParameters.o: $(INCLUDEDIRHIGGS)/src/JetCorrectorParameters.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectorParameters.o $<

$(OUTLIB)SimpleJetCorrectionUncertainty.o: $(INCLUDEDIRHIGGS)/src/SimpleJetCorrectionUncertainty.cc \
	$(OUTLIB)JetCorrectorParameters.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)SimpleJetCorrectionUncertainty.o $<

$(OUTLIB)JetCorrectionUncertainty.o: $(INCLUDEDIRHIGGS)/src/JetCorrectionUncertainty.cc \
	$(OUTLIB)SimpleJetCorrectionUncertainty.o
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)JetCorrectionUncertainty.o $<

$(OUTLIB)RochCor2012.o: $(INCLUDEDIR)/src/RochCor2012.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RochCor2012.o $<

$(OUTLIB)HiggsSelection.o: $(INCLUDEDIRHIGGS)/src/HiggsSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsSelection.o $<

$(OUTLIB)HiggsMLSelection.o: $(INCLUDEDIRHIGGS)/src/HiggsMLSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsMLSelection.o $<

$(OUTLIB)HiggsMLSelection.o: $(INCLUDEDIRHIGGS)/src/tHSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)tHSelection.o $<

$(OUTLIB)HiggsEleIdOptimToyMC.o: $(INCLUDEDIRHIGGS)/src/HiggsEleIdOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsEleIdOptimToyMC.o $<

$(OUTLIB)RedEleIDOptimTree.o: $(INCLUDEDIRHIGGS)/src/RedEleIDOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDOptimTree.o $<

$(OUTLIB)RedLikeOptimTree.o: $(INCLUDEDIRHIGGS)/src/RedLikeOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedLikeOptimTree.o $<

$(OUTLIB)HiggsIsolationOptimToyMC.o: $(INCLUDEDIRHIGGS)/src/HiggsIsolationOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsIsolationOptimToyMC.o $<

$(OUTLIB)HiggsKinematicsOptimToyMC.o: $(INCLUDEDIRHIGGS)/src/HiggsKinematicsOptimToyMC.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsKinematicsOptimToyMC.o $<

$(OUTLIB)RedIsolationOptimTree.o: $(INCLUDEDIRHIGGS)/src/RedIsolationOptimTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedIsolationOptimTree.o $<

$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIRHIGGS)/src/RedEleIDTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<

$(OUTLIB)kFactorEvaluator.o: $(INCLUDEDIRHIGGS)/src/kFactorEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)kFactorEvaluator.o $<

$(OUTLIB)ElectronIDMVA.o: $(INCLUDEDIRHIGGS)/src/ElectronIDMVA.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronIDMVA.o $<

$(OUTLIB)IDForBsMVA.o: $(INCLUDEDIRHIGGS)/src/IDForBsMVA.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)IDForBsMVA.o $<

$(OUTLIB)tH_leptonMcCorrections.o: $(INCLUDEDIRHIGGS)/src/tH_leptonMcCorrections.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)tH_leptonMcCorrections.o $<

$(OUTLIB)GetDYMVA.o: $(INCLUDEDIRHIGGS)/src/GetDYMVA.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GetDYMVA.o $<

$(OUTLIB)RedHiggsTree.o: $(INCLUDEDIRHIGGS)/src/RedHiggsTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedHiggsTree.o $<

$(OUTLIB)RedTopHiggsTree.o: $(INCLUDEDIR)/src/RedTopHiggsTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedTopHiggsTree.o $<

$(OUTLIB)RedTriggerTree.o: $(INCLUDEDIRHIGGS)/src/RedTriggerTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedTriggerTree.o $<

$(OUTLIB)CommonHiggsPreselector.o: $(INCLUDEDIRHIGGS)/src/CommonHiggsPreselector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CommonHiggsPreselector.o $<

$(OUTLIB)CutBasedHiggsSelector.o: $(INCLUDEDIRHIGGS)/src/CutBasedHiggsSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedHiggsSelector.o $<

$(OUTLIB)CutBasedTopHiggsSelector.o: $(INCLUDEDIR)/src/CutBasedTopHiggsSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)CutBasedTopHiggsSelector.o $<

$(OUTLIB)ZplusJetsSelection.o: $(INCLUDEDIRHIGGS)/src/ZplusJetsSelection.C
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ZplusJetsSelection.o $<

$(OUTLIB)LeptonPlusFakeSelection.o: $(INCLUDEDIRHIGGS)/src/LeptonPlusFakeSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeSelection.o $<

$(OUTLIB)LeptonPlusFakeMLSelection.o: $(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeMLSelection.o $<

$(OUTLIB)LeptonPlusFakeMLSelection_fullEE.o: $(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_fullEE.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeMLSelection_fullEE.o $<

$(OUTLIB)LeptonPlusFakeMLSelection_fullME.o: $(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_fullME.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeMLSelection_fullME.o $<

$(OUTLIB)LeptonPlusFakeMLSelection_ME.o: $(INCLUDEDIRHIGGS)/src/LeptonPlusFakeMLSelection_ME.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)LeptonPlusFakeMLSelection_ME.o $<

$(OUTLIB)HiggsVertexing.o: $(INCLUDEDIRHIGGS)/src/HiggsVertexing.cpp
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)HiggsVertexing.o $<

$(OUTLIB)VertexTree.o: $(INCLUDEDIRHIGGS)/src/VertexTree.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)VertexTree.o $<

$(OUTLIB)ElectronBestCandidateSelector.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronBestCandidateSelector.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)ElectronBestCandidateSelector.o $<

$(OUTLIB)BestLeptonSelectorWjets.o: $(INCLUDEDIRHIGGS)/src/BestLeptonSelectorWjets.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BestLeptonSelectorWjets.o $<

$(OUTLIB)LikelihoodPdf.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdf.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)LikelihoodPdf.o $<

$(OUTLIB)LikelihoodSpecies.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodSpecies.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)LikelihoodSpecies.o $<

$(OUTLIB)LikelihoodPdfProduct.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/LikelihoodPdfProduct.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)LikelihoodPdfProduct.o $<

$(OUTLIB)ElectronLikelihood.o: $(INCLUDEDIRCOMMON)/EgammaAnalysisTools/src/ElectronLikelihood.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIRCOMMON) -o $(OUTLIB)ElectronLikelihood.o $<

# no more used ------------------------------------
#$(OUTLIB)ElectronID.o: $(INCLUDEDIR)/src/ElectronID.C
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)ElectronID.o $<
#$(OUTLIB)RedEWKTree.o: $(INCLUDEDIR)/src/RedEWKTree.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEWKTree.o $<
#$(OUTLIB)RedEleIDTree.o: $(INCLUDEDIR)/src/RedEleIDTree.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)RedEleIDTree.o $<
#$(OUTLIB)eleID_Higgs_Studies.o: $(INCLUDEDIR)/src/eleID_Higgs_Studies.cc
#	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)eleID_Higgs_Studies.o $<

#----------------------------------------------------#

# ==================== HiggsApp =============================================
TopHiggsApp:  $(INCLUDEDIR)/src/TopHiggsApp.C \
	$(OUTLIB)TopHiggsBase.o \
	$(OUTLIB)TopHiggs.o \
	$(OUTLIB)RochCor2012.o \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)PUWeight.o \
	$(OUTLIBCOMMON)Selection.o \
	$(OUTLIBCOMMON)EfficiencyEvaluator.o \
	$(OUTLIBCOMMON)Counters.o \
	$(OUTLIBCOMMON)Monitor.o \
	$(OUTLIBCOMMON)SprDataFiller.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIB)kFactorEvaluator.o \
	$(OUTLIB)ElectronIDMVA.o \
	$(OUTLIB)IDForBsMVA.o \
	$(OUTLIB)tH_leptonMcCorrections.o \
	$(OUTLIB)GetDYMVA.o \
	$(OUTLIB)RedHiggsTree.o \
	$(OUTLIB)RedTopHiggsTree.o \
	$(OUTLIB)RedTriggerTree.o \
	$(OUTLIB)RedEleIDOptimTree.o \
	$(OUTLIB)RedLikeOptimTree.o \
	$(OUTLIB)RedIsolationOptimTree.o \
	$(OUTLIB)RedEleIDTree.o \
	$(OUTLIB)CutBasedEleIDSelector.o \
	$(OUTLIB)CiCBasedEleSelector.o \
	$(OUTLIB)EcalCleaner.o \
	$(OUTLIB)CommonHiggsPreselector.o \
	$(OUTLIB)CutBasedHiggsSelector.o \
	$(OUTLIB)CutBasedTopHiggsSelector.o \
	$(OUTLIB)VertexTree.o \
	$(OUTLIB)LikelihoodPdf.o \
	$(OUTLIB)LikelihoodSpecies.o \
	$(OUTLIB)LikelihoodPdfProduct.o \
	$(OUTLIB)ElectronLikelihood.o
	$(CXX) $(CXXFLAGS) -ldl -o TopHiggsApp $(OUTLIB)/*.o $(OUTLIBCOMMON)/*o $(GLIBS) $(LDFLAGS) $ $<
TopHiggsApp.clean:
	rm -f TopHiggsApp

# ==================== eleID =============================================
eleID_Higgs_Studies:  $(INCLUDEDIR)eleID_Higgs_Studies.cpp $(OUTLIB)RedEleIDTree.o 
	$(CXX) $(CXXFLAGS) -o eleID_Higgs_Studies $(OUTLIB)/*.o $(GLIBS) $ $<
#eleID_Higgs_Studies.clean:
#	rm -f eleID_Higgs_Studies

eleIDtableToy:  $(INCLUDEDIRHIGGS)/src/eleIDtableToy.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtableToy $(GLIBS) $ $<
eleIDtoyPlot_input:  $(INCLUDEDIRHIGGS)/src/eleIDtoyPlot_input.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtoyPlot_input $(GLIBS) $ $<
eleIDtableToy.clean:
	rm -f eleIDtableToy
eleIDtoyPlot_input.clean:
	rm -f eleIDtoyPlot_input
eleIDtableLike:  $(INCLUDEDIRHIGGS)/src/eleIDtableLike.cpp
	$(CXX) $(CXXFLAGS) -o eleIDtableLike $(GLIBS) $ $<
eleIDtableLike.clean:
	rm -f eleIDtableLike
isolationTableToy:  $(INCLUDEDIRHIGGS)/src/isolationTableToy.cpp
	$(CXX) $(CXXFLAGS) -o isolationTableToy $(GLIBS) $ $<
isolationTableToy.clean:
	rm -f isolationTableToy
kinematicsTableToy:  $(INCLUDEDIRHIGGS)/src/kinematicsTableToy.cpp
	$(CXX) $(CXXFLAGS) -o kinematicsTableToy $(GLIBS) $ $<
kinematicsTableToy.clean:
	rm -f kinematicsTableToy
vtxAndIsoOptim:  $(INCLUDEDIRHIGGS)/src/vtxAndIsoOptim.cpp
	$(CXX) $(CXXFLAGS) -o vtxAndIsoOptim $(GLIBS) $ $<
vtxAndIsoOptim.clean:
	rm -f vtxAndIsoOptim

# ================= other ===================
BestLeptonApp: $(INCLUDEDIRHIGGS)/src/BestLeptonApp.C \
	$(OUTLIBCOMMON)Conditions.o \
	$(OUTLIBCOMMON)TriggerMask.o \
	$(OUTLIBCOMMON)Utils.o \
	$(OUTLIB)HiggsBase.o \
	$(OUTLIB)ElectronBestCandidateSelector.o \
	$(OUTLIB)BestLeptonSelectorWjets.o
	$(CXX) $(CXXFLAGS) -o BestLeptonApp $(OUTLIB)/*o $(OUTLIBCOMMON)/*o $(GLIBS) $ $<
BestLeptonApp.clean:
	rm -f $(OUTLIB)BestLeptonSelectorWjets.o
	rm -f BestLeptonApp

# ==================== reduced trees =============================================
#ReducedTree_HwwEleId:  $(INCLUDEDIR)ReducedTree_HwwEleId.cpp $(OUTLIB)RedEleIDTree.o 
#	$(CXX) $(CXXFLAGS) -o ReducedTree_HwwEleId $(OUTLIB)/*.o $(GLIBS) $ $<
#ReducedTree_HwwEleId.clean:
#	rm -f ReducedTree_HwwEleId

clean:
	rm -f $(OUTLIB)*.o $(OUTLIBCOMMON)*.o
	rm -f TopHiggsApp
# rm -f ReducedTree_HwwEleId
# rm -f eleID_Higgs_Studies
	rm -f eleIDtableToy
	rm -f eleIDtoyPlot_input
	rm -f eleIDtableLike
	rm -f isolationTableToy
	rm -f kinematicsTableToy
	rm -f vtxAndIsoOptim

all:  TopHiggsApp
