import ROOT
import math, os,re,copy
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator

class jetRecalib(Module):
    def __init__(self,  globalTag, jetType = "AK4PFchs", METBranchName="MET", unclEnThreshold=15):

        self.metBranchName = METBranchName
        self.unclEnThreshold = unclEnThreshold
        if "AK4" in jetType : 
            self.jetBranchName = "Jet"
        elif "AK8" in jetType :
            self.jetBranchName = "FatJet"
            self.subJetBranchName = "SubJet"
        else:
            raise ValueError("ERROR: Invalid jet type = '%s'!" % jetType)
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        self.lenVar = "n" + self.jetBranchName
        # To do : change to real values
        self.jmsVals = [1.00, 0.99, 1.01]
        

        self.jesInputFilePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/jme/"

        self.jetReCalibrator = JetReCalibrator(globalTag, jetType , True, self.jesInputFilePath, calculateSeparateCorrections = False, calculateType1METCorrection  = False)
	
        # load libraries for accessing JES scale factors and uncertainties from txt files
        for library in [ "libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools" ]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '%s'" % library.replace("lib", ""))
                ROOT.gSystem.Load(library)

    def beginJob(self):
	pass

    def endJob(self):
	pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("%s_pt_nom" % self.jetBranchName, "F", lenVar=self.lenVar)
        self.out.branch("%s_pt_nom"%self.metBranchName , "F")
        self.out.branch("%s_phi_nom"%self.metBranchName, "F")
            
                        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets = Collection(event, self.jetBranchName )
        met = Object(event, self.metBranchName) 
        rawmet = Object(event, "RawMET") 

        jets_pt_nom = []
        ( met_px,         met_py        ) = ( met.pt*math.cos(met.phi),         met.pt*math.sin(met.phi) )
        ( raw_met_px,     raw_met_py    ) = ( rawmet.pt*math.cos(rawmet.phi),   rawmet.pt*math.sin(rawmet.phi) )

        rho = getattr(event, self.rhoBranchName)

        met_shift_x, met_shift_y = 0., 0.

        for jet in jets:
            # copy the jet_pt from the tuple
            jet_pt          = jet.pt
            # get the raw jet pt
            jet_pt_raw      = jet.pt*(1 - jet.rawFactor)
            # get the corrected jet pt
            jet_pt_nom      = self.jetReCalibrator.correct(jet,rho)
            jet.pt          = jet_pt_raw * ( 1 - jet.muEF )
            jet.rawFactor   = 0 # set raw factor to zero for jetReCalibrator tool which otherwise applies the rawFactor again
            jet_pt_noMu     = self.jetReCalibrator.correct(jet,rho) if jet.pt > 15 else jet_pt_raw*( 1 - jet.muEF ) # only correct the non-mu fraction of the jet if it's above 15 GeV, otherwise take the raw pt

            # only correct the non-muon fraction of the jet for T1 MET
            jet_pt_T1   = jet_pt_noMu + jet.muEF*jet_pt_raw

            if jet_pt_nom < 0.0:
                jet_pt_nom *= -1.0
            jets_pt_nom    .append(jet_pt_nom)

            print jet_pt, jet_pt_raw, jet_pt_noMu, jet_pt_T1, jet.muEF

            # only use jets with pt>15 GeV and EMF < 0.9 for T1 MET
            if jet_pt_T1 > self.unclEnThreshold and (jet.neEmEF+jet.chEmEF) < 0.9:
                jet_cosPhi = math.cos(jet.phi)
                jet_sinPhi = math.sin(jet.phi)
                if not ( self.metBranchName == 'METFixEE2017' and 2.65<abs(jet.eta)<3.14 and jet.pt*(1 - jet.rawFactor) < 50):

                    met_shift_x += (jet_pt_raw - jet_pt_T1) * jet_cosPhi
                    met_shift_y += (jet_pt_raw - jet_pt_T1) * jet_sinPhi

        met_px_nom = raw_met_px + met_shift_x
        met_py_nom = raw_met_py + met_shift_y

        self.out.fillBranch("%s_pt_nom" % self.jetBranchName, jets_pt_nom)
        self.out.fillBranch("%s_pt_nom" % self.metBranchName, math.sqrt(met_px_nom**2 + met_py_nom**2))
        self.out.fillBranch("%s_phi_nom" % self.metBranchName, math.atan2(met_py_nom, met_px_nom))        

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

jetRecalib2017B = lambda : jetRecalib("Fall17_17Nov2017B_V6_DATA")
jetRecalib2017C = lambda : jetRecalib("Fall17_17Nov2017C_V6_DATA")
jetRecalib2017D = lambda : jetRecalib("Fall17_17Nov2017D_V6_DATA")
jetRecalib2017E = lambda : jetRecalib("Fall17_17Nov2017E_V6_DATA")
jetRecalib2017F = lambda : jetRecalib("Fall17_17Nov2017F_V6_DATA")
