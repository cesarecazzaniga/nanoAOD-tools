import ROOT
import os
import numpy as np
import math
import itertools
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module



class METminProducer(Module):

    def __init__(self, isData=False, calcVariations=False):
        self.isData         = isData
        self.calcVariations = calcVariations
        if not isData and calcVariations:
            self.varList = ['', '_jesTotalUp', '_jesTotalDown', '_unclustEnUp', '_unclustEnDown']
        else:
            self.varList = ['']

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("MET_pt_min", "F")
        if self.calcVariations:
            self.out.branch("MET_pt_min_jesTotalUp", "F")
            self.out.branch("MET_pt_min_jesTotalDown", "F")
            self.out.branch("MET_pt_min_unclustEnUp", "F")
            self.out.branch("MET_pt_min_unclustEnDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        jets        = Collection(event, "Jet")
        for var in self.varList:
            MET_pt_var  = "MET_pt%s"%var
            MET_phi_var = "MET_phi%s"%var
            if not var.startswith('_unclust'):
                Jet_pt_var  = "Jet_pt%s"%var
            else:
                Jet_pt_var  = "Jet_pt"

            MET_min = []
            EE_jets = []
            MET_x = getattr(event, MET_pt_var) * math.cos(getattr(event, MET_phi_var))
            MET_y = getattr(event, MET_pt_var) * math.sin(getattr(event, MET_phi_var))
            
            if event.PV_npvsGood > 150:
                # there are problems with some weird events with super high nPV and nJet
                self.out.fillBranch("MET_pt_min%s"%var, float('NaN'))
                return False

            for i in range(event.nJet):
                # should we cut on pt>15?
                if 2.65<abs(event.Jet_eta[i])<3.14:
                    EE_jet = {'pt': getattr(event, Jet_pt_var)[i], 'phi':event.Jet_phi[i], 'eta':event.Jet_eta[i], 'neEmEF':event.Jet_neEmEF[i]}
                    EE_jets.append(EE_jet)
                    Jet_x = EE_jet['pt'] * math.cos(EE_jet['phi'])
                    Jet_y = EE_jet['pt'] * math.sin(EE_jet['phi'])
                    
                    if EE_jet['neEmEF'] > 0:
                        alpha_j = 1 + (MET_x*Jet_x + MET_y*Jet_y)/(EE_jet['neEmEF'] * (Jet_x**2 + Jet_y**2))
                    else:
                        alpha_j = 1

                    alpha_j     = max(min(alpha_j,1),0)
                    MET_min_j   = math.sqrt((MET_x - (alpha_j-1) * EE_jet['neEmEF'] * Jet_x)**2 + (MET_y - (alpha_j-1) * EE_jet['neEmEF'] * Jet_y)**2)
                    MET_min.append( MET_min_j )

            nEE_jets = min(len(EE_jets), 6) # build pseudo-jets out of combinations of maximum 6 jets
            for i in range(2,nEE_jets+1):
                combinations = itertools.combinations(EE_jets, i)
                for comb in combinations:
                    pseudoJet_x = 0.
                    pseudoJet_y = 0.
                    pseudoJet_neEmEF = 0.
                    pseudoJet_sumPt  = 0.

                    for jet in comb:
                        Jet_x = jet['pt'] * math.cos(jet['phi'])
                        Jet_y = jet['pt'] * math.sin(jet['phi'])

                        pseudoJet_x += Jet_x
                        pseudoJet_y += Jet_y

                        pseudoJet_sumPt     += jet['pt']
                        pseudoJet_neEmEF    += jet['neEmEF'] * jet['pt']
                    
                    if pseudoJet_neEmEF>0:
                        pseudoJet_neEmEF    = pseudoJet_neEmEF/pseudoJet_sumPt
                        alpha_j             = 1 + (MET_x*pseudoJet_x + MET_y*pseudoJet_y)/(pseudoJet_neEmEF*(pseudoJet_x**2 + pseudoJet_y**2))
                    else:
                        alpha_j = 1

                    alpha_j     = max(min(alpha_j,1),0)
                    MET_min_j   = math.sqrt((MET_x - (alpha_j-1) * pseudoJet_neEmEF * pseudoJet_x)**2 + (MET_y - (alpha_j-1) * pseudoJet_neEmEF * pseudoJet_y)**2)
                    MET_min.append( MET_min_j )

            if len(MET_min)>0:
                MET_pt_min = min(MET_min)
            else:
                MET_pt_min = getattr(event, MET_pt_var)
            del EE_jets, MET_min
            self.out.fillBranch("MET_pt_min%s"%var, MET_pt_min)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tools = lambda : METminProducer( )

