import ROOT
import os
import numpy as np
import math
import itertools
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module



class METSigTools(Module):

    def __init__(self):
        pass

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Muon_isGoodMuon", "I", lenVar="nMuon")
        self.out.branch("dl_mass", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodMuon(self, muon):
        if muon.pt > 20 and muon.mediumId and abs(muon.eta)<2.4 and muon.pfRelIso03_all < 0.2:
            return True
        else:
            return False

    def findBestZcandidate(self, leptons):
        mZ = 91.2
        inds = range(len(leptons))
        vecs = [ ROOT.TLorentzVector() for i in inds ]
        for i, v in enumerate(vecs):
            v.SetPtEtaPhiM(leptons[i].pt, leptons[i].eta, leptons[i].phi, 0.)
        dlMasses = [((vecs[comb[0]] + vecs[comb[1]]).M(), comb[0], comb[1])  for comb in itertools.combinations(inds, 2) if leptons[comb[0]].pdgId*leptons[comb[1]].pdgId < 0 and abs(leptons[comb[0]].pdgId) == abs(leptons[comb[1]].pdgId) ]
        return min(dlMasses, key=lambda (m,i1,i2):abs(m-mZ)) if len(dlMasses)>0 else (float('nan'), -1, -1)

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = Collection(event, "Muon")

        goodMuons = []
        isGoodMuon = []
        nGoodElectrons = 0
        
        for m in muons:
            if self.isGoodMuon(m):
                isGoodMuon += [1]
                goodMuons += [m]
            else:
                isGoodMuon += [0]

        dl_mass = -1.
        if len(goodMuons) > 1:
            dl_mass = self.findBestZcandidate(goodMuons)[0]
        

        self.out.fillBranch("Muon_isGoodMuon", isGoodMuon)
        self.out.fillBranch("dl_mass", dl_mass)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tools = lambda : METSigTools( )

