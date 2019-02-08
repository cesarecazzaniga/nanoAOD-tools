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
        self.out.branch("Electron_isGoodElectron", "I", lenVar="nElectron")
        self.out.branch("Jet_cleanmaskMETSigRec", "I", lenVar="nJet")
        self.out.branch("Jet_cleanmaskMETSig", "I", lenVar="nJet")
        self.out.branch("dl_mass", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def isGoodMuon(self, muon):
        if muon.pt > 20 and muon.mediumId and abs(muon.eta)<2.4 and muon.pfRelIso03_all < 0.2:
            return True
        else:
            return False

    def isGoodElectron(self, electron):
        if electron.pt > 20 and electron.cutBased >= 4 and abs(electron.eta) < 2.4 and electron.pfRelIso03_all < 0.15 and electron.sip3d < 4.0 and abs(electron.dxy) < 0.05 and abs(electron.dz) < 0.1:
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

    def deltaPhi(self, phi1, phi2):
        dphi = phi2-phi1
        if  dphi > math.pi:
            dphi -= 2.0*math.pi
        if dphi <= -math.pi:
            dphi += 2.0*math.pi
        return abs(dphi)

    def deltaR2(self, l1, l2):
        return self.deltaPhi(l1.phi, l2.phi)**2 + (l1.eta - l2.eta)**2

    def deltaR(self, l1, l2):
        return math.sqrt(self.deltaR2(l1,l2))

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons       = Collection(event, "Muon")
        electrons   = Collection(event, "Electron")
        photons     = Collection(event, "Photon")
        jets        = Collection(event, "Jet")

        cleanMaskV = []
        for j in jets:
            cleanMask = 1
            for p in photons:
                if p.pt > 10:
                    if self.deltaR(j, p) < 0.4:
                        cleanMask = 0
            for coll in [electrons, muons]:
                for p in coll:
                    if p.pt > 10:
                        if self.deltaR(j, p) < 0.4 and j.cleanmask==0: #recover some jets. need to check
                            cleanMask = 0
            
            cleanMaskV.append(cleanMask)

        cleanMaskV2 = []
        for j in jets:
            cleanMask = 1
            for coll in [electrons, muons, photons]:
                for p in coll:
                    if p.pt > 10:
                        if self.deltaR(j, p) < 0.4: #classical approach
                            cleanMask = 0

            cleanMaskV2.append(cleanMask)

        goodMuons       = []
        goodElectrons   = []
        isGoodMuon      = []
        isGoodElectron  = []
        nGoodElectrons  = 0
        
        for m in muons:
            if self.isGoodMuon(m):
                isGoodMuon += [1]
                goodMuons += [m]
            else:
                isGoodMuon += [0]


        for e in electrons:
            if self.isGoodElectron(e):
                isGoodElectron += [1]
                goodElectrons += [e]
            else:
                isGoodElectron += [0]

        dl_mass_muons = -1.
        if len(goodMuons) > 1:
            dl_mass_muons = self.findBestZcandidate(goodMuons)[0]
        dl_mass_electrons = -1
        if len(goodElectrons) > 1:
            dl_mass_electrons = self.findBestZcandidate(goodElectrons)[0]

        if abs(dl_mass_muons-91.2)<abs(dl_mass_electrons-91.2):
            dl_mass = dl_mass_muons
        else:
            dl_mass = dl_mass_electrons

        self.out.fillBranch("Muon_isGoodMuon", isGoodMuon)
        self.out.fillBranch("Electron_isGoodElectron", isGoodElectron)
        self.out.fillBranch("Jet_cleanmaskMETSigRec",  cleanMaskV)
        self.out.fillBranch("Jet_cleanmaskMETSig",  cleanMaskV2)
        self.out.fillBranch("dl_mass", dl_mass)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

tools = lambda : METSigTools( )

