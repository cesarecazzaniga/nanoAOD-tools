from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetReCalibrator import JetReCalibrator
from PhysicsTools.NanoAODTools.postprocessing.tools import matchObjectCollection, matchObjectCollectionMultiple
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
import ROOT
import math, os,re,copy
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True


class jetRecalib(Module):
    def __init__(self,  globalTag, jetType = "AK4PFchs", METBranchName="MET", unclEnThreshold=15):

        self.metBranchName = METBranchName
        self.unclEnThreshold = unclEnThreshold
        if "AK4" in jetType : 
            self.jetBranchName = "Jet"
        elif "AK8" in jetType:
            self.jetBranchName = "FatJet"
            self.subJetBranchName = "SubJet"
        else:
            raise ValueError("ERROR: Invalid jet type = '%s'!" % jetType)
        self.rhoBranchName = "fixedGridRhoFastjetAll"
        self.lenVar = "n" + self.jetBranchName

        self.jesInputArchivePath = os.environ['CMSSW_BASE'] + \
            "/src/PhysicsTools/NanoAODTools/data/jme/"
        # Text files are now tarred so must extract first into temporary
        # directory (gets deleted during python memory management at script exit)
        self.jesArchive = tarfile.open(
            self.jesInputArchivePath + archive + ".tgz", "r:gz")
        self.jesInputFilePath = tempfile.mkdtemp()
        self.jesArchive.extractall(self.jesInputFilePath)

        self.jetReCalibrator = JetReCalibrator(
            globalTag,
            jetType,
            True,
            self.jesInputFilePath,
            calculateSeparateCorrections=False,
            calculateType1METCorrection=False)

        # load libraries for accessing JES scale factors and
        # uncertainties from txt files
        for library in [
                "libCondFormatsJetMETObjects", "libPhysicsToolsNanoAODTools"
        ]:
            if library not in ROOT.gSystem.GetLibraries():
                print("Load Library '%s'" % library.replace("lib", ""))
                ROOT.gSystem.Load(library)

        self.puppiCorrFile = ROOT.TFile.Open(
            os.environ['CMSSW_BASE'] +
            "/src/PhysicsTools/NanoAODTools/data/jme/puppiCorr.root")
        self.puppisd_corrGEN = self.puppiCorrFile.Get("puppiJECcorr_gen")
        self.puppisd_corrRECO_cen = self.puppiCorrFile.Get(
            "puppiJECcorr_reco_0eta1v3")
        self.puppisd_corrRECO_for = self.puppiCorrFile.Get(
            "puppiJECcorr_reco_1v3eta2v5")

    def beginJob(self):
        pass

    def endJob(self):
        shutil.rmtree(self.jesInputFilePath)

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

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail,
        go to next event)"""
        jets = Collection(event, self.jetBranchName)
        subJets = Collection(event, self.subJetBranchName)
        met = Object(event, "MET")

        jets_pt_raw = []
        jets_pt_nom = []
        ( met_px,         met_py        ) = ( met.pt*math.cos(met.phi),         met.pt*math.sin(met.phi) )
        ( raw_met_px,     raw_met_py    ) = ( rawmet.pt*math.cos(rawmet.phi),   rawmet.pt*math.sin(rawmet.phi) )

        rho = getattr(event, self.rhoBranchName)

        met_shift_x, met_shift_y = 0., 0.

        for jet in jets:
            # get the raw jet pt
            jet_pt_raw  = jet.pt*(1 - jet.rawFactor)
            # get the corrected jet pt
            jet_pt_nom  = self.jetReCalibrator.correct(jet,rho)

            # get the correction factor
            jec         = jet_pt_nom/jet_pt_raw

            # only correct the non-muon fraction of the jet for T1 MET
            jet_pt_T1   = jec*(1-jet.muEF)*jet_pt_raw + jet.muEF*jet_pt_raw

            if jet_pt_nom < 0.0:
                jet_pt_nom *= -1.0
            jets_pt_nom    .append(jet_pt_nom)

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


# define modules using the syntax 'name = lambda : constructor' to avoid
# having them loaded when not needed
jetRecalib2016BCD = lambda: jetRecalib("Summer16_07Aug2017BCD_V11_DATA",
                                       "Summer16_07Aug2017_V11_DATA")
jetRecalib2016EF = lambda: jetRecalib("Summer16_07Aug2017EF_V11_DATA",
                                      "Summer16_07Aug2017_V11_DATA")
jetRecalib2016GH = lambda: jetRecalib("Summer16_07Aug2017GH_V11_DATA",
                                      "Summer16_07Aug2017_V11_DATA")

jetRecalib2016BCDAK8Puppi = lambda: jetRecalib(
    "Summer16_07Aug2017BCD_V11_DATA",
    "Summer16_07Aug2017_V11_DATA",
    jetType="AK8PFPuppi")
jetRecalib2016EFAK8Puppi = lambda: jetRecalib("Summer16_07Aug2017EF_V11_DATA",
                                              "Summer16_07Aug2017_V11_DATA",
                                              jetType="AK8PFPuppi")
jetRecalib2016GHAK8Puppi = lambda: jetRecalib("Summer16_07Aug2017GH_V11_DATA",
                                              "Summer16_07Aug2017_V11_DATA",
                                              jetType="AK8PFPuppi")

jetRecalib2017B = lambda: jetRecalib("Fall17_17Nov2017B_V32_DATA",
                                     "Fall17_17Nov2017_V32_DATA")
jetRecalib2017C = lambda: jetRecalib("Fall17_17Nov2017C_V32_DATA",
                                     "Fall17_17Nov2017_V32_DATA")
jetRecalib2017DE = lambda: jetRecalib("Fall17_17Nov2017DE_V32_DATA",
                                      "Fall17_17Nov2017_V32_DATA")
jetRecalib2017F = lambda: jetRecalib("Fall17_17Nov2017F_V32_DATA",
                                     "Fall17_17Nov2017_V32_DATA")

jetRecalib2017BAK8Puppi = lambda: jetRecalib("Fall17_17Nov2017B_V32_DATA",
                                             "Fall17_17Nov2017_V32_DATA",
                                             jetType="AK8PFPuppi")
jetRecalib2017CAK8Puppi = lambda: jetRecalib("Fall17_17Nov2017C_V32_DATA",
                                             "Fall17_17Nov2017_V32_DATA",
                                             jetType="AK8PFPuppi")
jetRecalib2017DEAK8Puppi = lambda: jetRecalib("Fall17_17Nov2017DE_V32_DATA",
                                              "Fall17_17Nov2017_V32_DATA",
                                              jetType="AK8PFPuppi")
jetRecalib2017FAK8Puppi = lambda: jetRecalib("Fall17_17Nov2017F_V32_DATA",
                                             "Fall17_17Nov2017_V32_DATA",
                                             jetType="AK8PFPuppi")

jetRecalib2018A = lambda: jetRecalib(
    "Autumn18_RunA_V8_DATA", "Autumn18_V8_DATA", redoJEC=True)
jetRecalib2018B = lambda: jetRecalib(
    "Autumn18_RunB_V8_DATA", "Autumn18_V8_DATA", redoJEC=True)
jetRecalib2018C = lambda: jetRecalib(
    "Autumn18_RunC_V8_DATA", "Autumn18_V8_DATA", redoJEC=True)
jetRecalib2018D = lambda: jetRecalib(
    "Autumn18_RunD_V8_DATA", "Autumn18_V8_DATA", redoJEC=True)

jetRecalib2018AAK8Puppi = lambda: jetRecalib("Autumn18_RunA_V8_DATA",
                                             "Autumn18_V8_DATA",
                                             jetType="AK8PFPuppi",
                                             redoJEC=True)
jetRecalib2018BAK8Puppi = lambda: jetRecalib("Autumn18_RunB_V8_DATA",
                                             "Autumn18_V8_DATA",
                                             jetType="AK8PFPuppi",
                                             redoJEC=True)
jetRecalib2018CAK8Puppi = lambda: jetRecalib("Autumn18_RunC_V8_DATA",
                                             "Autumn18_V8_DATA",
                                             jetType="AK8PFPuppi",
                                             redoJEC=True)
jetRecalib2018DAK8Puppi = lambda: jetRecalib("Autumn18_RunD_V8_DATA",
                                             "Autumn18_V8_DATA",
                                             jetType="AK8PFPuppi",
                                             redoJEC=True)
