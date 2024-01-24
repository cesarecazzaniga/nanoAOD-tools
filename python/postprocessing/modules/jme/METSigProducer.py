import ROOT
import os
import numpy as np
import math, tarfile, tempfile
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module



class METSigProducer(Module):

    def __init__(self, JERera, parameters, METCollection="MET", useRecorr=True, calcVariations=False, jetThreshold=15., vetoEtaRegion=(10,10), pTdependent=True):
        jetCorrParam = ROOT.JetCorrectorParameters()

        self.pars               = parameters
        self.pTdependent        = pTdependent
        
        #here hack for UL2017, since the JER files are not yet available take the ones from 2018
        if JERera == "Summer19UL17_JRV2_MC" or JERera == "Summer19UL17_JRV2_DATA":
            if "MC" in JERera:
                JERera = "Summer19UL18_JRV2_MC"
            else:
                JERera = "Summer19UL18_JRV2_DATA"

        self.JERera             = JERera
        self.METCollection      = METCollection + "_T1"     #this is a hack to get the correct MET collection for UL
        self.useRecorr          = useRecorr
        self.calcVariations     = calcVariations
        self.jetThreshold       = jetThreshold
        self.JERdirectory       = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/jme/")

        # read jet energy scale (JES) uncertainties
        # (downloaded from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC )
        #self.jesInputArchivePath = os.environ['CMSSW_BASE'] + "/src/PhysicsTools/NanoAODTools/data/jme/"
        # Text files are now tarred so must extract first into temporary directory (gets deleted during python memory management at script exit)
        try:
            self.jerArchive = tarfile.open(self.JERdirectory+JERera+".tgz", "r:gz")# if not archive else tarfile.open(self.jesInputArchivePath+archive+".tgz", "r:gz")
        except:
            try:
                #try to open tar
                self.jerArchive = tarfile.open(self.JERdirectory+JERera+".tar", "r")
            except:
                #raise exception if neither works
                #print path 
                print(self.JERdirectory+JERera)
                raise Exception("Could not open JER archive")
        
        self.jerInputFilePath = tempfile.mkdtemp()
        self.jerArchive.extractall(self.jerInputFilePath)

        self.vetoEtaRegion      = vetoEtaRegion     


    def beginJob(self):
        self.JERdirectory   = os.path.expandvars(self.JERdirectory)
        self.res_pt         = ROOT.JME.JetResolution("%s/%s_PtResolution_AK4PFchs.txt"%(self.jerInputFilePath, self.JERera))
        self.res_phi        = ROOT.JME.JetResolution("%s/%s_PhiResolution_AK4PFchs.txt"%(self.jerInputFilePath, self.JERera))

    def endJob(self):
        del self.res_pt
        del self.res_phi
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("MET_significance", "F")
        if self.calcVariations:
            for var in ['_jesTotalUp', '_jesTotalDown', '_jer', '_jerUp', '_jerDown', '_unclustEnUp', '_unclustEnDown']:
                self.out.branch("MET_significance"+var, "F")
        #self.out.branch("MET_significance_nom", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def getBin(self, abseta):
        etabins = [0.8,1.3,1.9,2.5,100]
        for i, a in enumerate(etabins):
            if abseta < a:
                return int(i)
                break

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
        jets        = Collection(event, "Jet")
        electrons   = Collection(event, "Electron")
        muons       = Collection(event, "Muon")
        photons     = Collection(event, "Photon")
        met         = Object(event, self.METCollection)        #here there are the corrected MET objects with up and down variations
        metStd      = Object(event, "MET")                     #here there are the standard MET objects
        rho         = getattr(event, "fixedGridRhoFastjetAll")


        ## calculate MET significance for all variants of MET/jets
        if self.useRecorr: variations = ['_nom']
        else: variations = ['']

        if self.calcVariations:
            #variations += ['_jesTotalUp', '_jesTotalDown', '_jerUp', '_jer', '_jerDown', '_unclustEnUp', '_unclustEnDown']
            variations += ['_jesTotalUp', '_jesTotalDown', '_jerUp', '_jerDown', '_unclustEnUp', '_unclustEnDown']   #CZZ: FIX - _jer is not yet available for UL

        for var in variations:
            jetPtVar = 'pt' if not self.useRecorr else 'pt_nom'
            # no unclustered energy uncertainty for jets
            if var in ['_jesTotalUp', '_jesTotalDown', '_jerUp', '_jerDown']:
                jetPtVar = 'pt'+var
            metPtVar = 'pt'+var
            phiVar = 'phi'+var

            # clean against electrons, muons and photons with pt>10 GeV
            sumPtFromJets = 0
            cleanJets = []
            for j in jets:
                correctJER = j.corr_JER if var in ['_jer'] else 1
                clean = True
                for coll in [electrons,muons,photons]:
                    for l in coll:
                        if l.pt > 10 and self.deltaR(j, l) < 0.4:
                            clean = False
                if clean:
                    if not (self.vetoEtaRegion[0] < abs(j.eta) < self.vetoEtaRegion[1]):

                        if getattr(j, jetPtVar)*correctJER > self.jetThreshold:
                            cleanJets += [j]
                        else:
                            sumPtFromJets += getattr(j, jetPtVar)*correctJER

            # get the JER
            jet = ROOT.JME.JetParameters()
            for j in cleanJets:
                correctJER = j.corr_JER if var in ['_jer'] else 1
                jet.setJetEta(j.eta).setJetPt(getattr(j, jetPtVar)/correctJER).setRho(rho)
                j.dpt   = self.res_pt.getResolution(jet)
                j.dphi  = self.res_phi.getResolution(jet)

            cov_xx  = 0
            cov_xy  = 0
            cov_yy  = 0
            i = 0
            for j in cleanJets:
                index       = self.getBin(abs(j.eta))
                correctJER = j.corr_JER if var in ['_jer'] else 1
                jet_index   = 0 if getattr(j, jetPtVar)*correctJER < 40 else 1 # split into high/low pt jets

                cj = math.cos(j.phi)
                sj = math.sin(j.phi)
                
                dpt = self.pars[2*index + jet_index if self.pTdependent else index] * getattr(j, jetPtVar)*correctJER * j.dpt
                dph =                                                                 getattr(j, jetPtVar)*correctJER * j.dphi

                dpt *= dpt
                dph *= dph

                cov_xx += dpt*cj*cj + dph*sj*sj
                cov_xy += (dpt-dph)*cj*sj
                cov_yy += dph*cj*cj + dpt*sj*sj


                i += 1

            # unclustered energy - here applies the correction for unclustered energy (here is taking sumPtUnclustered from MET collection, since it is not available in MET_T1 collection)                                                                                                                                                 
            if hasattr(metStd, 'sumPt'):                                                                                                                                 
                if 'unclust' in var:
                    if 'Up' in var:
                        totalSumPt = metStd.sumPt + metStd.sumPtUnclustEnDeltaUp + sumPtFromJets
                    else:
                        totalSumPt = metStd.sumPt - metStd.sumPtUnclustEnDeltaUp + sumPtFromJets
                else:
                    totalSumPt = metStd.sumPt + sumPtFromJets

            elif hasattr(metStd, 'sumPtUnclustered') :                                                                                                                                             
                if 'unclust' in var:
                    sumPtUnclustEnDeltaUp = math.sqrt(metStd.MetUnclustEnUpDeltaX*metStd.MetUnclustEnUpDeltaX  + metStd.MetUnclustEnUpDeltaY*metStd.MetUnclustEnUpDeltaY )
                    if 'Up' in var:
                        totalSumPt = metStd.sumPtUnclustered + sumPtUnclustEnDeltaUp + sumPtFromJets
                    else:
                        totalSumPt = metStd.sumPtUnclustered - sumPtUnclustEnDeltaUp + sumPtFromJets
                else:
                    totalSumPt = metStd.sumPtUnclustered + sumPtFromJets

            else:
                print("Cannot compute MET significance !")
                #throw exception
                return False
              

            ind = 10 if self.pTdependent else 5
            cov_tt = self.pars[ind]**2 + self.pars[ind+1]**2*totalSumPt
            cov_xx += cov_tt
            cov_yy += cov_tt

            det = cov_xx*cov_yy - cov_xy*cov_xy

            #Invert the covariance matrix
            if det>0:
                ncov_xx =  cov_yy / det
                ncov_xy = -cov_xy / det
                ncov_yy =  cov_xx / det
            else:
                #print cov_xx, cov_yy, cov_xy
                ncov_xx = cov_xx if cov_xx > 0 else 1
                ncov_yy = cov_yy if cov_yy > 0 else 1
                ncov_xy = cov_xy if cov_xy > 0 else 1

            met_pt  = getattr(met, metPtVar)
            met_phi = getattr(met, phiVar)

            met_x = met_pt * math.cos(met_phi)
            met_y = met_pt * math.sin(met_phi)

            #Definition of the significance
            MET_sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy


            if var == '' or var == '_nom':
                self.out.fillBranch("MET_significance", MET_sig)
            else:
                self.out.fillBranch("MET_significance"+var, MET_sig)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
#metSig = lambda : METSigProducer( "Summer16_25nsV1_MC", [1.0,1.0,1.0,1.0,1.0,0.0,0.5] )

