import ROOT
import os
import numpy as np
import math
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

## 
from FWCore.PythonUtilities.LumiList import LumiList

class applyJSON(Module):

    def __init__(self, json_file):
        if json_file:
            self.lumiList = LumiList(os.path.expandvars(json_file))
        else:
            self.lumiList = None

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("jsonPassed", "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if self.lumiList:
            jsonPassed = self.lumiList.contains(event.run, event.luminosityBlock)
        else:
            jsonPassed = 1

        self.out.fillBranch("jsonPassed", jsonPassed)
        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

#json = lambda : applyJSON( 'filepath' )
