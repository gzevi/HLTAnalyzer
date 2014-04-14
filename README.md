HLTAnalyzer
===========

Test code for miniAOD

To run this, first checkout the miniAOD code: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MiniAOD

Then do:
cd $CMSSW_BASE/src
mkdir HLTtest
cd HLTtest/
git clone https://github.com/gzevi/HLTAnalyzer HLTAnalyzer
cd HLTAnalyzer
scramv1 b
cmsRun python/ConfFile_cfg.py
