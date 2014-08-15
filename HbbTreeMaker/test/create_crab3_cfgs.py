#!/usr/bin/env python

import os
import sys
import shutil

datasets = {
    "hbb_spring14_qcd_Pt-0015to0030" : "/QCD_Pt-15to30_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0030to0050" : "/QCD_Pt-30to50_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0050to0080" : "/QCD_Pt-50to80_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0080to0120" : "/QCD_Pt-80to120_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0120to0170" : "/QCD_Pt-120to170_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0170to0300" : "/QCD_Pt-170to300_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0300to0470" : "/QCD_Pt-300to470_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0470to0600" : "/QCD_Pt-470to600_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0600to0800" : "/QCD_Pt-600to800_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-0800to1000" : "/QCD_Pt-800to1000_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-1000to1400" : "/QCD_Pt-1000to1400_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
    "hbb_spring14_qcd_Pt-1400to1800" : "/QCD_Pt-1400to1800_TuneZ2star_13TeV_pythia6/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM",
}

for label,dataset in datasets.iteritems():
    cfg_name = "crab3_"+label+"_cfg.py"
    with open(cfg_name,"w") as cfg:
        cfg.write("from WMCore.Configuration import Configuration\n")
        cfg.write("config = Configuration()\n")
        cfg.write("\n\n")
        cfg.write("config.section_('General')\n")
        cfg.write("config.General.requestName = '"+label+"'\n")
        cfg.write("config.General.workArea = '/afs/desy.de/user/m/matsch/xxl/crab'\n")
        cfg.write("\n\n")
        cfg.write("config.section_('JobType')\n")
        cfg.write("config.JobType.pluginName = 'Analysis'\n")
        cfg.write("config.JobType.psetName = 'makeHbbTree-v2_cfg.py'\n")
        cfg.write("\n\n")
        cfg.write("config.section_('Data')\n")
        cfg.write("config.Data.inputDataset = '"+dataset+"'\n")
        cfg.write("config.Data.dbsUrl = 'global'\n")
        cfg.write("config.Data.splitting = 'FileBased'\n")
        cfg.write("config.Data.unitsPerJob = 10\n")
        cfg.write("config.Data.publication = False\n")
        cfg.write("config.Data.outlfn = '/store/user/mschrode/spring14/"+label+"'\n")
        cfg.write("\n\n")
        cfg.write("config.section_('Site')\n")
        cfg.write("config.Site.storageSite = 'T2_DE_DESY'\n")
    print("Wrote '"+cfg_name+"'")
