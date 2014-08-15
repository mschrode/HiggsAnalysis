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
    cfg_name = "crab2_"+label+".cfg"
    with open(cfg_name,"w") as cfg:
        cfg.write("[CRAB]\n")
        cfg.write("jobtype    = cmssw\n")
        cfg.write("scheduler  = remoteGlidein\n")
        cfg.write("use_server = 0\n")
        cfg.write("\n")
        cfg.write("[CMSSW]\n")
        cfg.write("datasetpath            = "+dataset+"\n")
        cfg.write("pset                   = makeHbbTree-v2_cfg.py\n")
        cfg.write("total_number_of_events = -1\n")
        cfg.write("events_per_job         = 25000\n")
        cfg.write("output_file            = HbbTree.root\n")
        cfg.write("\n")
        cfg.write("[USER]\n")
        cfg.write("return_data         = 0\n")
        cfg.write("copy_data           = 1\n")
        cfg.write("storage_element     = dcache-se-cms.desy.de\n")
        cfg.write("storage_path        = /srm/managerv2?SFN=/pnfs/desy.de/cms/tier2/store/user/mschrode/hbbspring14/v2\n")
        cfg.write("user_remote_dir     = "+label+"\n")
        cfg.write("ui_working_dir      = /afs/desy.de/user/m/matsch/xxl/crab/"+label+"\n")
        cfg.write("dontCheckSpaceLeft  = 1\n")
        cfg.write("\n")
        cfg.write("\n")
        cfg.write("[GRID]\n")
        cfg.write("se_black_list = T0\n")
        cfg.write("virtual_organization = cms\n")
    print("Wrote '"+cfg_name+"'")
