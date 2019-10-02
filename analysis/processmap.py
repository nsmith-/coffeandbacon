from collections import OrderedDict
from coffea import hist

process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_map = OrderedDict()

process_map["tthqq125"] = [
    "ttHTobb_M125_TuneCP5_13TeV_powheg_pythia8",
]
process_map["whqq125"] = [
    "WminusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8",
    "WplusH_HToBB_WToQQ_M125_13TeV_powheg_pythia8",
]
process_map["zhqq125"] = [
    "ZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8",
    "ZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8",
    "ggZH_HToBB_ZToQQ_M125_13TeV_powheg_pythia8",
    "ggZH_HToBB_ZToNuNu_M125_13TeV_powheg_pythia8",
]
process_map["vbfhqq125"] = [
    "VBFHToBB_M_125_13TeV_powheg_pythia8_weightfix",
]
process_map["hqq125"] = [
    "GluGluHToBB_M125_13TeV_powheg_pythia8",
]
process_map["hcc125"] = [
    "GluGluHToCC_M125_LHEHpT_250_Inf_13TeV_amcatnloFXFX_pythia8",
]

process_map["zll"] = [
    "DYJetsToLL_M_50_HT_2500toInf_TuneCP5_13TeV",
    "DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV",
    "DYJetsToLL_M_50_HT_800to1200_TuneCP5_13TeV",
    "DYJetsToLL_M_50_HT_1200to2500_TuneCP5_13TeV",
    "DYJetsToLL_M_50_HT_600to800_TuneCP5_13TeV",
]
process_map["wlnu"] = [
    "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8",
    "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8",
]

process_map["vvqq"] = [
    "WW_TuneCP5_13TeV-pythia8",
    "ZZ_TuneCP5_13TeV-pythia8",
    "WZ_TuneCP5_13TeV-pythia8",
]
process_map["stqq"] = [
    "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
    "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8",
]
process_map["tqq"] = [
    "TTToHadronic_TuneCP5_13TeV_powheg_pythia8",
    "TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8",
    "TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8",
]
process_map["zqq"] = [
    "ZJetsToQQ_HT600to800_qc19_4j_TuneCP5_13TeV",
    "ZJetsToQQ_HT-800toInf_qc19_4j_TuneCP5_13TeV",
    "ZJetsToQQ_HT400to600_qc19_4j_TuneCP5_13TeV",
]

process_map["wqq"] = [
    "WJetsToQQ_HT400to600_qc19_3j_TuneCP5_13TeV",
    "WJetsToQQ_HT600to800_qc19_3j_TuneCP5_13TeV",
    "WJetsToQQ_HT-800toInf_qc19_3j_TuneCP5_13TeV",
]
process_map["qcd"] = [
    "QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8",
    "QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8",
    "QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8",
    "QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8",
    "QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8",
]
process_map["data_obs"] = [
    "JetHT",
    "SingleMuon",
]

def apply(h):
    return h.group(process_cat, process, process_map)
