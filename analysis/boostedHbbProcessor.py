#!/usr/bin/env python
import pprint
import copy
import numpy as np
from coffea import hist, processor
from coffea.util import load, save
import argparse
import warnings


def deltaphi(a, b):
    return (a - b + np.pi)%(2*np.pi) - np.pi


class BoostedHbbProcessor(processor.ProcessorABC):
    def __init__(self, corrections, columns=[], debug=False, year='2017', skipPileup=False, skipTrigger=False):
        self._columns = columns
        self._corrections = corrections
        self._debug = debug
        self._year = year
        self._skipPileup = skipPileup
        self._skipTrigger = skipTrigger

        dataset_axis = hist.Cat("dataset", "Primary dataset")
        gencat_axis = hist.Bin("AK8Puppijet0_isHadronicV", "V matching index", [0,1,2,3,9,10,11])
        jetpt_axis = hist.Bin("AK8Puppijet0_pt", r"Jet $p_T$", [450, 500, 550, 600, 675, 800, 1200])
        jetmass_axis = hist.Bin("AK8Puppijet0_msd", r"Jet $m_{sd}$", 23, 40, 201)
        jetpt_coarse_axis = hist.Bin("AK8Puppijet0_pt", r"Jet $p_T$", [450, 600, 1200])
        jetmass_coarse_axis = hist.Bin("AK8Puppijet0_msd", r"Jet $m_{sd}$", [40, 82, 103, 152, 201])
        jetrho_axis = hist.Bin("ak8jet_rho", r"Jet $\rho$", 13, -6, -2.1)
        doubleb_axis = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", 20, 0., 1)
        doublec_axis = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", 20, 0., 1.)
        doublecvb_axis = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", 20, 0., 1.)
        doubleb_wps = [1., 0.9, 0.89, 0.85, 0.7]
        doubleb_coarse_axis = hist.Bin("AK8Puppijet0_deepdoubleb", "Double-b", doubleb_wps[::-1])
        doublec_wps = [0.87, 0.84, 0.83, 0.79, 0.69]
        doublec_coarse_axis = hist.Bin("AK8Puppijet0_deepdoublec", "Double-c", doublec_wps[::-1])
        doublecvb_wps = [0.93, 0.91, 0.6, 0.2, 0.17]
        doublecvb_coarse_axis = hist.Bin("AK8Puppijet0_deepdoublecvb", "Double-cvb", doublecvb_wps[::-1])

        hists = processor.dict_accumulator()
        hist.Hist.DEFAULT_DTYPE = 'f'  # save some space by keeping float bin counts instead of double
        hists['sumw'] = processor.defaultdict_accumulator(int)
        hists['genVpt_noselection'] = hist.Hist("Events / 20 GeV",
                                                dataset_axis,
                                                gencat_axis,
                                                hist.Bin("genVPt", "Gen. V $p_T$", 60, 0, 1200),
                                                )
        hists['jetpt_preselection'] = hist.Hist("Events",
                                                dataset_axis,
                                                gencat_axis,
                                                hist.Bin("runNum", "Run number", [294927, 297046, 299368, 302030, 303824, 305040, 306546]),
                                                hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 50, 300, 1300),
                                                )
        hists['jeteta_preselection'] = hist.Hist("Events",
                                                 dataset_axis,
                                                 gencat_axis,
                                                 hist.Bin("AK8Puppijet0_eta", r"Jet $\eta$", 50, -3, 3),
                                                 )
        hists['jetpt_muoncontrol'] = hist.Hist("Events",
                                               dataset_axis,
                                               gencat_axis,
                                               hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 100, 300, 1300),
                                               )
        hists['muonpt_muoncontrol'] = hist.Hist("Events",
                                                dataset_axis,
                                                gencat_axis,
                                                hist.Bin("vmuoLoose0_pt", "Leading muon $p_T$", 100, 0, 1000),
                                                )
        hists['muoneta_muoncontrol'] = hist.Hist("Events",
                                                 dataset_axis,
                                                 gencat_axis,
                                                 hist.Bin("vmuoLoose0_eta", r"Leading muon $\eta$", 50, -3, 3),
                                                 )
        hists['jetpt_signalregion'] = hist.Hist("Events",
                                                dataset_axis, 
                                                gencat_axis,
                                                hist.Bin("AK8Puppijet0_pt", "Jet $p_T$", 50, 300, 1300),
                                                )
        hists['sculpt_signalregion'] = hist.Hist("Events",
                                                 dataset_axis,
                                                 gencat_axis,
                                                 jetpt_axis,
                                                 jetmass_axis,
                                                 doubleb_coarse_axis,
                                                 doublec_coarse_axis,
                                                 doublecvb_coarse_axis
                                                 )
        hists['ddb_noselection'] = hist.Hist("Events",
                                              dataset_axis,
                                              gencat_axis,
                                              jetpt_coarse_axis,
                                              jetmass_coarse_axis,
                                              doubleb_axis,
                                             )
        hists['ddb_signalregion'] = hist.Hist("Events",
                                              dataset_axis,
                                              gencat_axis,
                                              jetpt_coarse_axis,
                                              jetmass_coarse_axis,
                                              doubleb_axis,
                                             )
        hists['opposite_ak8_n3sdb1_signalregion'] = hist.Hist("Events",
                                                              dataset_axis,
                                                              gencat_axis,
                                                              jetpt_coarse_axis,
                                                              jetmass_coarse_axis,
                                                              hist.Bin("opposite_ak8_n3sdb1", r"Jet $N_{3,sd}^{\beta=1}$", 40, 0.5, 3)
                                                              )
        hists['opposite_ak8_tau32_signalregion'] = hist.Hist("Events",
                                                             dataset_axis,
                                                             gencat_axis,
                                                             jetpt_coarse_axis,
                                                             jetmass_coarse_axis,
                                                             hist.Bin("opposite_ak8_tau32", r"Jet $\tau_{32}$", 40, 0, 1)
                                                             )
        hists['opposite_ak8_msd_signalregion'] = hist.Hist("Events",
                                                           dataset_axis,
                                                           gencat_axis,
                                                           jetpt_coarse_axis,
                                                           jetmass_coarse_axis,
                                                           hist.Bin("opposite_ak8_msd", r"Jet $\m_{sd}$", 40, 50, 200)
                                                           )
        hists['njets_ak4_signalregion'] = hist.Hist("Events",
                                                    dataset_axis,
                                                    gencat_axis,
                                                    jetpt_coarse_axis,
                                                    jetmass_coarse_axis,
                                                    hist.Bin("nAK4PuppijetsPt30", "Number AK4 Jets", 8, 0, 8)
                                                    )

        hists['nminus1_antiak4btagMediumOppHem_signalregion'] = hist.Hist("Events",
                                                                          dataset_axis,
                                                                          gencat_axis,
                                                                          jetpt_coarse_axis,
                                                                          jetmass_coarse_axis,
                                                                          hist.Bin("opposite_ak4_leadingDeepCSV", r"Max(DeepCSV) (of $\leq4$ leading)", 40, 0, 1)
                                                                          )
        hists['nminus1_pfmet_signalregion'] = hist.Hist("Events",
                                                        dataset_axis,
                                                        gencat_axis,
                                                        jetpt_coarse_axis,
                                                        jetmass_coarse_axis,
                                                        doubleb_coarse_axis,
                                                        hist.Bin("pfmet", r"PF $p_{T}^{miss}$", 40, 0, 200)
                                                        )
        hists['nminus1_n2ddtPass_signalregion'] = hist.Hist("Events",
                                                            dataset_axis,
                                                            gencat_axis,
                                                            jetmass_coarse_axis,
                                                            doubleb_coarse_axis,
                                                            hist.Bin("ak8jet_n2ddt", r"Jet $N_{2,DDT}^{\beta=1}$", 40, -.25, .25)
                                                            )
        hists['nminus1_ak4btagMediumDR08_muoncontrol'] = hist.Hist("Events",
                                                                   dataset_axis,
                                                                   gencat_axis,
                                                                   jetmass_coarse_axis,
                                                                   doubleb_coarse_axis,
                                                                   hist.Bin("ak4_leadingDeepCSV_dR08", r"Max(DeepCSV) ($\DeltaR(ak4, ak8)>0.8$)", 40, 0, 1)
                                                                   )
        hists['nminus1_muonDphiAK8_muoncontrol'] = hist.Hist("Events",
                                                             dataset_axis,
                                                             gencat_axis,
                                                             jetmass_coarse_axis,
                                                             doubleb_coarse_axis,
                                                             hist.Bin("muon_dphi", r"$\Delta\phi(\mu, j)$", 40, 0, np.pi)
                                                             )
        hists['templates_signalregion'] = hist.Hist("Events",
                                                    dataset_axis,
                                                    gencat_axis,
                                                    hist.Cat("systematic", "Systematic"),
                                                    jetpt_axis,
                                                    jetmass_axis,
                                                    doubleb_coarse_axis
                                                    )
        hists['templates_muoncontrol'] = hist.Hist("Events",
                                                   dataset_axis,
                                                   gencat_axis,
                                                   hist.Cat("systematic", "Systematic"),
                                                   jetpt_axis,
                                                   jetmass_axis,
                                                   doubleb_coarse_axis
                                                   )
        hists['templates_hCCsignalregion'] = hist.Hist("Events",
                                                    dataset_axis,
                                                    gencat_axis,
                                                    hist.Cat("systematic", "Systematic"),
                                                    jetpt_axis,
                                                    jetmass_axis,
                                                    doublec_coarse_axis
                                                    )
        hists['templates_hCCmuoncontrol'] = hist.Hist("Events",
                                                   dataset_axis,
                                                   gencat_axis,
                                                   hist.Cat("systematic", "Systematic"),
                                                   jetpt_axis,
                                                   jetmass_axis,
                                                   doublec_coarse_axis
                                                   )
        #                                           hist.Bin("ak8jet_n2ddt", r"Jet N_{2,ddt}^{\beta=1}", 8, -.2, .2),
        #                                           hist.Bin("ak8jet_p_bb", r"Jet p_{bb}", 50, 0, 1),
        #                                           hist.Bin("ak8jet_p_cc", r"Jet p_{cc}", 50, 0, 1),
        self._accumulator = hists

    @property
    def columns(self):
        return self._columns

    @property
    def accumulator(self):
        return self._accumulator

    def clean(self, df, val, default, positive=False):
        temp = df[val].copy()
        if positive:
            temp[np.isnan(df[val])|(df[val]<=0.)] = default
        else:
            temp[np.isnan(df[val])|(df[val]==-999.)] = default
        df[val] = temp

    def build_leading_ak8_variables(self, df):
        # jet |eta|<2.5 sometimes gives no events
        # or other cuts in: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/VJetLoader.cc#L270-L272
        # set safe dummy values to avoid domain errors (FPU exceptions slow things down!)
        self.clean(df, 'AK8Puppijet0_pt', 100., positive=True) # msdweight goes negative for pt < 13
        self.clean(df, 'AK8Puppijet0_msd', 1e-7, positive=True)
        self.clean(df, 'AK8Puppijet0_N2sdb1', np.inf)
        self.clean(df, 'AK8Puppijet0_pt_JESUp', 1e-3)
        self.clean(df, 'AK8Puppijet0_pt_JESDown', 1e-3)
        self.clean(df, 'AK8Puppijet0_pt_JERUp', 1e-3)
        self.clean(df, 'AK8Puppijet0_pt_JERDown', 1e-3)
        self.clean(df, 'AK8Puppijet0_deepdoubleb', -1.)
        self.clean(df, 'AK8Puppijet0_deepdoublec', -1.)
        self.clean(df, 'AK8Puppijet0_deepdoublecvb', -1.)
        df['AK8Puppijet0_msd_raw'] = df['AK8Puppijet0_msd']
        # for very large pt values, correction can become negative
        df['AK8Puppijet0_msd'] = np.maximum(1e-7, df['AK8Puppijet0_msd']*self._corrections['msdweight'](df['AK8Puppijet0_pt'], df['AK8Puppijet0_eta']))
        df['ak8jet_rho'] = 2*np.log(df['AK8Puppijet0_msd']/df['AK8Puppijet0_pt'])
        df['ak8jet_n2ddt'] = df['AK8Puppijet0_N2sdb1'] - self._corrections[f'{self._year}_n2ddt_rho_pt'](df['ak8jet_rho'], df['AK8Puppijet0_pt'])
        if 'multiclass' in self._corrections:
            multiclass = self._corrections['multiclass'](df['AK8Puppijet0_deepdoubleb'], df['AK8Puppijet0_deepdoublec'], df['AK8Puppijet0_deepdoublecvb'])
            df['ak8jet_p_bb'] = multiclass[:, 2]
            df['ak8jet_p_cc'] = multiclass[:, 1]
        else:
            df['ak8jet_p_bb'] = (df['AK8Puppijet0_deepdoubleb'] + 1 - df['AK8Puppijet0_deepdoublecvb']) / 3.
            df['ak8jet_p_cc'] = (df['AK8Puppijet0_deepdoublec'] + df['AK8Puppijet0_deepdoublecvb']) / 3.

    def subleading_n3(self, df):
        e4_v2_jet1 = self.clean(df, 'AK8Puppijet1_e4_v2_sdb1', 1.)
        e3_v1_jet1 = self.clean(df, 'AK8Puppijet1_e3_v1_sdb1', 1e-4, positive=True)
        return df['AK8Puppijet1_e4_v2_sdb1']/df['AK8Puppijet1_e3_v1_sdb1']**2

    def build_subleading_ak8_variables(self, df):
        dphi = np.abs(deltaphi(df['AK8Puppijet1_phi'], df['AK8Puppijet0_phi']))
        df['opposite_ak8_n3sdb1'] = np.where(dphi > np.pi/2., self.subleading_n3(df), np.inf)
        df['opposite_ak8_tau32'] = np.where(dphi > np.pi/2., df['AK8Puppijet1_tau32'], np.inf)
        df['opposite_ak8_msd'] = np.where(dphi > np.pi/2., df['AK8Puppijet1_msd'], np.inf)

    def build_ak4_variables(self, df):
        # dR08, dPhi08 with respect to leading ak8 jet: https://github.com/DAZSLE/BaconAnalyzer/blob/102x/Analyzer/src/JetLoader.cc#L478-L479
        n_ak4 = 4
        def stack(var): return np.column_stack([df['AK4Puppijet%d_%s' % (i, var)] for i in range(n_ak4)])
        dR = stack('dR08')
        dphi = stack('dPhi08')
        btag = stack('deepcsvb') + stack('deepcsvbb')
        pt = stack('pt')
        # seems |eta|<2.5 already in tuple
        require = (np.abs(dphi) > np.pi/2) & (pt > 30.)
        btag_ttrej = np.where(require, btag, -np.inf)
        df['opposite_ak4_leadingDeepCSV'] = np.max(btag_ttrej, axis=1)
        require = (dR > 0.8) & (pt > 50.)
        btag_muCR = np.where(require, btag, -np.inf)
        df['ak4_leadingDeepCSV_dR08'] = np.max(btag_muCR, axis=1)

    def build_met_systematics(self, df):
        metx = df['pfmet']*np.sin(df['pfmetphi'])
        mety = df['pfmet']*np.cos(df['pfmetphi'])
        df['pfmet_JESUp'] = np.hypot(metx + df['MetXCorrjesUp'], mety + df['MetYCorrjesUp'])
        df['pfmet_JESDown'] = np.hypot(metx + df['MetXCorrjesDown'], mety + df['MetYCorrjesDown'])
        df['pfmet_JERUp'] = np.hypot(metx + df['MetXCorrjerUp'], mety + df['MetYCorrjerUp'])
        df['pfmet_JERDown'] = np.hypot(metx + df['MetXCorrjerDown'], mety + df['MetYCorrjerDown'])

    def process(self, df):
        dataset = df['dataset']
        if self._debug:
            print("Processing dataframe from", dataset)
        isRealData = dataset in ["JetHT", "SingleMuon", "data_obs_mu", "data_obs_jet"]

        self.build_leading_ak8_variables(df)
        self.build_subleading_ak8_variables(df)
        self.build_ak4_variables(df)
        self.build_met_systematics(df)
        df['muon_dphi'] = np.abs(deltaphi(df['vmuoLoose0_phi'], df['AK8Puppijet0_phi']))

        selection = processor.PackedSelection()
        if isRealData:
            # Only take jet triggers from JetHT, single muon triggers from SingleMuon dataset
            # necessary but not sufficient condition to prevent double-counting
            # (this plus mutually exclusive offline selections are sufficient)
            selection.add('trigger', (df['triggerBits'] & self._corrections[f'{self._year}_triggerMask']).astype('bool') & (dataset=="JetHT"))
            selection.add('mutrigger', ((df['triggerBits']&1) & df['passJson']).astype('bool') & (dataset=="SingleMuon"))
            if self._debug:
                print("Trigger pass/all", selection.all('trigger').sum(), df.size)
                print("Muon trigger pass/all", selection.all('mutrigger').sum(), df.size)
        else:
            selection.add('trigger', np.ones(df.size, dtype='bool'))
            selection.add('mutrigger', np.ones(df.size, dtype='bool'))

        btagLooseWPs = {
            '2016': 0.6321,
            '2017': 0.4941,
            '2018': 0.4184,
        }

        selection.add('noLeptons', (df['neleLoose']==0) & (df['nmuLoose']==0) & (df['ntau']==0))
        selection.add('oneMuon', (df['neleLoose']==0) & (df['nmuLoose']==1) & (df['ntau']==0))
        selection.add('muonAcceptance', (df['vmuoLoose0_pt'] > 55.) & (np.abs(df['vmuoLoose0_eta']) < 2.1))
        selection.add('muonDphiAK8', df['muon_dphi'] > 2*np.pi/3)
        selection.add('ak4btagMediumDR08', df['ak4_leadingDeepCSV_dR08'] > btagLooseWPs[self._year])  # at least one passes medium cut
        selection.add('antiak4btagMediumOppHem', df['opposite_ak4_leadingDeepCSV'] < btagLooseWPs[self._year])  # none pass
        selection.add('tightVjet', df['AK8Puppijet0_isTightVJet'] != 0)
        selection.add('n2ddtPass', df['ak8jet_n2ddt'] < 0)
        selection.add('jetMass', df['AK8Puppijet0_msd'] > 40.)
        selection.add('deepcvb', df['AK8Puppijet0_deepdoublecvb'] > 0.2)

        selection.add('jetKinematics', df['AK8Puppijet0_pt'] > 450.)
        selection.add('jetKinematicsMuonCR', df['AK8Puppijet0_pt'] > 400.)
        selection.add('pfmet', df['pfmet'] < 140.)

        regions = {}
        regions['noselection'] = {}
        regions['preselection'] = {'trigger', 'noLeptons'}
        regions['signalregion'] = {'trigger', 'noLeptons', 'jetKinematics', 'pfmet', 'n2ddtPass', 'tightVjet', 'antiak4btagMediumOppHem'}
        regions['muoncontrol'] = {'mutrigger', 'oneMuon', 'muonAcceptance', 'jetKinematicsMuonCR', 'n2ddtPass', 'tightVjet', 'ak4btagMediumDR08', 'muonDphiAK8'}
        regions['hCCsignalregion'] = {'trigger', 'noLeptons', 'jetKinematics', 'pfmet', 'n2ddtPass', 'tightVjet', 'antiak4btagMediumOppHem', 'deepcvb'}
        regions['hCCmuoncontrol'] = {'mutrigger', 'oneMuon', 'muonAcceptance', 'jetKinematicsMuonCR', 'n2ddtPass', 'tightVjet', 'ak4btagMediumDR08', 'muonDphiAK8', 'deepcvb'}

        shiftSystematics = ['JESUp', 'JESDown', 'JERUp', 'JERDown']
        shiftedQuantities = {'AK8Puppijet0_pt', 'pfmet'}
        shiftedSelections = {'jetKinematics', 'jetKinematicsMuonCR', 'pfmet'}
        for syst in shiftSystematics:
            selection.add('jetKinematics'+syst, df['AK8Puppijet0_pt_'+syst] > 450)
            selection.add('jetKinematicsMuonCR'+syst, df['AK8Puppijet0_pt_'+syst] > 400.)
            selection.add('pfmet'+syst, df['pfmet_'+syst] < 140.)

        # mass shift applied only to V-matched data
        # https://github.com/kakwok/ZPrimePlusJet/blob/PerBinEff/fitting/PbbJet/buildRhalphabetHbb.py#L30
        if not isRealData:
            shiftSystematics.append('matchedUp')
            shiftedQuantities.add('AK8Puppijet0_msd')
            msdshifts = {'2016': 1.001, '2017': 0.979, '2018': 0.970}
            df['AK8Puppijet0_msd_matchedUp'] = msdshifts[self._year] * df['AK8Puppijet0_msd']

        weights = processor.Weights(df.size)

        if not isRealData:
            # SumWeights is sum(scale1fb), so we need to use full value here
            weights.add('genweight', df['scale1fb'])

        if not self._skipPileup:
            if self._year == '2017' and dataset in self._corrections['2017_pileupweight_dataset']:
                weights.add('pileupweight',
                            self._corrections['2017_pileupweight_dataset'][dataset](df['npu']),
                            self._corrections['2017_pileupweight_dataset_puUp'][dataset](df['npu']),
                            self._corrections['2017_pileupweight_dataset_puDown'][dataset](df['npu']),
                            )
            elif self._year != '2017':
                weights.add('pileupweight',
                            self._corrections[f'{self._year}_pileupweight'](df['npu']),
                            self._corrections[f'{self._year}_pileupweight_puUp'](df['npu']),
                            self._corrections[f'{self._year}_pileupweight_puDown'](df['npu']),
                            )

        # TODO unc.
        if self._year == '2017' and 'ZJetsToQQ_HT' in dataset:
            nlo_over_lo_qcd = self._corrections['2017_Z_nlo_qcd'](df['genVPt'])
            nlo_over_lo_ewk = self._corrections['Z_nlo_over_lo_ewk'](df['genVPt'])
            weights.add('kfactor', nlo_over_lo_qcd * nlo_over_lo_ewk)
        elif self._year == '2017' and 'WJetsToQQ_HT' in dataset:
            nlo_over_lo_qcd = self._corrections['2017_W_nlo_qcd'](df['genVPt'])
            nlo_over_lo_ewk = self._corrections['W_nlo_over_lo_ewk'](df['genVPt'])
            weights.add('kfactor', nlo_over_lo_qcd * nlo_over_lo_ewk)
        elif self._year == '2016' and 'DYJetsToQQ' in dataset:
            nlo_over_lo_qcd = self._corrections['2016_Z_nlo_qcd'](df['genVPt'])
            nlo_over_lo_ewk = self._corrections['Z_nlo_over_lo_ewk'](df['genVPt'])
            weights.add('kfactor', nlo_over_lo_qcd * nlo_over_lo_ewk)
        elif self._year == '2016' and 'WJetsToQQ' in dataset:
            nlo_over_lo_qcd = self._corrections['2016_W_nlo_qcd'](df['genVPt'])
            nlo_over_lo_ewk = self._corrections['W_nlo_over_lo_ewk'](df['genVPt'])
            weights.add('kfactor', nlo_over_lo_qcd * nlo_over_lo_ewk)

        if not isRealData:
            # handle weight systematics for signal region
            def regionMask(w):
                if self._skipTrigger:
                    return np.ones(df.size)
                return np.where(selection.all('noLeptons'), w, 1.)
            weights.add('trigweight',
                        regionMask(self._corrections[f'{self._year}_trigweight_msd_pt'](df['AK8Puppijet0_msd_raw'], df['AK8Puppijet0_pt'])),
                        regionMask(self._corrections[f'{self._year}_trigweight_msd_pt_trigweightUp'](df['AK8Puppijet0_msd_raw'], df['AK8Puppijet0_pt'])),
                        regionMask(self._corrections[f'{self._year}_trigweight_msd_pt_trigweightDown'](df['AK8Puppijet0_msd_raw'], df['AK8Puppijet0_pt'])),
                        )
            vmatch = (np.abs(deltaphi(df['AK8Puppijet0_phi'], df['genVPhi'])) < 0.8) & (np.abs(df['AK8Puppijet0_pt']-df['genVPt'])/df['genVPt'] < 0.5) & (np.abs(df['AK8Puppijet0_msd']-df['genVMass'])/df['genVMass'] < 0.3)
            weights.add('matched', np.ones(df.size, dtype='f'), vmatch.astype('f'), 1.-vmatch)

            # handle weight systematics for muon CR
            def regionMask(w):
                if self._skipTrigger:
                    return np.ones(df.size)
                return np.where(selection.all('oneMuon'), w, 1.)
            mu_abseta = np.abs(df['vmuoLoose0_eta'])
            weights.add('mutrigweight',
                        regionMask(self._corrections[f'{self._year}_mutrigweight_pt_abseta'](df['vmuoLoose0_pt'], mu_abseta)),
                        regionMask(self._corrections[f'{self._year}_mutrigweight_pt_abseta_mutrigweightShift'](df['vmuoLoose0_pt'], mu_abseta)),
                        shift=True
                        )
            weights.add('muidweight',
                        regionMask(self._corrections[f'{self._year}_muidweight_abseta_pt'](mu_abseta, df['vmuoLoose0_pt'])),
                        regionMask(self._corrections[f'{self._year}_muidweight_abseta_pt_muidweightShift'](mu_abseta, df['vmuoLoose0_pt'])),
                        shift=True
                        )
            weights.add('muisoweight',
                        regionMask(self._corrections[f'{self._year}_muisoweight_abseta_pt'](mu_abseta, df['vmuoLoose0_pt'])),
                        regionMask(self._corrections[f'{self._year}_muisoweight_abseta_pt_muisoweightShift'](mu_abseta, df['vmuoLoose0_pt'])),
                        shift=True
                        )

        if self._debug:
            print("Weight statistics:")
            pprint.pprint(weights._weightStats, indent=4)

        hout = self.accumulator.identity()
        for histname, h in hout.items():
            if not isinstance(h, hist.Hist):
                continue
            if not all(k in df or k == 'systematic' for k in h.fields):
                # Cannot fill this histogram due to missing fields
                # is this an error, warning, or ignorable?
                if self._debug:
                    print("Missing fields %r from %r" % (set(h.fields) - set(df.keys()), h))
                continue
            fields = {k: df[k] for k in h.fields if k in df}
            region = [r for r in regions.keys() if r in histname.split('_')]

            if 'nminus1' in histname:
                _, sel, region = histname.split('_')
                cut = regions[region] - {sel}
                weight = weights.weight() * selection.all(*cut)
                h.fill(**fields, weight=weight)
            elif len(region) == 1:
                region = region[0]
                weight = weights.weight()
                cut = selection.all(*regions[region])
                if 'systematic' not in h.fields:
                    h.fill(**fields, weight=weight*cut)
                else:
                    h.fill(systematic="", **fields, weight=weight*cut)
                    if self._debug:
                        print("Filling systematics for %s" % histname)
                    systs = set(weights.variations)
                    systs.update(shiftSystematics)
                    for syst in systs:
                        if self._debug:
                            print("  Filling systematic %s" % syst)
                        fields_syst = copy.copy(fields)
                        for val in shiftedQuantities:
                            if val+'_'+syst in df and val in fields:
                                fields_syst[val] = df[val+'_'+syst]
                                if self._debug:
                                    print("    Replacing field %s with %s" % (val, val+'_'+syst))
                        if syst in weights.variations:
                            weight_syst = weights.weight(syst)
                            if self._debug:
                                print("    Using modified weight")
                        else:
                            weight_syst = weight
                        if syst in set(shiftSystematics):
                            cut_syst = set()
                            for sel in regions[region]:
                                if sel in shiftedSelections and sel+syst in selection.names:
                                    cut_syst.add(sel+syst)
                                    if self._debug:
                                        print("    Replacing cut %s with systematic-shifted %s" % (sel, sel+syst))
                                else:
                                    cut_syst.add(sel)
                            cut_syst = selection.all(*cut_syst)
                        else:
                            cut_syst = cut
                        h.fill(systematic=syst, **fields_syst, weight=weight_syst*cut_syst)
            elif len(region) > 1:
                raise ValueError("Histogram '%s' has a name matching multiple region definitions: %r" % (histname, region))
            else:
                raise ValueError("Histogram '%s' does not fall into any region definitions." % (histname, ))

        if not isRealData:
            if 'skim_sumw' in df:
                # hacky way to only accumulate file-level information once
                if df['skim_sumw'] is not None:
                    hout['sumw'][dataset] += df['skim_sumw']
            else:
                hout['sumw'][dataset] += np.sum(df['scale1fb'])
        return hout

    def postprocess(self, accumulator):
        # set everything to 1/fb scale
        lumi = 1000  # [1/pb]
        
        if 'sumw_external' in self._corrections:
            normlist = self._corrections['sumw_external']
            for key in accumulator['sumw'].keys():
                accumulator['sumw'][key] = normlist[key].value

        scale = {}
        for dataset, dataset_sumw in accumulator['sumw'].items():
            if dataset in self._corrections['xsections']:
                scale[dataset] = lumi*self._corrections['xsections'][dataset]/dataset_sumw
            else:
                warnings.warn("Missing cross section for dataset %s.  Normalizing to 1 pb" % dataset, RuntimeWarning)
                scale[dataset] = lumi / dataset_sumw
            
        for h in accumulator.values():
            if isinstance(h, hist.Hist):
                h.scale(scale, axis="dataset")

        return accumulator


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Boosted Hbb processor')
    parser.add_argument('--year', choices=['2016', '2017', '2018'], default='2017', help='Which data taking year to correct MC to.')
    parser.add_argument('--debug', action='store_true', help='Enable debug printouts')
    parser.add_argument('--skipPileup', action='store_true', help='Do not apply pileup reweight corrections to MC')
    parser.add_argument('--skipTrigger', action='store_true', help='Do not apply trigger weight corrections to MC')
    parser.add_argument('--externalSumW', help='Path to external sum weights file (if provided, will be used in place of self-determined sumw)')
    parser.add_argument('--multiclass', help='Path to multiclass transform NN model')
    args = parser.parse_args()

    corrections = load('corrections.coffea')

    if args.externalSumW is not None:
        corrections['sumw_external'] = load(args.externalSumW)

    if args.multiclass is not None:
        corrections['multiclass'] = load(args.multiclass)

    from columns import gghbbcolumns, gghbbcolumns_mc
    allcolumns = gghbbcolumns + gghbbcolumns_mc
        
    processor_instance = BoostedHbbProcessor(corrections=corrections,
                                             columns=allcolumns,
                                             debug=args.debug,
                                             year=args.year,
                                             skipPileup=args.skipPileup,
                                             skipTrigger=args.skipTrigger,
                                             )

    save(processor_instance, 'boostedHbbProcessor.coffea')
