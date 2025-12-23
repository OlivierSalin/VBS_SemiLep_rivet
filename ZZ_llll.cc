// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/RivetHepMC.hh"

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/Cutflow.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "TFile.h"
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "EventWeights.h"

namespace Rivet
{

    /// @brief Add a short analysis description here
    class ZZ_llll : public Analysis
    {
    public:
        /// Constructor
        RIVET_DEFAULT_ANALYSIS_CTOR(ZZ_llll);



        const double Centrality(const Jets &tagjets, const FourMomentum &jet_Vlep, const FourMomentum &jet_Vhad)
        {
            double min_eta_tag_jet = std::min(tagjets[0].eta(), tagjets[1].eta());
            double max_eta_tag_jet = std::max(tagjets[0].eta(), tagjets[1].eta());
            double min_eta_Vjet = std::min(jet_Vlep.eta(), jet_Vhad.eta());
            double max_eta_Vjet = std::max(jet_Vlep.eta(), jet_Vhad.eta());
            double delta_eta_pos = max_eta_tag_jet - min_eta_Vjet;
            double delta_eta_neg = max_eta_Vjet - min_eta_tag_jet;
            return std::min(delta_eta_pos, delta_eta_neg);
        }

        double cosCollinsSoper(const FourMomentum &l1, const FourMomentum &l2)
        {
            const FourMomentum ll = l1 + l2;
            const double nom = (l1.E() + l1.pz()) * (l2.E() - l2.pz()) - (l1.E() - l1.pz()) * (l2.E() + l2.pz());
            const double denom = ll.mass() * sqrt(sqr(ll.mass()) + sqr(ll.pt()));
            return sign(ll.pz()) * safediv(nom, denom); // protect against division by zero, you never know...
        }

        FourMomentum RotateZ(double phi, double costheta, FourMomentum &p_)
        {
            FourMomentum p = p_;
            double sintheta = sqrt(1. - pow(costheta, 2.));
            double RotE = p.E();
            double RotPx = cos(phi) * costheta * p.px() - sin(phi) * p.py() + cos(phi) * sintheta * p.pz();
            double RotPy = sin(phi) * costheta * p.px() + cos(phi) * p.py() + sin(phi) * sintheta * p.pz();
            double RotPz = -sintheta * p.px() + costheta * p.pz();
            p.setXYZE(RotPx, RotPy, RotPz, RotE);
            return p;
        }

        double pair_m_dist_m_z(const Particle &lep_1, const Particle &lep_2) {
          double i_m_ll = (lep_1.mom() + lep_2.mom()).mass() / GeV;
          double i_m_ll_dist_m_z = fabs(i_m_ll - 91.18);
          return i_m_ll_dist_m_z;
        }

        /// @name Analysis methods
        /// @{

        /// Book histograms and initialise projections before the run
        void init()
        {
            std::string txt_dir = "/exp/atlas/salin/ATLAS/VBS_mc/VBS_Pol_Rivet/VBS_rivet/";

            std::string out_dir = getOption("OUTDIR");
            // std::string out_dir = histoDir();
            // std::string out_dir = "";
            // auto cross_section_fb = crossSection()/femtobarn;

            std::cout << "out_dir Rivet: " << out_dir << std::endl;

            // cross_section_fb = crossSection()/femtobarn;

            if (out_dir.find("SM") != std::string::npos)
                _label = 0;
            else
            {
                if (out_dir.find("FS") != std::string::npos)
                    _label = 1;
                if (out_dir.find("FM") != std::string::npos && out_dir.find("odd") == std::string::npos)
                    _label = 2;
                if (out_dir.find("FT") != std::string::npos)
                    _label = 3;
                if (out_dir.find("FM") != std::string::npos && out_dir.find("odd") != std::string::npos)
                    _label = 4;
                if (out_dir.find("FT") != std::string::npos && out_dir.find("odd") != std::string::npos)
                    _label = 5;
            }

            if (out_dir.find("SM") != std::string::npos)
                _label_binary = 0;
            if (out_dir.find("QUAD") != std::string::npos)
                _label_binary = 1;

            std::string ntuple_dir = out_dir;

            _docut = 0; // most cuts on number of particles are always applied to avoid segfault
            if (out_dir.find("DOCUT_YES") != string::npos)
                _docut = 1;
            std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

            std::string jsonfilestr = txt_dir + "Cuts/ZZ_llll_cuts.json";
            std::cout << "++++++assume .json for this ZZ_llll" << " is " << jsonfilestr << "\n";
            std::ifstream json_file(jsonfilestr);
            _jcuts = json::parse(json_file);

            std::cout << "++++++ to check json 1 var got photon pt min " << _jcuts["m_tagjets"] << "\n";
            // CUT Lepton

            _electron_eta_cut =
                    (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) ||
                    (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
            _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
            _electron_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_electron"]) * GeV;
            _muon_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_muon"]) * GeV;
            // The basic final-state projection:
            // all final-state particles within
            // the given eta acceptance
            const FinalState fs;

            // FinalState of direct photons and bare muons and electrons in the event - ignore taus but if want to include use TauFinder
            DirectFinalState bare_e(Cuts::abspid == PID::ELECTRON);
            DirectFinalState bare_mu(Cuts::abspid == PID::MUON);
            DirectFinalState photons_for_dressing(Cuts::abspid == PID::PHOTON);
            // Dress the bare direct leptons with direct photons within dR < 0.1,
            // and apply some fiducial cuts on the dressed leptons depending on param passed
            DressedLeptons dressed_e(photons_for_dressing, bare_e, 0.1);
            DressedLeptons dressed_mu(photons_for_dressing, bare_mu, 0.1);
            // declare(dressed_leps, "leptons_stable");
            declare(dressed_e, "e_stable");
            declare(dressed_mu, "mu_stable");

            // The final-state particles declared above are clustered using FastJet with
            // the anti-kT algorithm and a jet-radius parameter 0.4
            // muons and neutrinos are excluded from the clustering, also veto electrons(+muons but this is redundant) there
            VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, _jcuts["eta_tagjets"])));
            hadrons.addVetoOnThisFinalState(dressed_e);
            hadrons.addVetoOnThisFinalState(dressed_mu);
            declare(hadrons, "hadrons");
            FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
            declare(jetsfs, "jets");


            // Merged histograms

            // plots common with others


            // plots that are not in other ana
            std::ifstream ana_hist_min_file(txt_dir + "/Hists_Lep/2lepton_hists_min2.json");
            json ana_hist_min = json::parse(ana_hist_min_file);
            for (json::iterator it = ana_hist_min.begin(); it != ana_hist_min.end(); ++it)
            {
                book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
                _hist_names.push_back(it.key());
            }


            _tf = make_unique<TFile>(getOption("ROOTFILE", ntuple_dir + "ntuple_rivet.root").c_str(), "RECREATE");
            _tt_SR = make_unique<TTree>("SR", "Rivet_physics");
            _tt_SR->Branch("EventNumber", &EventNumber);
            _tt_SR->Branch("EventWeight", &EventWeight);
      
    
            _tt_SR->Branch("Label", &_label);
            _tt_SR->Branch("Label_binary", &_label_binary);
            for (auto &var_ : varMap)
            {
                _tt_SR->Branch(var_.first.c_str(), var_.second);
            }
            for (auto &var_ : varMapInt)
            {
                _tt_SR->Branch(var_.first.c_str(), var_.second);
            }

            for (auto &var_ : weightMap)
            {
                _tt_SR->Branch(var_.first.c_str(), &var_.second);
            }
            for (auto &var_ : weightMap_cross)
            {
                _tt_SR->Branch(var_.first.c_str(), &var_.second);
            }
            for (auto &var_ : weightMap_int)
            {
                _tt_SR->Branch(var_.first.c_str(), &var_.second);
            }


            _tt_bef_cut = make_unique<TTree>("Bef_cut", "Rivet_physics");
            _tt_bef_cut->Branch("EventNumber", &EventNumber);
            _tt_bef_cut->Branch("EventWeight", &EventWeight);
            _tt_bef_cut->Branch("VBS_event", &VBS_event);
            _tt_bef_cut->Branch("Label", &_label);
            for (auto &var_ : weightMap)
            {
                _tt_bef_cut->Branch(var_.first.c_str(), &var_.second);
            }
            for (auto &var_ : weightMap_cross)
            {
                _tt_bef_cut->Branch(var_.first.c_str(), &var_.second);
            }
            for (auto &var_ : weightMap_int)
            {
                _tt_bef_cut->Branch(var_.first.c_str(), &var_.second);
            }
            



            // counter for efficiency
            book(_c["pos_w_initial"], "pos_w_initial");
            book(_c["pos_w_final"], "pos_w_final");

            book(_c["neg_w_initial"], "neg_w_initial");
            book(_c["neg_w_final"], "neg_w_final");


            // Cut-flows merged region
            _cutflows.addCutflow("ZZ_llll_selections", {"have_four_lep", "pt_lep1_2", "dR_all_pairs", "SFOC_2pairs_min",
                                                        "m_llll", "n_jets", "pt_tagjet1_2", "m_tagjets", "dy_tagjets",
                                                        "centrality_quadjj"});
            // Cut-flows resolved region


            // Read systWeights from JSON file
            std::string weight_json_file_path = out_dir + "/systWeights.json";
            std::ifstream weight_json_file(weight_json_file_path);
            json systWeights_json;
            weight_json_file >> systWeights_json;



            for (auto& [key, value] : systWeights_json.items()) {
                weightNameToIndex[key] = value;
                indexToWeightName[value] = key;
            }
        }

        void analyze(const Event &event)
        {
            // save weights before cuts
            EventNumber = event.genEvent()->event_number();
            double ev_nominal_weight = event.weights()[0];
            EventWeight = ev_nominal_weight;

            std::vector<string> weight_names = Rivet::HepMCUtils::weightNames(*event.genEvent());

 
            std::vector<double> weights_mc = event.genEvent()->weights();
            const Particles all_particles = event.allParticles();

            if (ev_nominal_weight >= 0)
            { _c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
            else {_c["neg_w_initial"]->fill();}

            _tt_bef_cut->Fill();
            _cutflows.fillinit();
            // Set the weight values
            for (const auto& [key, value] : weightNameToIndex) {
                if (key.find("quad") != std::string::npos || key.find("cross") != std::string::npos
                    || key.find("QUAD") != std::string::npos || key.find("CROSS") != std::string::npos
                    || key.find("int") != std::string::npos || key.find("INT") != std::string::npos) {
                    std::string key_lower = key;
                    std::transform(key_lower.begin(), key_lower.end(), key_lower.begin(), ::tolower);
                    std::string weightName = "EventWeight_" + key_lower;
                    double weight = weights_mc[value];
                    //std::cout << "Weight name: " << weightName << " Weight: " << weight << std::endl;
                    if (key.find("quad")!= std::string::npos ||key.find("QUAD") != std::string::npos){
                        weightMap[weightName] = weight;
                        // Fill polarisation map for _quad_ll, _quad_lt, _quad_tl, _quad_tt (case-insensitive)
                        if (
                            key_lower.find("_quad_ll") != std::string::npos || key.find("_QUAD_LL") != std::string::npos ||
                            key_lower.find("_quad_lt") != std::string::npos || key.find("_QUAD_LT") != std::string::npos ||
                            key_lower.find("_quad_tl") != std::string::npos || key.find("_QUAD_TL") != std::string::npos ||
                            key_lower.find("_quad_tt") != std::string::npos || key.find("_QUAD_TT") != std::string::npos
                        ) {
                            weightMap_Polarisation[weightName] = weight;
                        }
                    }
                    else if (key.find("cross")!= std::string::npos || key.find("CROSS") != std::string::npos){
                        weightMap_cross[weightName] = weight;
                    }
                    else if (key.find("int")!= std::string::npos || key.find("INT") != std::string::npos){
                        weightMap_int[weightName] = weight;
                        //std::cout << "Weight name: " << weightName << " Weight: " << weightMap_int[weightName] << std::endl;
                    }
                    //std::cout << "Weight name: " << weightName << " Weight: " << weightMap[weightName] << std::endl;
                }
            }

            // Retrieve dressed leptons, sorted by pT
            Particles e_stable;
            Particles mu_stable;
            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _electron_pt_cut);
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _muon_pt_cut);

            Particles leptons = e_stable + mu_stable;
            // sort by pt as not clear if e+m is in order invidually not but necessary together
            std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
                          return a.pT() > b.pT(); // biggest pT will be first in array
                      });

            int nlep = leptons.size();
            // CHANGING THE TEMPLATE FROM HERE
            if (nlep < _jcuts["n_lepton_stable"]) vetoEvent;
            _cutflows.fillnext();

            if ((leptons[0].pT() < _jcuts["pt_lepton1"] || leptons[1].pT() < _jcuts["pt_lepton2"]))
                vetoEvent;
            _cutflows.fillnext();

            // deltaR between leptons for all pairs
            int num_pair_bad_dR = 0;
            for (int i = 0; i < nlep; i++) {
                const Particle &i_lep = leptons[i];
                for (int j = 0; j < nlep; j++) {
                // to avoid comparing to itself and double couting like (i,j),(j,i) as two separate pairs
                if (i == j || j > i) continue;
                const Particle &j_lep = leptons[j];
                double i_dR = deltaR(i_lep.rap(), i_lep.phi(), j_lep.rap(), j_lep.phi());
                if (i_dR < _jcuts["dR_all_pairs"]) num_pair_bad_dR += 1;
                }
            }
            if ( num_pair_bad_dR > 0) vetoEvent;
            _cutflows.fillnext();

            // check SFOC and m_ll and save possible pairs
            std::vector<std::vector<int>> pairs_ind;
            for (int i = 0; i < nlep; i++) {
                const Particle &i_lep = leptons[i];
                for (int j = 0; j < nlep; j++) {
                // to avoid comparing to itself and double couting like (i,j),(j,i) as two separate pairs
                if (i == j || j > i) continue;
                const Particle &j_lep = leptons[j];
                int i_sum_pids = i_lep.pid() + j_lep.pid(); // to have SFOC will be 0
                double i_m_ll = (i_lep.mom() + j_lep.mom()).mass() / GeV;
                if (i_sum_pids == 0 && i_m_ll > _jcuts["m_ll_all_pairs"]) pairs_ind.push_back({i, j});
                }
            }
            if (pairs_ind.size() < 2) vetoEvent;
            _cutflows.fillnext();


            // order found pairs by how close they are m_z
            std::sort(pairs_ind.begin(), pairs_ind.end(),
                        [this, leptons](std::vector<int> &ind_pair_1, std::vector<int> &ind_pair_2) {
                            double dist_pair_1 = pair_m_dist_m_z(leptons[ind_pair_1[0]], leptons[ind_pair_1[1]]);
                            double dist_pair_2 = pair_m_dist_m_z(leptons[ind_pair_2[0]], leptons[ind_pair_2[1]]);
                            return dist_pair_1 < dist_pair_2; // closest to m_z will be first
                        });
            // define two pairs as the ones with smallest dist
            const Particle &pair_1_lep_1 = leptons[pairs_ind[0][0]];
            const Particle &pair_1_lep_2 = leptons[pairs_ind[0][1]];
            const Particle &pair_2_lep_1 = leptons[pairs_ind[1][0]];
            const Particle &pair_2_lep_2 = leptons[pairs_ind[1][1]];

            const Particle Z1_lepton_minus = (pair_1_lep_1.charge() < 0) ? pair_1_lep_1 : pair_1_lep_2;
            const Particle Z2_lepton_minus = (pair_2_lep_1.charge() < 0) ? pair_2_lep_1 : pair_2_lep_2;

            int sum_abs_pids_quadruplet = fabs(pair_1_lep_1.pid()) + fabs(pair_1_lep_2.pid()) + fabs(pair_2_lep_1.pid()) +
                                            fabs(pair_2_lep_2.pid());

            double m_ll_pair_1 = (pair_1_lep_1.mom() + pair_1_lep_2.mom()).mass() / GeV;
            double m_ll_pair_2 = (pair_2_lep_1.mom() + pair_2_lep_2.mom()).mass() / GeV;

            const FourMomentum fourvec_Z1 = pair_1_lep_1.mom() + pair_1_lep_2.mom();
            const FourMomentum fourvec_Z2 = pair_2_lep_1.mom() + pair_2_lep_2.mom();
            const FourMomentum fourvec_ZZ = fourvec_Z1 + fourvec_Z2;

            const FourMomentum fourvec_llll =
                    pair_1_lep_1.mom() + pair_1_lep_2.mom() + pair_2_lep_1.mom() + pair_2_lep_2.mom();
            double m_llll = fourvec_llll.mass() / GeV;
            if ( m_llll < _jcuts["m_llll"]) vetoEvent;
            _cutflows.fillnext();

            // // Retrieve clustered jets, sorted by pT, with a minimum pT cut
            Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_jet"]) * GeV);
            idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

            int n_jets = jets.size();
            if (n_jets < _jcuts["n_jets"]) vetoEvent;
            _cutflows.fillnext();

            const FourMomentum tag1_jet = jets[0].mom();
            const FourMomentum tag2_jet = jets[1].mom();
            if (
                (tag1_jet.pT() < dbl(_jcuts["pt_tagjet1"]) || tag2_jet.pT() < dbl(_jcuts["pt_tagjet2"])))
                vetoEvent;
            _cutflows.fillnext();

            const double m_tagjets = (tag1_jet + tag2_jet).mass() / GeV;
            if ( m_tagjets < _jcuts["m_tagjets"]) vetoEvent;
            _cutflows.fillnext();

            const FourMomentum fourvec_jj = tag1_jet + tag2_jet;
            const FourMomentum fourvec_lllljj = fourvec_llll + fourvec_jj;

            const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
            if ( dy_tagjets < _jcuts["dy_tagjets"]) vetoEvent;
            if ( tag1_jet.eta() * tag2_jet.eta() > 0) vetoEvent;
            _cutflows.fillnext();

            const double centrality_quadjj = fabs(0.5 * (fourvec_llll.rap() - (tag1_jet.rap() + tag2_jet.rap()) / 2) /
                                                    (tag1_jet.rap() - tag2_jet.rap()));
            if ( centrality_quadjj > _jcuts["centrality_quadjj"]) vetoEvent;
            _cutflows.fillnext();

            int n_gap_jets = 0;
            for (int i = 0; i < n_jets; i++) {
                const double i_jet_rap = jets[i].rap();
                if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) ||
                    (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))
                ++n_gap_jets;
            }
            // do clipping - sometimes there are two Z and one gamma - in this case to avoid much work take Z with highest pt and gamma
            std::vector<FourMomentum> hs_bosons_z = {};
            for (const Particle &p_rivet: event.allParticles()) {
                ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
                int status = p_hepmc->status();
                if (abs(status) == 23 or abs(status) == 22) {
                int i_pid = p_hepmc->pid();
                FourMomentum i_mom = p_hepmc->momentum();
                if (abs(i_pid) == 23) { hs_bosons_z.push_back(i_mom); }
                }
            }
            std::sort(hs_bosons_z.begin(), hs_bosons_z.end(), [](FourMomentum const &a, FourMomentum const &b) {
                return a.pT() > b.pT();
            }); // biggest pT will be first in array
            bool have_two_hs_bosons = false;
            double hs_diboson_mass = 0.0;
            if (hs_bosons_z.size() > 1) {
                hs_diboson_mass = (hs_bosons_z[0] + hs_bosons_z[1]).mass() / GeV;
                have_two_hs_bosons = true;
            }
            if (!have_two_hs_bosons) vetoEvent;


            FourMomentum beam_lab; beam_lab.setXYZE(0.0, 0.0, 1.0, 1.0);

            LorentzTransform boost_ZZ_rf;
            boost_ZZ_rf.setBetaVec(-fourvec_ZZ.betaVec());
            FourMomentum fourvec_ZZ_ZZrf = boost_ZZ_rf.transform(fourvec_ZZ);
            FourMomentum fourvec_Z1_ZZrf = boost_ZZ_rf.transform(fourvec_Z1);
            FourMomentum fourvec_lepminus_Z1_ZZrf = boost_ZZ_rf.transform(Z1_lepton_minus.mom());
            FourMomentum fourvec_lepminus_Z2_ZZrf = boost_ZZ_rf.transform(Z2_lepton_minus.mom());
            FourMomentum fourvec_Z2_ZZrf = boost_ZZ_rf.transform(fourvec_Z2);
            FourMomentum beam_ZZrf = boost_ZZ_rf.transform(beam_lab);

            cos_theta_star_ZZrf = cos(fourvec_Z1_ZZrf.p3().angle(beam_ZZrf.p3())); 
            cos_theta_star_VV = cos(fourvec_Z1_ZZrf.p3().angle(beam_ZZrf.p3())); 


            LorentzTransform boost_Z1_rf;
            boost_Z1_rf.setBetaVec(-fourvec_Z1.betaVec());
            FourMomentum fourvec_Z1_Z1rf = boost_Z1_rf.transform(fourvec_Z1);
            FourMomentum fourvec_Z2_Z1rf = boost_Z1_rf.transform(fourvec_Z2);
            FourMomentum fourvec_lepminus_Z1_rf = boost_Z1_rf.transform(Z1_lepton_minus.mom());
            FourMomentum beam_Z1rf = boost_Z1_rf.transform(beam_lab);
            cos_theta_star_Z1rf = cos(fourvec_lepminus_Z1_rf.p3().angle(beam_Z1rf.p3()));           
            cos_theta_star_V1 = cos(fourvec_lepminus_Z1_rf.p3().angle(beam_Z1rf.p3()));           

            LorentzTransform boost_Z2_rf;
            boost_Z2_rf.setBetaVec(-fourvec_Z2.betaVec());
            FourMomentum fourvec_Z1_Z2rf = boost_Z2_rf.transform(fourvec_Z1);
            FourMomentum fourvec_Z2_Z2rf = boost_Z2_rf.transform(fourvec_Z2);
            FourMomentum fourvec_lepminus_Z2_rf = boost_Z2_rf.transform(Z2_lepton_minus.mom());
            FourMomentum beam_Z2rf = boost_Z2_rf.transform(beam_lab);
            cos_theta_star_Z2rf = cos(fourvec_lepminus_Z2_rf.p3().angle(beam_Z2rf.p3()));   
            cos_theta_star_V2 = cos(fourvec_lepminus_Z2_rf.p3().angle(beam_Z2rf.p3()));   

            LorentzTransform boost_Z1_ZZrf;
            boost_Z1_ZZrf.setBetaVec(-fourvec_Z1_ZZrf.betaVec());
            FourMomentum ModHel_fourvec_lepminus_Z1_rf = boost_Z1_ZZrf.transform(fourvec_lepminus_Z1_ZZrf);
            FourMomentum ModHel_beam_Z1rf = boost_Z1_ZZrf.transform(beam_ZZrf);
            ModHel_cos_theta_star_Z1rf = cos(ModHel_fourvec_lepminus_Z1_rf.p3().angle(ModHel_beam_Z1rf.p3()));

            LorentzTransform boost_Z2_ZZrf;
            boost_Z2_ZZrf.setBetaVec(-fourvec_Z2_ZZrf.betaVec());
            FourMomentum ModHel_fourvec_lepminus_Z2_rf = boost_Z2_ZZrf.transform(fourvec_lepminus_Z2_ZZrf);
            FourMomentum ModHel_beam_Z2rf = boost_Z2_ZZrf.transform(beam_ZZrf);
            ModHel_cos_theta_star_Z2rf = cos(ModHel_fourvec_lepminus_Z2_rf.p3().angle(ModHel_beam_Z2rf.p3()));
  

            n_lep = leptons.size();
            n_jets = jets.size();


            tagjet1_pt = tag1_jet.pt();
            tagjet2_pt = tag2_jet.pt();
            
            tagjets_pt = (tag1_jet + tag2_jet).pT();
            
            tagjets_delta_pt = abs(tag1_jet.pt() - tag2_jet.pt());

            tagjets_delta_eta = abs(tag1_jet.eta() - tag2_jet.eta());
            
            tagjet1_eta = tag1_jet.eta();
            tagjet2_eta = tag2_jet.eta();
            tagjet1_phi = tag1_jet.phi();
            tagjet2_phi = tag2_jet.phi();
            tagjets_m = m_tagjets;

            pt_Z1_lep1 = pair_1_lep_1.pT();
            pt_Z1_lep2 = pair_1_lep_2.pT();
            pt_Z2_lep1 = pair_2_lep_1.pT();
            pt_Z2_lep2 = pair_2_lep_2.pT();
            m_ll_Z1 = m_ll_pair_1;
            m_ll_Z2 = m_ll_pair_2;
            m_VV = m_llll;
            m_diboson_ZZ = hs_diboson_mass;

            centrality_ZZjj = centrality_quadjj;
            pt_full_system = fourvec_lllljj.pT();



            // STOP REPLACING TEMPLATE

            _tt_SR->Fill();

            // save weights after cuts
            if (ev_nominal_weight >= 0){_c["pos_w_final"]->fill();}
            else {_c["neg_w_final"]->fill();}
        }

        /// Normalise histograms etc., after the run
        void finalize()
        {
            if (nsys < 1)
            {
                _tf->Write();

                // cout << "Number of systematics: " << nsys << endl;
                std::string cut_str = _cutflows.str();
                std::string cutflow_SR_file = getOption("OUTDIR") + "/cutflow_SR.txt";
                std::ofstream ofs_merged(cutflow_SR_file, std::ofstream::out);
                ofs_merged << cut_str;
                ofs_merged.close();


                const double sumw = sumOfWeights();
                const double cross_section_fb = crossSection() / femtobarn;
                std::cout << "Sum of weights: " << sumw << std::endl;
                std::cout << "Cross section in fb: " << crossSection() / femtobarn << std::endl;
                
                std::string SumW_str = std::to_string(sumOfWeights());
                std::string SumW_file = getOption("OUTDIR") + "/SumW.txt";
                std::ofstream ofs_sumW(SumW_file, std::ofstream::out);
                ofs_sumW << SumW_str;
                ofs_sumW.close();

                std::string X_section_str = std::to_string(cross_section_fb);
                std::string X_section_fb_file = getOption("OUTDIR") + "/X_section_fb.txt";
                std::ofstream ofs_xsec(X_section_fb_file, std::ofstream::out);
                ofs_xsec << X_section_str;
                ofs_xsec.close();



                for (auto &i_name : _hist_names)
                {
                    // std::cout << "normalizeing hist " << i_name <<" to 1; " ;
                    normalize(_h[i_name], 1.0);
                }
            }
            nsys = nsys + 1;
        }

        /// @}

        /// @name Histograms
        /// @{
        fastjet::Filter _trimmer;
        map<string, Histo1DPtr> _h;
        map<string, Histo2DPtr> _h2;
        // map<string, Profile1DPtr> _p;
        map<string, CounterPtr> _c;
        int _docut;
        Cut _electron_eta_cut;
        Cut _muon_eta_cut;
        Cut _el_eta_cut;
        Cut _electron_pt_cut;
        Cut _muon_pt_cut;

        json _jcuts;
        Cutflows _cutflows;
        std::vector<std::string> _hist_names;
        int EventNumber;
        double EventWeight;
        int nsys = 0;
        std::string operator_strings = "";
        int _label = -1;
        int _label_binary = -1;
        double cross_section_fb;

        unique_ptr<TFile> _tf;
        unique_ptr<TTree> _tt_SR;
        unique_ptr<TTree> _tt_bef_cut;


        int VBS_event;
        int n_lep, n_jets;


        double CS_V_cos_theta, cos_theta_star;

        double tagjet1_pt, tagjet2_pt, tagjets_pt, tagjets_delta_pt, tagjets_delta_eta;
        double tagjet1_eta, tagjet2_eta, tagjet1_phi, tagjet2_phi, tagjets_m;
        double pt_Z1_lep1, pt_Z1_lep2, pt_Z2_lep1, pt_Z2_lep2;
        double m_ll_Z1, m_ll_Z2, m_VV, m_diboson_ZZ;
        double centrality_ZZjj, pt_full_system;

        double cos_theta_star_Z1rf, cos_theta_star_Z2rf, cos_theta_star_ZZrf;
        double cos_theta_star_V1, cos_theta_star_V2, cos_theta_star_VV;
        double ModHel_cos_theta_star_Z1rf, ModHel_cos_theta_star_Z2rf;



        /// @}

        std::map<std::string, int> weightNameToIndex;
        std::map<std::string, int> weightNameToIndex_DynScale;
        std::map<int, std::string> indexToWeightName;
        std::map<int, std::string> indexToWeightName_DynScale;

        //double eventWeight;
    
        std::map<std::string, double *> varMap = {
            {"CS_V_cos_theta", &CS_V_cos_theta},
            {"cos_theta_star", &cos_theta_star},
            {"tagjet1_pt", &tagjet1_pt},
            {"tagjet2_pt", &tagjet2_pt},
            {"tagjets_pt", &tagjets_pt},
            {"tagjets_delta_pt", &tagjets_delta_pt},
            {"tagjets_delta_eta", &tagjets_delta_eta},
            {"tagjet1_eta", &tagjet1_eta},
            {"tagjet2_eta", &tagjet2_eta},
            {"tagjet1_phi", &tagjet1_phi},
            {"tagjet2_phi", &tagjet2_phi},
            {"tagjets_m", &tagjets_m},
            {"pt_Z1_lep1", &pt_Z1_lep1},
            {"pt_Z1_lep2", &pt_Z1_lep2},
            {"pt_Z2_lep1", &pt_Z2_lep1},
            {"pt_Z2_lep2", &pt_Z2_lep2},
            {"m_ll_Z1", &m_ll_Z1},
            {"m_ll_Z2", &m_ll_Z2},
            {"m_VV", &m_VV},
            {"m_diboson_ZZ", &m_diboson_ZZ},
            {"centrality_ZZjj", &centrality_ZZjj},
            {"pt_full_system", &pt_full_system},
            {"cos_theta_star_Z1rf", &cos_theta_star_Z1rf},
            {"cos_theta_star_Z2rf", &cos_theta_star_Z2rf},
            {"cos_theta_star_ZZrf", &cos_theta_star_ZZrf},
            {"cos_theta_star_V1", &cos_theta_star_V1},
            {"cos_theta_star_V2", &cos_theta_star_V2},
            {"cos_theta_star_VV", &cos_theta_star_VV},
            {"ModHel_cos_theta_star_Z1rf", &ModHel_cos_theta_star_Z1rf},
            {"ModHel_cos_theta_star_Z2rf", &ModHel_cos_theta_star_Z2rf},
        };

        std::map<std::string, int *> varMapInt = {
            {"n_jets", &n_jets},
            {"n_lep", &n_lep},
        };
    };


    RIVET_DECLARE_PLUGIN(ZZ_llll);

}

