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
    class WpWp_lvlv : public Analysis
    {
    public:
        /// Constructor
        RIVET_DEFAULT_ANALYSIS_CTOR(WpWp_lvlv);



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

            std::ifstream json_file("Cuts/ssWW_lvlv_cuts.json");
            _jcuts = json::parse(json_file);
            std::cout << "++++++ to check json 1 var got photon pt min" << _jcuts["abs_diff_m_z"] << "\n";

            _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_electron"][0][0], _jcuts["eta_electron"][0][1])) ||
                                (Cuts::absetaIn(_jcuts["eta_electron"][1][0], _jcuts["eta_electron"][1][1]));
            _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_muon"]);

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


            declare(MissingMomentum(), "METFinder");

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
        // Cut-flows
            _cutflows.addCutflow("ssWW_lvlv_selections", {"pt_MET","n_lep","lep_charge_dR", "m_ll","n_jets_bjets",
                                                        "m_tagjets","dy_tagjets","two_HS_bosons"});
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
            const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
            const double scalar_MET = METfinder.missingPt()/GeV;
            if (scalar_MET<_jcuts["pt_MET"]) vetoEvent;
            _cutflows.fillnext();

            const FourMomentum fourvec_MET = METfinder.missingMomentum();

            // Retrieve clustered jets, sorted by pT, with a minimum pT cut
            Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 4.5);
            Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);

            // Retrieve dressed leptons, sorted by pT
            Particles e_stable;
            Particles mu_stable;
            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut);
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut);

            // Remove jet if lepton too close
            idiscardIfAny(jets, e_stable, [](const Jet& jet, const Particle& lep) {
                return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() > 0.5);
            });
            // otherwise remove lepton
            idiscardIfAny(e_stable, jets, [](const Particle& lep, const Jet& jet) {
                return deltaR(lep, jet) < 0.4 && (lep.pT() / jet.pT() < 0.5);
            });

            Particles leptons = e_stable + mu_stable;
            idiscardIfAnyDeltaRLess(btagging_jets, leptons, 0.2);
            // sort by pt as not clear if e+m is in order invidually not but necessary together
            std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
                return a.pT() > b.pT(); // biggest pT will be first in array
            });

            int nlep_stable = leptons.size();
            int cut_nlep = _jcuts["n_lepton_stable"];
            if (nlep_stable!=cut_nlep)  vetoEvent;
            _cutflows.fillnext();

            const Particle& lep1 = leptons[0];
            const Particle& lep2 = leptons[1];
            if (lep1.pT()<dbl(_jcuts["pt_lepton12"])*GeV) vetoEvent;
            if (lep2.pT()<dbl(_jcuts["pt_lepton12"])*GeV) vetoEvent;
            double eta_ee = _jcuts["eta_electron_ee"];
            if (fabs(lep1.pid())==11 && fabs(lep2.pid())==11 && (lep1.eta()>eta_ee || lep2.eta()>eta_ee)) vetoEvent;


            if (lep1.charge() * lep2.charge() < 0) vetoEvent; // want same charge leptons
            if (deltaR(lep1, lep2) < 0.2)  vetoEvent;
            _cutflows.fillnext();

            const FourMomentum fourvec_ll = lep1.mom() + lep2.mom();
            const double m_ll = fourvec_ll.mass()/GeV;
            if (m_ll<_jcuts["m_ll"]) vetoEvent;
            _cutflows.fillnext();
        
            const double lep1_abs_pid = lep1.pid();
            const double lep2_abs_pid = lep2.pid();
            double m_z = 91.18;
            if ((lep1_abs_pid == 11 && lep2_abs_pid == 11 && fabs(m_ll - m_z) < _jcuts["abs_diff_m_z"])) vetoEvent;

            int njets = jets.size();
            if (njets < _jcuts["n_jets"])  vetoEvent;
            int nbtags = count(btagging_jets, hasBTag());
            if (nbtags>_jcuts["n_b_jets"]) vetoEvent;
            _cutflows.fillnext();

            const FourMomentum tag1_jet = jets[0].mom();
            const FourMomentum tag2_jet = jets[1].mom();
            if ((tag1_jet.pT()<_jcuts["pt_tagjet1"] || tag2_jet.pT()<_jcuts["pt_tagjet2"])) vetoEvent;

            const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
            if (m_tagjets<_jcuts["m_tagjetsSR"]) vetoEvent;
            _cutflows.fillnext();

            const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
            if (dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
            _cutflows.fillnext();
        

            const double m_T = (fourvec_MET + fourvec_ll).mass()/GeV;

            // do clipping - in case with 2+ w to avoid much work take two w with highest pt
            std::vector<FourMomentum> hs_bosons = {}; // also can be WZ for WZ CR
            const Particles all_particles = event.allParticles();
            for(const Particle& p_rivet : all_particles){
                ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
                int status = p_hepmc->status();
                if (abs(status)==23 or abs(status)==22){
                    int i_pid = p_hepmc->pid();
                    FourMomentum i_mom = p_hepmc->momentum();
                    if (abs(i_pid) == 24 or abs(i_pid) == 23){
                        hs_bosons.push_back(i_mom);
                    }
                }
            }
            std::sort(hs_bosons.begin(), hs_bosons.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); }); // biggest pT will be first in array
            bool have_two_hs_bosons = false;
            double hs_diboson_mass = 0.0;
            int n_z_hs = hs_bosons.size();
            if (n_z_hs > 1){
                hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
                have_two_hs_bosons = true;
            }
            if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have ww somehow
            _cutflows.fillnext();

            costhetastar = leptons.size()>1 ? fabs(tanh((leptons[0].eta() - leptons[1].eta()) / 2)) : -0.2;

            // STOP REPLACING TEMPLATE

            n_lep = leptons.size();
            n_jets = jets.size();


            tagjet1_pt = tag1_jet.pt();
            tagjet2_pt = tag2_jet.pt();
            
            tagjets_pt = (tag1_jet + tag2_jet).pT();
            
            tagjets_delta_pt = abs(tag1_jet.pt() - tag2_jet.pt());

            tagjets_delta_eta = abs(tag1_jet.eta() - tag2_jet.eta());
            tagjets_delta_phi = deltaPhi(tag1_jet, tag2_jet);
            
            tagjet1_eta = tag1_jet.eta();
            tagjet2_eta = tag2_jet.eta();
            tagjet1_phi = tag1_jet.phi();
            tagjet2_phi = tag2_jet.phi();

            tagjets_m = m_tagjets;

            met= scalar_MET;
            n_b_tags = nbtags;

            pt_lep1 = lep1.pT();
            pt_lep2 = lep2.pT();
            delta_eta_lepton = abs(lep1.eta() - lep2.eta());
            delta_phi_lepton = deltaPhi(lep1, lep2);
            m_VV = m_T;

            m_diboson = hs_diboson_mass;
            
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


        int VBS_event, n_jets, n_lep;
        int n_b_tags;


        double CS_V_cos_theta, cos_theta_star;

        double tagjet1_pt, tagjet2_pt, tagjets_pt, tagjets_delta_pt, tagjets_delta_eta, tagjets_delta_phi;
        double tagjet1_eta, tagjet2_eta, tagjet1_phi, tagjet2_phi, tagjets_m;
        double met, m_VV;
        double m_diboson;
        double pt_lep1, pt_lep2;
        double delta_eta_lepton, delta_phi_lepton;
        double costhetastar;


        /// @}

        std::map<std::string, int> weightNameToIndex;
        std::map<std::string, int> weightNameToIndex_DynScale;
        std::map<int, std::string> indexToWeightName;
        std::map<int, std::string> indexToWeightName_DynScale;

        //double eventWeight;
    
        std::map<std::string, double *> varMap = {
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
            {"met", &met},
            {"pt_lep1", &pt_lep1},
            {"pt_lep2", &pt_lep2},
            {"delta_eta_lepton", &delta_eta_lepton},
            {"delta_phi_lepton", &delta_phi_lepton},
            {"costhetastar", &costhetastar},
            {"m_diboson", &m_diboson},
            {"m_VV", &m_VV},

        };

        std::map<std::string, int *> varMapInt = {
            {"n_jets", &n_jets},
            {"n_lep", &n_lep},
        };
    };


    RIVET_DECLARE_PLUGIN(WpWp_lvlv);

}

