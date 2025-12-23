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
    class WpZ_lllv : public Analysis
    {
    public:
        /// Constructor
        RIVET_DEFAULT_ANALYSIS_CTOR(WpZ_lllv);



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

            std::string jsonfilestr = txt_dir + "Cuts/WZ_lllv_cuts.json";
            std::cout << "++++++assume .json for this WpZ_lllv" << " is " << jsonfilestr << "\n";
            std::ifstream json_file(jsonfilestr);

            _jcuts = json::parse(json_file);
            std::cout << "++++++ to check json 1 var got photon pt min " << _jcuts["m_tagjets"] << "\n";

            // CUT Lepton

            _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                        (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
            _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
            _lepton_stage1_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton"])*GeV;      


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

            IdentifiedFinalState nu_id;
            nu_id.acceptNeutrinos();
            PromptFinalState neutrinos(nu_id);
            neutrinos.acceptTauDecays(false);
            declare(neutrinos, "Neutrinos");


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


            _cutflows.addCutflow("WpZ_lllv_selections", {"three_leptons_one_neutrino","res_shape_step_1", "mZ_mismatch",
                                                        "pt_W_lepton", "m_T_W", "dR_leptons", "n_jets","n_b_jets",
                                                        "tagj_opposite_hemisph_and_pt","m_tagjets" ,"two_WZ_HS"});

            


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
            Particles e_stable;
            Particles mu_stable;

            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton_stage1_pt_cut);
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton_stage1_pt_cut);
            
            Particles leptons = e_stable + mu_stable;
            // sort by pt as not clear if e+m is in order invidually not but necessary together
            std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
                          return a.pT() > b.pT(); // biggest pT will be first in array
                      });

            int nlep = leptons.size();
            // CHANGING THE TEMPLATE FROM HERE
            if (nlep < _jcuts["n_lepton_stable"])  vetoEvent; // meaning both are e,mu and not tau

            const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
            if (neutrinos.size() < _jcuts["n_neutrinos"]) vetoEvent;
            _cutflows.fillnext();

            /////
            // resonant shape algo copied from their previous routine
            // https://gitlab.cern.ch/atlas-physics/pmg/rivet-routines/-/blob/master/STDM-2017-23_WZ_VBS_36ifb/ATLAS_2018_I1711223.cc
            ////
            int i, j, k;
            double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
            double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
            double WeightZ1, WeightZ2, WeightZ3;
            double WeightW1, WeightW2, WeightW3;
            double M1, M2, M3;
            double WeightTotal1, WeightTotal2, WeightTotal3;

            int icomb=0;
            // try Z pair  of leptons 01

            if ( (leptons[0].pid() ==-(leptons[1].pid()))  && (leptons[2].pid()*neutrinos[0].pid()< 0) && (leptons[2].abspid()==neutrinos[0].abspid()-1)) {
                MassZ01 = (leptons[0].momentum() + leptons[1].momentum()).mass();
                MassW2 = (leptons[2].momentum() + neutrinos[0].momentum()).mass();
                icomb = 1;
            }
            // try Z pair of leptons 02
            if ( (leptons[0].pid()==-(leptons[2].pid()))  && (leptons[1].pid()*neutrinos[0].pid()< 0) && (leptons[1].abspid()==neutrinos[0].abspid()-1)) {
                MassZ02 = (leptons[0].momentum() + leptons[2].momentum()).mass();
                MassW1 = (leptons[1].momentum() + neutrinos[0].momentum()).mass();
                icomb = 2;
            }
            // try Z pair of leptons 12
            if ( (leptons[1].pid()==-(leptons[2].pid())) && (leptons[0].pid()*neutrinos[0].pid()< 0) && (leptons[0].abspid()==neutrinos[0].abspid()-1)) {
                MassZ12 = (leptons[1].momentum() + leptons[2].momentum()).mass();
                MassW0 = (leptons[0].momentum() + neutrinos[0].momentum()).mass();
                icomb = 3;
            }

            if (icomb<=0)  vetoEvent;
            _cutflows.fillnext();

            WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
            WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
            WeightTotal1 = WeightZ1*WeightW1;
            M1 = -1*WeightTotal1;

            WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
            WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
            WeightTotal2 = WeightZ2*WeightW2;
            M2 = -1*WeightTotal2;

            WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
            WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
            WeightTotal3 = WeightZ3*WeightW3;
            M3 = -1*WeightTotal3;

            if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ) {
                i = 0; j = 1; k = 2;
            }
            if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ) {
                i = 0; j = 2; k = 1;
            }
            if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ) {
                i = 1; j = 2; k = 0;
            }

            const Particle Z_lepton1 = leptons[i];
            const Particle Z_lepton2 = leptons[j];

            const FourMomentum& Z_lep_1 = leptons[i].momentum();
            const FourMomentum& Z_lep_2 = leptons[j].momentum();

            const Particle Z_lepton_minus = (Z_lepton1.charge() < 0) ? Z_lepton1 : Z_lepton2;
            const Particle W_lepton = leptons[k];

            const FourMomentum& W_lep  = leptons[k].momentum();
            const FourMomentum Wboson   = leptons[k].momentum()+neutrinos[0].momentum();
            const FourMomentum& Z_boson   = leptons[i].momentum()+leptons[j].momentum();

            const FourMomentum fourvec_Z = Z_lep_1 + Z_lep_2;
            const FourMomentum fourvec_W = W_lep + neutrinos[0].momentum();
            const FourMomentum fourvec_WZ = fourvec_Z + fourvec_W;



            if (fabs(Z_boson.mass()/GeV - MZ_PDG)>_jcuts["abs_diff_m_z"]) vetoEvent;
            _cutflows.fillnext();

            if (W_lep.pT() < _jcuts["pt_W_lepton"]) vetoEvent;
            _cutflows.fillnext();
            
            const double pt_neutrino = neutrinos[0].pt();
            const double d_phi_MET_lep = deltaPhi(neutrinos[0].phi(), W_lep.phi());
            const double m_T_W = sqrt( 2*W_lep.pT()*pt_neutrino*(1 - cos(d_phi_MET_lep)) );
            if (m_T_W<_jcuts["m_T_W"]) vetoEvent;
            _cutflows.fillnext();

            // Retrieve clustered jets, sorted by pT, with a minimum pT cut
            Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); 
            Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
            // Remove all jets within certain dR of a dressed lepton
            idiscardIfAnyDeltaRLess(jets, leptons, _jcuts["dR_lepton_jet"]);
            idiscardIfAnyDeltaRLess(btagging_jets, leptons, _jcuts["dR_lepton_jet"]);

            double dR_Z_lep_1_lep_2 = deltaR(Z_lep_1, Z_lep_2);  
            double dR_Z_lep_1_lep_W = deltaR(Z_lep_1, W_lep);  
            double dR_Z_lep_2_lep_W = deltaR(Z_lep_2, W_lep);  
            if (dR_Z_lep_1_lep_2 < _jcuts["dR_lepton1_lepton2_Z"]) vetoEvent; 
            if (dR_Z_lep_1_lep_W < _jcuts["dR_Z_leptons_W_lepton"]) vetoEvent; 
            if (dR_Z_lep_2_lep_W < _jcuts["dR_Z_leptons_W_lepton"]) vetoEvent; 
            _cutflows.fillnext();

            int n_jets = jets.size();
            if (n_jets < _jcuts["n_jets"]) vetoEvent;
            _cutflows.fillnext();

            int n_b_jets = count(btagging_jets, hasBTag());
            if (n_b_jets>_jcuts["n_b_jets"]) vetoEvent;
            _cutflows.fillnext();


            FourMomentum tag1_jet = jets[0].mom();
            FourMomentum tag2_jet;
            bool foundVBSJetPair = false;
            for (const Jet& i_jet : jets) {
                if(i_jet.pT() > dbl(_jcuts["pt_tagjet1"])*GeV && i_jet.eta()*tag1_jet.eta() < 0.) {
                tag2_jet = i_jet.mom();
                foundVBSJetPair = true;
                break;
                }
            }
            if (!foundVBSJetPair)  vetoEvent;
            _cutflows.fillnext();

            if ((tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"])*GeV)) vetoEvent;

            const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
            if (m_tagjets<dbl(_jcuts["m_tagjets"])*GeV) vetoEvent;
            _cutflows.fillnext();

            const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));

            const double m_T_WZ_term1 =  pow(Z_lep_1.pT() + Z_lep_2.pT() + W_lep.pT() + pt_neutrino, 2);
            const double m_T_WZ_term2 =  pow(Z_lep_1.px() + Z_lep_2.px() + W_lep.px() + neutrinos[0].px(), 2);
            const double m_T_WZ_term3 =  pow(Z_lep_1.py() + Z_lep_2.py() + W_lep.py() + neutrinos[0].py(), 2);
            const double m_T_WZ = sqrt(m_T_WZ_term1 - m_T_WZ_term2 - m_T_WZ_term3);

            // do clipping - in case with 2+ w/z to avoid much work take two w with highest pt
            std::vector<FourMomentum> hs_bosons = {};
            const Particles all_particles = event.allParticles();
            for(const Particle& p_rivet : all_particles){
                ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
                int status = p_hepmc->status();
                if (abs(status)==23 or abs(status)==22){
                int i_pid = p_hepmc->pid();
                FourMomentum i_mom = p_hepmc->momentum();
                if (abs(i_pid)==24 or abs(i_pid)==23){hs_bosons.push_back(i_mom);
                }
                }
            }
            // biggest pT will be first in array
            std::sort(hs_bosons.begin(), hs_bosons.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); });
            bool have_two_hs_bosons = false;
            double hs_diboson_mass = 0.0;
            int n_hs = hs_bosons.size();
            if (n_hs > 1){
                hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
                have_two_hs_bosons = true;
            }
            if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have wz somehow
            _cutflows.fillnext();

            // STOP REPLACING TEMPLATE
            FourMomentum beam_lab; beam_lab.setXYZE(0.0, 0.0, 1.0, 1.0);

            LorentzTransform boost_WZ_rf;
            boost_WZ_rf.setBetaVec(-fourvec_WZ.betaVec());
            FourMomentum fourvec_WZ_WZrf = boost_WZ_rf.transform(fourvec_WZ);
            FourMomentum fourvec_Z_WZrf = boost_WZ_rf.transform(fourvec_Z);
            FourMomentum fourvec_W_WZrf = boost_WZ_rf.transform(fourvec_W);
            FourMomentum fourvec_lepminus_Z_WZrf = boost_WZ_rf.transform(Z_lepton_minus.mom());
            FourMomentum fourvec_lep_W_WZrf = boost_WZ_rf.transform(W_lepton.mom());
            FourMomentum beam_WZrf = boost_WZ_rf.transform(beam_lab);

            cos_theta_star_WZrf = cos(fourvec_Z_WZrf.p3().angle(beam_WZrf.p3())); 
            cos_theta_star_VV = cos(fourvec_Z_WZrf.p3().angle(beam_WZrf.p3())); 


            LorentzTransform boost_Z_rf;
            boost_Z_rf.setBetaVec(-fourvec_Z.betaVec());
            FourMomentum fourvec_Z_Zrf = boost_Z_rf.transform(fourvec_Z);
            FourMomentum fourvec_W_Zrf = boost_Z_rf.transform(fourvec_W);
            FourMomentum fourvec_lepminus_Z_rf = boost_Z_rf.transform(Z_lepton_minus.mom());
            FourMomentum beam_Zrf = boost_Z_rf.transform(beam_lab);
            cos_theta_star_Zrf = cos(fourvec_lepminus_Z_rf.p3().angle(beam_Zrf.p3())); 
            cos_theta_star_V2 = cos(fourvec_lepminus_Z_rf.p3().angle(beam_Zrf.p3())); 

            LorentzTransform boost_W_rf;
            boost_W_rf.setBetaVec(-fourvec_W.betaVec());
            FourMomentum fourvec_Z_Wrf = boost_W_rf.transform(fourvec_Z);
            FourMomentum fourvec_W_Wrf = boost_W_rf.transform(fourvec_W);
            FourMomentum fourvec_lep_W_rf = boost_W_rf.transform(W_lepton.mom());
            FourMomentum beam_Wrf = boost_W_rf.transform(beam_lab);
            cos_theta_star_Wrf = cos(fourvec_lep_W_rf.p3().angle(beam_Wrf.p3()));   
            cos_theta_star_V1 = cos(fourvec_lep_W_rf.p3().angle(beam_Wrf.p3()));   

            LorentzTransform boost_Z_WZrf;
            boost_Z_WZrf.setBetaVec(-fourvec_Z_WZrf.betaVec());
            FourMomentum ModHel_fourvec_lepminus_Z_rf = boost_Z_WZrf.transform(fourvec_lepminus_Z_WZrf);
            FourMomentum ModHel_beam_Zrf = boost_Z_WZrf.transform(beam_WZrf);
            ModHel_cos_theta_star_Zrf = cos(ModHel_fourvec_lepminus_Z_rf.p3().angle(ModHel_beam_Zrf.p3()));

            LorentzTransform boost_W_WZrf;
            boost_W_WZrf.setBetaVec(-fourvec_W_WZrf.betaVec());
            FourMomentum ModHel_fourvec_lep_W_rf = boost_W_WZrf.transform(fourvec_lep_W_WZrf);
            FourMomentum ModHel_beam_Wrf = boost_W_WZrf.transform(beam_WZrf);
            ModHel_cos_theta_star_Wrf = cos(ModHel_fourvec_lep_W_rf.p3().angle(ModHel_beam_Wrf.p3()));
  

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

            pt_Z_lep1 = Z_lep_1.pT();
            pt_Z_lep2 = Z_lep_2.pT();
            pt_W_lep = W_lep.pT();
            pt_v = pt_neutrino;
            eta_Z_lep1 = Z_lep_1.eta();
            eta_Z_lep2 = Z_lep_2.eta();
            eta_W_lep = W_lep.eta();
            delta_eta_Z_leps = abs(Z_lep_1.eta() -  Z_lep_2.eta());
            m_VV = (Z_boson + Wboson).mass()/GeV;
            n_bjets = n_b_jets;
            mT_W = m_T_W;
            m_diboson_truth = hs_diboson_mass;

            
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
        Cut _lepton_stage1_pt_cut;

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

        double MZ_PDG = 91.1876;
        double MW_PDG = 80.385;
        double GammaZ_PDG = 2.4952;
        double GammaW_PDG = 2.085;

        unique_ptr<TFile> _tf;
        unique_ptr<TTree> _tt_SR;
        unique_ptr<TTree> _tt_bef_cut;


        int VBS_event;
        int n_lep, n_jets;

        double tagjet1_pt, tagjet2_pt, tagjets_pt, tagjets_delta_pt, tagjets_delta_eta;
        double tagjet1_eta, tagjet2_eta, tagjet1_phi, tagjet2_phi, tagjets_m;
        double pt_Z_lep1, pt_Z_lep2, pt_W_lep, pt_v;
        double eta_Z_lep1, eta_Z_lep2, eta_W_lep, delta_eta_Z_leps;
        double m_VV, n_bjets, mT_W, m_diboson_truth;

        double cos_theta_star_WZrf, cos_theta_star_Zrf, cos_theta_star_Wrf;
        double cos_theta_star_VV, cos_theta_star_V1, cos_theta_star_V2;
        double ModHel_cos_theta_star_Zrf, ModHel_cos_theta_star_Wrf;



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
            {"pt_Z_lep1", &pt_Z_lep1},
            {"pt_Z_lep2", &pt_Z_lep2},
            {"pt_W_lep", &pt_W_lep},
            {"pt_v", &pt_v},
            {"eta_Z_lep1", &eta_Z_lep1},
            {"eta_Z_lep2", &eta_Z_lep2},
            {"eta_W_lep", &eta_W_lep},
            {"delta_eta_Z_leps", &delta_eta_Z_leps},
            {"m_VV", &m_VV},
            {"n_bjets", &n_bjets},
            {"mT_W", &mT_W},
            {"m_diboson_truth", &m_diboson_truth},
            {"cos_theta_star_WZrf", &cos_theta_star_WZrf},
            {"cos_theta_star_Zrf", &cos_theta_star_Zrf},
            {"cos_theta_star_Wrf", &cos_theta_star_Wrf},
            {"cos_theta_star_VV", &cos_theta_star_VV},
            {"cos_theta_star_V1", &cos_theta_star_V1},
            {"cos_theta_star_V2", &cos_theta_star_V2},
            {"ModHel_cos_theta_star_Zrf", &ModHel_cos_theta_star_Zrf},
            {"ModHel_cos_theta_star_Wrf", &ModHel_cos_theta_star_Wrf},
        };

        std::map<std::string, int *> varMapInt = {
            {"n_jets", &n_jets},
            {"n_lep", &n_lep},
        };
    };


    RIVET_DECLARE_PLUGIN(WpZ_lllv);

}

