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
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/Cutflow.hh"
#include "Rivet/Tools/RivetHepMC.hh"
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
    class Zy_lly : public Analysis
    {
    public:
        /// Constructor
        RIVET_DEFAULT_ANALYSIS_CTOR(Zy_lly);



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

            std::string jsonfilestr = "Cuts/Zy_lly_cuts.json";
            std::cout << "++++++assume .json for this Zy_lly" << "is " << jsonfilestr << "\n";
            std::ifstream json_file(jsonfilestr);
            
            _jcuts = json::parse(json_file);
            std::cout << "++++++ to check 1 var got photon pt min" << _jcuts["pt_photon"] << "\n";
            _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                        (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
            _photon_eta_cut = (Cuts::absetaIn(_jcuts["eta_photon"][0][0], _jcuts["eta_photon"][0][1])) || 
                                        (Cuts::absetaIn(_jcuts["eta_photon"][1][0], _jcuts["eta_photon"][1][1]));
            _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
            _lepton2_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton2"])*GeV;      

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

            // FS excluding the leading big pt photons, muons and neutrinos to calculate cone energy of bit pt photons
            // so basically electrons + jets left

            // photons as separate particles for final state, not dressing
            DirectFinalState photons(Cuts::abspid == PID::PHOTON);
            declare(photons, "photons");

            DirectFinalState bare_muon_for_ph_iso(Cuts::abspid == PID::MUON);
            DressedLeptons dressed_muons_for_ph_iso(photons_for_dressing, bare_muon_for_ph_iso, 0.1, _muon_eta_cut &&_lepton2_pt_cut);
            declare(dressed_muons_for_ph_iso, "muons_for_ph_iso");
            VetoedFinalState vfs;
            vfs.addVetoOnThisFinalState(photons);
            vfs.addVetoOnThisFinalState(dressed_muons_for_ph_iso);
            vfs.addVetoOnThisFinalState(InvisibleFinalState());
            declare(vfs, "ElAndJetsForPhotonIsoCalc");
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
            _cutflows.addCutflow("Zy_lly_selections", {"n_lep_ok_pt2_eta", "lep_pid_charge", "lep_pt1", "m_ll", "have_iso_photons_ok_pt_eta",
                                        "m_ll_plus_m_lly", "n_jets", "pt_tagjet1_2","m_tagjets","dy_tagjets","centrality_lly",
                                        "n_gap_jets"});
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
            Particles e_stable;
            Particles mu_stable;

            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton2_pt_cut);
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton2_pt_cut);
            
            Particles leptons = e_stable + mu_stable;
            // will ahve either e or mu pair so don't need to sort combined e+mu array
            Particles leptons_stable = e_stable + mu_stable; 


            // sort by pt as not clear if e+m is in order invidually not but necessary together
            std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
                          return a.pT() > b.pT(); // biggest pT will be first in array
                      });

            int nlep = leptons.size();
            // CHANGING THE TEMPLATE FROM HERE
            int nlep_stable = leptons_stable.size();
            if (nlep_stable!=_jcuts["n_lepton_stable"])  vetoEvent; // meaning both are e,mu and not tau
            _cutflows.fillnext();

            const Particle& lep1 = leptons_stable[0];
            const Particle& lep2 = leptons_stable[1];

            const Particle Z_lepton_minus = (lep1.charge() < 0) ? lep1 : lep2;          

            if (lep1.pid()+lep2.pid()!=0) vetoEvent; // want opposite charge leptons of same fravour
            _cutflows.fillnext();
            if ((lep1.pT() < dbl(_jcuts["pt_lepton1"])*GeV)) vetoEvent;
            _cutflows.fillnext();
            const FourMomentum fourvec_ll = lep1.mom() + lep2.mom(); 
            const FourMomentum fourvec_Z = lep1.mom() + lep2.mom(); 
            const double m_ll = fourvec_ll.mass()/GeV;
            if (m_ll < dbl(_jcuts["m_ll"])*GeV) vetoEvent;
            _cutflows.fillnext();

            //photons
            Particles photons = apply<FinalState>(event, "photons").particlesByPt(_photon_eta_cut);
            if (photons.empty())  vetoEvent;
            //photon cone calculation and photon OR with leptons 
            Particles isolated_photons;
            std::vector<double> cone_to_photon_fracs = {};
            Particles cone_sum_particles = apply<VetoedFinalState>(event, "ElAndJetsForPhotonIsoCalc").particles();
            for (const Particle &i_high_pt_photon : photons) {
                // check photon isolation
                double i_high_pt_photon_cone_E = 0.0;
                for (const Particle &i_p_cone : cone_sum_particles) {
                if (deltaR(i_high_pt_photon, i_p_cone) < 0.2) { // etcone20
                    i_high_pt_photon_cone_E += i_p_cone.Et();
                }
                }
                double i_cone_to_photon_frac = i_high_pt_photon_cone_E / i_high_pt_photon.pT(); 
                if (i_cone_to_photon_frac > _jcuts["cone_frac_photon"])  continue;
                if (any(leptons_stable, deltaRLess(i_high_pt_photon, 0.4))) continue;
                isolated_photons += i_high_pt_photon;
                cone_to_photon_fracs.push_back(i_cone_to_photon_frac);
            }
            if (isolated_photons.empty())  vetoEvent;
            const Particle& lead_iso_photon = isolated_photons[0];
            if (lead_iso_photon.pT() < dbl(_jcuts["pt_photon"])*GeV) vetoEvent;
            _cutflows.fillnext();
            
            const FourMomentum fourvec_lly = lep1.mom() + lep2.mom() + lead_iso_photon.mom();
            const double m_lly = fourvec_lly.mass();
            if ((m_ll + m_lly)<=dbl(_jcuts["m_lly"])*GeV)  vetoEvent;
            _cutflows.fillnext();

            // Retrieve clustered jets, sorted by pT, with a minimum pT cut
            Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); // then will do cut on two leading pt>50
            // Remove all jets within certain dR of a dressed lepton
            idiscardIfAnyDeltaRLess(jets, leptons_stable, 0.3);
            idiscardIfAnyDeltaRLess(jets, isolated_photons, 0.4);

            int n_jets = jets.size();
            if (n_jets < _jcuts["n_jets"]) vetoEvent;
            _cutflows.fillnext();

            const FourMomentum tag1_jet = jets[0].mom();
            const FourMomentum tag2_jet = jets[1].mom();
            if ((tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV || tag2_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV)) vetoEvent; 
            _cutflows.fillnext();

            const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
            if (m_tagjets<dbl(_jcuts["m_tagjets"])*GeV) vetoEvent;
            _cutflows.fillnext();

            const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
            if (dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
            _cutflows.fillnext();

            const double centrality_lly = fabs(0.5 * (fourvec_lly.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
            if (centrality_lly > _jcuts["centrality_lly"])  vetoEvent;
            _cutflows.fillnext();

            int n_gap_jets = 0;
            for (int i = 0; i < n_jets; i++) {
                const double i_jet_rap = jets[i].rap();
                if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
            }
            if (n_gap_jets > _jcuts["n_gap_jets"])  vetoEvent;
            _cutflows.fillnext();
            // STOP REPLACING TEMPLATE

            FourMomentum beam_lab; beam_lab.setXYZE(0.0, 0.0, 1.0, 1.0);

            LorentzTransform boost_Z_rf;
            boost_Z_rf.setBetaVec(-fourvec_Z.betaVec());
            FourMomentum fourvec_Z_Zrf = boost_Z_rf.transform(fourvec_Z);
            FourMomentum fourvec_lepminus_Z_rf = boost_Z_rf.transform(Z_lepton_minus.mom());
            FourMomentum beam_Zrf = boost_Z_rf.transform(beam_lab);
            cos_theta_star_Zrf = cos(fourvec_lepminus_Z_rf.p3().angle(beam_Zrf.p3()));           


            n_lep = leptons.size();
            n_jets = jets.size();


            tagjet1_pt = tag1_jet.pt();
            tagjet2_pt = tag2_jet.pt();
            
            tagjets_pt = (tag1_jet + tag2_jet).pT();
            
            tagjets_delta_pt = abs(tag1_jet.pt() - tag2_jet.pt());

            tagjets_delta_eta = abs(tag1_jet.eta() - tag2_jet.eta());
            tagjets_delta_phi = abs(deltaPhi(tag1_jet, tag2_jet));
            
            tagjet1_eta = tag1_jet.eta();
            tagjet2_eta = tag2_jet.eta();
            tagjet1_phi = tag1_jet.phi();
            tagjet2_phi = tag2_jet.phi();
            tagjets_m = m_tagjets;

            pt_lepton1 = lep1.pT();
            eta_lepton1 = lep1.eta();
            pt_lepton2 = lep2.pT();
            eta_lepton2 = lep2.eta();


            pt_photon = lead_iso_photon.pT();
            eta_photon = lead_iso_photon.eta();
            phi_photon = lead_iso_photon.phi();
            n_photons_iso = isolated_photons.size();

            mass_ll = m_ll;
            mass_lly = m_lly;
            m_VV = m_lly; 
            
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
        Cut _photon_eta_cut;
        Cut _muon_eta_cut;
        Cut _lepton2_pt_cut;

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

        double tagjet1_pt, tagjet2_pt, tagjets_pt, tagjets_delta_pt, tagjets_delta_eta;
        double tagjets_delta_phi, tagjet1_eta, tagjet2_eta, tagjet1_phi, tagjet2_phi, tagjets_m;
        double pt_lepton1, eta_lepton1, pt_lepton2, eta_lepton2;
        double pt_photon, eta_photon, phi_photon;
        int n_photons_iso;
        double mass_ll, mass_lly, m_VV;
        double cos_theta_star_Zrf;



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
            {"tagjets_delta_phi", &tagjets_delta_phi},
            {"tagjet1_eta", &tagjet1_eta},
            {"tagjet2_eta", &tagjet2_eta},
            {"tagjet1_phi", &tagjet1_phi},
            {"tagjet2_phi", &tagjet2_phi},
            {"tagjets_m", &tagjets_m},
            {"pt_lepton1", &pt_lepton1},
            {"eta_lepton1", &eta_lepton1},
            {"pt_lepton2", &pt_lepton2},
            {"eta_lepton2", &eta_lepton2},
            {"pt_photon", &pt_photon},
            {"eta_photon", &eta_photon},
            {"phi_photon", &phi_photon},
            {"mass_ll", &mass_ll},
            {"mass_lly", &mass_lly},
            {"m_VV", &m_VV},
            {"cos_theta_star_Zrf", &cos_theta_star_Zrf}
        };

        std::map<std::string, int *> varMapInt = {
            {"n_jets", &n_jets},
            {"n_lep", &n_lep},
            {"n_photons_iso", &n_photons_iso}
        };
    };


    RIVET_DECLARE_PLUGIN(Zy_lly);

}

