// -*- C++ -*-
#include "Rivet/Analysis.hh"
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

namespace Rivet {


    /// @brief Add a short analysis description here
    class WmZ_llqq : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WmZ_llqq);

    // Calculate the number of charged tracks in a jet
    int CountChargedTracks(Jet& jet, double pTcut = 0.5) 
    {
        int n_trk = 0;
        for (const Particle& p : jet.particles()) {
            if (p.pT() > pTcut && p.isCharged()) {
                ++n_trk;
            }
        }
        return n_trk;
    }


    const double Centrality(const Jets& tagjets, const FourMomentum& jet_Vlep, const FourMomentum& jet_Vhad){
        double min_eta_tag_jet=std::min(tagjets[0].eta(), tagjets[1].eta());
        double max_eta_tag_jet=std::max(tagjets[0].eta(), tagjets[1].eta());
        double min_eta_Vjet=std::min(jet_Vlep.eta(), jet_Vhad.eta());
        double max_eta_Vjet=std::max(jet_Vlep.eta(), jet_Vhad.eta());
        double delta_eta_pos= max_eta_tag_jet - min_eta_Vjet;   
        double delta_eta_neg= max_eta_Vjet - min_eta_tag_jet;   
        return std::min(delta_eta_pos, delta_eta_neg);
    }

    double cosCollinsSoper(const FourMomentum& l1, const FourMomentum& l2) {
        const FourMomentum ll = l1 + l2;
        const double nom  = (l1.E() + l1.pz()) * (l2.E() - l2.pz()) - (l1.E() - l1.pz()) * (l2.E() + l2.pz());
        const double denom = ll.mass() * sqrt( sqr(ll.mass()) + sqr(ll.pt()) );
        return sign(ll.pz()) * safediv(nom, denom); // protect against division by zero, you never know...
    }

    
    FourMomentum RotateZ(double phi, double costheta, FourMomentum &p_)
    {
        FourMomentum p = p_;
        double sintheta = sqrt(1. - pow(costheta, 2.));
        double RotE     = p.E();
        double RotPx    = cos(phi) * costheta * p.px() - sin(phi) * p.py() + cos(phi) * sintheta * p.pz(); 
        double RotPy    = sin(phi) * costheta * p.px() + cos(phi) * p.py() + sin(phi) * sintheta * p.pz();
        double RotPz    = -sintheta * p.px() + costheta * p.pz();
        p.setXYZE(RotPx, RotPy, RotPz, RotE);
        return p;

    }


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
        std::string txt_dir = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/";
        
        std::string out_dir = getOption("OUTDIR");

        std::cout << "out_dir Rivet: " << getOption("OUTDIR") << std::endl;

        //cross_section_fb = crossSection()/femtobarn;

        if (out_dir.find("SM") != std::string::npos) _label = 0;
        else{
            if (out_dir.find("FM0") != std::string::npos) _label = 1;
            if (out_dir.find("FM2") != std::string::npos) _label = 2;
            if (out_dir.find("FS1") != std::string::npos) _label = 3;
            if (out_dir.find("FT0") != std::string::npos) _label = 4;
            if (out_dir.find("FT1") != std::string::npos) _label = 4;
            if (out_dir.find("FT5") != std::string::npos) _label = 5;
            if (out_dir.find("FT8") != std::string::npos) _label = 6;
        }

        if (out_dir.find("SM") != std::string::npos) _label_binary = 0;
        if (out_dir.find("QUAD") != std::string::npos) _label_binary = 1;

        std::string ntuple_dir = out_dir;
        
        _docut = 0; // most cuts on number of particles are always applied to avoid segfault
        if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
        std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

        std::string jsonfilestr =  txt_dir + "Cuts_def.json"; 
        std::cout << "++++++assume .json for this WmZ_llqq" << " is " << jsonfilestr << "\n";
        std::ifstream json_file(jsonfilestr);
        
        _jcuts = json::parse(json_file);
        std::cout << "++++++ to check json 1 var got photon pt min " << _jcuts["m_tagjets"] << "\n";
        _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));



        _el_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
        _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
        _electron_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_electron"])*GeV; 
        _muon_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_muon"])*GeV; 

        _jet_pt20_eta_cut_1 = ((Cuts::absetaIn(0.0, _jcuts["eta_jets"][0])) && (Cuts::pt > dbl(_jcuts["pt_jet_"][0])*GeV));
        _jet_pt30_eta_cut_2 = ((Cuts::absetaIn(_jcuts["eta_jets"][0], _jcuts["eta_jets"][1])) && (Cuts::pt > dbl(_jcuts["pt_jet_"][1])*GeV));

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
        VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, _jcuts["eta_jets_max"])));
        hadrons.addVetoOnThisFinalState(dressed_e);
        hadrons.addVetoOnThisFinalState(dressed_mu);
        declare(hadrons, "hadrons");
        FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetsfs, "jets");

        // Fat jet (merged regime) --> the anti-kT algorithm and a jet-radius parameter 1.0       
        FastJets fatjetsfs(hadrons, FastJets::ANTIKT, 1.0, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        // Define AntiKt10TruthTrimmedPtFrac5SmallR20Jets
        _trimmer = fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));
        declare(fatjetsfs, "fjets");
        declare(MissingMomentum(), "METFinder");



        // Merged histograms
        
        // plots common with others
        std::ifstream jet_hist_merged_file(txt_dir + "/Hists_bis/jet_hists_merged.json");      
        json jet_hist_merged = json::parse(jet_hist_merged_file);
        for (json::iterator it = jet_hist_merged.begin(); it != jet_hist_merged.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }

        // plots that are not in other ana
        std::ifstream ana_hist_merged_file(txt_dir + "/Hists_bis/2lepton_hists_merged.json");      
        json ana_hist_merged = json::parse(ana_hist_merged_file);
        for (json::iterator it = ana_hist_merged.begin(); it != ana_hist_merged.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }



        // Resolved histograms
        // plots common with others
        std::ifstream jet_hist_resolved_file(txt_dir + "/Hists_bis/jet_hists_resolved.json");      
        json jet_hist_resolved = json::parse(jet_hist_resolved_file);
        for (json::iterator it = jet_hist_resolved.begin(); it != jet_hist_resolved.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }

        // plots that are not in other ana
        std::ifstream ana_hist_resolved_file(txt_dir + "/Hists_bis/2lepton_hists_resolved.json");      
        json ana_hist_resolved = json::parse(ana_hist_resolved_file);
        for (json::iterator it = ana_hist_resolved.begin(); it != ana_hist_resolved.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }


        _tf = make_unique<TFile>(getOption("ROOTFILE", ntuple_dir+ "ntuple_rivet.root").c_str(), "RECREATE");
        _tt = make_unique<TTree>("Merged", "Rivet_physics");
        _tt->Branch("EventNumber", &merged_EventNumber);
        _tt->Branch("EventWeight", &merged_EventWeight);
        _tt->Branch("Label", &_label);
        _tt->Branch("Label_binary", &_label_binary);
        for (auto& var_ : varMap) {
            _tt->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        for (auto& var_ : varMapInt) {
            _tt->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        
        _tt_angle = make_unique<TTree>("Angle", "Rivet_physics");
        _tt_angle->Branch("EventWeight", &Angle_EventWeight);
        _tt_angle->Branch("Label", &_label);
        _tt_angle->Branch("cos_theta_star", &cos_theta_star);

        //counter for efficiency
        book(_c["pos_w_initial"],"pos_w_initial");
        book(_c["pos_w_final"],"pos_w_final");
        book(_c["pos_w_final_merged"],"pos_w_final_merged");
        book(_c["pos_w_final_resolved"],"pos_w_final_resolved");
        book(_c["neg_w_initial"],"neg_w_initial");
        book(_c["neg_w_final"],"neg_w_final");
        book(_c["neg_w_final_merged"],"neg_w_final_merged");
        book(_c["neg_w_final_resolved"],"neg_w_final_resolved");


        // Cut-flows merged region
        _cutflows_merged.addCutflow("WmZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "At_least_one_fjets","fjets_is_W/Z","Total_Merged_selec",});
        // Cut-flows resolved region
        _cutflows_resolved.addCutflow("WmZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "Failed_Merged_selection","At_least_two_signal_jets","signal_jets_pT","signal_mjj","M_jjj","Total_Resolved_selec",});

    }


    /// Perform the per-event analysis

    void analyze(const Event& event) {
        // save weights before cuts
        EventNumber= event.genEvent()->event_number();       
        double ev_nominal_weight =  event.weights()[0];


        //std::cout << "Type of event.genEvent()->event_number(): " << typeid(event.genEvent()->event_number()).name() << std::endl;
        //std::cout << "Value of event.genEvent()->event_number(): " << EvntNumber << std::endl;

        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
        else {_c["neg_w_initial"]->fill();}

        const Particles all_particles = event.allParticles();
        // Fill the truth histograms


        _cutflows_merged.fillinit();
        _cutflows_resolved.fillinit();
        // Retrieve dressed leptons, sorted by pT
        Particles e_stable;
        Particles mu_stable;
        if (_docut==1){
            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_el_eta_cut && _electron_pt_cut);
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _muon_pt_cut);
        }
        else{
            e_stable = apply<FinalState>(event, "e_stable").particlesByPt();
            mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt();
        }
        Particles leptons = e_stable + mu_stable; 
        // sort by pt as not clear if e+m is in order invidually not but necessary together
        std::sort(leptons.begin(), leptons.end(), [](Particle const &a, Particle const &b) {
            return a.pT() > b.pT(); // biggest pT will be first in array
            });



        int nlep = leptons.size();
        if (nlep !=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];
        FourMomentum lepton1 = lep1.mom();
        FourMomentum lepton2 = lep2.mom();

        const Particle lepton_minus = (lep1.charge() < 0) ? lep1 : lep2;
        const Particle lepton_plus = (lep1.charge() > 0) ? lep1 : lep2;


        // Cuts on the pT of the leptons
        if (_docut==1 && (leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();
        

        //const FourMomentum fourvec_Vlep = lep1.mom() + lep2.mom();
        FourMomentum fourvec_Vlep = lep1.mom() + lep2.mom();
        double m_Vlep = fourvec_Vlep.mass()/GeV;



        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(_jet_pt20_eta_cut_1 || _jet_pt30_eta_cut_2);
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);


        int n_jets = jets.size();
        if (n_jets < _jcuts["n_jets"])  vetoEvent;  
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext(); 

        // Retrive VBS tagging jets : look in opposite hemispheres and pair should have highest mjj

        Jets tag_jets;

        bool foundVBSJetPair = false; 
        double max_m_tag_jj = 0;
        int tag1_jet_index = -1 ,tag2_jet_index = -1;
        for (int i = 0; i < n_jets; i++) {
            const Jet& i_jet = jets[i];
            for (int j = i + 1; j < n_jets; j++) {
                const Jet& j_jet = jets[j];
                const double m_tag_jj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
                const double eta_prod = i_jet.eta()*j_jet.eta();
                if  (eta_prod < 0.0 && m_tag_jj>max_m_tag_jj){
                    max_m_tag_jj = m_tag_jj;
                    foundVBSJetPair = true;
                    tag1_jet_index = i;
                    tag2_jet_index = j;
                }
            }
        }
        if (tag2_jet_index < tag1_jet_index) swap(tag1_jet_index, tag2_jet_index); // organize tag jets by pt  
        if (!foundVBSJetPair)  vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();    

        tag_jets.push_back(jets[tag1_jet_index]);
        tag_jets.push_back(jets[tag2_jet_index]);


        const FourMomentum tag1_jet = jets[tag1_jet_index].mom();
        const FourMomentum tag2_jet = jets[tag2_jet_index].mom();

        if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"]) || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"]))) vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();



        const FourMomentum fourvec_tag_jj = tag1_jet+tag2_jet;

        const double m_tagjets = (fourvec_tag_jj).mass()/GeV;
        if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());

         


        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut

        //printf("Size n_fjets after: %d\n", fjets.size());

        Jets fjets_ = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > 200*GeV && Cuts::abseta < 2.0);
        //Jets fjets_ = apply<FastJets>(event, "fjets").jetsByPt();
        
        PseudoJets ljets;
        for (const Jet& fjet : fjets_) { ljets += _trimmer(fjet); }
        sort(ljets.begin(), ljets.end(), [](PseudoJet const &l, PseudoJet const &r) { return l.pt() > r.pt(); });
        Jets fjets;
        PseudoJets tr_ljets;
        for (const PseudoJet &lj : ljets) { 
            if((Jet(lj).pt() > 200*GeV && Jet(lj).abseta() < 2.0)) {
                fjets.push_back(Jet(lj));
                tr_ljets.push_back(lj); 
            }
        }

        //printf("\nSize n_fjets: %d\n", fjets.size());
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4);
        idiscardIfAnyDeltaRLess(fjets, leptons, 1.0);
        //printf("Size n_fjets after: %d\n", fjets.size());
        int n_fjets = fjets.size();

        // Create a new Jets object for the signal jets
        Jets sjets_sig;
                

        // Fill sjets_sig with the jets that are not the tagging jets
        for (int i = 0; i < n_jets; i++) {
            if (i != tag1_jet_index && i != tag2_jet_index) {
                sjets_sig.push_back(jets[i]);
            }
        }

        // Verify that sjets_sig are well sorted by momentum
        std::sort(sjets_sig.begin(), sjets_sig.end(), [](const Jet& a, const Jet& b) {
            return a.mom().pT() > b.mom().pT();
        });
        

        // Mean eta tagging jets for Zeppenfeld variable
        const double eta_tag_jet_mean = (tag1_jet.eta() + tag2_jet.eta())/2;


        LorentzTransform boost_Vlep_rf__;
        boost_Vlep_rf__.setBetaVec(-fourvec_Vlep.betaVec());
        FourMomentum fourvec_Vlep_Vlep_rf__ = boost_Vlep_rf__.transform(fourvec_Vlep);
        FourMomentum fourvec_lep_minus_Vlep_rf__ = boost_Vlep_rf__.transform(lepton_minus.mom());
        cos_theta_star = cos(fourvec_lep_minus_Vlep_rf__.p3().angle(fourvec_Vlep.p3()));

        Angle_EventWeight = ev_nominal_weight;
        _tt_angle->Fill();

            
        // Check if we are in the Merged signal region
        if (n_fjets > 0) {
            _cutflows_merged.fillnext();

            const double beta = 1;
            const fastjet::PseudoJet &LJet= tr_ljets[0];
            fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

            double recf3 = ECF3(LJet);
            double recf2 = ECF2(LJet);
            double recf1 = ECF1(LJet);
            double d2_fjets = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);

            const FourMomentum fourvec_fjets = fjets[0].mom();
            if ((_docut==1 && (fourvec_fjets.mass() >= _jcuts["m_fjet_WZ"][0] && fourvec_fjets.mass() <= _jcuts["m_fjet_WZ"][1]))||_docut==0) {
                _cutflows_merged.fillnext();
                
                // Total cutflow of the merged region
                _cutflows_merged.fillnext();

                Vector3 Axis_lab(1,1,1);
                Vector3 Axis_x_lab(1,0,0);Vector3 Axis_y_lab(0,1,0);Vector3 Axis_Vlep_lab(0,0,1);
                // Four vector of the QGC system in resolved region llJ
                const FourMomentum fourvec_fjets_Vlep = fourvec_Vlep + fjets[0].mom();
                FourMomentum fourvec_V = fjets[0].mom();
                // Polarization variable VZ system
                FourMomentum fourvec_VlepVhad = fourvec_Vlep + fjets[0].mom();

                LorentzTransform boost_VlepVhad;
                FourMomentum fourvec_VlepVhad_rotZ = fourvec_VlepVhad;
                FourMomentum taggjet1 = tag1_jet;
                FourMomentum taggjet2 = tag2_jet;

                double alpha_ = atan2(fourvec_VlepVhad_rotZ.py(), fourvec_VlepVhad_rotZ.px()); 
                fourvec_VlepVhad_rotZ = RotateZ(-alpha_, 1, fourvec_VlepVhad_rotZ);

                boost_VlepVhad.setBetaVec(-fourvec_VlepVhad_rotZ.betaVec());
                Vector3 BoostVZ_vec= -fourvec_VlepVhad.betaVec();
                const Vector3 unitboostVZvec = BoostVZ_vec.unit();

                FourMomentum fourvec_Vlep_CS = boost_VlepVhad.transform(RotateZ(-alpha_, 1,fourvec_Vlep));
                FourMomentum fourvec_VlepVhad_CS = boost_VlepVhad.transform(RotateZ(-alpha_, 1,fourvec_VlepVhad));

                // Calculate the theta and phi angles for each variable in the CS frame and fill the histograms
                merged_CS_V_cos_theta = abs(cos(fourvec_Vlep_CS.p3().theta()));
                //_h["merged_CS_V_cos_theta"]->fill(merged_CS_V_cos_theta);

                LorentzTransform boost_Vlep_rf;
                boost_Vlep_rf.setBetaVec(-fourvec_Vlep.betaVec());
                FourMomentum fourvec_Vlep_Vlep_rf = boost_Vlep_rf.transform(fourvec_Vlep);
                FourMomentum fourvec_V_Vlep_rf = boost_Vlep_rf.transform(fourvec_V);
                FourMomentum fourvec_lep_Vlep_rf = boost_Vlep_rf.transform(lep1.mom());
                FourMomentum fourvec_lepm_Vlep_rf = boost_Vlep_rf.transform(lepton_minus.mom());

                merged_cos_theta_star = cos(fourvec_lepm_Vlep_rf.p3().angle(fourvec_Vlep.p3()));
                //_h["merged_cos_theta_star"]->fill(merged_cos_theta_star);

                // Four vector of the Full system in resolved region 
                const FourMomentum fourvec_fjets_full = fourvec_Vlep + fjets[0].mom() + tag1_jet + tag2_jet;


                double ZeppMerged = 0.0;
                if (n_fjets > 1) {
                    const FourMomentum fourvec_fjets2 = fjets[1].mom();
                    ZeppMerged = abs(fourvec_fjets2.eta() - eta_tag_jet_mean);
                    //_h["merged_ZeppMerged"]->fill(ZeppMerged);
                }

                // We are in the Merged signal region
                // Fill in the histogram of the merged region
                merged_n_jets = n_jets;
                _h["merged_n_jets"]->fill(merged_n_jets);
                merged_tagjet1_pt = tag1_jet.pt();
                _h["merged_tagjet1_pt"]->fill(merged_tagjet1_pt);
                merged_tagjet2_pt = tag2_jet.pt();
                _h["merged_tagjet2_pt"]->fill(merged_tagjet2_pt);
                merged_tagjets_pt = fourvec_tag_jj.pT();
                _h["merged_tagjets_pt"]->fill(merged_tagjets_pt);
                merged_tagjets_delta_pt = abs(tag1_jet.pt()-tag2_jet.pt());
                _h["merged_tagjets_delta_pt"]->fill(merged_tagjets_delta_pt);

                _h["merged_tagjets_eta"]->fill(tag1_jet.eta()); _h["merged_tagjets_eta"]->fill(tag2_jet.eta());
                merged_tagjets_delta_eta = abs(tag1_jet.eta()-tag2_jet.eta());
                _h["merged_tagjets_delta_eta"]->fill(merged_tagjets_delta_eta);
                merged_tagjet1_eta = tag1_jet.eta();
                _h["merged_tagjet1_eta"]->fill(merged_tagjet1_eta);
                merged_tagjet2_eta = tag2_jet.eta();
                _h["merged_tagjet2_eta"]->fill(merged_tagjet2_eta);                
                _h["merged_tagjets_phi"]->fill(tag1_jet.phi()); _h["merged_tagjets_phi"]->fill(tag2_jet.phi());
                merged_tagjet1_phi = tag1_jet.phi(); merged_tagjet2_phi = tag2_jet.phi();
                _h["merged_tagjet1_phi"]->fill(tag1_jet.phi()); _h["merged_tagjet2_phi"]->fill(tag2_jet.phi());
                merged_tagjets_m = m_tagjets;
                _h["merged_tagjets_m"]->fill(merged_tagjets_m);
                merged_tagjets_eta = fourvec_tag_jj.eta();
                _h["merged_tagjets_eta"]->fill(merged_tagjets_eta);
                merged_tagjets_dy = dy_tagjets;
                _h["merged_tagjets_dy"]->fill(merged_tagjets_dy);
                merged_tagjets_dphi = deltaPhi(tag1_jet,tag2_jet);
                _h["merged_tagjets_dphi"]->fill(merged_tagjets_dphi);
                //lepton plots
                merged_n_lepton_stable = nlep;
                _h["merged_n_lepton_stable"]->fill(merged_n_lepton_stable);
                _h["merged_lepton_pt"]->fill(lep1.pT()); _h["merged_lepton_pt"]->fill(lep2.pT()); 
                merged_lepton_delta_pt = abs(lep1.pT()-lep2.pT());
                _h["merged_lepton_delta_pt"]->fill(merged_lepton_delta_pt);
                merged_lepton1_pt = lep1.pT();
                _h["merged_lepton1_pt"]->fill(merged_lepton1_pt);
                merged_lepton2_pt = lep2.pT();
                _h["merged_lepton2_pt"]->fill(merged_lepton2_pt); 
                _h["merged_lepton_eta"]->fill(lep1.eta()); _h["merged_lepton_eta"]->fill(lep2.eta()); 
                merged_lepton1_eta = lep1.eta();
                merged_lepton2_eta = lep2.eta();
                _h["merged_lepton1_eta"]->fill(merged_lepton1_eta); _h["merged_lepton2_eta"]->fill(merged_lepton2_eta);
                    
                merged_lepton_delta_eta = abs(lep1.eta()-lep2.eta());
                _h["merged_lepton_delta_eta"]->fill(merged_lepton_delta_eta);  

                merged_lepton_delta_phi = deltaPhi(lep1, lep2);
                _h["merged_lepton_delta_phi"]->fill(merged_lepton_delta_phi);   
                //ana-specific
                merged_Vlep_mass = m_Vlep;
                _h["merged_Vlep_mass"]->fill(merged_Vlep_mass);
                merged_Vlep_pt = fourvec_Vlep.pT();
                _h["merged_Vlep_pt"]->fill(merged_Vlep_pt);
                merged_Vlep_eta = fourvec_Vlep.eta();
                _h["merged_Vlep_eta"]->fill(merged_Vlep_eta);
                merged_Vlep_phi = fourvec_Vlep.phi();
                _h["merged_Vlep_phi"]->fill(merged_Vlep_phi);
                merged_Vlep_DR_tagjet1 = deltaR(fourvec_Vlep, tag1_jet);
                _h["merged_Vlep_DR_tagjet"]->fill(merged_Vlep_DR_tagjet1);
                merged_Vlep_DR_tagjet2 = deltaR(fourvec_Vlep, tag2_jet);
                _h["merged_Vlep_DR_tagjet"]->fill(merged_Vlep_DR_tagjet2);
                merged_Vlep_Dphi_tagjet1 = deltaPhi(fourvec_Vlep, tag1_jet);
                _h["merged_Vlep_Dphi_tagjet"]->fill(merged_Vlep_Dphi_tagjet1);
                merged_Vlep_Dphi_tagjet2 = deltaPhi(fourvec_Vlep, tag2_jet);
                _h["merged_Vlep_Dphi_tagjet"]->fill(merged_Vlep_Dphi_tagjet2);
                merged_Vlep_Deta_tagjet1 = abs(fourvec_Vlep.eta() - tag1_jet.eta());
                _h["merged_Vlep_Deta_tagjet"]->fill(merged_Vlep_Deta_tagjet1);
                merged_Vlep_Deta_tagjet2 = abs(fourvec_Vlep.eta() - tag2_jet.eta());
                _h["merged_Vlep_Deta_tagjet"]->fill(merged_Vlep_Deta_tagjet2);
                merged_lepton1_pids = lep1.pid();
                _h["merged_leptons_pids"]->fill(merged_lepton1_pids);
                merged_lepton2_pids = lep2.pid();
                _h["merged_leptons_pids"]->fill(merged_lepton2_pids);


                merged_fjet_eta = fourvec_fjets.eta();
                _h["merged_fjet_eta"]->fill(merged_fjet_eta);
                merged_fjet_pt = fourvec_fjets.pT();
                _h["merged_fjet_pt"]->fill(merged_fjet_pt);
                merged_fjet_D2 = d2_fjets;
                _h["merged_fjet_D2"]->fill(merged_fjet_D2);
                merged_fjet_n = n_fjets;
                _h["merged_fjet_n"]->fill(merged_fjet_n);
                merged_fjet_DR_lepton1 = deltaR(fourvec_fjets, lep1);
                merged_fjet_DR_lepton2 = deltaR(fourvec_fjets, lep2);
                _h["merged_fjet_DR_lepton1"]->fill(merged_fjet_DR_lepton1);
                _h["merged_fjet_DR_lepton2"]->fill(merged_fjet_DR_lepton2);
                merged_fjet_DR_tagjet1 = deltaR(fourvec_fjets, tag1_jet);
                _h["merged_fjet_DR_tagjet1"]->fill(merged_fjet_DR_tagjet1);
                merged_fjet_DR_tagjet2 = deltaR(fourvec_fjets, tag2_jet);
                _h["merged_fjet_DR_tagjet2"]->fill(merged_fjet_DR_tagjet2);
                merged_fjet_Dphi_tagjet1 = deltaPhi(fourvec_fjets, tag1_jet);
                _h["merged_fjet_Dphi_tagjet1"]->fill(merged_fjet_Dphi_tagjet1);
                merged_fjet_Dphi_tagjet2 = deltaPhi(fourvec_fjets, tag2_jet);
                _h["merged_fjet_Dphi_tagjet2"]->fill(merged_fjet_Dphi_tagjet2);
                merged_fjet_Deta_tagjet1 = abs(fourvec_fjets.eta() - tag1_jet.eta());
                _h["merged_fjet_Deta_tagjet1"]->fill(merged_fjet_Deta_tagjet1);
                merged_fjet_Deta_tagjet2 = abs(fourvec_fjets.eta() - tag2_jet.eta());
                _h["merged_fjet_Deta_tagjet2"]->fill(merged_fjet_Deta_tagjet2);
                merged_Vhad_DeltaR_Vlep = deltaR(fourvec_fjets, fourvec_Vlep);
                _h["merged_Vhad_DeltaR_Vlep"]->fill(merged_Vhad_DeltaR_Vlep);
                merged_Vhad_DeltaPhi_Vlep = deltaPhi(fourvec_fjets, fourvec_Vlep);
                _h["merged_Vhad_DeltaPhi_Vlep"]->fill(merged_Vhad_DeltaPhi_Vlep);
                merged_Vhad_DeltaEta_Vlep = abs(fourvec_fjets.eta() - fourvec_Vlep.eta());
                _h["merged_Vhad_DeltaEta_Vlep"]->fill(merged_Vhad_DeltaEta_Vlep);
                merged_fjet_mass = fourvec_fjets.mass();
                _h["merged_fjet_mass"]->fill(merged_fjet_mass);
                merged_VlepVhad_mass = fourvec_fjets_Vlep.mass();
                _h["merged_VlepVhad_mass"]->fill(merged_VlepVhad_mass);
                merged_VlepVhad_pt = fourvec_fjets_Vlep.pt();
                _h["merged_VlepVhad_pt"]->fill(merged_VlepVhad_pt);
                merged_VlepVhad_eta = fourvec_fjets_Vlep.eta();
                _h["merged_VlepVhad_eta"]->fill(merged_VlepVhad_eta);
                merged_Full_mass = fourvec_fjets_full.mass();
                _h["merged_Full_mass"]->fill(merged_Full_mass);
                merged_Full_pt = fourvec_fjets_full.pt();
                _h["merged_Full_pt"]->fill(merged_Full_pt);
 
                // Centrality variables
                merged_Centrality = Centrality(tag_jets, fourvec_Vlep, fourvec_fjets);
                merged_CentralityVhad = Centrality(tag_jets, fourvec_fjets, fourvec_fjets);
                merged_CentralityVlep = Centrality(tag_jets, fourvec_Vlep,fourvec_Vlep);
                merged_CentralityVlepVhad = Centrality(tag_jets, fourvec_fjets_Vlep, fourvec_fjets_Vlep);
                _h["merged_Centrality"]->fill(merged_Centrality);
                _h["merged_CentralityVhad"]->fill(merged_CentralityVhad);
                _h["merged_CentralityVlep"]->fill(merged_CentralityVlep);
                _h["merged_CentralityVlepVhad"]->fill(merged_CentralityVlepVhad);

                // Zeppenfeld variables
                merged_ZeppVlep = abs(fourvec_Vlep.eta() - eta_tag_jet_mean);
                merged_ZeppVhad = abs(fourvec_fjets.eta() - eta_tag_jet_mean);
                merged_ZeppVhadVlep = abs(fourvec_fjets_Vlep.eta() - eta_tag_jet_mean);
                _h["merged_ZeppVlep"]->fill(merged_ZeppVlep);
                _h["merged_ZeppVhad"]->fill(merged_ZeppVhad);
                _h["merged_ZeppVhadVlep"]->fill(merged_ZeppVhadVlep);
                
                //_h["merged_ZeppMerged"]->fill(ZeppMerged);
                merged_Ntrk_tagjets1 = CountChargedTracks(tag_jets[0]);
                _h["merged_Ntrk_tagjets"]->fill(merged_Ntrk_tagjets1);
                merged_Ntrk_tagjets2 = CountChargedTracks(tag_jets[1]);
                _h["merged_Ntrk_tagjets"]->fill(merged_Ntrk_tagjets2);
                merged_Ntrk_fjets = CountChargedTracks(fjets_[0]);
                _h["merged_Ntrk_fjets"]->fill(merged_Ntrk_fjets);


                merged_EventNumber = EventNumber;
                merged_EventWeight = ev_nominal_weight;

                _tt->Fill();

                if (ev_nominal_weight>=0){_c["pos_w_final_merged"]->fill();}
                else {_c["neg_w_final_merged"]->fill();}             

            } else {
                goto resolvedRegion;
            }
        // Check if we are in the Resolved signal region
        } else {
            resolvedRegion:
            // We are in the Resolved signal region
            // Check if we have at least two signal jets
            _cutflows_resolved.fillnext();
            if (sjets_sig.size() < 2) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            const FourMomentum signal_jet1 = sjets_sig[0].mom();
            const FourMomentum signal_jet2 = sjets_sig[1].mom(); 

            // Make a cut on the leading pT and subleading pT
            if ((_docut==1 && (signal_jet1.pT() < 40 || signal_jet2.pT() < 20)) ) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            

            
            // Make a cut on signal_mjj
            const FourMomentum fourvec_signal_jets = signal_jet1 + signal_jet2;
            const FourMomentum fourvec_Vhad = signal_jet1 + signal_jet2;
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if ((_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1]) )) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            
            const FourMomentum fourvec_signal_jets_Vlep = fourvec_Vlep + fourvec_signal_jets;

            // Four vector of the Full system in resolved region 
            const FourMomentum fourvec_signal_jets_full = fourvec_Vlep + fourvec_signal_jets + tag1_jet + tag2_jet;
            // More than 2 signal jets candidates
            FourMomentum fourvec_signal_jjj;

            // Centrality variables
            const double Centrality_resolved = Centrality(tag_jets, fourvec_Vlep, fourvec_signal_jets);
            const double CentralityVhad_resolved = Centrality(tag_jets, fourvec_signal_jets, fourvec_signal_jets);
            const double CentralityVlep_resolved = Centrality(tag_jets, fourvec_Vlep,fourvec_Vlep);
            const double CentralityVlepVhad_resolved = Centrality(tag_jets, fourvec_signal_jets_Vlep, fourvec_signal_jets_Vlep);

            // Zeppenfeld variables
            const double ZeppVlep_resolved = abs(fourvec_Vlep.eta() - eta_tag_jet_mean);
            const double ZeppVhad_resolved = abs(fourvec_signal_jets.eta() - eta_tag_jet_mean);
            const double ZeppVhadVlep_resolved = abs(fourvec_signal_jets_Vlep.eta() - eta_tag_jet_mean);

            double ZeppRes = 0.0;


            double topMass = 173.0; // Top quark mass in GeV
            double sdM = std::numeric_limits<double>::max(); // Initialize to maximum possible value
            FourMomentum signal_jet3;
            double signal_mjjj = 0.0;

            for (const Jet& jetcand : jets) {
                if (jetcand.mom() == signal_jet1 || jetcand.mom() == signal_jet2) continue;

                double mjjj = (signal_jet1 + signal_jet2 + jetcand.mom()).mass();
                if (abs(mjjj - topMass) < sdM) {
                    sdM = abs(mjjj - topMass);
                    signal_jet3 = jetcand.mom();
                }
            }

            if (_docut==1 && signal_jet3.isZero()) {
                // No suitable third jet found
                vetoEvent;
            } else {
                // Update fourvec_signal_jjj and signal_mjjj with the found third jet
                fourvec_signal_jjj = signal_jet1 + signal_jet2 + signal_jet3;
                signal_mjjj = fourvec_signal_jjj.mass();

                // Add condition for signal_mjjj to be above 220 GeV
                if (_docut==1 && signal_mjjj < 220.0) {
                    vetoEvent;
                }
            }



            _cutflows_resolved.fillnext();
            _cutflows_resolved.fillnext();
            // Fill in the histograms for the resolved region
            //jet plots
            _h["resolved_n_jets"]->fill(n_jets);
            _h["resolved_tagjet1_pt"]->fill(tag1_jet.pt());
            _h["resolved_tagjet2_pt"]->fill(tag2_jet.pt());
            _h["resolved_tagjets_pt"]->fill(fourvec_tag_jj.pT());
            _h["resolved_tagjets_delta_pt"]->fill(abs(tag1_jet.pt()-tag2_jet.pt()));
            _h["resolved_tagjets_eta"]->fill(tag1_jet.eta()); _h["resolved_tagjets_eta"]->fill(tag2_jet.eta());
            _h["resolved_tagjets_delta_eta"]->fill(abs(tag1_jet.eta()-tag2_jet.eta()));
            _h["resolved_tagjet1_eta"]->fill(tag1_jet.eta()); _h["resolved_tagjet2_eta"]->fill(tag2_jet.eta());
            _h["resolved_tagjets_phi"]->fill(tag1_jet.phi()); _h["resolved_tagjets_phi"]->fill(tag2_jet.phi());
            _h["resolved_tagjets_m"]->fill(m_tagjets);

            _h["resolved_tagjets_eta"]->fill(fourvec_tag_jj.eta());
            _h["resolved_tagjets_dy"]->fill(dy_tagjets);
            _h["resolved_tagjets_dphi"]->fill(deltaPhi(tag1_jet,tag2_jet));
            //lepton plots
            _h["resolved_n_lepton_stable"]->fill(nlep);
            _h["resolved_lepton_pt"]->fill(lep1.pT()); _h["resolved_lepton_pt"]->fill(lep2.pT()); 
            _h["resolved_lepton_delta_pt"]->fill(abs(lep1.pT()-lep2.pT()));
            _h["resolved_lepton1_pt"]->fill(lep1.pT()); _h["resolved_lepton2_pt"]->fill(lep2.pT()); 
            _h["resolved_lepton_eta"]->fill(lep1.eta()); _h["resolved_lepton_eta"]->fill(lep2.eta());
            _h["resolved_lepton_delta_eta"]->fill(abs(lep1.eta()-lep2.eta()));    
            //ana-specific
            _h["resolved_Vlep_mass"]->fill(m_Vlep);
            _h["resolved_Vlep_pt"]->fill(fourvec_Vlep.pT());
            _h["resolved_Vlep_eta"]->fill(fourvec_Vlep.eta());
            _h["resolved_Vlep_phi"]->fill(fourvec_Vlep.phi());
            _h["resolved_Vlep_DR_tagjet"]->fill(deltaR(fourvec_Vlep, tag1_jet));_h["resolved_Vlep_DR_tagjet"]->fill(deltaR(fourvec_Vlep, tag2_jet));
            _h["resolved_Vlep_Dphi_tagjet"]->fill(deltaPhi(fourvec_Vlep, tag1_jet));_h["resolved_Vlep_Dphi_tagjet"]->fill(deltaPhi(fourvec_Vlep, tag2_jet));
            _h["resolved_Vlep_Deta_tagjet"]->fill(abs(fourvec_Vlep.eta() - tag1_jet.eta()));_h["resolved_Vlep_Deta_tagjet"]->fill(abs(fourvec_Vlep.eta() - tag2_jet.eta()));
            _h["resolved_leptons_pids"]->fill(lep1.pid());_h["resolved_leptons_pids"]->fill(lep2.pid());

            _h["resolved_DR_min_lepton_tagjets1"]->fill(std::min(deltaR(lep1, tag1_jet), deltaR(lep2, tag1_jet)));
            _h["resolved_DR_min_lepton_tagjets2"]->fill(std::min(deltaR(lep1, tag2_jet), deltaR(lep2, tag2_jet)));
            _h["resolved_DR_min_lepton_sigjets1"]->fill(std::min(deltaR(signal_jet1, lep1),deltaR(signal_jet1, lep2)));
            _h["resolved_DR_min_lepton_sigjets2"]->fill(std::min(deltaR(signal_jet2, lep1),deltaR(signal_jet2, lep2)));

            
            _h["resolved_signal_jets_pt"]->fill(signal_jet1.pT());_h["resolved_signal_jets_pt"]->fill(signal_jet2.pT());
            _h["resolved_signal_jets_eta"]->fill(signal_jet1.eta());_h["resolved_signal_jets_eta"]->fill(signal_jet2.eta());
            _h["resolved_signal_jets_DeltaEta"]->fill(abs(signal_jet1.eta()-signal_jet2.eta()));
            _h["resolved_signal_jets_DeltaR"]->fill(deltaR(signal_jet2, signal_jet1));
            _h["resolved_signal_jets_DeltaPhi"]->fill(deltaPhi(signal_jet2, signal_jet1));
            _h["resolved_signal_jets_mass"]->fill(signal_mjj);

            _h["resolved_Vhad_DR_tagjet"]->fill(deltaR(fourvec_Vhad, tag1_jet));_h["resolved_Vhad_DR_tagjet"]->fill(deltaR(fourvec_Vhad, tag2_jet));
            _h["resolved_Vhad_Dphi_tagjet"]->fill(deltaPhi(fourvec_Vhad, tag1_jet));_h["resolved_Vhad_Deta_tagjet"]->fill(deltaPhi(fourvec_Vhad, tag2_jet));
            _h["resolved_Vhad_Deta_tagjet"]->fill(abs(fourvec_Vhad.eta() - tag1_jet.eta()));_h["resolved_Vhad_DR_tagjet"]->fill(abs(fourvec_Vhad.eta() - tag2_jet.eta()));


            _h["resolved_Vhad_DeltaR_Vlep"]->fill(deltaR(fourvec_Vhad, fourvec_Vlep));
            _h["resolved_Vhad_DeltaPhi_Vlep"]->fill(deltaPhi(fourvec_Vhad, fourvec_Vlep));
            _h["resolved_Vhad_DeltaEta_Vlep"]->fill(abs(fourvec_Vhad.eta() - fourvec_Vlep.eta()));

            _h["resolved_VlepVhad_mass"]->fill(fourvec_signal_jets_Vlep.mass());                
            _h["resolved_VlepVhad_pt"]->fill(fourvec_signal_jets_Vlep.pt());                
            _h["resolved_VlepVhad_eta"]->fill(fourvec_signal_jets_Vlep.eta());                
            _h["resolved_Full_mass"]->fill(fourvec_signal_jets_full.mass());                
            _h["resolved_Full_pt"]->fill(fourvec_signal_jets_full.pt());  
            _h["resolved_Centrality"]->fill(Centrality_resolved);
            _h["resolved_CentralityVhad"]->fill(CentralityVhad_resolved);
            _h["resolved_CentralityVlep"]->fill(CentralityVlep_resolved);
            _h["resolved_CentralityVlepVhad"]->fill(CentralityVlepVhad_resolved);
            _h["resolved_ZeppVlep"]->fill(ZeppVlep_resolved);
            _h["resolved_ZeppVhad"]->fill(ZeppVhad_resolved);
            _h["resolved_ZeppVhadVlep"]->fill(ZeppVhadVlep_resolved);
            
            _h["resolved_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[0])); _h["resolved_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[1]));
            _h["resolved_Ntrk_signal_jets"]->fill(CountChargedTracks(sjets_sig[0])); _h["resolved_Ntrk_signal_jets"]->fill(CountChargedTracks(sjets_sig[1]));
            _h["resolved_mjjj"]->fill(fourvec_signal_jjj.mass());
        
            // More than 2 signal jets candidates
            if (sjets_sig.size() > 2){
                _h["resolved_DeltaPt_sig_jet2_jet3"]-> fill(abs(signal_jet2.pT() - sjets_sig[2].mom().pT()));   
                
                _h["resolved_ZeppRes"]->fill(ZeppRes);
            }

            if (ev_nominal_weight>=0){_c["pos_w_final_resolved"]->fill();}
            else {_c["neg_w_final_resolved"]->fill();}

        } 
 
        // save weights after cuts
        if (ev_nominal_weight>=0){_c["pos_w_final"]->fill();}
        else {_c["neg_w_final"]->fill();}

    }

    /// Normalise histograms etc., after the run
    void finalize() {
        if(nsys<1) {
            _tf->Write();
        
            //cout << "Number of systematics: " << nsys << endl;
            std::string cut_merged_str = _cutflows_merged.str();
            std::string cutflow_merged_file = getOption("OUTDIR") + "/cutflow_merged.txt";

            std::ofstream ofs_merged (cutflow_merged_file, std::ofstream::out); 
            ofs_merged << cut_merged_str;
            ofs_merged.close();

            std::string cut_resolved_str = _cutflows_resolved.str();
            std::string cutflow_resolved_file = getOption("OUTDIR") + "/cutflow_resolved.txt";

            std::ofstream ofs_resolved (cutflow_resolved_file, std::ofstream::out); 
            ofs_resolved << cut_resolved_str;
            ofs_resolved.close();

 /*            double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
            double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);

            double pos_w_sum_final_merged = dbl(*_c["pos_w_final_merged"]);
            double neg_w_sum_final_merged = dbl(*_c["neg_w_final_merged"]);
            double pos_w_sum_final_resolved = dbl(*_c["pos_w_final_resolved"]);
            double neg_w_sum_final_resolved = dbl(*_c["neg_w_final_resolved"]); */


            // normalize all to 1 since in case of mostly negative weights not clear what it will do
/*             const double total_weight_initial = pos_w_sum_initial + neg_w_sum_initial;
            const double weight_merged = pos_w_sum_final_merged + neg_w_sum_final_merged;
            const double weight_resolved = pos_w_sum_final_resolved + neg_w_sum_final_resolved;


            const double xs = crossSection()/femtobarn;
            const double lumi= 139.0; // fb-1 */

            for (auto & i_name : _hist_names){ 
                //std::cout << "normalizeing hist " << i_name <<" to 1; " ;
                normalize(_h[i_name], 1.0);
            }
        }
        nsys=nsys+1;
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
    Cut _jet_pt20_eta_cut_1;
    Cut _jet_pt30_eta_cut_2;
    json _jcuts;
    Cutflows _cutflows_merged;
    Cutflows _cutflows_resolved;
    std::vector<std::string> _hist_names;
    int EventNumber; 
    int merged_EventNumber; 
    double merged_EventWeight; 
    int nsys=0;
    std::string operator_strings="";
    int _label=-1;
    int _label_binary=-1;
    double cross_section_fb;

    unique_ptr<TFile> _tf;
    unique_ptr<TTree> _tt;
    unique_ptr<TTree> _tt_angle;

    double cos_theta_star, Angle_EventWeight;

    double merged_CS_V_cos_theta,merged_cos_theta_star;

    double merged_tagjet1_pt,merged_tagjet2_pt,merged_tagjets_pt,merged_tagjets_delta_pt,merged_tagjets_delta_eta;
    double merged_tagjet1_eta,merged_tagjet2_eta,merged_tagjet1_phi,merged_tagjet2_phi,merged_tagjets_m;
    double merged_tagjets_eta,merged_tagjets_dy,merged_tagjets_dphi;
    double merged_lepton_delta_pt,merged_lepton1_pt,merged_lepton2_pt,merged_lepton1_eta,merged_lepton2_eta,merged_lepton_delta_eta,merged_lepton_delta_phi,merged_Vlep_mass;
    double merged_Vlep_pt,merged_Vlep_eta,merged_Vlep_phi,merged_Vlep_DR_tagjet1,merged_Vlep_DR_tagjet2;
    double merged_Vlep_Dphi_tagjet1,merged_Vlep_Dphi_tagjet2,merged_Vlep_Deta_tagjet1,merged_Vlep_Deta_tagjet2;
    double merged_fjet_eta,merged_fjet_pt,merged_fjet_D2,merged_fjet_DR_lepton1;
    double merged_fjet_DR_lepton2,merged_fjet_DR_tagjet1,merged_fjet_DR_tagjet2,merged_fjet_Dphi_tagjet1,merged_fjet_Dphi_tagjet2;
    double merged_fjet_Deta_tagjet1,merged_fjet_Deta_tagjet2,merged_Vhad_DeltaR_Vlep,merged_Vhad_DeltaPhi_Vlep,merged_Vhad_DeltaEta_Vlep;
    double merged_fjet_mass,merged_VlepVhad_mass,merged_VlepVhad_pt,merged_VlepVhad_eta,merged_Full_mass;
    double merged_Full_pt,merged_Centrality,merged_CentralityVhad,merged_CentralityVlep,merged_CentralityVlepVhad,merged_ZeppVlep,merged_ZeppVhad,merged_ZeppVhadVlep;
    int merged_n_jets,merged_n_lepton_stable,merged_fjet_n,merged_lepton1_pids,merged_lepton2_pids;
    int merged_Ntrk_tagjets1,merged_Ntrk_tagjets2,merged_Ntrk_tagjets,merged_Ntrk_fjets;
    /// @}
    std::map<std::string, double*> varMap = {
        {"merged_CS_V_cos_theta", &merged_CS_V_cos_theta},
        {"merged_cos_theta_star", &merged_cos_theta_star},
        {"merged_tagjet1_pt", &merged_tagjet1_pt},
        {"merged_tagjet2_pt", &merged_tagjet2_pt},
        {"merged_tagjets_pt", &merged_tagjets_pt},
        {"merged_tagjets_delta_pt", &merged_tagjets_delta_pt},
        {"merged_tagjets_delta_eta", &merged_tagjets_delta_eta},
        {"merged_tagjet1_eta", &merged_tagjet1_eta},
        {"merged_tagjet2_eta", &merged_tagjet2_eta},
        {"merged_tagjet1_phi", &merged_tagjet1_phi},
        {"merged_tagjet2_phi", &merged_tagjet2_phi},
        {"merged_tagjets_m", &merged_tagjets_m},
        {"merged_tagjets_eta", &merged_tagjets_eta},
        {"merged_tagjets_dy", &merged_tagjets_dy},
        {"merged_tagjets_dphi", &merged_tagjets_dphi},
        {"merged_lepton_delta_pt", &merged_lepton_delta_pt},
        {"merged_lepton1_pt", &merged_lepton1_pt},
        {"merged_lepton2_pt", &merged_lepton2_pt},
        {"merged_lepton1_eta", &merged_lepton1_eta},
        {"merged_lepton2_eta", &merged_lepton2_eta},
        {"merged_lepton_delta_eta", &merged_lepton_delta_eta},
        {"merged_lepton_delta_phi", &merged_lepton_delta_phi},
        {"merged_Vlep_mass", &merged_Vlep_mass},
        {"merged_Vlep_pt", &merged_Vlep_pt},
        {"merged_Vlep_eta", &merged_Vlep_eta},
        {"merged_Vlep_phi", &merged_Vlep_phi},
        {"merged_Vlep_DR_tagjet1", &merged_Vlep_DR_tagjet1},
        {"merged_Vlep_DR_tagjet2", &merged_Vlep_DR_tagjet2},
        {"merged_Vlep_Dphi_tagjet1", &merged_Vlep_Dphi_tagjet1},
        {"merged_Vlep_Dphi_tagjet2", &merged_Vlep_Dphi_tagjet2},
        {"merged_Vlep_Deta_tagjet1", &merged_Vlep_Deta_tagjet1},
        {"merged_Vlep_Deta_tagjet2", &merged_Vlep_Deta_tagjet2},
        {"merged_fjet_eta", &merged_fjet_eta},
        {"merged_fjet_pt", &merged_fjet_pt},
        {"merged_fjet_D2", &merged_fjet_D2},
        {"merged_fjet_DR_lepton1", &merged_fjet_DR_lepton1},
        {"merged_fjet_DR_lepton2", &merged_fjet_DR_lepton2},
        {"merged_fjet_DR_tagjet1", &merged_fjet_DR_tagjet1},
        {"merged_fjet_DR_tagjet2", &merged_fjet_DR_tagjet2},
        {"merged_fjet_Dphi_tagjet1", &merged_fjet_Dphi_tagjet1},
        {"merged_fjet_Dphi_tagjet2", &merged_fjet_Dphi_tagjet2},
        {"merged_fjet_Deta_tagjet1", &merged_fjet_Deta_tagjet1},
        {"merged_fjet_Deta_tagjet2", &merged_fjet_Deta_tagjet2},
        {"merged_Vhad_DeltaR_Vlep", &merged_Vhad_DeltaR_Vlep},
        {"merged_Vhad_DeltaPhi_Vlep", &merged_Vhad_DeltaPhi_Vlep},
        {"merged_Vhad_DeltaEta_Vlep", &merged_Vhad_DeltaEta_Vlep},
        {"merged_fjet_mass", &merged_fjet_mass},
        {"merged_VlepVhad_mass", &merged_VlepVhad_mass},
        {"merged_VlepVhad_pt", &merged_VlepVhad_pt},
        {"merged_VlepVhad_eta", &merged_VlepVhad_eta},
        {"merged_Full_mass", &merged_Full_mass},
        {"merged_Full_pt", &merged_Full_pt},
        {"merged_Centrality", &merged_Centrality},
        {"merged_CentralityVhad", &merged_CentralityVhad},
        {"merged_CentralityVlep", &merged_CentralityVlep},
        {"merged_CentralityVlepVhad", &merged_CentralityVlepVhad},
        {"merged_ZeppVlep", &merged_ZeppVlep},
        {"merged_ZeppVhad", &merged_ZeppVhad},
        {"merged_ZeppVhadVlep", &merged_ZeppVhadVlep}
    };
    std::map<std::string, int*> varMapInt = {
        {"merged_n_jets", &merged_n_jets},
        {"merged_n_lepton_stable", &merged_n_lepton_stable},
        {"merged_fjet_n", &merged_fjet_n},
        {"merged_lepton1_pids", &merged_lepton1_pids},
        {"merged_lepton2_pids", &merged_lepton2_pids},
        {"merged_Ntrk_tagjets1", &merged_Ntrk_tagjets1},
        {"merged_Ntrk_tagjets2", &merged_Ntrk_tagjets2},
        {"merged_Ntrk_tagjets", &merged_Ntrk_tagjets},
        {"merged_Ntrk_fjets", &merged_Ntrk_fjets},
    };

    };


    RIVET_DECLARE_PLUGIN(WmZ_llqq);

}