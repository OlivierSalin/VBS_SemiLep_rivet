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
    double CountChargedTracks(Jet& jet, double pTcut = 0.5) 
    {
        int n_trk = 0;
        for (const Particle& p : jet.particles()) {
            if (p.pT() > pTcut && p.isCharged()) {
                ++n_trk;
            }
        }
        return static_cast<double>(n_trk);
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


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
        std::string txt_dir = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/";
        
        std::string out_dir = getOption("OUTDIR");
        
        _docut = 0; // most cuts on number of particles are always applied to avoid segfault
        if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
        std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

        std::string jsonfilestr =  txt_dir + "WpZ_llqq_cuts.json"; 
        std::cout << "++++++assume .json for this WpZ_llqq" << " is " << jsonfilestr << "\n";
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


        // Before cut histograms

        std::ifstream bef_cut_file(txt_dir + "/Hists/WpZ_llqq_hists_bef_cuts.json");      
        json bef_cut = json::parse(bef_cut_file);
        for (json::iterator it = bef_cut.begin(); it != bef_cut.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }



        // Merged histograms
        
        // plots common with others
        std::ifstream jet_hist_merged_file(txt_dir + "/Hists/jet_hists_merged.json");      
        json jet_hist_merged = json::parse(jet_hist_merged_file);
        for (json::iterator it = jet_hist_merged.begin(); it != jet_hist_merged.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        std::ifstream lep_hist_merged_file(txt_dir + "/Hists/lepton_hists_merged.json");      
        json lep_hist_merged = json::parse(lep_hist_merged_file);
        for (json::iterator it = lep_hist_merged.begin(); it != lep_hist_merged.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_merged_file(txt_dir + "/Hists/WpZ_llqq_hists_merged.json");      
        json ana_hist_merged = json::parse(ana_hist_merged_file);
        for (json::iterator it = ana_hist_merged.begin(); it != ana_hist_merged.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }

        // Resolved histograms
        // plots common with others
        std::ifstream jet_hist_resolved_file(txt_dir + "/Hists/jet_hists_resolved.json");      
        json jet_hist_resolved = json::parse(jet_hist_resolved_file);
        for (json::iterator it = jet_hist_resolved.begin(); it != jet_hist_resolved.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        std::ifstream lep_hist_resolved_file(txt_dir + "/Hists/lepton_hists_resolved.json");      
        json lep_hist_resolved = json::parse(lep_hist_resolved_file);
        for (json::iterator it = lep_hist_resolved.begin(); it != lep_hist_resolved.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_resolved_file(txt_dir + "/Hists/WpZ_llqq_hists_resolved.json");      
        json ana_hist_resolved = json::parse(ana_hist_resolved_file);
        for (json::iterator it = ana_hist_resolved.begin(); it != ana_hist_resolved.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }


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
        double ev_nominal_weight =  event.weights()[0];
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

        // Cuts on the number of the leptons

        int nlep = leptons.size();
        _h["bef_cut_nlep"]->fill(nlep);
        if (nlep !=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];

/*         const Particle parent_lep1 = GetParent(lep1);
        const Particle parent_lep2 = GetParent(lep2); */


        // Cuts on the pT of the leptons
        if (_docut==1 && (leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();
        

        const FourMomentum fourvec_ll = lep1.mom() + lep2.mom();
        double m_ll = fourvec_ll.mass()/GeV;

        // Cuts on the invariant mass of the lepton pair   
        //if (_docut==1 && (m_ll>_jcuts["m_ll"][1] || m_ll<_jcuts["m_ll"][0])) vetoEvent;
        //_cutflows_merged.fillnext();
        //_cutflows_resolved.fillnext();

        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(_jet_pt20_eta_cut_1 || _jet_pt30_eta_cut_2);
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);


        int n_jets = jets.size();
        _h["bef_cut_n_jets"]->fill(n_jets);
        if (n_jets < _jcuts["n_jets"])  vetoEvent;  
        _h["bef_cutTag_n_jets"]->fill(n_jets);
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
        _h["aft_cutTag_n_jets"]->fill(n_jets);          

        tag_jets.push_back(jets[tag1_jet_index]);
        tag_jets.push_back(jets[tag2_jet_index]);


        const FourMomentum tag1_jet = jets[tag1_jet_index].mom();
        const FourMomentum tag2_jet = jets[tag2_jet_index].mom();

        if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"]) || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"]))) vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();



        const FourMomentum fourvec_tag_jj = tag1_jet+tag2_jet;

        const double m_tagjets = (fourvec_tag_jj).mass()/GeV;
        _h["bef_cut_m_tag_jets"]->fill(m_tagjets);
        if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());

         


        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut
/*         Jets fjets_ = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > 100*GeV && Cuts::abseta < 2.0);
        //printf("\nSize n_fjets: %d\n", fjets.size());
        Jets fjets;
        for (const Jet &LargeRjet : fjets_) {
            if (LargeRjet.pt() > 200*GeV) {
                fjets.push_back(LargeRjet);
            }
        }
        _h["bef_cutDR_n_fjets"]->fill(fjets.size());
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4); */
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
        _h["bef_cutDR_n_fjets"]->fill(fjets.size());
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4);
        idiscardIfAnyDeltaRLess(fjets, leptons, 1.0);
        //printf("Size n_fjets after: %d\n", fjets.size());
        int n_fjets = fjets.size();
        _h["bef_cut_n_fjets"]->fill(fjets.size());

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
        _h["resolved_DR_min_lepton_tagjets1"]->fill(std::min(deltaR(lep1, tag1_jet), deltaR(lep2, tag1_jet)));
        _h["resolved_DR_min_lepton_tagjets2"]->fill(std::min(deltaR(lep1, tag2_jet), deltaR(lep2, tag2_jet)));



            
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

            
                // Four vector of the QGC system in resolved region llJ
                const FourMomentum fourvec_fjets_ll = lep1.mom() + lep2.mom() + fjets[0].mom();
                // Four vector of the Full system in resolved region 
                const FourMomentum fourvec_fjets_full = lep1.mom() + lep2.mom() + fjets[0].mom() + tag1_jet + tag2_jet;


                // Centrality variables
                const double Centrality_merged = Centrality(tag_jets, fourvec_ll, fourvec_fjets);
                const double CentralityV_merged = Centrality(tag_jets, fourvec_fjets, fourvec_fjets);
                const double CentralityZ_merged = Centrality(tag_jets, fourvec_ll,fourvec_ll);
                const double CentralityZV_merged = Centrality(tag_jets, fourvec_fjets_ll, fourvec_fjets_ll);

                // Zeppenfeld variables
                const double ZeppZ_merged = abs(fourvec_ll.eta() - eta_tag_jet_mean);
                const double ZeppV_merged = abs(fourvec_fjets.eta() - eta_tag_jet_mean);
                const double ZeppZV_merged = abs(fourvec_fjets_ll.eta() - eta_tag_jet_mean);

                double ZeppMerged = 0.0;
                if (n_fjets > 1) {
                    const FourMomentum fourvec_fjets2 = fjets[1].mom();
                    ZeppMerged = abs(fourvec_fjets2.eta() - eta_tag_jet_mean);
                }

                // We are in the Merged signal region
                // Fill in the histogram of the merged region
                _h["merged_n_jets"]->fill(n_jets);
                _h["merged_pt_tagjet1"]->fill(tag1_jet.pt());
                _h["merged_pt_tagjet2"]->fill(tag2_jet.pt());
                _h["merged_eta_tagjets"]->fill(tag1_jet.eta()); _h["merged_eta_tagjets"]->fill(tag2_jet.eta());
                _h["merged_eta_tagjet1"]->fill(tag1_jet.eta()); _h["merged_eta_tagjet2"]->fill(tag2_jet.eta());
                _h["merged_eta_tagjet1_2"]->fill(tag1_jet.eta()+ tag2_jet.eta());
                _h["merged_phi_tagjets"]->fill(tag1_jet.phi()); _h["merged_phi_tagjets"]->fill(tag2_jet.phi());
                _h["merged_m_tagjets"]->fill(m_tagjets);
                _h["merged_dy_tagjets"]->fill(dy_tagjets);
                _h["merged_dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
                //lepton plots
                _h["merged_n_lepton_stable"]->fill(nlep);
                _h["merged_pt_lepton"]->fill(lep1.pT()); _h["merged_pt_lepton"]->fill(lep2.pT()); 
                _h["merged_pt_lepton1"]->fill(lep1.pT()); _h["merged_pt_lepton2"]->fill(lep2.pT()); 
                _h["merged_eta_lepton"]->fill(lep1.eta()); _h["merged_eta_lepton"]->fill(lep2.eta());    
                //ana-specific
                _h["merged_m_ll"]->fill(m_ll);
                _h["merged_pt_ll"]->fill(fourvec_ll.pT());
                _h["merged_leptons_pids"]->fill(lep1.pid());_h["merged_leptons_pids"]->fill(lep2.pid());
                _h["merged_pt_tagjets"]->fill(fourvec_tag_jj.pT());
                _h["merged_eta_fjet"]->fill(fourvec_fjets.eta());
                _h["merged_pt_fjet"]->fill(fourvec_fjets.pT());
                _h["merged_D2_fjet"]->fill(d2_fjets);
                _h["merged_n_fjets"]->fill(n_fjets);
                _h["merged_DR_lepton_fjets"]->fill(std::min(deltaR(fourvec_fjets, lep1),deltaR(fourvec_fjets, lep2)));
                _h["merged_DR_ll_fjets"]->fill(deltaR(fourvec_fjets, fourvec_ll));
                _h["merged_DR_tagJet_fjets"]->fill(deltaR(fourvec_fjets, tag1_jet));_h["merged_DR_tagJet_fjets"]->fill(deltaR(fourvec_fjets, tag2_jet));
                _h["merged_mass_fjets"]->fill(fourvec_fjets.mass());
                _h["merged_mass_WZ"]->fill(fourvec_fjets_ll.mass());                
                _h["merged_pt_WZ"]->fill(fourvec_fjets_ll.pt());                
                _h["merged_mass_Full"]->fill(fourvec_fjets_full.mass());                
                _h["merged_pt_Full"]->fill(fourvec_fjets_full.pt());
                _h["merged_Centrality"]->fill(Centrality_merged);
                _h["merged_CentralityVhad"]->fill(CentralityV_merged);
                _h["merged_CentralityZ"]->fill(CentralityZ_merged);
                _h["merged_CentralityZV"]->fill(CentralityZV_merged);
                _h["merged_ZeppZ"]->fill(ZeppZ_merged);
                _h["merged_ZeppV"]->fill(ZeppV_merged);
                _h["merged_ZeppZV"]->fill(ZeppZV_merged);
                //_h["merged_ZeppMerged"]->fill(ZeppMerged);
                _h["merged_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[0])); _h["merged_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[1]));  
                _h["merged_Ntrk_fjets"]->fill(CountChargedTracks(fjets[0]));

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
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if ((_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1]) )) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            
            const FourMomentum fourvec_signal_jets_ll = lep1.mom() + lep2.mom() + fourvec_signal_jets;

            // Four vector of the Full system in resolved region 
            const FourMomentum fourvec_signal_jets_full = lep1.mom() + lep2.mom() + fourvec_signal_jets + tag1_jet + tag2_jet;
            // More than 2 signal jets candidates
            FourMomentum fourvec_signal_jjj;

            // Centrality variables
            const double Centrality_resolved = Centrality(tag_jets, fourvec_ll, fourvec_signal_jets);
            const double CentralityV_resolved = Centrality(tag_jets, fourvec_signal_jets, fourvec_signal_jets);
            const double CentralityZ_resolved = Centrality(tag_jets, fourvec_ll,fourvec_ll);
            const double CentralityZV_resolved = Centrality(tag_jets, fourvec_signal_jets_ll, fourvec_signal_jets_ll);

            // Zeppenfeld variables
            const double ZeppZ_resolved = abs(fourvec_ll.eta() - eta_tag_jet_mean);
            const double ZeppV_resolved = abs(fourvec_signal_jets.eta() - eta_tag_jet_mean);
            const double ZeppZV_resolved = abs(fourvec_signal_jets_ll.eta() - eta_tag_jet_mean);

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
            _h["resolved_pt_tagjet1"]->fill(tag1_jet.pt());
            _h["resolved_pt_tagjet2"]->fill(tag2_jet.pt());
            _h["resolved_delta_pt_tagjet"]->fill(abs(tag1_jet.pt()-tag2_jet.pt()));
            _h["resolved_eta_tagjets"]->fill(tag1_jet.eta()); _h["resolved_eta_tagjets"]->fill(tag2_jet.eta());
            _h["resolved_phi_tagjets"]->fill(tag1_jet.phi()); _h["resolved_phi_tagjets"]->fill(tag2_jet.phi());
            _h["resolved_m_tagjets"]->fill(m_tagjets);
            _h["resolved_dy_tagjets"]->fill(dy_tagjets);
            _h["resolved_dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
            //lepton plots
            _h["resolved_n_lepton_stable"]->fill(nlep);
            _h["resolved_pt_lepton"]->fill(lep1.pT()); _h["resolved_pt_lepton"]->fill(lep2.pT()); 
            _h["resolved_pt_lepton1"]->fill(lep1.pT()); _h["resolved_pt_lepton2"]->fill(lep2.pT()); 
            _h["resolved_eta_lepton"]->fill(lep1.eta()); _h["resolved_eta_lepton"]->fill(lep2.eta());    
            //ana-specific
            _h["resolved_m_ll"]->fill(m_ll);
            _h["resolved_pt_ll"]->fill(fourvec_ll.pT());
            _h["resolved_leptons_pids"]->fill(lep1.pid());_h["resolved_leptons_pids"]->fill(lep2.pid());
            _h["resolved_DR_min_lepton_sigjets1"]->fill(std::min(deltaR(signal_jet1, lep1),deltaR(signal_jet1, lep2)));
            _h["resolved_DR_min_lepton_sigjets2"]->fill(std::min(deltaR(signal_jet2, lep1),deltaR(signal_jet2, lep2)));
            _h["resolved_pt_tagjets"]->fill(fourvec_tag_jj.pT());
            _h["resolved_n_fjets"]->fill(n_fjets);
            _h["resolved_pt_signal_jets"]->fill(signal_jet1.pT());_h["resolved_pt_signal_jets"]->fill(signal_jet2.pT());
            _h["resolved_mass_signal_jets"]->fill(signal_mjj);
            _h["resolved_mass_WZ"]->fill(fourvec_signal_jets_ll.mass());                
            _h["resolved_pt_WZ"]->fill(fourvec_signal_jets_ll.pt());                
            _h["resolved_mass_Full"]->fill(fourvec_signal_jets_full.mass());                
            _h["resolved_pt_Full"]->fill(fourvec_signal_jets_full.pt());  
            _h["resolved_Centrality"]->fill(Centrality_resolved);
            _h["resolved_CentralityVhad"]->fill(CentralityV_resolved);
            _h["resolved_CentralityZ"]->fill(CentralityZ_resolved);
            _h["resolved_CentralityZV"]->fill(CentralityZV_resolved);
            _h["resolved_ZeppZ"]->fill(ZeppZ_resolved);
            _h["resolved_ZeppV"]->fill(ZeppV_resolved);
            _h["resolved_ZeppZV"]->fill(ZeppZV_resolved);
            
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

        double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
        double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
        double pos_w_sum_final = dbl(*_c["pos_w_final"]);
        double neg_w_sum_final = dbl(*_c["neg_w_final"]);
        double pos_w_sum_final_merged = dbl(*_c["pos_w_final_merged"]);
        double neg_w_sum_final_merged = dbl(*_c["neg_w_final_merged"]);
        double pos_w_sum_final_resolved = dbl(*_c["pos_w_final_resolved"]);
        double neg_w_sum_final_resolved = dbl(*_c["neg_w_final_resolved"]);
        MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

/*         MSG_INFO("\n pos weights initial final ratio Merged " << pos_w_sum_initial <<" " << pos_w_sum_final_merged <<" "<< pos_w_sum_final_merged/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio Merged " << neg_w_sum_initial <<" " << neg_w_sum_final_merged <<" "<< neg_w_sum_final_merged/neg_w_sum_initial << "\n" );

        MSG_INFO("\n pos weights initial final ratio Resolved  " << pos_w_sum_initial <<" " << pos_w_sum_final_resolved <<" "<< pos_w_sum_final_resolved/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio Resolved " << neg_w_sum_initial <<" " << neg_w_sum_final_resolved <<" "<< neg_w_sum_final_resolved/neg_w_sum_initial << "\n" );
 */

        // normalize all to 1 since in case of mostly negative weights not clear what it will do
        for (auto & i_name : _hist_names){ 
            std::cout << "normalizeing hist " << i_name <<" to 1; " ;
            normalize(_h[i_name], 1.0);
        }
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
    /// @}


    };


    RIVET_DECLARE_PLUGIN(WmZ_llqq);

}