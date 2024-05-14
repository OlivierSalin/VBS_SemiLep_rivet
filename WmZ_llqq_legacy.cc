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
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


    /// @brief Add a short analysis description here
    class WmZ_llqq : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WmZ_llqq);

    const Particle GetParent(const Particle& p) 
    {
	const Particles parents = p.parents();
	if (parents.size() == 0) return p;
        for (const Particle& p2 : parents) {
	    if (p2.pid() != p.pid()) return p2;
	    else return GetParent(p2);
	}
	return p;
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

        std::string jsonfilestr =  txt_dir + "WmZ_llqq_cuts.json"; 
        std::cout << "++++++assume .json for this WmZ_llqq" << " is " << jsonfilestr << "\n";
        std::ifstream json_file(jsonfilestr);
        
        _jcuts = json::parse(json_file);
        std::cout << "++++++ to check json 1 var got photon pt min " << _jcuts["m_tagjets"] << "\n";
        _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
        _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
        _electron_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_electron"])*GeV; 
        _muon_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton_muon"])*GeV; 

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
        declare(fatjetsfs, "fjets");
        declare(MissingMomentum(), "METFinder");


        // Merged histograms
        
        // plots common with others
        std::ifstream jet_hist_merged_file(txt_dir + "/jet_hists_merged.json");      
        json jet_hist_merged = json::parse(jet_hist_merged_file);
        for (json::iterator it = jet_hist_merged.begin(); it != jet_hist_merged.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        std::ifstream lep_hist_merged_file(txt_dir + "/lepton_hists_merged.json");      
        json lep_hist_merged = json::parse(lep_hist_merged_file);
        for (json::iterator it = lep_hist_merged.begin(); it != lep_hist_merged.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_merged_file(txt_dir + "/WmZ_llqq_hists_merged.json");      
        json ana_hist_merged = json::parse(ana_hist_merged_file);
        for (json::iterator it = ana_hist_merged.begin(); it != ana_hist_merged.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }

        // Resolved histograms
        // plots common with others
        std::ifstream jet_hist_resolved_file(txt_dir + "/jet_hists_resolved.json");      
        json jet_hist_resolved = json::parse(jet_hist_resolved_file);
        for (json::iterator it = jet_hist_resolved.begin(); it != jet_hist_resolved.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        std::ifstream lep_hist_resolved_file(txt_dir + "/lepton_hists_resolved.json");      
        json lep_hist_resolved = json::parse(lep_hist_resolved_file);
        for (json::iterator it = lep_hist_resolved.begin(); it != lep_hist_resolved.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_resolved_file(txt_dir + "/WmZ_llqq_hists_resolved.json");      
        json ana_hist_resolved = json::parse(ana_hist_resolved_file);
        for (json::iterator it = ana_hist_resolved.begin(); it != ana_hist_resolved.end(); ++it) {
            book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
            _hist_names.push_back(it.key());
        }


        //counter for efficiency
        book(_c["pos_w_initial"],"pos_w_initial");
        book(_c["pos_w_final"],"pos_w_final");
        book(_c["neg_w_initial"],"neg_w_initial");
        book(_c["neg_w_final"],"neg_w_final");

        // Cut-flows merged region
        _cutflows_merged.addCutflow("WmZ_llqq_selections", {"have_two_lep","SFOC","pt_lep1_2","m_ll",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "At least one fjets","fjets is W/Z (mJ cut)","Total Merged selec",});
        // Cut-flows resolved region
        _cutflows_resolved.addCutflow("WmZ_llqq_selections", {"have_two_lep","SFOC","pt_lep1_2","m_ll",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "Failed Merged selection","At least two signal jets","signal jets pT","signal_mjj","Total Resolved selec",});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
        // save weights before cuts
        double ev_nominal_weight =  event.weights()[0];
        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
        else {_c["neg_w_initial"]->fill();}

        _cutflows_merged.fillinit();
        _cutflows_resolved.fillinit();
        // Retrieve dressed leptons, sorted by pT
        Particles e_stable;
        Particles mu_stable;
        if (_docut==1){
            e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _electron_pt_cut);
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
        if (nlep!=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];

        const Particle parent_lep1 = GetParent(lep1);
        const Particle parent_lep2 = GetParent(lep2);

        printf("Parent of lep1 PID: %d\n", parent_lep1.pid());
        printf("Parent of lep2 PID: %d\n", parent_lep2.pid());

        // opposite charge and same flavour
        if (lep1.pid()+lep2.pid()!=0) vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();


        // Cuts on the pT of the leptons
        if (_docut==1 && (leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();
        

        const FourMomentum fourvec_ll = lep1.mom() + lep2.mom();
        double m_ll = fourvec_ll.mass()/GeV;

        // Cuts on the invariant mass of the lepton pair   
        if (_docut==1 && (m_ll>_jcuts["m_ll"][1] || m_ll<_jcuts["m_ll"][0])) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_jet"])*GeV);
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
            for (int j = 0; j < n_jets; j++) {
                if (i!=j){
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
        Jets fjets = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_fjet"])*GeV);
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4);

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
        
        // Check if we are in the Merged signal region
        if (n_fjets > 0) {
            _cutflows_merged.fillnext();
            const FourMomentum fourvec_fjets = fjets[0].mom();
            
            if ((_docut==1 && (fourvec_fjets.mass() >= _jcuts["m_fjet_WZ"][0] && fourvec_fjets.mass() <= _jcuts["m_fjet_WZ"][1])) || _docut==0) {
                _cutflows_merged.fillnext();
                _cutflows_merged.fillnext();

                // The four vector of the vertex
                const FourMomentum fourvec_fjets_ll = lep1.mom() + lep2.mom() + fjets[0].mom();
                // We are in the Merged signal region
                // Fill in the histogram of the merged region
                _h["merged_n_jets"]->fill(n_jets);
                _h["merged_pt_tagjet1"]->fill(tag1_jet.pt());
                _h["merged_pt_tagjet2"]->fill(tag2_jet.pt());
                _h["merged_eta_tagjets"]->fill(tag1_jet.eta()); _h["merged_eta_tagjets"]->fill(tag2_jet.eta());
                _h["merged_phi_tagjets"]->fill(tag1_jet.phi()); _h["merged_phi_tagjets"]->fill(tag2_jet.phi());
                _h["merged_m_tagjets"]->fill(m_tagjets);
                _h["merged_dy_tagjets"]->fill(dy_tagjets);
                _h["merged_dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
                //lepton plots
                _h["merged_n_lepton_stable"]->fill(nlep);
                _h["merged_pt_lepton"]->fill(lep1.pT()); _h["merged_pt_lepton"]->fill(lep2.pT()); 
                _h["merged_pt_lepton"]->fill(lep1.pT()); _h["merged_pt_lepton"]->fill(lep2.pT()); 
                _h["merged_eta_lepton"]->fill(lep1.eta()); _h["merged_eta_lepton"]->fill(lep2.eta());    
                _h["merged_eta_lepton"]->fill(lep1.eta()); _h["merged_eta_lepton"]->fill(lep2.eta());    
                //ana-specific
                _h["merged_m_ll"]->fill(m_ll);
                _h["merged_pt_ll"]->fill(fourvec_ll.pT());
                _h["merged_leptons_pids"]->fill(lep1.pid());_h["merged_leptons_pids"]->fill(lep2.pid());
                _h["merged_pt_tagjets"]->fill(fourvec_tag_jj.pT());
                _h["merged_pt_fjet"]->fill(fourvec_fjets.pT());
                _h["merged_n_fjets"]->fill(n_fjets);
                _h["merged_mass_fjets"]->fill(fourvec_fjets.mass());
                _h["merged_mass_WZ"]->fill(fourvec_fjets_ll.mass());                

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
            if (_docut==1 && (signal_jet1.pT() < 40 || signal_jet2.pT() < 20)) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();

            // Make a cut on signal_mjj
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if (_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1])) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            const FourMomentum fourvec_signaljets_ll = lep1.mom() + lep2.mom() + signal_jet1 + signal_jet2;

            _cutflows_resolved.fillnext();
            // Fill in the histograms for the resolved region
            //jet plots
            _h["resolved_n_jets"]->fill(n_jets);
            _h["resolved_pt_tagjet1"]->fill(tag1_jet.pt());
            _h["resolved_pt_tagjet2"]->fill(tag2_jet.pt());
            _h["resolved_eta_tagjets"]->fill(tag1_jet.eta()); _h["resolved_eta_tagjets"]->fill(tag2_jet.eta());
            _h["resolved_phi_tagjets"]->fill(tag1_jet.phi()); _h["resolved_phi_tagjets"]->fill(tag2_jet.phi());
            _h["resolved_m_tagjets"]->fill(m_tagjets);
            _h["resolved_dy_tagjets"]->fill(dy_tagjets);
            _h["resolved_dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
            //lepton plots
            _h["resolved_n_lepton_stable"]->fill(nlep);
            _h["resolved_pt_lepton"]->fill(lep1.pT()); _h["resolved_pt_lepton"]->fill(lep2.pT()); 
            _h["resolved_pt_lepton"]->fill(lep1.pT()); _h["resolved_pt_lepton"]->fill(lep2.pT()); 
            _h["resolved_eta_lepton"]->fill(lep1.eta()); _h["resolved_eta_lepton"]->fill(lep2.eta());    
            _h["resolved_eta_lepton"]->fill(lep1.eta()); _h["resolved_eta_lepton"]->fill(lep2.eta());    
            //ana-specific
            _h["resolved_m_ll"]->fill(m_ll);
            _h["resolved_pt_ll"]->fill(fourvec_ll.pT());
            _h["resolved_leptons_pids"]->fill(lep1.pid());_h["resolved_leptons_pids"]->fill(lep2.pid());
            _h["resolved_pt_tagjets"]->fill(fourvec_tag_jj.pT());
            _h["resolved_n_fjets"]->fill(n_fjets);
            _h["resolved_pt_signal_jets"]->fill(signal_jet1.pT());_h["resolved_pt_signal_jets"]->fill(signal_jet2.pT());
            _h["resolved_mass_signal_jets"]->fill(signal_mjj);
            _h["resolved_mass_WZ"]->fill(fourvec_signaljets_ll.mass());

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
        MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
        MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

        // normalize all to 1 since in case of mostly negative weights not clear what it will do
        for (auto & i_name : _hist_names){ 
            std::cout << "normalizeing hist " << i_name <<" to 1; " ;
            normalize(_h[i_name], 1.0);
        }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Histo2DPtr> _h2;
    // map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    int _docut;
    Cut _electron_eta_cut;
    Cut _muon_eta_cut;
    Cut _electron_pt_cut;
    Cut _muon_pt_cut;
    json _jcuts;
    Cutflows _cutflows_merged;
    Cutflows _cutflows_resolved;
    std::vector<std::string> _hist_names;
    /// @}


    };


    RIVET_DECLARE_PLUGIN(WmZ_llqq);

}