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
#include <algorithm>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


    /// @brief Add a short analysis description here
    class WpZ_llqq : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WpZ_llqq);

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


    const bool DescendantsWithParticle(const Particle& p, const Particle& q) 
    {
        for (const Particle& descendant : p.allDescendants()) {
            //printf("Func Descendant PID: %d and Descendant Status: %d\n", descendant.pid(), descendant.genParticle()->status());
            if (descendant.genParticle() == q.genParticle()) {
                return true;
            }
        }
        return false;
    }

    const bool AncestorWithParticle(const Particle& p, const Particle& q) 
    {
        //printf("Ancestor size: %d\n", p.ancestors(Cuts::OPEN, false).size());
        for (const Particle& ancestor : p.ancestors(Cuts::OPEN, false)) {
            
            //printf("Ancestor PID: %d and Ancestor Status: %d\n", ancestor.pid(), ancestor.genParticle()->status());
            if (ancestor.genParticle() == q.genParticle()) {
                return true;
            }
        }
        return false;
    }

    double MatchingJet(const Jet& jet, const Particle& particle) {
        int matchingParticles = 0;
        int totalParticles = 0;

        // Loop over all particles in the jet
        for(const Particle& jet_particle : jet.particles()){
            totalParticles++;

            // Check if the particle in the jet is a descendant of the given particle
            if(DescendantsWithParticle(particle, jet_particle)){
                matchingParticles++;
            }
        }
        double fraction = static_cast<double>(matchingParticles) / totalParticles;
        // Calculate and return the fraction of particles in the jet that are descendants of the given particle
        return fraction;
    }

    std::vector<Particle> ValidQuark_W(const Particle& W_boson){
        std::vector<Particle> valid_descendants;
        //printf("W boson PID: %d and Status: %d\n", W_boson.pid(), W_boson.genParticle()->status());
        for(const Particle& descendant : W_boson.allDescendants()){
            int descendant_status = descendant.genParticle()->status();
            if(descendant.abspid() < 7 && descendant_status > 30){
                //printf("W descendant quarks PID: %d and Status: %d\n", descendant.pid(), descendant_status);
                valid_descendants.push_back(descendant);
            }
        }
        while(valid_descendants.size() > 2){
            bool erased = false;
            for(int i = 0; i < valid_descendants.size(); i++){
                for(int j = i + 1; j < valid_descendants.size(); j++){
                    if(valid_descendants[i].pid() == -valid_descendants[j].pid()){
                        bool found = false;
                        for(const Particle& ancestor_i : valid_descendants[i].ancestors(Cuts::OPEN, false)){
                            if(ancestor_i.pid()==21 && ancestor_i.genParticle()->status() == 51){
                                if(AncestorWithParticle(valid_descendants[j],ancestor_i)){
                                    found = true;
                                    break;
                                }
                            }
                        }
                        if(found){
                            valid_descendants.erase(std::remove_if(valid_descendants.begin(), valid_descendants.end(), [&](const Particle& p){
                                return p.genParticle() == valid_descendants[i].genParticle() || p.genParticle() == valid_descendants[j].genParticle();
                            }), valid_descendants.end());
                            erased = true;
                            break;
                        }
                    } 
                }
                if(erased) break;
            }
            if(!erased){
                valid_descendants.clear();
                break;
            }  // If no elements were erased in this iteration, break the loop to avoid infinite loop reutrn empty vector
        }

        return valid_descendants;
    }

    std::vector<Particle> Tagging_quarks(const Particles& all_particles){
        std::vector<Particle> tagging_quarks;
        bool found_ = false;
        for(const Particle& p : all_particles){
            if (found_) break;
            ConstGenParticlePtr p_ = p.genParticle();
            int status = p_->status();
            if (abs(status) == 21){
                found_ = true;
                for(const Particle& child : p.children()){
                    if(child.abspid() < 7 && child.genParticle()->status() == 23){
                        if(child.children().size() == 1){
                            for(const Particle& grandchild : child.children()) {
                                tagging_quarks.push_back(grandchild);
                            }} }    
                }}}
        if(tagging_quarks.size() != 2){
            tagging_quarks.clear(); // Clear the vector if it doesn't contain exactly two particles
        }
        return tagging_quarks;
    }

    bool Check_VBS_event(const std::vector<Particle>& all_particles){
        int W_boson_count = 0;
        int Z_boson_count = 0;
        int H_boson_count = 0;

        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            if (status == 21){
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    if(child.abspid() == 24){
                        W_boson_count++;
                    }
                    else if(child.abspid() == 23){
                        Z_boson_count++;
                    }
                    else if(child.abspid() == 25){
                        H_boson_count++;
                    }
                }
                break;
            }
        }

        if(W_boson_count != 1 || Z_boson_count != 1 || H_boson_count != 0 || Tagging_quarks(all_particles).size() != 2){
            return false;
        }

        return true;
    }

    std::vector<double> Truth_q_minDR_jets(const std::vector<Particle>& all_particles, const std::vector<Jet>& jets, bool merged = false){
        std::vector<double> minDeltaRs;
        int W_boson_count = 0;
        Particle W_boson;

        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            if (status == 21){
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    if(child.abspid() == 24){
                        W_boson = child;
                        W_boson_count++;
                    }
                }
                break;
            }
        }

        if(W_boson_count != 1) return std::vector<double>(); // Return an empty vector if there isn't exactly one W boson

        if(merged){
            double minDeltaR = std::numeric_limits<double>::max();
            for(const Jet& jet : jets){
                const double dR = deltaR(W_boson, jet);
                if(dR < minDeltaR) minDeltaR = dR;
            }
            minDeltaRs.push_back(minDeltaR);
        } else {
            std::vector<Particle> valid_quarks = ValidQuark_W(W_boson);
            for(const Particle& descendant : ValidQuark_W(W_boson)){
                double minDeltaR = std::numeric_limits<double>::max();
                for(const Jet& jet : jets){
                    const double dR = deltaR(descendant, jet);
                    if(dR < minDeltaR) minDeltaR = dR;
                }
                minDeltaRs.push_back(minDeltaR);
            }
        }                   

        return minDeltaRs;
    }

    bool IsTruthTagJet(const Particles& taggingQuarks, const Jet& tag_jet, double dR_cut = 0.1){        
        for(const Particle& taggingQuark : taggingQuarks){
            const double dR = deltaR(taggingQuark, tag_jet);
            if(dR < dR_cut) return true;
        }
        
        return false;
    }

    int NbTagJetMisID(const std::vector<Particle>& all_particles, const std::vector<Jet>& tag_jets, double dR_cut = 0.1){
        int tag_jet_misID = 0;
        std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() != 2){
            tag_jet_misID = -1; // Set to -1 to indicate a bad event
        } else {
            for(const Jet& jet : tag_jets){
                bool truth_DRmatch = false;
                for(const Particle& taggingQuark : taggingQuarks){
                    const double dR = deltaR(taggingQuark, jet);
                    if(dR < dR_cut) {
                        truth_DRmatch = true;
                        break;
                    }
                }
                if(!truth_DRmatch) tag_jet_misID++;
            }
        }
        return tag_jet_misID;
    }   
     

    int IsTrueWboson(const std::vector<Particle>& all_particles, const Jet& jet, double dR_cut = 0.1, bool merged = false){
        int W_boson_count = 0;
        Particle W_boson;

        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            if (status == 21){
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    if(child.abspid() == 24){
                        W_boson = child;
                        W_boson_count++;
                    }
                }
                break; // Break the loop to avoid looking at the 2nd parton
            }
        }

        if(W_boson_count != 1) return -1; // Return -1 if there isn't exactly one W boson

        if(merged){
            const double dR = deltaR(W_boson, jet);
            if(dR < dR_cut) return 1; // Return 1 if the jet matches the W boson
        } else {
            std::vector<Particle> valid_quarks = ValidQuark_W(W_boson);
            for(const Particle& descendant : valid_quarks){
                const double dR = deltaR(descendant, jet);
                if(dR < dR_cut) return 1; // Return 1 if the jet matches a valid quark from the W boson
            }
        }                   

        return 0; // Return 0 if the jet does not match the W boson or any valid quark from the W boson
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


        // Before cut histograms

        std::ifstream bef_cut_file(txt_dir + "/WpZ_llqq_hists_bef_cuts.json");      
        json bef_cut = json::parse(bef_cut_file);
        for (json::iterator it = bef_cut.begin(); it != bef_cut.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }

        // Truth information histograms

        std::ifstream truth_hist_file(txt_dir + "/WpZ_llqq_hists_truth.json");      
        json truth_hist = json::parse(truth_hist_file);
        for (json::iterator it = truth_hist.begin(); it != truth_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }


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
        std::ifstream ana_hist_merged_file(txt_dir + "/WpZ_llqq_hists_merged.json");      
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
        std::ifstream ana_hist_resolved_file(txt_dir + "/WpZ_llqq_hists_resolved.json");      
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
        _cutflows_merged.addCutflow("WpZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "At_least_one_fjets","fjets_is_W/Z","Total_Merged_selec",});
        // Cut-flows resolved region
        _cutflows_resolved.addCutflow("WpZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "Failed_Merged_selection","At_least_two_signal_jets","signal_jets_pT","signal_mjj","Total_Resolved_selec",});

    }


    /// Perform the per-event analysis

    void analyze(const Event& event) {
        // save weights before cuts
        double ev_nominal_weight =  event.weights()[0];
        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
        else {_c["neg_w_initial"]->fill();}

        const Particles all_particles = event.allParticles();
        // Fill the truth histograms
        if(Check_VBS_event(all_particles)){
            _h["bef_cut_Is_VBS_event"]->fill(1);
        } else {
            _h["bef_cut_Is_VBS_event"]->fill(0);
        }
        const double cutof_DR = 0.15;

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
        _h["bef_cut_nlep"]->fill(nlep);
        if (nlep !=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];

        const Particle parent_lep1 = GetParent(lep1);
        const Particle parent_lep2 = GetParent(lep2);

/*         printf("PID and Status of Lepton from FS\n");
        printf("Lepton PID: %d Status lepton %d\n",lep1.pid(),lep1.genParticle()->status());
        printf("Parent Lepton PID : %d boson %d\n\n\n",parent_lep1.pid(),parent_lep1.genParticle()->status()); */

/*         printf("Process \n");
        printf("Check VBS_event : %d\n", Check_VBS_event(all_particles));
        //const Particles all_particles = event.allParticles();
        bool found = false;
        for(const Particle& p : all_particles){
            if (found) break;
            ConstGenParticlePtr p_ = p.genParticle();
            int status = p_->status();
            if (abs(status) == 21){
                found = true;
                printf("\nParticle PID: %d, Status: %d, Eta: %f, Pt: %f\n", p.pid(), status, p.eta(), p.mom().pt());
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    printf("Child PID: %d, Status: %d, Eta: %f, Pt: %f\n", child.pid(), child_->status(), child.eta(), child.mom().pt());
                }
            }
        } */


        _h["bef_cut_pt_lepton1"]->fill(lep1.pT());_h["bef_cut_pt_lepton2"]->fill(lep2.pT());
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
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_jet"])*GeV);
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
        Jets fjets = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_fjet"])*GeV);
        _h["bef_cutDR_n_fjets"]->fill(fjets.size());
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4);

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



/*         if(sjets_sig.size()>0){
        // Loop over all particles
            printf("\nTruth information: matching Jet particles with W Boson quark \n");
            for(const Particle& p : all_particles){
                ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
                int status = p_->status(); // Get the status of the particle
                
                // Check if the particle is a W boson with a status of 23, 22, or 21
                if(p.pid() == 24 && (abs(status) == 23 || abs(status) == 22 || abs(status) == 21)){
                    const Particle W_boson= p;
                    std::vector<Particle> valid_quarks = ValidQuark_W(W_boson);
                    for(const Particle& descendant : ValidQuark_W(W_boson)){
                        int descendant_status = descendant.genParticle()->status();
                        //printf("Descendant PID: %d and Descendant Status: %d\n", descendant.pid(), descendant_status);
                        int jetCounter = 0;
                        double minDeltaR = std::numeric_limits<double>::max();
                        for(const Jet& jet : jets){
                            const double dR = deltaR(descendant, jet);
                            if(dR < minDeltaR) minDeltaR = dR;
                            printf("\nQuark pid : %d\n ", descendant.pid());
                            std::cout << "Jet " << jetCounter << ", DeltaR quark with jet: " << dR << std::endl;
                            printf("Matching Jet particles with W Boson: %f\n", MatchingJet(jet, W_boson));
                            std::cout << " DeltaR jet and W boson : " << deltaR(W_boson, jet) << std::endl;
                            jetCounter++;
                        }
                        printf("\nQuark pid : %d\n ", descendant.pid());
                        std::cout << "Min DeltaR quark with jets: " << minDeltaR << std::endl;
                        
                    
                    }                   
                }
            }

        } */
        if(sjets_sig.size()>0){
            std::vector<double> minDeltaRs = Truth_q_minDR_jets(all_particles, sjets_sig);
            //std::cout << "\nTruth min DRs: ";
            for(size_t i = 0; i < minDeltaRs.size(); ++i){
                _h["Min_DR_Wq_jets"]->fill(minDeltaRs[i]);
                //std::cout << minDeltaRs[i] << " ";
            }
            //std::cout << "\n";
        }

        std::vector<double> minDeltaRs_tag = Truth_q_minDR_jets(all_particles, tag_jets);
        //std::cout << "\nTruth Tag min DRs: ";
        for(size_t i = 0; i < minDeltaRs_tag.size(); ++i){
            _h["Min_DR_q_Tagjets"]->fill(minDeltaRs_tag[i]);
            //std::cout << minDeltaRs_tag[i] << " ";
        }
        //std::cout << "\n";


        //printf("\nTagging Quarks \n");
        std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() == 2){
            for(const Particle& taggingQuark : taggingQuarks){
                //printf("\nTagging Quark PID: %d, Status: %d\n", taggingQuark.pid(), taggingQuark.genParticle()->status());
                int jetCount_ = 0;
                double minDeltaR_ = std::numeric_limits<double>::max();
                for(const Jet& jet : tag_jets){
                    const double dR = deltaR(taggingQuark, jet);
                    if(dR < minDeltaR_) minDeltaR_ = dR;
                    //std::cout << "Jet " << jetCount_ << ", DeltaR quark with jet: " << dR << std::endl;
                    //printf("Matching Jet particles with Truth: %f\n", MatchingJet(jet, taggingQuark));
                    //printf("Is Tagging Jet: %d\n", IsTruthTagJet(taggingQuarks, jet));
                    jetCount_++;
                }
                //std::cout << "   min DeltaR Tag quark with jets: " << minDeltaR_ << std::endl;
                _h["Min_DR_VBSq_jets"]->fill(minDeltaR_);
            }
        }
        
        
        int tag_jet_misID = 0;
        double dR_cut = 0.1;
        //std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() != 2){
            tag_jet_misID = -1; // Set to -1 to indicate a bad event
        } else {
            for(const Jet& jet : tag_jets){
                bool truth_DRmatch = false;
                for(const Particle& taggingQuark : taggingQuarks){
                    const double dR = deltaR(taggingQuark, jet);
                    if(dR < dR_cut) {
                        truth_DRmatch = true;
                        break;
                    }
                }
                if(!truth_DRmatch) tag_jet_misID++;
            }
        }
        //printf("\nTagging Jets MisID: %d\n", tag_jet_misID);
        //printf("NbTagJetMisID: %d\n", NbTagJetMisID(all_particles, tag_jets, 0.1));  
        _h["misID_tag_jets"]->fill(tag_jet_misID);
        _h["FromTruth_W_q_tagjets"]->fill(IsTrueWboson(all_particles, tag_jets[0], cutof_DR, true)+ IsTrueWboson(all_particles, tag_jets[1], cutof_DR, true));
        

            
        // Check if we are in the Merged signal region
        if (n_fjets > 0) {
            _cutflows_merged.fillnext();
            std::vector<double> minDeltaRs_merged = Truth_q_minDR_jets(all_particles, fjets, true);
            //std::cout << "\nTruth merged min DRs: ";
            for(size_t i = 0; i < minDeltaRs_merged.size(); ++i){
                _h["merged_min_DR_Wq_fjets"]->fill(minDeltaRs_merged[i]);
                //std::cout << minDeltaRs_merged[i] << ",  ";
            }
            //std::cout << "\n";

            _h["merged_fjet_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, fjets[0], cutof_DR, true));


            const FourMomentum fourvec_fjets = fjets[0].mom();
            if ((_docut==1 && (fourvec_fjets.mass() >= _jcuts["m_fjet_WZ"][0] && fourvec_fjets.mass() <= _jcuts["m_fjet_WZ"][1]))||_docut==0) {
                _cutflows_merged.fillnext();
                
                // Total cutflow of the merged region
                _cutflows_merged.fillnext();

                

                _h["merged_misID_tag_jets"]->fill(NbTagJetMisID(all_particles, tag_jets, cutof_DR));
            

                // Four vector of the QGC system in resolved region llJ
                const FourMomentum fourvec_fjets_ll = lep1.mom() + lep2.mom() + fjets[0].mom();

                // Four vector of the Full system in resolved region 
                const FourMomentum fourvec_fjets_full = lep1.mom() + lep2.mom() + fjets[0].mom() + tag1_jet + tag2_jet;

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
                _h["merged_pt_WZ"]->fill(fourvec_fjets_ll.pt());                
                _h["merged_mass_Full"]->fill(fourvec_fjets_full.mass());                
                _h["merged_pt_Full"]->fill(fourvec_fjets_full.pt());
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

            std::vector<double> minDeltaRs_resolved = Truth_q_minDR_jets(all_particles, sjets_sig);
            //std::cout << "\nTruth resolved min DRs: ";
            for(size_t i = 0; i < minDeltaRs_resolved.size(); ++i){
                _h["resolved_min_DR_Wq_jets"]->fill(minDeltaRs_resolved[i]);
                //std::cout << minDeltaRs_resolved[i] << " ";
            }
            //std::cout << "\n";


            _h["resolved_sjets_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, sjets_sig[0], cutof_DR, false)+IsTrueWboson(all_particles, sjets_sig[1], cutof_DR, false));
            _h["resolved_sjet1_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, sjets_sig[0], cutof_DR, false));
            _h["resolved_sjet2_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, sjets_sig[1], cutof_DR, false));



            
            // Make a cut on signal_mjj
            const FourMomentum fourvec_signal_jets = signal_jet1 + signal_jet2;
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if ((_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1]) )) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();

            

            _h["resolved_misID_tag_jets"]->fill(NbTagJetMisID(all_particles, tag_jets, cutof_DR));
            

            const FourMomentum fourvec_signal_jets_ll = lep1.mom() + lep2.mom() + fourvec_signal_jets;

            // Four vector of the Full system in resolved region 
            const FourMomentum fourvec_signal_jets_full = lep1.mom() + lep2.mom() + fourvec_signal_jets + tag1_jet + tag2_jet;
            // More than 2 signal jets candidates
            FourMomentum fourvec_signal_jjj;

            // Zeppenfeld variables
            const double ZeppZ_resolved = abs(fourvec_ll.eta() - eta_tag_jet_mean);
            const double ZeppV_resolved = abs(fourvec_signal_jets.eta() - eta_tag_jet_mean);
            const double ZeppZV_resolved = abs(fourvec_signal_jets_ll.eta() - eta_tag_jet_mean);

            double ZeppRes = 0.0;
            if (sjets_sig.size() > 2){
                const FourMomentum signal_jet3 = sjets_sig[2].mom();

                ZeppRes = abs(signal_jet3.eta() - eta_tag_jet_mean);

                fourvec_signal_jjj = signal_jet1 + signal_jet2 + signal_jet3;
                double signal_mjjj = fourvec_signal_jjj.mass();                
                // if ((_docut==1 && (signal_mjjj < _jcuts["m_jjj"] )) ) {vetoEvent;}
                // _cutflows_resolved.fillnext();            
            }



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
            _h["resolved_mass_WZ"]->fill(fourvec_signal_jets_ll.mass());                
            _h["resolved_pt_WZ"]->fill(fourvec_signal_jets_ll.pt());                
            _h["resolved_mass_Full"]->fill(fourvec_signal_jets_full.mass());                
            _h["resolved_pt_Full"]->fill(fourvec_signal_jets_full.pt());  
            _h["resolved_ZeppZ"]->fill(ZeppZ_resolved);
            _h["resolved_ZeppV"]->fill(ZeppV_resolved);
            _h["resolved_ZeppZV"]->fill(ZeppZV_resolved);
            
            _h["resolved_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[0])); _h["resolved_Ntrk_tagjets"]->fill(CountChargedTracks(tag_jets[1]));
            _h["resolved_Ntrk_signal_jets"]->fill(CountChargedTracks(sjets_sig[0])); _h["resolved_Ntrk_signal_jets"]->fill(CountChargedTracks(sjets_sig[1]));

        
            // More than 2 signal jets candidates
            if (sjets_sig.size() > 2){
                _h["resolved_DeltaPt_sig_jet2_jet3"]-> fill(abs(signal_jet2.pT() - sjets_sig[2].mom().pT()));   
                _h["resolved_mjjj"]->fill(fourvec_signal_jjj.mass());
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


    RIVET_DECLARE_PLUGIN(WpZ_llqq);

}