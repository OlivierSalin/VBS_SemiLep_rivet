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
    class ZZ_llqq : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZZ_llqq);

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

        if (!W_boson.genParticle() || (abs(W_boson.pid()) != 24) || W_boson.pT() == 0) {
            // If W_boson is not valid, return an empty vector immediately
            return valid_descendants;
        }  

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

        if(W_boson_count != 0 || Z_boson_count != 2 || H_boson_count != 0 || Tagging_quarks(all_particles).size() != 2){
            return false;
        }

        return true;
    }

    const Particle GetWboson(const std::vector<Particle>& all_particles){
        Particle W_boson;
        int W_boson_count = 0;

        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            if (status == 21){
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    if(child.abspid() == 23){
                        W_boson = child;
                        W_boson_count++;
                    }
                }
                if(W_boson_count == 1) break; // Break the loop if exactly one W boson has been found
            }
        }

        return W_boson;
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

    int IsTruthTagJet(const std::vector<Particle>& all_particles, const Jet& tag_jet, double dR_cut = 0.4){
        std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() != 2) return -1; // Return -1 if the size of taggingQuarks is not 2

        for(const Particle& taggingQuark : taggingQuarks){
            const double dR = deltaR(taggingQuark, tag_jet);
            if(dR < dR_cut) return 1; // Return 1 if the jet matches a tagging quark
        }
        
        return 0; // Return 0 if the jet does not match any tagging quark
    }

    int NbTagJetMisID(const std::vector<Particle>& all_particles, const std::vector<Jet>& tag_jets, double dR_cut = 0.4){
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

    int NbTagJetMatched(const std::vector<Particle>& all_particles, const std::vector<Jet>& tag_jets, double dR_cut = 0.4){
        int tag_jet_matched = 0; // Initialize to 0 to count matched jets
        std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() != 2){
            return -1; // Return -1 to indicate a bad event
        } else {
            for(const Jet& jet : tag_jets){
                for(const Particle& taggingQuark : taggingQuarks){
                    const double dR = deltaR(taggingQuark, jet);
                    if(dR < dR_cut) {
                        tag_jet_matched++; // Increment if a jet is matched
                        break; // Break to avoid double counting the same jet
                    }
                }
            }
        }
        return tag_jet_matched;
    }

    int IsTrueWboson(const std::vector<Particle>& all_particles, const Jet& jet, double dR_cut = 0.4, bool merged = false){
        int W_boson_count = 0;
        Particle W_boson;

        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            if (status == 21){
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    if(child.abspid() == 23){
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
        std::cout << "++++++assume .json for this ZZ_llqq" << " is " << jsonfilestr << "\n";
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
        _tt_merged = make_unique<TTree>("Merged", "Rivet_physics");
        _tt_merged->Branch("EventNumber", &merged_EventNumber);
        _tt_merged->Branch("EventWeight", &merged_EventWeight);
        _tt_merged->Branch("Label", &_label);
        _tt_merged->Branch("Label_binary", &_label_binary);
        for (auto& var_ : varMap) {
            _tt_merged->Branch(var_.first.c_str(), var_.second);
        }
        for (auto& var_ : varMapInt) {
            _tt_merged->Branch(var_.first.c_str(), var_.second);
        }


        _tt_bef_cut = make_unique<TTree>("Bef_cut", "Rivet_physics");
        _tt_bef_cut->Branch("EventNumber", &EventNumber);
        _tt_bef_cut->Branch("VBS_event", &VBS_event);
        _tt_bef_cut->Branch("Label", &_label);

        _tt_aft_lep_cut = make_unique<TTree>("Lep_cut", "Rivet_physics");
        _tt_aft_lep_cut->Branch("EventNumber", &EventNumber);
        _tt_aft_lep_cut->Branch("VBS_event", &VBS_event);
        _tt_aft_lep_cut->Branch("Label", &_label);
        _tt_aft_lep_cut->Branch("n_jets", &n_jets);
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_aft_lep_cut->Branch(var_.first.c_str(), var_.second);
        }

        _tt_bef_tag_jet = make_unique<TTree>("Tagjet_bef_cut", "Rivet_physics");
        _tt_bef_tag_jet->Branch("EventNumber", &EventNumber);
        _tt_bef_tag_jet->Branch("VBS_event", &VBS_event);
        _tt_bef_tag_jet->Branch("Label", &_label);
        _tt_bef_tag_jet->Branch("n_jets", &n_jets);
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_bef_tag_jet->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_bef_tag_jet->Branch(var_.first.c_str(), var_.second);
        }

        // TRUTH TREEE

        _tt_truth = make_unique<TTree>("Truth", "Rivet_physics");
        _tt_truth->Branch("EventNumber", &EventNumber);
        _tt_truth->Branch("VBS_event", &VBS_event);
        _tt_truth->Branch("Label", &_label);
        _tt_truth->Branch("n_jets", &n_jets);
        for (auto& var_ : varMapDouble_Tagjets_truth) {
            _tt_truth->Branch(var_.first.c_str(), var_.second);
        }   
        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_truth->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_truth->Branch(var_.first.c_str(), var_.second);
        }

        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth->Branch(var_.first.c_str(), var_.second);
        }


        _tt_truth_mismatch = make_unique<TTree>("Truth_mismatch", "Rivet_physics");
        _tt_truth_mismatch->Branch("EventNumber", &EventNumber);
        _tt_truth_mismatch->Branch("VBS_event", &VBS_event);
        _tt_truth_mismatch->Branch("Label", &_label);
        _tt_truth_mismatch->Branch("n_jets", &n_jets);
        for (auto& var_ : varMapDouble_Tagjets_truth) {
            _tt_truth_mismatch->Branch(var_.first.c_str(), var_.second);
        }   

        _tt_truth_mismatch->Branch("Mismatch_Tagjet_pT", &Mismatch_Tagjet_pT);    
        _tt_truth_mismatch->Branch("Mismatch_Tagjet_eta", &Mismatch_Tagjet_eta);   
        _tt_truth_mismatch->Branch("Mismatch_Tagjet_mass", &Mismatch_Tagjet_mass);   
        _tt_truth_mismatch->Branch("Delta_eta_Mismatch_Tagjet_VBS_q", &Delta_eta_Mismatch_Tagjet_VBS_q);   
        _tt_truth_mismatch->Branch("Delta_R_Mismatch_Tagjet_VBS_q", &Delta_R_Mismatch_Tagjet_VBS_q);   
        _tt_truth_mismatch->Branch("Delta_Phi_Mismatch_Tagjet_VBS_q", &Delta_Phi_Mismatch_Tagjet_VBS_q);   
        _tt_truth_mismatch->Branch("Delta_eta_Mismatch_Tagjet_Wboson", &Delta_eta_Mismatch_Tagjet_Wboson);   
        _tt_truth_mismatch->Branch("tag_jet_mismatched_fromVhad", &tag_jet_mismatched_fromVhad);   

        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_truth_mismatch->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_truth_mismatch->Branch(var_.first.c_str(), var_.second);
        }

        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_mismatch->Branch(var_.first.c_str(), var_.second);
        }


        _tt_truth_mismatch_noVhad = make_unique<TTree>("Truth_mismatch_noVhad", "Rivet_physics");
        _tt_truth_mismatch_noVhad->Branch("EventNumber", &EventNumber);
        _tt_truth_mismatch_noVhad->Branch("VBS_event", &VBS_event);
        _tt_truth_mismatch_noVhad->Branch("Label", &_label);
        _tt_truth_mismatch_noVhad->Branch("n_jets", &n_jets);
        for (auto& var_ : varMapDouble_Tagjets_truth) {
            _tt_truth_mismatch_noVhad->Branch(var_.first.c_str(), var_.second);
        }   

        _tt_truth_mismatch_noVhad->Branch("Mismatch_Tagjet_pT", &Mismatch_Tagjet_pT);    
        _tt_truth_mismatch_noVhad->Branch("Mismatch_Tagjet_eta", &Mismatch_Tagjet_eta);   
        _tt_truth_mismatch_noVhad->Branch("Mismatch_Tagjet_mass", &Mismatch_Tagjet_mass);   
        _tt_truth_mismatch_noVhad->Branch("Delta_eta_Mismatch_Tagjet_VBS_q", &Delta_eta_Mismatch_Tagjet_VBS_q);   
        _tt_truth_mismatch_noVhad->Branch("Delta_Phi_Mismatch_Tagjet_VBS_q", &Delta_Phi_Mismatch_Tagjet_VBS_q);   
        _tt_truth_mismatch_noVhad->Branch("Delta_eta_Mismatch_Tagjet_Wboson", &Delta_eta_Mismatch_Tagjet_Wboson);   
        _tt_truth_mismatch_noVhad->Branch("tag_jet_mismatched_fromVhad", &tag_jet_mismatched_fromVhad);   

        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_truth_mismatch_noVhad->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_truth_mismatch_noVhad->Branch(var_.first.c_str(), var_.second);
        }

        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_mismatch_noVhad->Branch(var_.first.c_str(), var_.second);
        }

        _tt_truth_matched = make_unique<TTree>("Truth_matched", "Rivet_physics");
        _tt_truth_matched->Branch("EventNumber", &EventNumber);
        _tt_truth_matched->Branch("VBS_event", &VBS_event);
        _tt_truth_matched->Branch("Label", &_label);
        _tt_truth_matched->Branch("n_jets", &n_jets);
        for (auto& var_ : varMapDouble_Tagjets_truth) {
            _tt_truth_matched->Branch(var_.first.c_str(), var_.second);
        }   

        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_truth_matched->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_truth_matched->Branch(var_.first.c_str(), var_.second);
        }

        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_matched->Branch(var_.first.c_str(), var_.second);
        }


        _tt_truth_merged = make_unique<TTree>("Truth_Merged", "Rivet_physics");
        _tt_truth_merged->Branch("EventNumber", &EventNumber);
        _tt_truth_merged->Branch("VBS_event", &VBS_event);
        _tt_truth_merged->Branch("Label", &_label);
        _tt_truth_merged->Branch("n_jets", &n_jets);
        for (auto& var_ : varMapDouble_Tagjets_truth) {
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
        }   

        for(auto& var_ : varMapDouble_VBS_q_truth){
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
        }
        for(auto& var_ : varMapDouble_Vhadboson_truth){
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
        }

        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
        }



        _tt_angle = make_unique<TTree>("Angle", "Rivet_physics");
        _tt_angle->Branch("EventWeight", &Angle_EventWeight);
        _tt_angle->Branch("Label", &_label);
        _tt_angle->Branch("cos_theta_star", &cos_theta_star);

        _tt_fjet = make_unique<TTree>("Fjet_def", "Rivet_physics");
        _tt_fjet->Branch("EventWeight", &Angle_EventWeight);
        _tt_fjet->Branch("Label", &_label);
        _tt_fjet->Branch("Leading_fjet_pt", &Leading_fjet_pt);
        _tt_fjet->Branch("Leading_fjet_eta", &Leading_fjet_eta);
        _tt_fjet->Branch("Leading_fjet_mass", &Leading_fjet_mass);
        _tt_fjet->Branch("Leading_fjet_TrueWboson", &Leading_fjet_TrueWboson);

        _tt_fjet->Branch("SubLeading_fjet_pt", &SubLeading_fjet_pt);
        _tt_fjet->Branch("SubLeading_fjet_eta", &SubLeading_fjet_eta);
        _tt_fjet->Branch("SubLeading_fjet_mass", &SubLeading_fjet_mass);
        _tt_fjet->Branch("SubLeading_fjet_TrueWboson", &SubLeading_fjet_TrueWboson);
        _tt_fjet->Branch("Delta_R_fjets", &Delta_R_fjets);
        _tt_fjet->Branch("Delta_M_fjets", &Delta_M_fjets);


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
        _cutflows_merged.addCutflow("ZZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","m_tagjets",
                            "At_least_one_fjets","fjets_is_W/Z","Total_Merged_selec",});
        // Cut-flows resolved region
        _cutflows_resolved.addCutflow("ZZ_llqq_selections", {"have_two_lep","pt_lep1_2",
                            "n_jets","found_tag_jets","m_tagjets",
                            "Failed_Merged_selection","At_least_two_signal_jets","signal_jets_pT","signal_mjj","M_jjj","Total_Resolved_selec",});

    }

    int Nb_Events_pre_VBS_jet = 0;
    int Nb_Events_post_VBS_jet = 0;
    int Nb_Events_VBS_jet_pure = 0;
    int Nb_pure_tag_quarks = 0;
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
        double Vhad_mass= 80.4;

        if(Check_VBS_event(all_particles)){
            VBS_event = 1;
        } else {
            VBS_event = 0;
        }

        _tt_bef_cut->Fill();


        const double cutof_DR = 0.4;


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

        const Particle W_bson_ = GetWboson(all_particles);
        W_pT = W_bson_.pT();
        FourMomentum W_fourvec_truth = W_bson_.mom();
        W_eta = W_bson_.eta();
        W_phi = W_bson_.phi();
        if(W_pT > 0) {
            std::vector<Particle> W_quarks_= ValidQuark_W(W_bson_);
            if (W_quarks_.size() >= 2) {
                W_Quark1_pT = W_quarks_[0].pT();
                W_Quark1_eta = W_quarks_[0].eta();
                W_Quark1_mass = W_quarks_[0].mom().mass();

                W_Quark2_pT = W_quarks_[1].pT();
                W_Quark2_eta = W_quarks_[1].eta();
                W_Quark2_mass = W_quarks_[1].mom().mass();

                W_Quarks_delta_Eta = abs(W_quarks_[0].eta() - W_quarks_[1].eta());
            }
            _tt_aft_lep_cut->Fill();
        }
        

        //const FourMomentum fourvec_Vlep = lep1.mom() + lep2.mom();
        FourMomentum fourvec_Vlep = lep1.mom() + lep2.mom();
        double m_Vlep = fourvec_Vlep.mass()/GeV;




        // // Retrieve clustered small R jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(_jet_pt20_eta_cut_1 || _jet_pt30_eta_cut_2);
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);


        n_jets = jets.size();
        if (n_jets < _jcuts["n_jets"])  vetoEvent;  
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext(); 


        FourMomentum fourvec_truth_VZ = fourvec_Vlep + W_fourvec_truth;
        std::vector<Particle> taggingQuarks_ = Tagging_quarks(all_particles);
        if(taggingQuarks_.size() == 2 && W_pT > 0){

            VBS_Quarks_truth_Delta_Eta = abs(taggingQuarks_[0].eta() - taggingQuarks_[1].eta());
            VBS_Quarks_truth_Delta_Phi = abs(taggingQuarks_[0].phi() - taggingQuarks_[1].phi());
            VBS_Quarks_truth_mass= (taggingQuarks_[0].mom() + taggingQuarks_[1].mom()).mass();

            VBS_Quark1_truth_pT = taggingQuarks_[0].pT();
            VBS_Quark1_truth_eta = taggingQuarks_[0].eta();
            VBS_Quark1_truth_mass = taggingQuarks_[0].mom().mass();

            VBS_Quark2_truth_pT = taggingQuarks_[1].pT();
            VBS_Quark2_truth_eta = taggingQuarks_[1].eta();
            VBS_Quark2_truth_mass = taggingQuarks_[1].mom().mass();

            const FourMomentum VBS_q_mom = taggingQuarks_[0].mom() + taggingQuarks_[1].mom();

            VBS_Quarks_truth_pT= (VBS_q_mom).pT();
            VBS_Quarks_truth_eta= (VBS_q_mom).eta();
            VBS_Quarks_truth_mass= (VBS_q_mom).mass();

            Truth_Delta_Eta_VBS_q_VZ = abs((VBS_q_mom).eta() - fourvec_truth_VZ.eta());
            Truth_Delta_Phi_VBS_q_VZ = abs((VBS_q_mom).phi() - fourvec_truth_VZ.phi());

            Delta_Phi_VBS_q_Vlep = abs((VBS_q_mom).phi() - fourvec_Vlep.phi());
            _tt_bef_tag_jet->Fill();
        }



        // SELECTION TAG JETS FIRST

        Jets tag_jets;
        Nb_Events_pre_VBS_jet++;
        bool foundVBSJetPair = false; 
        double max_m_tag_jj = 0;
        double max_eta_jj = 0;
        double max_eta_mTag = 0;
        int tag1_jet_index = -1 ,tag2_jet_index = -1;


/*         for (int i = 0; i < n_jets; i++) {
            const Jet& i_jet = jets[i];
            for (int j = i + 1; j < n_jets; j++) {
                const Jet& j_jet = jets[j];
                if(i_jet.pT() < _jcuts["pt_tagjet1"] || j_jet.pT() < _jcuts["pt_tagjet2"]) continue;
                const double m_tag_jj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
                const double eta_jj = abs(i_jet.eta()-j_jet.eta());
                const double eta_prod = i_jet.eta()*j_jet.eta();
                if  (eta_prod < 0.0 && m_tag_jj>max_m_tag_jj){
                    max_m_tag_jj = m_tag_jj;
                    foundVBSJetPair = true;
                    tag1_jet_index = i;
                    tag2_jet_index = j;
                }
            }
        } */


        // MAX ETA
        for (int i = 0; i < n_jets; i++) {
            const Jet& i_jet = jets[i];
            for (int j = i + 1; j < n_jets; j++) {
                const Jet& j_jet = jets[j];
                if(i_jet.pT() < _jcuts["pt_tagjet1"] || j_jet.pT() < _jcuts["pt_tagjet2"]) continue;
                const double m_tag_jj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
                const double eta_jj = abs(i_jet.eta()-j_jet.eta());
                const double eta_prod = i_jet.eta()*j_jet.eta();
                if  (eta_prod < 0.0 && eta_jj>max_eta_jj && m_tag_jj > 300 ){
                    max_eta_jj = eta_jj;
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
        const FourMomentum fourvec_tag_jj = tag1_jet+tag2_jet;

        const double m_tagjets = (fourvec_tag_jj).mass()/GeV;
        if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
        Nb_Events_post_VBS_jet++;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());


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
        idiscardIfAnyDeltaRLess(fjets, leptons, 1.0);
        idiscardIfAnyDeltaRLess(fjets, tag_jets, 1.4);
        int n_fjets = fjets.size(); 

        // Define truth variables to study the VBS signal region
        //std::vector<Particle> taggingQuarks_ = Tagging_quarks(all_particles);
        const FourMomentum tag_jets_mom = tag_jets[0].mom() + tag_jets[1].mom();
        
        if(Check_VBS_event(all_particles)){
            if(taggingQuarks_.size() == 2 && W_pT > 0 ){
                Tagging_Jets_VZ_truth_Delta_Phi = abs(tag_jets_mom.phi() - fourvec_truth_VZ.phi());
                VBS_Quarks_truth_Delta_Eta = abs(taggingQuarks_[0].eta() - taggingQuarks_[1].eta());
                VBS_Quarks_truth_Delta_Phi = abs(taggingQuarks_[0].phi() - taggingQuarks_[1].phi());
                VBS_Quarks_truth_mass= (taggingQuarks_[0].mom() + taggingQuarks_[1].mom()).mass();
                std::vector<Particle> W_quarks_= ValidQuark_W(W_bson_);
                if(VBS_Quarks_truth_Delta_Eta > 1.0 && VBS_Quarks_truth_mass > 250 && W_quarks_.size() >= 2){
                    Nb_pure_tag_quarks++;

                    Tagging_Jet1_pT = tag_jets[0].pT();
                    Tagging_Jet1_eta = tag_jets[0].eta();
                    Tagging_Jet1_mass = tag_jets[0].mom().mass();

                    Tagging_Jet2_pT = tag_jets[1].pT();
                    Tagging_Jet2_eta = tag_jets[1].eta();
                    Tagging_Jet2_mass = tag_jets[1].mom().mass();

                    
                    Tagging_Jets_mass = (tag1_jet + tag2_jet).mass();
                    Tagging_Jets_delta_Eta = abs(tag1_jet.eta() - tag2_jet.eta());

                    FourMomentum tag_quarks=taggingQuarks_[0].mom() + taggingQuarks_[1].mom();

                    absDelta_Eta_Tagging_jets_VBS_q_truth = abs(tag_quarks.eta() - tag_jets_mom.eta());
                    absDelta_pT_Tagging_jets_VBS_q_truth = abs(tag_quarks.pT() - tag_jets_mom.pT());
                    absDelta_Mass_Tagging_jets_VBS_q_truth = abs(tag_quarks.mass() - tag_jets_mom.mass());

                    Delta_Eta_Tagging_jets_VBS_q_truth = tag_jets_mom.eta()- tag_quarks.eta();
                    Delta_Mass_Tagging_jets_VBS_q_truth = tag_jets_mom.mass() - tag_quarks.mass();
                    Delta_pT_Tagging_jets_VBS_q_truth =  tag_jets_mom.pT()- tag_quarks.pT();



                    
                    //W_bson_= GetWboson(all_particles);
                    
                    Delta_Eta_TaggingQuark1_Wquarks.clear();
                    Delta_Eta_TaggingQuark2_Wquarks.clear();
                    for (auto& W_quark : W_quarks_) {
                        Delta_Eta_TaggingQuark1_Wquarks.push_back(abs(taggingQuarks_[0].eta() - W_quark.eta()));
                        Delta_Eta_TaggingQuark2_Wquarks.push_back(abs(taggingQuarks_[1].eta() - W_quark.eta()));
                    }

                    Delta_eta_Tagging_Jet1_VBS_q = std::min(abs(tag_jets[0].eta() - taggingQuarks_[0].eta()), abs(tag_jets[0].eta() - taggingQuarks_[1].eta()));
                    Delta_eta_Tagging_Jet2_VBS_q = std::min(abs(tag_jets[1].eta() - taggingQuarks_[0].eta()), abs(tag_jets[1].eta() - taggingQuarks_[1].eta()));

                    Delta_Phi_Tagging_Jet1_VBS_q = std::min(abs(tag_jets[0].phi() - taggingQuarks_[0].phi()), abs(tag_jets[0].phi() - taggingQuarks_[1].phi()));
                    Delta_Phi_Tagging_Jet2_VBS_q = std::min(abs(tag_jets[1].phi() - taggingQuarks_[0].phi()), abs(tag_jets[1].phi() - taggingQuarks_[1].phi()));

                    


                    tag_jets_matched = NbTagJetMatched(all_particles, tag_jets, cutof_DR);
                    tag_jet1_matched = IsTruthTagJet(all_particles, tag_jets[0], cutof_DR);
                    tag_jet2_matched = IsTruthTagJet(all_particles, tag_jets[1], cutof_DR);


                    if(tag_jets_matched ==1) {
                        FourMomentum mismatch_tag_jets;
                        
                        if(tag_jet1_matched == 0) {
                            mismatch_tag_jets=tag_jets[0].mom();
                        }
                        else {
                            mismatch_tag_jets=tag_jets[1].mom();
                        }
                        Mismatch_Tagjet_pT = mismatch_tag_jets.pT();
                        Mismatch_Tagjet_eta = mismatch_tag_jets.eta();
                        Mismatch_Tagjet_mass = mismatch_tag_jets.mass();
                        Delta_eta_Mismatch_Tagjet_VBS_q = std::min(abs(mismatch_tag_jets.eta() - taggingQuarks_[0].eta()), abs(mismatch_tag_jets.eta() - taggingQuarks_[1].eta()));                            
                        Delta_R_Mismatch_Tagjet_VBS_q = std::min(deltaR(mismatch_tag_jets, taggingQuarks_[0]), deltaR(mismatch_tag_jets, taggingQuarks_[1]));
                        Delta_Phi_Mismatch_Tagjet_VBS_q = std::min(abs(mismatch_tag_jets.phi() - taggingQuarks_[0].phi()), abs(mismatch_tag_jets.phi() - taggingQuarks_[1].phi()));
                        Delta_eta_Mismatch_Tagjet_Wboson= std::min(abs(mismatch_tag_jets.eta() - W_quarks_[0].eta()), abs(mismatch_tag_jets.eta() - W_quarks_[1].eta()));

                        


                        tag_jet_mismatched_fromVhad = 0;
                        if (tag_jet1_matched == 0) {
                            tag_jet_mismatched_fromVhad=IsTrueWboson(all_particles, tag_jets[0], cutof_DR);
                        }
                        if (tag_jet2_matched == 0) {
                            tag_jet_mismatched_fromVhad=IsTrueWboson(all_particles, tag_jets[1], cutof_DR);
                        }
                        if(tag_jet_mismatched_fromVhad == 0) {
                            _tt_truth_mismatch_noVhad->Fill();
                        }
                        _tt_truth_mismatch->Fill();
                    }
                    if(tag_jets_matched == 2) {
                        _tt_truth_matched->Fill();
                        Nb_Events_VBS_jet_pure++;
                    }

                    _tt_truth->Fill();

                }

            }
        }
            

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
            Leading_fjet_eta = fjets[0].eta(); Leading_fjet_mass = fjets[0].mass(); Leading_fjet_pt = fjets[0].pt();
            Leading_fjet_TrueWboson = IsTrueWboson(all_particles, fjets[0], 0.4, true);
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
                    SubLeading_fjet_eta = fjets[1].eta(); SubLeading_fjet_mass = fjets[1].mass(); SubLeading_fjet_pt = fjets[1].pt();
                    SubLeading_fjet_TrueWboson = IsTrueWboson(all_particles, fjets[1], 0.4, true);
                    Delta_R_fjets = deltaR(fjets[0], fjets[1]);
                    Delta_M_fjets = abs(fjets[0].mass()- fjets[1].mass());

                    _tt_fjet->Fill();
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

                merged_VlepVhad_Deta_tagjets = abs(fourvec_fjets_Vlep.eta() - fourvec_tag_jj.eta());
                merged_VlepVhad_Dphi_tagjets = deltaPhi(fourvec_fjets_Vlep, fourvec_tag_jj);

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

                _tt_merged->Fill();
                _tt_truth_merged->Fill();

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

            double VBS_jet_efficiency = static_cast<double>(Nb_Events_post_VBS_jet) / Nb_Events_pre_VBS_jet;
            std::string VBS_jet_efficiency_str = std::to_string(VBS_jet_efficiency);
            std::string VBS_jet_efficiency_file = getOption("OUTDIR") + "/VBS_jet_efficiency.txt";
            std::ofstream ofs_VBS_jet_efficiency (VBS_jet_efficiency_file, std::ofstream::out); 
            ofs_VBS_jet_efficiency << VBS_jet_efficiency_str;
            ofs_VBS_jet_efficiency.close();

            if (Nb_Events_post_VBS_jet > 0){
                double VBS_jet_purity = static_cast<double>(Nb_Events_VBS_jet_pure) / Nb_Events_post_VBS_jet;
                std::string VBS_jet_purity_str = std::to_string(VBS_jet_purity);
                std::string VBS_jet_purity_file = getOption("OUTDIR") + "/VBS_jet_purity.txt";
                std::ofstream ofs_VBS_jet_purity (VBS_jet_purity_file, std::ofstream::out); 
                ofs_VBS_jet_purity << VBS_jet_purity_str;
                ofs_VBS_jet_purity.close();
            }

            if (Nb_pure_tag_quarks > 0){
                double VBS_jet_purity_bis = static_cast<double>(Nb_Events_VBS_jet_pure) / Nb_pure_tag_quarks;
                std::string VBS_jet_purity_bis_str = std::to_string(VBS_jet_purity_bis);
                std::string VBS_jet_purity_bis_file = getOption("OUTDIR") + "/VBS_jet_purity_bis.txt";
                std::ofstream ofs_VBS_jet_purity_bis (VBS_jet_purity_bis_file, std::ofstream::out); 
                ofs_VBS_jet_purity_bis << VBS_jet_purity_bis_str;
                ofs_VBS_jet_purity_bis.close();
            }

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
    unique_ptr<TTree> _tt_merged;
    unique_ptr<TTree> _tt_bef_cut;
    unique_ptr<TTree> _tt_aft_lep_cut;
    unique_ptr<TTree> _tt_bef_tag_jet;
    unique_ptr<TTree> _tt_truth;
    unique_ptr<TTree> _tt_truth_merged;
    unique_ptr<TTree> _tt_truth_mismatch;
    unique_ptr<TTree> _tt_truth_mismatch_noVhad;
    unique_ptr<TTree> _tt_truth_matched;
    unique_ptr<TTree> _tt_angle;
    unique_ptr<TTree> _tt_fjet;

    double cos_theta_star, Angle_EventWeight;

    int VBS_event,n_jets;

    double Leading_fjet_pt,Leading_fjet_eta,Leading_fjet_mass,Leading_fjet_TrueWboson;
    double SubLeading_fjet_pt,SubLeading_fjet_eta,SubLeading_fjet_mass,SubLeading_fjet_TrueWboson,Delta_R_fjets,Delta_M_fjets;


    double Pass_VBS_jet;

    double VBS_Quark1_truth_pT, VBS_Quark1_truth_eta, VBS_Quark1_truth_mass;
    double VBS_Quark2_truth_pT, VBS_Quark2_truth_eta, VBS_Quark2_truth_mass;
    double VBS_Quarks_truth_pT, VBS_Quarks_truth_eta, VBS_Quarks_truth_mass, VBS_Quarks_truth_Delta_Eta, VBS_Quarks_truth_Delta_Phi;
    double Truth_Delta_Eta_VBS_q_VZ, Truth_Delta_Phi_VBS_q_VZ,Delta_Phi_VBS_q_Vlep;

    double W_pT,W_eta,W_phi,W_Quark1_pT, W_Quark1_eta, W_Quark1_mass;
    double W_Quark2_pT, W_Quark2_eta, W_Quark2_mass;
    double W_Quarks_delta_Eta;

    double Tagging_Jet1_pT, Tagging_Jet1_eta, Tagging_Jet1_mass;
    double Mismatch_Tagjet_pT, Mismatch_Tagjet_eta, Mismatch_Tagjet_mass;
    double Tagging_Jet2_pT, Tagging_Jet2_eta, Tagging_Jet2_mass;
    double Tagging_Jets_mass,Tagging_Jets_delta_Eta;
    double Tagging_Jets_VZ_truth_Delta_Phi;

    double Delta_Eta_Tagging_jets_VBS_q_truth, Delta_Mass_Tagging_jets_VBS_q_truth, Delta_pT_Tagging_jets_VBS_q_truth;
    double absDelta_Eta_Tagging_jets_VBS_q_truth, absDelta_Mass_Tagging_jets_VBS_q_truth, absDelta_pT_Tagging_jets_VBS_q_truth;
    double Delta_eta_Tagging_Jet1_VBS_q,Delta_eta_Tagging_Jet2_VBS_q,Delta_eta_Mismatch_Tagjet_VBS_q,Delta_R_Mismatch_Tagjet_VBS_q;
    double Delta_Phi_Tagging_Jet1_VBS_q,Delta_Phi_Tagging_Jet2_VBS_q,Delta_Phi_Mismatch_Tagjet_VBS_q,Delta_eta_Mismatch_Tagjet_Wboson;

    double tag_jets_matched, tag_jet1_matched, tag_jet2_matched;
    std::vector<double> Min_DR_Wq_jets, Min_DR_q_Tagjets,Delta_Eta_TaggingQuark1_Wquarks,Delta_Eta_TaggingQuark2_Wquarks;

    int tag_jet_mismatched_fromVhad;

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
    double merged_fjet_mass,merged_VlepVhad_mass,merged_VlepVhad_pt,merged_VlepVhad_eta,merged_VlepVhad_Deta_tagjets,merged_VlepVhad_Dphi_tagjets,merged_Full_mass;
    double merged_Full_pt,merged_Centrality,merged_CentralityVhad,merged_CentralityVlep,merged_CentralityVlepVhad,merged_ZeppVlep,merged_ZeppVhad,merged_ZeppVhadVlep;
    int merged_n_jets,merged_n_lepton_stable,merged_fjet_n,merged_lepton1_pids,merged_lepton2_pids;
    int merged_Ntrk_tagjets1,merged_Ntrk_tagjets2,merged_Ntrk_tagjets,merged_Ntrk_fjets;
    /// @}


    std::map<std::string, double*> varMapDouble_Tagjets_truth = {


        {"Tagging_Jet1_pT", &Tagging_Jet1_pT},
        {"Tagging_Jet1_eta", &Tagging_Jet1_eta},
        {"Tagging_Jet1_mass", &Tagging_Jet1_mass},

        {"Tagging_Jet2_pT", &Tagging_Jet2_pT},
        {"Tagging_Jet2_eta", &Tagging_Jet2_eta},
        {"Tagging_Jet2_mass", &Tagging_Jet2_mass},
        {"Tagging_Jets_mass", &Tagging_Jets_mass},
        {"Tagging_Jets_delta_Eta", &Tagging_Jets_delta_Eta},

        {"Tagging_Jets_VZ_truth_Delta_Phi", &Tagging_Jets_VZ_truth_Delta_Phi},


        {"Delta_Eta_Tagging_jets_VBS_q_truth", &Delta_Eta_Tagging_jets_VBS_q_truth},
        {"Delta_Mass_Tagging_jets_VBS_q_truth", &Delta_Mass_Tagging_jets_VBS_q_truth},
        {"Delta_pT_Tagging_jets_VBS_q_truth", &Delta_pT_Tagging_jets_VBS_q_truth},

        {"absDelta_Eta_Tagging_jets_VBS_q_truth", &absDelta_Eta_Tagging_jets_VBS_q_truth},
        {"absDelta_Mass_Tagging_jets_VBS_q_truth", &absDelta_Mass_Tagging_jets_VBS_q_truth},
        {"absDelta_pT_Tagging_jets_VBS_q_truth", &absDelta_pT_Tagging_jets_VBS_q_truth},

        {"Delta_eta_Tagging_Jet1_VBS_q", &Delta_eta_Tagging_Jet1_VBS_q},
        {"Delta_eta_Tagging_Jet2_VBS_q", &Delta_eta_Tagging_Jet2_VBS_q},

        {"tag_jets_matched", &tag_jets_matched},
        {"tag_jet1_matched", &tag_jet1_matched},
        {"tag_jet2_matched", &tag_jet2_matched},
        {"Pass_VBS_jet", &Pass_VBS_jet},


    };



    std::map<std::string, double*> varMapDouble_VBS_q_truth = {
        {"VBS_Quarks_truth_Delta_Eta", &VBS_Quarks_truth_Delta_Eta},
        {"VBS_Quarks_truth_Delta_Phi", &VBS_Quarks_truth_Delta_Phi},
        {"VBS_Quark1_truth_pT", &VBS_Quark1_truth_pT},
        {"VBS_Quark1_truth_eta", &VBS_Quark1_truth_eta},
        {"VBS_Quark1_truth_mass", &VBS_Quark1_truth_mass},
        {"VBS_Quark2_truth_pT", &VBS_Quark2_truth_pT},
        {"VBS_Quark2_truth_eta", &VBS_Quark2_truth_eta},
        {"VBS_Quark2_truth_mass", &VBS_Quark2_truth_mass},

        {"VBS_Quarks_truth_pT", &VBS_Quarks_truth_pT},
        {"VBS_Quarks_truth_eta", &VBS_Quarks_truth_eta},

        {"VBS_Quarks_truth_mass", &VBS_Quarks_truth_mass},
        {"Truth_Delta_Eta_VBS_q_VZ", &Truth_Delta_Eta_VBS_q_VZ},
        {"Truth_Delta_Phi_VBS_q_VZ", &Truth_Delta_Phi_VBS_q_VZ},

        {"Delta_Phi_VBS_q_Vlep", &Delta_Phi_VBS_q_Vlep},

    };



    std::map<std::string, double*> varMapDouble_Vhadboson_truth = {
        {"W_Quark1_pT", &W_Quark1_pT},
        {"W_pT", &W_pT},
        {"W_eta", &W_eta},
        {"W_phi", &W_phi},
        {"W_Quark1_eta", &W_Quark1_eta},
        {"W_Quark1_mass", &W_Quark1_mass},
        {"W_Quark2_pT", &W_Quark2_pT},
        {"W_Quark2_eta", &W_Quark2_eta},
        {"W_Quark2_mass", &W_Quark2_mass},
        {"W_Quarks_delta_Eta", &W_Quarks_delta_Eta},

    };

    std::map<std::string, std::vector<double>*> varMapVectorDouble_truth = {
        {"Delta_Eta_TaggingQuark1_Wquarks", &Delta_Eta_TaggingQuark1_Wquarks},
        {"Delta_Eta_TaggingQuark2_Wquarks", &Delta_Eta_TaggingQuark2_Wquarks},
    };


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
        {"merged_VlepVhad_Deta_tagjets", &merged_VlepVhad_Deta_tagjets},
        {"merged_VlepVhad_Dphi_tagjets", &merged_VlepVhad_Dphi_tagjets},
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


    RIVET_DECLARE_PLUGIN(ZZ_llqq);

}