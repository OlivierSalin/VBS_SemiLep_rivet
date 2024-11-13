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

    const Particle GetWboson(const std::vector<Particle>& all_particles){
        Particle W_boson;
        int W_boson_count = 0;

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
     

    int IsTrueWboson(const std::vector<Particle>& all_particles, const Jet& jet, double dR_cut = 0.4, bool merged = false){
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


      // Initialise ROOT file & tree



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

        // Truth information histograms

        std::ifstream truth_hist_file(txt_dir + "/Hists/WpZ_llqq_hists_truth.json");      
        json truth_hist = json::parse(truth_hist_file);
        for (json::iterator it = truth_hist.begin(); it != truth_hist.end(); ++it) {
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

        // plots that are not in other ana
        std::ifstream ana_hist_polar_merged_file(txt_dir + "/Hists/WpZ_llqq_hists_polar_merged.json");      
        json ana_hist_polar_merged = json::parse(ana_hist_polar_merged_file);
        for (json::iterator it = ana_hist_polar_merged.begin(); it != ana_hist_polar_merged.end(); ++it) {
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

        _tf_truth = make_unique<TFile>(getOption("ROOTFILE", ntuple_dir+ "ntuple_truth.root").c_str(), "RECREATE");
        _tt_truth = make_unique<TTree>("Truth", "Rivet_physics");
        for (auto& var_ : varMapDouble_truth) {
            _tt_truth->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }   
        _tt_truth->Branch("Delta_R_tagging_quark", &Delta_R_tagging_quark);
        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        _tt_truth->Branch("FromTruth_W_q_misIDtagjets", &FromTruth_W_q_misIDtagjets);
        _tt_truth->Branch("misID_tag_jets", &misID_tag_jets);
        _tt_truth->Branch("FromTruth_W_q_tag_jets", &FromTruth_W_q_tag_jets); 

        _tt_truth_merged_befcutV = make_unique<TTree>("Truth_merged_befcutV", "Rivet_physics");
        for (auto& var_ : varMapDouble_truth) {
            _tt_truth_merged_befcutV->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }   
        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_merged_befcutV->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        _tt_truth_merged_befcutV->Branch("FromTruth_W_q_misIDtagjets", &FromTruth_W_q_misIDtagjets);
        _tt_truth_merged_befcutV->Branch("misID_tag_jets", &misID_tag_jets);
        _tt_truth_merged_befcutV->Branch("FromTruth_W_q_tag_jets", &FromTruth_W_q_tag_jets);   
        _tt_truth_merged_befcutV->Branch("merged_fjet_FromTruth_W_q", &merged_fjet_FromTruth_W_q),
        _tt_truth_merged_befcutV->Branch("merged_fjetmisID_FromTruth_Tagjet_q", &merged_fjetmisID_FromTruth_Tagjet_q),
        _tt_truth_merged_befcutV->Branch("merged_misID_tag_jets", &merged_misID_tag_jets), 

        _tt_truth_merged = make_unique<TTree>("Truth_merged", "Rivet_physics");
        for (auto& var_ : varMapDouble_truth) {
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }   
        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_merged->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        _tt_truth_merged->Branch("FromTruth_W_q_misIDtagjets", &FromTruth_W_q_misIDtagjets);
        _tt_truth_merged->Branch("misID_tag_jets", &misID_tag_jets);
        _tt_truth_merged->Branch("FromTruth_W_q_tag_jets", &FromTruth_W_q_tag_jets);   
        _tt_truth_merged->Branch("merged_fjet_FromTruth_W_q", &merged_fjet_FromTruth_W_q),
        _tt_truth_merged->Branch("merged_fjetmisID_FromTruth_Tagjet_q", &merged_fjetmisID_FromTruth_Tagjet_q),
        _tt_truth_merged->Branch("merged_misID_tag_jets", &merged_misID_tag_jets), 

        _tt_truth_resolved = make_unique<TTree>("Truth_resolved", "Rivet_physics");
        for (auto& var_ : varMapDouble_truth) {
            _tt_truth_resolved->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }   
        for (auto& var_ : varMapVectorDouble_truth) {
            _tt_truth_resolved->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }
        _tt_truth_resolved->Branch("FromTruth_W_q_misIDtagjets", &FromTruth_W_q_misIDtagjets);
        _tt_truth_resolved->Branch("misID_tag_jets", &misID_tag_jets);
        _tt_truth_resolved->Branch("FromTruth_W_q_tag_jets", &FromTruth_W_q_tag_jets);    
        for (auto& var_ : varMapInt_truth_resolved) {
            _tt_truth_resolved->Branch(var_.first.c_str(), var_.second);
            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        }   
        //_tf = make_unique<TFile>(getOption("ROOTFILE", ntuple_dir+ "ntuple_rivet_truth.root").c_str(), "RECREATE");

            //std::cout << "Variable name: " << var_.first << ", Variable value: " << var_.second << std::endl;
        

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
                            "n_jets","At_least_one_fjets","fjets_is_W/Z","found_tag_jets","pt_tagjet1_2","m_tagjets",
                            "Total_Merged_selec",});
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

        const Particles all_particles = event.allParticles();
        
        // Fill the truth histograms
        if(Check_VBS_event(all_particles)){
            _h["bef_cut_Is_VBS_event"]->fill(1);
        } else {
            _h["bef_cut_Is_VBS_event"]->fill(0);
        }
        const double cutof_DR = 0.4;

        //std::cout << "Type of event.genEvent()->event_number(): " << typeid(event.genEvent()->event_number()).name() << std::endl;
        //std::cout << "Value of event.genEvent()->event_number(): " << EvntNumber << std::endl;

        if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
        else {_c["neg_w_initial"]->fill();}

        // const Particles all_particles = event.allParticles();
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
        _h["bef_cut_nlep"]->fill(nlep);
        if (nlep !=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];
        FourMomentum lepton1 = lep1.mom();
        FourMomentum lepton2 = lep2.mom();

        const Particle lepton_minus = (lep1.charge() < 0) ? lep1 : lep2;
        const Particle lepton_plus = (lep1.charge() > 0) ? lep1 : lep2;

/*         const Particle parent_lep1 = GetParent(lep1);
        const Particle parent_lep2 = GetParent(lep2); */


        // Cuts on the pT of the leptons
        if (_docut==1 && (leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();
        

        const FourMomentum fourvec_ll = lep1.mom() + lep2.mom();
        FourMomentum fourvec_Z = lep1.mom() + lep2.mom();
        double m_ll = fourvec_ll.mass()/GeV;



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
        Jets fjets_ = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > 200*GeV && Cuts::abseta < 2.0);
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
        idiscardIfAnyDeltaRLess(fjets, leptons, 1.0);
        int n_fjets = fjets.size();

        bool merged = false, resolved = true;
        if (n_fjets > 0) {
            const FourMomentum fourvec_fjets = fjets[0].mom();
            _cutflows_merged.fillnext();

            std::vector<double> minDeltaRs_merged = Truth_q_minDR_jets(all_particles, fjets, true);
            //std::cout << "\nTruth merged min DRs: ";
            if(Check_VBS_event(all_particles)){
                for(size_t i = 0; i < minDeltaRs_merged.size(); ++i){
                    _h["merged_min_DR_Wq_fjets"]->fill(minDeltaRs_merged[i]);
                    //std::cout << minDeltaRs_merged[i] << ",  ";
                }
                //std::cout << "\n";

                merged_fjet_FromTruth_W_q = IsTrueWboson(all_particles, fjets[0], cutof_DR, true);
                _h["merged_fjet_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, fjets[0], cutof_DR, true));
                if(IsTrueWboson(all_particles, fjets[0], cutof_DR, true)==0){
                    merged_fjetmisID_FromTruth_Tagjet_q= IsTruthTagJet(all_particles, fjets[0], cutof_DR);
                    _h["merged_fjetmisID_FromTruth_Tagjet_q"]->fill(IsTruthTagJet(all_particles, fjets[0], cutof_DR));
                }

            }

            _tt_truth_merged_befcutV->Fill();

            if ((_docut==1 && (fourvec_fjets.mass() >= _jcuts["m_fjet_WZ"][0] && fourvec_fjets.mass() <= _jcuts["m_fjet_WZ"][1])) || _docut==0) {
                merged = true;
                resolved = false;
                _cutflows_merged.fillnext();
            }
        }
        bool closeToFjets;
        Jets eligibleJets;
        if (merged) {
            for (const Jet& jet : jets) {
                closeToFjets = false;
                if (deltaR(jet, fjets[0]) <= 1.4) {
                    closeToFjets = true;
                }

                if (!closeToFjets) {
                    eligibleJets.push_back(jet);
                }
            }
            }
        else {
            eligibleJets = jets;
        }
        if(merged){
            for (const Jet& jet : eligibleJets) {
                //std::cout << "Delta R eligible jet and fjet: " << deltaR(jet, fjets[0]) << std::endl;  

            }
        }

        Jets tag_jets;
        bool foundVBSJetPair = false; 
        // Now use eligibleJets instead of jets for finding VBS tagging jets
        double max_m_tag_jj = 0;
        int tag1_jet_index = -1, tag2_jet_index = -1;
        for (int i = 0; i < eligibleJets.size(); i++) {
            const Jet& i_jet = eligibleJets[i];
            for (int j = i + 1; j < eligibleJets.size(); j++) {   
                const Jet& j_jet = eligibleJets[j];
                const double m_tag_jj = (i_jet.mom() + j_jet.mom()).mass()/GeV;
                const double eta_prod = i_jet.eta()*j_jet.eta();
                if (eta_prod < 0.0 && m_tag_jj > max_m_tag_jj) {
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

        tag_jets.push_back(eligibleJets[tag1_jet_index]);
        tag_jets.push_back(eligibleJets[tag2_jet_index]);

        bool foundTagJet1 = false;
        bool foundTagJet2 = false;



        const FourMomentum tag1_jet = tag_jets[0].mom();
        const FourMomentum tag2_jet = tag_jets[1].mom();
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

        ///// TRUTH INFORMATION /////

        if(Check_VBS_event(all_particles)){
            const Particle W_bson = GetWboson(all_particles);
            std::vector<Particle> W_quarks= ValidQuark_W(W_bson);
            
            //printf("w_boson pt: %f",W_bson.mom().pt());
            
            std::vector<double> minDeltaRs = Truth_q_minDR_jets(all_particles, jets);
            //std::cout << "\nTruth min DRs: ";
            Min_DR_Wq_jets.clear();
            for(size_t i = 0; i < minDeltaRs.size(); ++i){
                _h["Min_DR_Wq_jets"]->fill(minDeltaRs[i]);
                Min_DR_Wq_jets.push_back(minDeltaRs[i]);
                //std::cout << minDeltaRs[i] << " ";
            }
            //std::cout << "\n";
            

            std::vector<double> minDeltaRs_tag = Truth_q_minDR_jets(all_particles, tag_jets);
            //std::cout << "\nTruth Tag min DRs: ";
            Min_DR_q_Tagjets.clear();
            for(size_t i = 0; i < minDeltaRs_tag.size(); ++i){
                _h["Min_DR_q_Tagjets"]->fill(minDeltaRs_tag[i]);
                Min_DR_q_Tagjets.push_back(minDeltaRs_tag[i]);
                //std::cout << minDeltaRs_tag[i] << " ";
            }
            //std::cout << "\n";


            //printf("\nTagging Quarks \n");
            std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
            if(taggingQuarks.size() == 2){
                int quarkCount = 0;
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
                    double minDR_VBSq_jets = minDeltaR_;
                    _h["Min_DR_VBSq_jets"]->fill(minDR_VBSq_jets);
                    if(quarkCount == 0){
                        minDR_VBSq1_jets = minDeltaR_;
                    }
                    else if(quarkCount == 1){
                        minDR_VBSq2_jets = minDeltaR_;
                    }
                    quarkCount++;

                }
            }
            
            
            int tag_jet_misID = 0;
            double dR_cut = 0.2;
            //std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
            if(taggingQuarks.size() != 2){
                tag_jet_misID = -1; // Set to -1 to indicate a bad event
            } else {
                FromTruth_W_q_misIDtagjets.clear();
                for(const Jet& jet : tag_jets){
                    bool truth_DRmatch = false;
                    for(const Particle& taggingQuark : taggingQuarks){
                        const double dR = deltaR(taggingQuark, jet);
                        if(dR < dR_cut) {
                            truth_DRmatch = true;
                            break;
                        }
                    }
                    if(!truth_DRmatch) {
                        tag_jet_misID++;
                        _h["FromTruth_W_q_misIDtagjets"]->fill(IsTrueWboson(all_particles, jet, cutof_DR, false));
                        FromTruth_W_q_misIDtagjets.push_back(IsTrueWboson(all_particles, jet, cutof_DR, false));
                    }
                }
            }
            //printf("\nTagging Jets MisID: %d\n", tag_jet_misID);
            //printf("NbTagJetMisID: %d\n", NbTagJetMisID(all_particles, tag_jets, 0.1));  
            misID_tag_jets= tag_jet_misID;
            _h["misID_tag_jets"]->fill(tag_jet_misID);
            FromTruth_W_q_tag_jets = IsTrueWboson(all_particles, tag_jets[0], cutof_DR, false)+ IsTrueWboson(all_particles, tag_jets[1], cutof_DR, false);
            _h["FromTruth_W_q_tagjets"]->fill(IsTrueWboson(all_particles, tag_jets[0], cutof_DR, false)+ IsTrueWboson(all_particles, tag_jets[1], cutof_DR, false));
            _tt_truth->Fill();
        }
        std::vector<Particle> taggingQuarks_ = Tagging_quarks(all_particles);
        if(taggingQuarks_.size() == 2){
            
            Delta_R_tagging_quark= deltaR(taggingQuarks_[0], taggingQuarks_[1]);
            _h["Delta_R_tagging_quark"]->fill(deltaR(taggingQuarks_[0], taggingQuarks_[1]));
            Delta_Eta_tagging_quark = abs(taggingQuarks_[0].eta() - taggingQuarks_[1].eta());
            Tagging_Quark1_pT = taggingQuarks_[0].pT();
            Tagging_Quark1_eta = taggingQuarks_[0].eta();
            Tagging_Quark1_mass = taggingQuarks_[0].mom().mass();

            Tagging_Quark2_pT = taggingQuarks_[1].pT();
            Tagging_Quark2_eta = taggingQuarks_[1].eta();
            Tagging_Quark2_mass = taggingQuarks_[1].mom().mass();

            Tagging_Quarks_mass= (taggingQuarks_[0].mom() + taggingQuarks_[1].mom()).mass();

            Tagging_Jet1_pT = tag_jets[0].pT();
            Tagging_Jet1_eta = tag_jets[0].eta();
            Tagging_Jet1_mass = tag_jets[0].mom().mass();

            Tagging_Jet2_pT = tag_jets[1].pT();
            Tagging_Jet2_eta = tag_jets[1].eta();
            Tagging_Jet2_mass = tag_jets[1].mom().mass();

            Tagging_Jets_mass = (tag1_jet + tag2_jet).mass();
            Tagging_Jets_delta_Eta = abs(tag1_jet.eta() - tag2_jet.eta());

            const Particle W_bson_ = GetWboson(all_particles);
            W_pT = W_bson_.pT();
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
            Delta_Eta_TaggingQuark1_Wquarks.clear();
            Delta_Eta_TaggingQuark2_Wquarks.clear();
            for (auto& W_quark : W_quarks_) {
                Delta_Eta_TaggingQuark1_Wquarks.push_back(abs(taggingQuarks_[0].eta() - W_quark.eta()));
                Delta_Eta_TaggingQuark2_Wquarks.push_back(abs(taggingQuarks_[1].eta() - W_quark.eta()));
            }
        }

            
        // Check if we are in the Merged signal region

        if (n_fjets>0 && merged) {
            const FourMomentum fourvec_fjets = fjets[0].mom();



            const double beta = 1;
            const fastjet::PseudoJet &LJet= tr_ljets[0];
            fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
            fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

            double recf3 = ECF3(LJet);
            double recf2 = ECF2(LJet);
            double recf1 = ECF1(LJet);
            double d2_fjets = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);            
            // Total cutflow of the merged region

            _cutflows_merged.fillnext();
            merged_misID_tag_jets = NbTagJetMisID(all_particles, tag_jets, cutof_DR);
            _h["merged_misID_tag_jets"]->fill(NbTagJetMisID(all_particles, tag_jets, cutof_DR));

            _tt_truth_merged->Fill();

            Vector3 Axis_lab(1,1,1);
            Vector3 Axis_x_lab(1,0,0);Vector3 Axis_y_lab(0,1,0);Vector3 Axis_z_lab(0,0,1);
            // Four vector of the QGC system in resolved region llJ
            const FourMomentum fourvec_fjets_ll = lep1.mom() + lep2.mom() + fjets[0].mom();
            FourMomentum fourvec_V = fjets[0].mom();
            // Polarization variable VZ system
            FourMomentum fourvec_VZ = lep1.mom() + lep2.mom() + fjets[0].mom();




            LorentzTransform boost_VZ;
            FourMomentum fourvec_VZ_rotZ = fourvec_VZ;
            FourMomentum taggjet1 = tag1_jet;
            FourMomentum taggjet2 = tag2_jet;

            double alpha_ = atan2(fourvec_VZ_rotZ.py(), fourvec_VZ_rotZ.px()); 
            fourvec_VZ_rotZ = RotateZ(-alpha_, 1, fourvec_VZ_rotZ);


            boost_VZ.setBetaVec(-fourvec_VZ_rotZ.betaVec());
            Vector3 BoostVZ_vec= -fourvec_VZ.betaVec();
            const Vector3 unitboostVZvec = BoostVZ_vec.unit();

            FourMomentum fourvec_Z_CS = boost_VZ.transform(RotateZ(-alpha_, 1,fourvec_Z));
            FourMomentum fourvec_VZ_CS = boost_VZ.transform(RotateZ(-alpha_, 1,fourvec_VZ));
            FourMomentum fourvec_V_CS = boost_VZ.transform(RotateZ(-alpha_, 1,fourvec_V));
            FourMomentum fourvec_lep1_CS = boost_VZ.transform(RotateZ(-alpha_, 1,lepton1));
            FourMomentum fourvec_lep2_CS = boost_VZ.transform(RotateZ(-alpha_, 1,lepton2));
            FourMomentum fourvec_Tagjet1_CS = boost_VZ.transform(RotateZ(-alpha_, 1,taggjet1));
            FourMomentum fourvec_Tagjet2_CS = boost_VZ.transform(RotateZ(-alpha_, 1,taggjet2));


            // Calculate the theta and phi angles for each variable in the CS frame and fill the histograms
            merged_CS_Z_cos_theta = cos(fourvec_Z_CS.p3().theta());
            _h["merged_CS_Z_cos_theta"]->fill(merged_CS_Z_cos_theta);
            merged_CS_Z_eta = fourvec_Z_CS.eta();
            _h["merged_CS_Z_eta"]->fill(merged_CS_Z_eta);
            merged_CS_Z_phi = fourvec_Z_CS.p3().phi();
            _h["merged_CS_Z_phi"]->fill(merged_CS_Z_phi);
            merged_CS_Z_pt = fourvec_Z_CS.pT();
            _h["merged_CS_Z_pt"]->fill(merged_CS_Z_pt);
            merged_CS_lepton_pt1 = fourvec_lep1_CS.pT();
            _h["merged_CS_lepton_pt"]->fill(merged_CS_lepton_pt1);
            merged_CS_lepton_pt2 = fourvec_lep2_CS.pT();
            _h["merged_CS_lepton_pt"]->fill(merged_CS_lepton_pt2);

            
            merged_CS_lepton_deltaEta = deltaR(fourvec_lep1_CS, fourvec_lep2_CS);
            _h["merged_CS_lepton_deltaEta"]->fill(merged_CS_lepton_deltaEta);
            merged_CS_lepton_deltaPhi = deltaPhi(fourvec_lep1_CS, fourvec_lep2_CS);
            _h["merged_CS_lepton_deltaPhi"]->fill(merged_CS_lepton_deltaPhi);
            merged_CS_lepton_deltaR_Vhad = deltaPhi(fourvec_lep1_CS, fourvec_V_CS);
            _h["merged_CS_lepton_deltaR_Vhad"]->fill(merged_CS_lepton_deltaR_Vhad);
            merged_CS_V_cos_theta = cos(fourvec_V_CS.p3().theta());
            _h["merged_CS_V_cos_theta"]->fill(merged_CS_V_cos_theta);
            merged_CS_V_eta = fourvec_V_CS.eta();
            _h["merged_CS_V_eta"]->fill(merged_CS_V_eta);
            merged_CS_V_phi = fourvec_V_CS.p3().phi();
            _h["merged_CS_V_phi"]->fill(merged_CS_V_phi);
            merged_CS_V_pt = fourvec_V_CS.pT();
            _h["merged_CS_V_pt"]->fill(merged_CS_V_pt);
            merged_CS_Tagjet1_theta = fourvec_Tagjet1_CS.eta();
            _h["merged_CS_Tagjet1_theta"]->fill(merged_CS_Tagjet1_theta);
            merged_CS_Tagjet1_phi = fourvec_Tagjet1_CS.p3().phi();
            _h["merged_CS_Tagjet1_phi"]->fill(merged_CS_Tagjet1_phi);
            merged_CS_Tagjet2_theta = fourvec_Tagjet2_CS.eta();
            _h["merged_CS_Tagjet2_theta"]->fill(merged_CS_Tagjet2_theta);
            merged_CS_Tagjet2_phi = fourvec_Tagjet2_CS.p3().phi();
            _h["merged_CS_Tagjet2_phi"]->fill(merged_CS_Tagjet2_phi);

            // Calculate the differences in phi, R, and eta betVeen the variables in the CS frame and fill the histograms
            merged_CS_deltaPhi_ZV = deltaPhi(fourvec_Z_CS, fourvec_V_CS);
            _h["merged_CS_deltaPhi_ZV"]->fill(merged_CS_deltaPhi_ZV);
            merged_CS_deltaR_ZV = deltaR(fourvec_Z_CS, fourvec_V_CS);
            _h["merged_CS_deltaR_ZV"]->fill(merged_CS_deltaR_ZV);
            merged_CS_deltaEta_ZV = deltaEta(fourvec_Z_CS, fourvec_V_CS);
            _h["merged_CS_deltaEta_ZV"]->fill(merged_CS_deltaEta_ZV);
            merged_CS_deltaPhi_ZTagjet1 = deltaPhi(fourvec_Z_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaPhi_ZTagjet1"]->fill(merged_CS_deltaPhi_ZTagjet1);
            merged_CS_deltaR_ZTagjet1 = deltaR(fourvec_Z_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaR_ZTagjet1"]->fill(merged_CS_deltaR_ZTagjet1);
            merged_CS_deltaEta_ZTagjet1 = deltaEta(fourvec_Z_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaEta_ZTagjet1"]->fill(merged_CS_deltaEta_ZTagjet1);
            merged_CS_deltaPhi_ZTagjet2 = deltaPhi(fourvec_Z_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaPhi_ZTagjet2"]->fill(merged_CS_deltaPhi_ZTagjet2);
            merged_CS_deltaR_ZTagjet2 = deltaR(fourvec_Z_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaR_ZTagjet2"]->fill(merged_CS_deltaR_ZTagjet2);
            merged_CS_deltaEta_ZTagjet2 = deltaEta(fourvec_Z_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaEta_ZTagjet2"]->fill(merged_CS_deltaEta_ZTagjet2);
            merged_CS_deltaPhi_VTagjet1 = deltaPhi(fourvec_V_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaPhi_VTagjet1"]->fill(merged_CS_deltaPhi_VTagjet1);
            merged_CS_deltaR_VTagjet1 = deltaR(fourvec_V_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaR_VTagjet1"]->fill(merged_CS_deltaR_VTagjet1);
            merged_CS_deltaEta_VTagjet1 = deltaEta(fourvec_V_CS, fourvec_Tagjet1_CS);
            _h["merged_CS_deltaEta_VTagjet1"]->fill(merged_CS_deltaEta_VTagjet1);
            merged_CS_deltaPhi_VTagjet2 = deltaPhi(fourvec_V_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaPhi_VTagjet2"]->fill(merged_CS_deltaPhi_VTagjet2);
            merged_CS_deltaR_VTagjet2 = deltaR(fourvec_V_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaR_VTagjet2"]->fill(merged_CS_deltaR_VTagjet2);
            merged_CS_deltaEta_VTagjet2 = deltaEta(fourvec_V_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaEta_VTagjet2"]->fill(merged_CS_deltaEta_VTagjet2);
            merged_CS_deltaPhi_Tagjet1Tagjet2 = deltaPhi(fourvec_Tagjet1_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaPhi_Tagjet1Tagjet2"]->fill(merged_CS_deltaPhi_Tagjet1Tagjet2);
            merged_CS_deltaR_Tagjet1Tagjet2 = deltaR(fourvec_Tagjet1_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaR_Tagjet1Tagjet2"]->fill(merged_CS_deltaR_Tagjet1Tagjet2);
            merged_CS_deltaEta_Tagjet1Tagjet2 = deltaEta(fourvec_Tagjet1_CS, fourvec_Tagjet2_CS);
            _h["merged_CS_deltaEta_Tagjet1Tagjet2"]->fill(merged_CS_deltaEta_Tagjet1Tagjet2);


            LorentzTransform boost_Zrf;
            boost_Zrf.setBetaVec(-fourvec_Z.betaVec());
            FourMomentum fourvec_Z_Zrf = boost_Zrf.transform(fourvec_Z);
            FourMomentum fourvec_V_Zrf = boost_Zrf.transform(fourvec_V);
            FourMomentum fourvec_lep1_Zrf = boost_Zrf.transform(lep1.mom());
            FourMomentum fourvec_lep2_Zrf = boost_Zrf.transform(lep2.mom());
            FourMomentum fourvec_lep_minus_Zrf = boost_Zrf.transform(lepton_minus.mom());
            FourMomentum fourvec_lep_plus_Zrf = boost_Zrf.transform(lepton_plus.mom());

            LorentzTransform boost_VZ_bis;
            boost_VZ_bis.setBetaVec(-fourvec_VZ.betaVec());
            FourMomentum fourvec_Z_VZ_bis = boost_VZ_bis.transform(fourvec_Z);
            FourMomentum fourvec_lep1_VZ_bis = boost_VZ_bis.transform(lep1.mom());
            FourMomentum fourvec_lep2_VZ_bis = boost_VZ_bis.transform(lep2.mom());
            FourMomentum fourvec_lep_minus_VZ_bis = boost_VZ_bis.transform(lepton_minus.mom());
            LorentzTransform boost_Zrf_bis;
            boost_Zrf_bis.setBetaVec(-fourvec_Z_VZ_bis.betaVec());
            FourMomentum fourvec_Z_Zrf_bis = boost_Zrf_bis.transform(fourvec_Z_VZ_bis);
            FourMomentum fourvec_V_Zrf_bis = boost_Zrf_bis.transform(fourvec_V);
            FourMomentum fourvec_lep1_Zrf_bis = boost_Zrf_bis.transform(fourvec_lep1_VZ_bis);
            FourMomentum fourvec_lep2_Zrf_bis = boost_Zrf_bis.transform(fourvec_lep2_VZ_bis);
            FourMomentum fourvec_lep_minus_Zrf_bis = boost_Zrf_bis.transform(fourvec_lep_minus_VZ_bis);

            merged_Zrf_Lepton1_cos_theta = cos(fourvec_lep1_Zrf.p3().theta());
            _h["merged_Zrf_Lepton1_cos_theta"]->fill(merged_Zrf_Lepton1_cos_theta);
            merged_Zrf_Lepton2_cos_theta = cos(fourvec_lep2_Zrf.p3().theta());
            _h["merged_Zrf_Lepton2_cos_theta"]->fill(merged_Zrf_Lepton2_cos_theta);
            merged_Zrf_Lepton_minus_cos_theta = cos(fourvec_lep_minus_Zrf.p3().theta());
            _h["merged_Zrf_Lepton_minus_cos_theta"]->fill(merged_Zrf_Lepton_minus_cos_theta);
            merged_Zrf_Lepton_minus_eta = fourvec_lep_minus_Zrf.eta();
            _h["merged_Zrf_Lepton_minus_eta"]->fill(merged_Zrf_Lepton_minus_eta);

            _h["merged_Zrf_Lepton_pt"]->fill(fourvec_lep1_Zrf.pt());_h["merged_Zrf_Lepton_pt"]->fill(fourvec_lep2_Zrf.pt());
            merged_Zrf_Lepton1_pt = fourvec_lep1_Zrf.pt();
            _h["merged_Zrf_Lepton1_pt"]->fill(merged_Zrf_Lepton1_pt);
            merged_Zrf_Lepton2_pt = fourvec_lep2_Zrf.pt();
            _h["merged_Zrf_Lepton2_pt"]->fill(merged_Zrf_Lepton2_pt);

            merged_Zrf_Lepton_minus_cos_theta_Zlab_axis = cos(fourvec_lep_minus_Zrf.p3().angle(fourvec_Z.p3()));
            _h["merged_Zrf_Lepton_minus_cos_theta_Zlab_axis"]->fill(merged_Zrf_Lepton_minus_cos_theta_Zlab_axis);

            merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis = cos(fourvec_lep_minus_Zrf_bis.p3().angle(fourvec_Z_VZ_bis));
            _h["merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis"]->fill(merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis);

            merged_Zrf_Lepton_deltaEta = deltaEta(fourvec_lep1_Zrf, fourvec_lep2_Zrf);
            _h["merged_Zrf_Lepton_deltaEta"]->fill(merged_Zrf_Lepton_deltaEta);

            merged_Zrf_Lepton_tanh_deltaEta= tanh(0.5*(fourvec_lep_minus_Zrf.eta()-fourvec_lep_plus_Zrf.eta()));
            _h["merged_Zrf_Lepton_deltaEta"]->fill(merged_Zrf_Lepton_tanh_deltaEta);

            merged_Zrf_Lepton_deltaPhi = deltaPhi(fourvec_lep1_Zrf, fourvec_lep2_Zrf);
            _h["merged_Zrf_Lepton_deltaPhi"]->fill(merged_Zrf_Lepton_deltaPhi);
            merged_Zrf_Lepton_deltaPhi_Vhad = deltaPhi(fourvec_lep1_Zrf, fourvec_V_Zrf);
            _h["merged_Zrf_Lepton_deltaPhi_Vhad"]->fill(merged_Zrf_Lepton_deltaPhi_Vhad);

            // Four vector of the Full system in resolved region 
            const FourMomentum fourvec_fjets_full = lep1.mom() + lep2.mom() + fjets[0].mom() + tag1_jet + tag2_jet;


            double ZeppMerged = 0.0;
            if (n_fjets > 1) {
                const FourMomentum fourvec_fjets2 = fjets[1].mom();
                ZeppMerged = abs(fourvec_fjets2.eta() - eta_tag_jet_mean);
                _h["merged_ZeppMerged"]->fill(ZeppMerged);
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
            //ana-specific
            merged_ll_mass = m_ll;
            _h["merged_ll_mass"]->fill(merged_ll_mass);
            merged_ll_pt = fourvec_ll.pT();
            _h["merged_ll_pt"]->fill(merged_ll_pt);
            merged_ll_eta = fourvec_ll.eta();
            _h["merged_ll_eta"]->fill(merged_ll_eta);
            merged_ll_phi = fourvec_ll.phi();
            _h["merged_ll_phi"]->fill(merged_ll_phi);
            merged_ll_DR_tagjet1 = deltaR(fourvec_ll, tag1_jet);
            _h["merged_ll_DR_tagjet"]->fill(merged_ll_DR_tagjet1);
            merged_ll_DR_tagjet2 = deltaR(fourvec_ll, tag2_jet);
            _h["merged_ll_DR_tagjet"]->fill(merged_ll_DR_tagjet2);
            merged_ll_Dphi_tagjet1 = deltaPhi(fourvec_ll, tag1_jet);
            _h["merged_ll_Dphi_tagjet"]->fill(merged_ll_Dphi_tagjet1);
            merged_ll_Dphi_tagjet2 = deltaPhi(fourvec_ll, tag2_jet);
            _h["merged_ll_Dphi_tagjet"]->fill(merged_ll_Dphi_tagjet2);
            merged_ll_Deta_tagjet1 = abs(fourvec_ll.eta() - tag1_jet.eta());
            _h["merged_ll_Deta_tagjet"]->fill(merged_ll_Deta_tagjet1);
            merged_ll_Deta_tagjet2 = abs(fourvec_ll.eta() - tag2_jet.eta());
            _h["merged_ll_Deta_tagjet"]->fill(merged_ll_Deta_tagjet2);
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
            merged_Vhad_DeltaR_Z = deltaR(fourvec_fjets, fourvec_ll);
            _h["merged_Vhad_DeltaR_Z"]->fill(merged_Vhad_DeltaR_Z);
            merged_Vhad_DeltaPhi_Z = deltaPhi(fourvec_fjets, fourvec_ll);
            _h["merged_Vhad_DeltaPhi_Z"]->fill(merged_Vhad_DeltaPhi_Z);
            merged_Vhad_DeltaEta_Z = abs(fourvec_fjets.eta() - fourvec_ll.eta());
            _h["merged_Vhad_DeltaEta_Z"]->fill(merged_Vhad_DeltaEta_Z);
            merged_fjet_mass = fourvec_fjets.mass();
            _h["merged_fjet_mass"]->fill(merged_fjet_mass);
            merged_VZ_mass = fourvec_fjets_ll.mass();
            _h["merged_VZ_mass"]->fill(merged_VZ_mass);
            merged_VZ_pt = fourvec_fjets_ll.pt();
            _h["merged_VZ_pt"]->fill(merged_VZ_pt);
            merged_VZ_eta = fourvec_fjets_ll.eta();
            _h["merged_VZ_eta"]->fill(merged_VZ_eta);
            merged_Full_mass = fourvec_fjets_full.mass();
            _h["merged_Full_mass"]->fill(merged_Full_mass);
            merged_Full_pt = fourvec_fjets_full.pt();
            _h["merged_Full_pt"]->fill(merged_Full_pt);

            // Centrality variables
            merged_Centrality = Centrality(tag_jets, fourvec_ll, fourvec_fjets);
            merged_CentralityV = Centrality(tag_jets, fourvec_fjets, fourvec_fjets);
            merged_CentralityZ = Centrality(tag_jets, fourvec_ll,fourvec_ll);
            merged_CentralityZV = Centrality(tag_jets, fourvec_fjets_ll, fourvec_fjets_ll);
            _h["merged_Centrality"]->fill(merged_Centrality);
            _h["merged_CentralityVhad"]->fill(merged_CentralityV);
            _h["merged_CentralityZ"]->fill(merged_CentralityZ);
            _h["merged_CentralityZV"]->fill(merged_CentralityZV);

            // Zeppenfeld variables
            merged_ZeppZ = abs(fourvec_ll.eta() - eta_tag_jet_mean);
            merged_ZeppV = abs(fourvec_fjets.eta() - eta_tag_jet_mean);
            merged_ZeppZV = abs(fourvec_fjets_ll.eta() - eta_tag_jet_mean);
            _h["merged_ZeppZ"]->fill(merged_ZeppZ);
            _h["merged_ZeppV"]->fill(merged_ZeppV);
            _h["merged_ZeppZV"]->fill(merged_ZeppZV);
            
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

        // Check if we are in the Resolved signal region
        } else {            
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

            resolved_sjets_FromTruth_W_q = IsTrueWboson(all_particles, sjets_sig[0], cutof_DR, false)+IsTrueWboson(all_particles, sjets_sig[1], cutof_DR, false);
            _h["resolved_sjets_FromTruth_W_q"]->fill(IsTrueWboson(all_particles, sjets_sig[0], cutof_DR, false)+IsTrueWboson(all_particles, sjets_sig[1], cutof_DR, false));
            int result1 = IsTrueWboson(all_particles, sjets_sig[0], cutof_DR, false);
            resolved_sjet1_FromTruth_W_q = result1;
            _h["resolved_sjet1_FromTruth_W_q"]->fill(result1);
            if(result1 == 0) {
                resolved_sjet1misID_FromTruth_Tagjet_q = IsTruthTagJet(all_particles, sjets_sig[0], cutof_DR);
                _h["resolved_sjetsmisID_FromTruth_Tagjet_q"]->fill(IsTruthTagJet(all_particles, sjets_sig[0], cutof_DR));
            }

            int result2 = IsTrueWboson(all_particles, sjets_sig[1], cutof_DR, false);
            resolved_sjet2_FromTruth_W_q = result2;
            _h["resolved_sjet2_FromTruth_W_q"]->fill(result2);
            if(result2 == 0) {
                resolved_sjet2misID_FromTruth_Tagjet_q = IsTruthTagJet(all_particles, sjets_sig[1], cutof_DR);
                _h["resolved_sjetsmisID_FromTruth_Tagjet_q"]->fill(IsTruthTagJet(all_particles, sjets_sig[1], cutof_DR));
            }
            
            
            // Make a cut on signal_mjj
            const FourMomentum fourvec_signal_jets = signal_jet1 + signal_jet2;
            const FourMomentum fourvec_Vhad = signal_jet1 + signal_jet2;
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if ((_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1]) )) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();
            _h["resolved_misID_tag_jets"]->fill(NbTagJetMisID(all_particles, tag_jets, cutof_DR));  
            _tt_truth_resolved->Fill();          
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
            _h["resolved_ll_mass"]->fill(m_ll);
            _h["resolved_ll_pt"]->fill(fourvec_ll.pT());
            _h["resolved_ll_eta"]->fill(fourvec_ll.eta());
            _h["resolved_ll_phi"]->fill(fourvec_ll.phi());
            _h["resolved_ll_DR_tagjet"]->fill(deltaR(fourvec_ll, tag1_jet));_h["resolved_ll_DR_tagjet"]->fill(deltaR(fourvec_ll, tag2_jet));
            _h["resolved_ll_Dphi_tagjet"]->fill(deltaPhi(fourvec_ll, tag1_jet));_h["resolved_ll_Dphi_tagjet"]->fill(deltaPhi(fourvec_ll, tag2_jet));
            _h["resolved_ll_Deta_tagjet"]->fill(abs(fourvec_ll.eta() - tag1_jet.eta()));_h["resolved_ll_Deta_tagjet"]->fill(abs(fourvec_ll.eta() - tag2_jet.eta()));
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


            _h["resolved_Vhad_DeltaR_Z"]->fill(deltaR(fourvec_Vhad, fourvec_ll));
            _h["resolved_Vhad_DeltaPhi_Z"]->fill(deltaPhi(fourvec_Vhad, fourvec_ll));
            _h["resolved_Vhad_DeltaEta_Z"]->fill(abs(fourvec_Vhad.eta() - fourvec_ll.eta()));

            _h["resolved_VZ_mass"]->fill(fourvec_signal_jets_ll.mass());                
            _h["resolved_VZ_pt"]->fill(fourvec_signal_jets_ll.pt());                
            _h["resolved_VZ_eta"]->fill(fourvec_signal_jets_ll.eta());                
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
        if(nsys<1) {
            _tf->Write();
            _tf_truth->Write();
        
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
            double pos_w_sum_final = dbl(*_c["pos_w_final"]);
            double neg_w_sum_final = dbl(*_c["neg_w_final"]);
            double pos_w_sum_final_merged = dbl(*_c["pos_w_final_merged"]);
            double neg_w_sum_final_merged = dbl(*_c["neg_w_final_merged"]);
            double pos_w_sum_final_resolved = dbl(*_c["pos_w_final_resolved"]);
            double neg_w_sum_final_resolved = dbl(*_c["neg_w_final_resolved"]); */
            //MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
            //MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

    /*         MSG_INFO("\n pos weights initial final ratio Merged " << pos_w_sum_initial <<" " << pos_w_sum_final_merged <<" "<< pos_w_sum_final_merged/pos_w_sum_initial << "\n" );
            MSG_INFO("\n neg weights initial final ratio Merged " << neg_w_sum_initial <<" " << neg_w_sum_final_merged <<" "<< neg_w_sum_final_merged/neg_w_sum_initial << "\n" );

            MSG_INFO("\n pos weights initial final ratio Resolved  " << pos_w_sum_initial <<" " << pos_w_sum_final_resolved <<" "<< pos_w_sum_final_resolved/pos_w_sum_initial << "\n" );
            MSG_INFO("\n neg weights initial final ratio Resolved " << neg_w_sum_initial <<" " << neg_w_sum_final_resolved <<" "<< neg_w_sum_final_resolved/neg_w_sum_initial << "\n" );
    */

            // normalize all to 1 since in case of mostly negative weights not clear what it will do
/*             const double total_weight_initial = pos_w_sum_initial + neg_w_sum_initial;
            const double weight_merged = pos_w_sum_final_merged + neg_w_sum_final_merged;
            const double weight_resolved = pos_w_sum_final_resolved + neg_w_sum_final_resolved;


            const double xs = crossSection()/femtobarn;
            const double lumi= 139.0; // fb-1 */
            //cout << "Cross section: " << xs << " fb" << endl;
            //cout << "Lumi scaled number of event: " << xs*lumi  << endl;
            //cout << "Scale fator : " << xs/sumOfWeights()  << endl;

            //cout << "Initial number of events: " << total_weight_initial << endl;

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
    double merged_CS_Z_cos_theta,merged_CS_Z_eta,merged_CS_Z_phi,merged_CS_Z_pt,merged_CS_lepton_pt1;
    double merged_CS_lepton_pt2,merged_CS_lepton_deltaEta,merged_CS_lepton_deltaPhi,merged_CS_lepton_deltaR_Vhad,merged_CS_V_cos_theta;
    double merged_CS_V_eta,merged_CS_V_phi,merged_CS_V_pt,merged_CS_Tagjet1_theta,merged_CS_Tagjet1_phi;
    double merged_CS_Tagjet2_theta,merged_CS_Tagjet2_phi,merged_CS_deltaPhi_ZV,merged_CS_deltaR_ZV,merged_CS_deltaEta_ZV;
    double merged_CS_deltaPhi_ZTagjet1,merged_CS_deltaR_ZTagjet1,merged_CS_deltaEta_ZTagjet1,merged_CS_deltaPhi_ZTagjet2,merged_CS_deltaR_ZTagjet2;
    double merged_CS_deltaEta_ZTagjet2,merged_CS_deltaPhi_VTagjet1,merged_CS_deltaR_VTagjet1,merged_CS_deltaEta_VTagjet1,merged_CS_deltaPhi_VTagjet2;
    double merged_CS_deltaR_VTagjet2,merged_CS_deltaEta_VTagjet2,merged_CS_deltaPhi_Tagjet1Tagjet2,merged_CS_deltaR_Tagjet1Tagjet2,merged_CS_deltaEta_Tagjet1Tagjet2;
    double merged_Zrf_Lepton1_cos_theta,merged_Zrf_Lepton2_cos_theta,merged_Zrf_Lepton_minus_cos_theta,merged_Zrf_Lepton_minus_eta,merged_Zrf_Lepton1_pt;
    double merged_Zrf_Lepton2_pt,merged_Zrf_Lepton_minus_cos_theta_Zlab_axis,merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis,merged_Zrf_Lepton_deltaEta,merged_Zrf_Lepton_tanh_deltaEta,merged_Zrf_Lepton_deltaPhi,merged_Zrf_Lepton_deltaPhi_Vhad;
    double merged_tagjet1_pt,merged_tagjet2_pt,merged_tagjets_pt,merged_tagjets_delta_pt,merged_tagjets_delta_eta;
    double merged_tagjet1_eta,merged_tagjet2_eta,merged_tagjet1_phi,merged_tagjet2_phi,merged_tagjets_m;
    double merged_tagjets_eta,merged_tagjets_dy,merged_tagjets_dphi,merged_lepton_delta_pt,merged_lepton1_pt;
    double merged_lepton2_pt,merged_lepton1_eta,merged_lepton2_eta,merged_lepton_delta_eta,merged_ll_mass;
    double merged_ll_pt,merged_ll_eta,merged_ll_phi,merged_ll_DR_tagjet1,merged_ll_DR_tagjet2;
    double merged_ll_Dphi_tagjet1,merged_ll_Dphi_tagjet2,merged_ll_Deta_tagjet1,merged_ll_Deta_tagjet2;
    double merged_fjet_eta,merged_fjet_pt,merged_fjet_D2,merged_fjet_DR_lepton1;
    double merged_fjet_DR_lepton2,merged_fjet_DR_tagjet1,merged_fjet_DR_tagjet2,merged_fjet_Dphi_tagjet1,merged_fjet_Dphi_tagjet2;
    double merged_fjet_Deta_tagjet1,merged_fjet_Deta_tagjet2,merged_Vhad_DeltaR_Z,merged_Vhad_DeltaPhi_Z,merged_Vhad_DeltaEta_Z;
    double merged_fjet_mass,merged_VZ_mass,merged_VZ_pt,merged_VZ_eta,merged_Full_mass;
    double merged_Full_pt,merged_Centrality,merged_CentralityV,merged_CentralityZ,merged_CentralityZV,merged_ZeppZ,merged_ZeppV,merged_ZeppZV;
    int merged_n_jets,merged_n_lepton_stable,merged_fjet_n,merged_lepton1_pids,merged_lepton2_pids;
    int merged_Ntrk_tagjets1,merged_Ntrk_tagjets2,merged_Ntrk_tagjets,merged_Ntrk_fjets;
    /// @}
    std::map<std::string, double*> varMap = {
        {"merged_CS_Z_cos_theta", &merged_CS_Z_cos_theta},
        {"merged_CS_Z_eta", &merged_CS_Z_eta},
        {"merged_CS_Z_phi", &merged_CS_Z_phi},
        {"merged_CS_Z_pt", &merged_CS_Z_pt},
        {"merged_CS_lepton_pt1", &merged_CS_lepton_pt1},
        {"merged_CS_lepton_pt2", &merged_CS_lepton_pt2},
        {"merged_CS_lepton_deltaEta", &merged_CS_lepton_deltaEta},
        {"merged_CS_lepton_deltaPhi", &merged_CS_lepton_deltaPhi},
        {"merged_CS_lepton_deltaR_Vhad", &merged_CS_lepton_deltaR_Vhad},
        {"merged_CS_V_cos_theta", &merged_CS_V_cos_theta},
        {"merged_CS_V_eta", &merged_CS_V_eta},
        {"merged_CS_V_phi", &merged_CS_V_phi},
        {"merged_CS_V_pt", &merged_CS_V_pt},
        {"merged_CS_Tagjet1_theta", &merged_CS_Tagjet1_theta},
        {"merged_CS_Tagjet1_phi", &merged_CS_Tagjet1_phi},
        {"merged_CS_Tagjet2_theta", &merged_CS_Tagjet2_theta},
        {"merged_CS_Tagjet2_phi", &merged_CS_Tagjet2_phi},
        {"merged_CS_deltaPhi_ZV", &merged_CS_deltaPhi_ZV},
        {"merged_CS_deltaR_ZV", &merged_CS_deltaR_ZV},
        {"merged_CS_deltaEta_ZV", &merged_CS_deltaEta_ZV},
        {"merged_CS_deltaPhi_ZTagjet1", &merged_CS_deltaPhi_ZTagjet1},
        {"merged_CS_deltaR_ZTagjet1", &merged_CS_deltaR_ZTagjet1},
        {"merged_CS_deltaEta_ZTagjet1", &merged_CS_deltaEta_ZTagjet1},
        {"merged_CS_deltaPhi_ZTagjet2", &merged_CS_deltaPhi_ZTagjet2},
        {"merged_CS_deltaR_ZTagjet2", &merged_CS_deltaR_ZTagjet2},
        {"merged_CS_deltaEta_ZTagjet2", &merged_CS_deltaEta_ZTagjet2},
        {"merged_CS_deltaPhi_VTagjet1", &merged_CS_deltaPhi_VTagjet1},
        {"merged_CS_deltaR_VTagjet1", &merged_CS_deltaR_VTagjet1},
        {"merged_CS_deltaEta_VTagjet1", &merged_CS_deltaEta_VTagjet1},
        {"merged_CS_deltaPhi_VTagjet2", &merged_CS_deltaPhi_VTagjet2},
        {"merged_CS_deltaR_VTagjet2", &merged_CS_deltaR_VTagjet2},
        {"merged_CS_deltaEta_VTagjet2", &merged_CS_deltaEta_VTagjet2},
        {"merged_CS_deltaPhi_Tagjet1Tagjet2", &merged_CS_deltaPhi_Tagjet1Tagjet2},
        {"merged_CS_deltaR_Tagjet1Tagjet2", &merged_CS_deltaR_Tagjet1Tagjet2},
        {"merged_CS_deltaEta_Tagjet1Tagjet2", &merged_CS_deltaEta_Tagjet1Tagjet2},
        {"merged_Zrf_Lepton1_cos_theta", &merged_Zrf_Lepton1_cos_theta},
        {"merged_Zrf_Lepton2_cos_theta", &merged_Zrf_Lepton2_cos_theta},
        {"merged_Zrf_Lepton_minus_cos_theta", &merged_Zrf_Lepton_minus_cos_theta},
        {"merged_Zrf_Lepton_minus_eta", &merged_Zrf_Lepton_minus_eta},
        {"merged_Zrf_Lepton1_pt", &merged_Zrf_Lepton1_pt},
        {"merged_Zrf_Lepton2_pt", &merged_Zrf_Lepton2_pt},
        {"merged_Zrf_Lepton_minus_cos_theta_Zlab_axis", &merged_Zrf_Lepton_minus_cos_theta_Zlab_axis},
        {"merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis", &merged_Zrf_Lepton_minus_cos_theta_VZrf_Zlab_axis},
        {"merged_Zrf_Lepton_deltaEta", &merged_Zrf_Lepton_deltaEta},
        {"merged_Zrf_Lepton_tanh_deltaEta", &merged_Zrf_Lepton_tanh_deltaEta},
        {"merged_Zrf_Lepton_deltaPhi", &merged_Zrf_Lepton_deltaPhi},
        {"merged_Zrf_Lepton_deltaPhi_Vhad", &merged_Zrf_Lepton_deltaPhi_Vhad},
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
        {"merged_ll_mass", &merged_ll_mass},
        {"merged_ll_pt", &merged_ll_pt},
        {"merged_ll_eta", &merged_ll_eta},
        {"merged_ll_phi", &merged_ll_phi},
        {"merged_ll_DR_tagjet1", &merged_ll_DR_tagjet1},
        {"merged_ll_DR_tagjet2", &merged_ll_DR_tagjet2},
        {"merged_ll_Dphi_tagjet1", &merged_ll_Dphi_tagjet1},
        {"merged_ll_Dphi_tagjet2", &merged_ll_Dphi_tagjet2},
        {"merged_ll_Deta_tagjet1", &merged_ll_Deta_tagjet1},
        {"merged_ll_Deta_tagjet2", &merged_ll_Deta_tagjet2},
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
        {"merged_Vhad_DeltaR_Z", &merged_Vhad_DeltaR_Z},
        {"merged_Vhad_DeltaPhi_Z", &merged_Vhad_DeltaPhi_Z},
        {"merged_Vhad_DeltaEta_Z", &merged_Vhad_DeltaEta_Z},
        {"merged_fjet_mass", &merged_fjet_mass},
        {"merged_VZ_mass", &merged_VZ_mass},
        {"merged_VZ_pt", &merged_VZ_pt},
        {"merged_VZ_eta", &merged_VZ_eta},
        {"merged_Full_mass", &merged_Full_mass},
        {"merged_Full_pt", &merged_Full_pt},
        {"merged_Centrality", &merged_Centrality},
        {"merged_CentralityV", &merged_CentralityV},
        {"merged_CentralityZ", &merged_CentralityZ},
        {"merged_CentralityZV", &merged_CentralityZV},
        {"merged_ZeppZ", &merged_ZeppZ},
        {"merged_ZeppV", &merged_ZeppV},
        {"merged_ZeppZV", &merged_ZeppZV}
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

    unique_ptr<TFile> _tf_truth;
    unique_ptr<TTree> _tt_truth;    
    unique_ptr<TTree> _tt_truth_merged;    
    unique_ptr<TTree> _tt_truth_merged_befcutV;    
    unique_ptr<TTree> _tt_truth_resolved; 

    double minDR_VBSq1_jets, minDR_VBSq2_jets,Delta_R_tagging_quark,Delta_Eta_tagging_quark;
    double Tagging_Quark1_pT, Tagging_Quark1_eta, Tagging_Quark1_mass;
    double Tagging_Quark2_pT, Tagging_Quark2_eta, Tagging_Quark2_mass;
    double Tagging_Quarks_mass;

    double W_pT,W_Quark1_pT, W_Quark1_eta, W_Quark1_mass;
    double W_Quark2_pT, W_Quark2_eta, W_Quark2_mass;
    double W_Quarks_delta_Eta;

    double Tagging_Jet1_pT, Tagging_Jet1_eta, Tagging_Jet1_mass;
    double Tagging_Jet2_pT, Tagging_Jet2_eta, Tagging_Jet2_mass;
    double Tagging_Jets_mass,Tagging_Jets_delta_Eta;

    std::vector<double> Min_DR_Wq_jets, Min_DR_q_Tagjets,Delta_Eta_TaggingQuark1_Wquarks,Delta_Eta_TaggingQuark2_Wquarks;

    int misID_tag_jets, FromTruth_W_q_tag_jets;
    int merged_fjet_FromTruth_W_q, merged_fjetmisID_FromTruth_Tagjet_q, merged_misID_tag_jets;
    int resolved_sjets_FromTruth_W_q, resolved_sjet1_FromTruth_W_q, resolved_sjet2_FromTruth_W_q, resolved_sjet1misID_FromTruth_Tagjet_q,resolved_sjet2misID_FromTruth_Tagjet_q;

    std::vector<int> FromTruth_W_q_misIDtagjets;   

    std::map<std::string, double*> varMapDouble_truth = {
        {"minDR_VBSq1_jets", &minDR_VBSq1_jets},
        {"minDR_VBSq2_jets", &minDR_VBSq2_jets},
        {"Delta_Eta_tagging_quark", &Delta_Eta_tagging_quark},
        {"Tagging_Quark1_pT", &Tagging_Quark1_pT},
        {"Tagging_Quark1_eta", &Tagging_Quark1_eta},
        {"Tagging_Quark1_mass", &Tagging_Quark1_mass},
        {"Tagging_Quark2_pT", &Tagging_Quark2_pT},
        {"Tagging_Quark2_eta", &Tagging_Quark2_eta},
        {"Tagging_Quark2_mass", &Tagging_Quark2_mass},
        {"Tagging_Quarks_mass", &Tagging_Quarks_mass},
        {"Tagging_Jet1_pT", &Tagging_Jet1_pT},
        {"Tagging_Jet1_eta", &Tagging_Jet1_eta},
        {"Tagging_Jet1_mass", &Tagging_Jet1_mass},
        {"Tagging_Jet2_pT", &Tagging_Jet2_pT},
        {"Tagging_Jet2_eta", &Tagging_Jet2_eta},
        {"Tagging_Jet2_mass", &Tagging_Jet2_mass},
        {"Tagging_Jets_mass", &Tagging_Jets_mass},
        {"Tagging_Jets_delta_Eta", &Tagging_Jets_delta_Eta},
        {"W_Quark1_pT", &W_Quark1_pT},
        {"W_pT", &W_pT},
        {"W_Quark1_eta", &W_Quark1_eta},
        {"W_Quark1_mass", &W_Quark1_mass},
        {"W_Quark2_pT", &W_Quark2_pT},
        {"W_Quark2_eta", &W_Quark2_eta},
        {"W_Quark2_mass", &W_Quark2_mass},
        {"W_Quarks_delta_Eta", &W_Quarks_delta_Eta}
    };

    std::map<std::string, std::vector<double>*> varMapVectorDouble_truth = {
        {"Min_DR_Wq_jets", &Min_DR_Wq_jets},
        {"Min_DR_q_Tagjets", &Min_DR_q_Tagjets},
        {"Delta_Eta_TaggingQuark1_Wquarks", &Delta_Eta_TaggingQuark1_Wquarks},
        {"Delta_Eta_TaggingQuark2_Wquarks", &Delta_Eta_TaggingQuark2_Wquarks},
    };

    std::map<std::string, int*> varMapInt_truth_resolved = {
        {"resolved_sjets_FromTruth_W_q", &resolved_sjets_FromTruth_W_q},
        {"resolved_sjet1_FromTruth_W_q", &resolved_sjet1_FromTruth_W_q},
        {"resolved_sjet2_FromTruth_W_q", &resolved_sjet2_FromTruth_W_q},
        {"resolved_sjet1misID_FromTruth_Tagjet_q", &resolved_sjet1misID_FromTruth_Tagjet_q},
        {"resolved_sjet2misID_FromTruth_Tagjet_q", &resolved_sjet2misID_FromTruth_Tagjet_q},
    };

    std::map<std::string, std::vector<int>*> varMapVectorInt = {
        {"FromTruth_W_q_misIDtagjets", &FromTruth_W_q_misIDtagjets},
    };

    };


    RIVET_DECLARE_PLUGIN(WmZ_llqq);

}