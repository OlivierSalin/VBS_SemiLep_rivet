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
    class WpZ_llqq_pid : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(WpZ_llqq_pid);

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

    const Particle GetWParent(const Particle& p) 
    {
        const Particles parents = p.parents();
        if (parents.size() == 0) return p;
        for (const Particle& p2 : parents) {
            if (abs(p2.pid()) == 24) return p2;  // Return if the parent is a W boson
            else return GetWParent(p2);  // Recursively find the parent
        }
        return p;
    }

    const bool GetWCandidate(const Jet& jet) 
    {
        Particles particles = jet.constituents();
        // Sort the particles by their pT in descending order
        std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) {
            return a.pt() > b.pt();
        });

        int WCount = 0;
        double totalWeight = 0.0;
        for (const Particle& p : particles) {
            if (p.pt() > 0.5) {  // Check if the particle's pT is greater than 0.5 GeV
                const Particle& WParent = GetWParent(p);
                if (abs(WParent.pid()) == 24) {  // Check if the parent is a W boson
                    WCount++;
                    totalWeight += p.pt();  // Add the particle's pT to the total weight
                }
            }
        }
        if (static_cast<double>(WCount) > particles.size() / 2 && totalWeight > particles.size() / 2) {  // Check if more than 50% of the particles have a W parent and the total weight is greater than 50%
            return true;  // Return true if the jet is a W candidate
        }
        else {
            return false;  // Return false if the jet is not a W candidate
        }
    }

/*     std::vector<int> GetAncestorPIDs(const Particle& p) 
    {
        std::vector<int> ancestorPIDs;
        const Particles ancestors = p.ancestors(ParticleSelector(), false);
        for (const Particle& ancestor : ancestors) {
            ancestorPIDs.push_back(ancestor.pid());
        }
        return ancestorPIDs;
    } */

    std::vector<int> GetAncestorPIDs(const Particle& p) 
    {
        std::vector<int> ancestorPIDs;
        const Particles parents = p.parents();
        //printf("Inside Parent size: %d\n", parents.size());       
        if (parents.empty()) {
            return ancestorPIDs;  // Base case: if the particle has no parents, return an empty vector
        }
        for (const Particle& parent : parents) {
            ancestorPIDs.push_back(parent.pid());  // Add the PID of the parent to the vector
            std::vector<int> grandparentPIDs = GetAncestorPIDs(parent);  // Recursive PIDs of the grandparents
            ancestorPIDs.insert(ancestorPIDs.end(), grandparentPIDs.begin(), grandparentPIDs.end());  // Add the PIDs of the grandparents
        }
        return ancestorPIDs;
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

    std::vector<Particle> ValidQuark_W(const Particle& W_boson){
        std::vector<Particle> valid_descendants;
        for(const Particle& descendant : W_boson.allDescendants()){
            int descendant_status = descendant.genParticle()->status();
            if(descendant.abspid() < 7 && descendant_status > 30){
                printf("\nDescendant PID: %d and Descendant Status: %d\n", descendant.pid(), descendant.genParticle()->status());
                for(const Particle& ancestor : descendant.ancestors(Cuts::OPEN, false)){
                    printf("Ancestor PID: %d and Ancestor Status: %d\n", ancestor.pid(), ancestor.genParticle()->status());
                }
                valid_descendants.push_back(descendant);
            }
        }
        if(valid_descendants.size() > 2){
            for(int i = 0; i < valid_descendants.size(); i++){
                for(int j = i + 1; j < valid_descendants.size(); j++){
                    if(valid_descendants[i].pid() == -valid_descendants[j].pid()){
                        bool found = false;
                        for(const Particle& ancestor_i : valid_descendants[i].ancestors(Cuts::OPEN, false)){
                            //printf("Ancestor PID: %d and Ancestor Status: %d\n", ancestor_i.pid(), ancestor_i.genParticle()->status());
                            if(ancestor_i.pid()==21 && ancestor_i.genParticle()->status() == 51){
                                if(AncestorWithParticle(valid_descendants[j],ancestor_i)){
                                    found = true;
                                    //break;
                                }
                            }
                        }
                        if(found){
                            printf("\nNon valid Quark PID: %d, Status: %d\n", valid_descendants[i].pid(), valid_descendants[i].genParticle()->status());
                            printf("Non valid Quark bis PID: %d, Status: %d\n", valid_descendants[j].pid(), valid_descendants[j].genParticle()->status());
                            valid_descendants.erase(std::remove_if(valid_descendants.begin(), valid_descendants.end(), [&](const Particle& p){
                                return p.genParticle() == valid_descendants[i].genParticle() || p.genParticle() == valid_descendants[j].genParticle();
                            }), valid_descendants.end());
                        }
                    } 
                }  
            }
        }

        return valid_descendants;
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


        // Merged histograms
        
        // plots common with others
        std::ifstream jet_hist_merged_file(txt_dir + "/jet_hists_merged.json");      
        json jet_hist_merged = json::parse(jet_hist_merged_file);
        for (json::iterator it = jet_hist_merged.begin(); it != jet_hist_merged.end(); ++it) {
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
  
        // plots that are not in other ana



        //counter for efficiency
        book(_c["pos_w_initial"],"pos_w_initial");
        book(_c["pos_w_final"],"pos_w_final");
        book(_c["neg_w_initial"],"neg_w_initial");
        book(_c["neg_w_final"],"neg_w_final");

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
        if (nlep !=_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext();

        const Particle& lep1 = leptons[0];
        const Particle& lep2 = leptons[1];

        const Particle parent_lep1 = GetParent(lep1);
        const Particle parent_lep2 = GetParent(lep2);
        printf("Parent of lep1 PID: %d\n", parent_lep1.pid());

        const Particles ancestor= lep1.ancestors();
        for (const Particle& p : ancestor) {
            std::cout << "Ancestor PID: " << p.pid() << std::endl;
        }
        std::vector<int> ancestor_bis= GetAncestorPIDs(lep1);

        std::cout << "Contents of ancestor_bis: " << std::endl;
        for (int pid : ancestor_bis) {
            std::cout << pid << std::endl;
        }

        const Particles ancestor= lep1.ancestors(Cuts::OPEN, false);
        if (!ancestor.empty()) {
            const Particle parton = ancestor.back();

            printf("\nParton PID: %d\n", parton.pid());
            

            const Particles descendants = parton.allDescendants();
            const Particles children_parton = parton.children();
            int counter = 0;

            // Print the size of the container
            printf("Number of children partons: %d\n", children_parton.size());

            // Check if the container is empty
            if (!children_parton.empty()) {
                // Print the type of the first element
                printf("Type of children parton: %s\n", typeid(children_parton[0]).name());

                // Print the properties of each particle
                for (const Particle& p : children_parton) {
                    printf("Particle ID: %d, pt: %f, eta: %f, phi: %f\n", p.pid(), p.pt(), p.eta(), p.phi());
                }
            } else {
                printf("children_parton is empty.\n");
            }
/*             printf("\nChildren parton PID: %d\n", children_parton[0].pid());
            int counter1 = 0;
            for (const Particle& p : children_parton[0].children()) {
                printf("Children num %d PID : %d\n", counter1,p.pid());
                counter1++;
            } */

/*             for (const Particle& p : descendants) {

                printf("Descendant num %d PID : %d\n", counter,p.pid());
                counter++;
            }
            printf("Total number of descendants: %d\n\n", counter); */
        }

        // Use of HepMC to get the parent of the lepton
  
        printf("Lep1 PID: %d\n", lep1.pid());
        printf("Parent of lep1 PID: %d\n", parent_lep1.pid());
        printf("Ancestor size: %d\n\n", ancestor.size());
        for (const Particle& p : ancestor) {
            std::cout << "Ancestor PID: " << p.pid() << std::endl;
        }
        printf("PID and Status of Lepton from FS\n");
        printf("Lepton PID: %d Status lepton %d\n",lep1.pid(),lep1.genParticle()->status());
        printf("Parent Lepton PID : %d boson %d\n\n\n",parent_lep1.pid(),parent_lep1.genParticle()->status());

        printf("Status Particles (only 22 or 23) \n");
        const Particles all_particles = event.allParticles();
        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle();
            int status = p_->status();
            if (abs(status) == 23 or abs(status) == 22 or abs(status) == 21){
                
                printf("Particle PID: %d, Status: %d\n", p.pid(), status);
            }
        }

        // More details about the descendant of event particles
/*         printf("Status Particles (only 22 or 23) \n");
        const Particles all_particles = event.allParticles();
        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle();
            int status = p_->status();
            if (abs(status) == 23 or abs(status) == 22 or abs(status) == 21){
                
                printf("\nParticle PID: %d, Status: %d, Eta: %f, Pt: %f\n", p.pid(), status, p.eta(), p.mom().pt());
                for(const Particle& child : p.children()){
                    ConstGenParticlePtr child_ = child.genParticle();
                    printf("Child PID: %d, Status: %d, Eta: %f, Pt: %f\n", child.pid(), child_->status(), child.eta(), child.mom().pt());
                    for(const Particle& grandchild : child.children()){
                        ConstGenParticlePtr grandchild_ = grandchild.genParticle();
                        printf("Grandchild PID: %d, Status: %d, Eta: %f, Pt: %f\n", grandchild.pid(), grandchild_->status(), grandchild.eta(), grandchild.mom().pt());
                        for(const Particle& great_grandchild : grandchild.children()){
                            ConstGenParticlePtr great_grandchild_ = great_grandchild.genParticle();
                            printf("Great Grandchild PID: %d, Status: %d, Eta: %f, Pt: %f\n", great_grandchild.pid(), great_grandchild_->status(), great_grandchild.eta(), great_grandchild.mom().pt());
                            for(const Particle& great_great_grandchild : great_grandchild.children()){
                                ConstGenParticlePtr great_great_grandchild_ = great_great_grandchild.genParticle();
                                printf("Great Great Grandchild PID: %d, Status: %d, Eta: %f, Pt: %f\n", great_great_grandchild.pid(), great_great_grandchild_->status(), great_great_grandchild.eta(), great_great_grandchild.mom().pt());
                            }}                    
                    }
                }
            }
        } */

        // Loop over all particles
        for(const Particle& p : all_particles){
            ConstGenParticlePtr p_ = p.genParticle(); // Get the underlying GenParticle
            int status = p_->status(); // Get the status of the particle

            // Check if the particle is a Z boson with a status of 23, 22, or 21
            if(p.pid() == 23 && (abs(status) == 23 || abs(status) == 22 || abs(status) == 21)){
                printf("Boson PID: %d, Status: %d\n", p.pid(), status);

                // Loop over all descendants of the Z boson
                for(const Particle& p_des : p.allDescendants()){
                    ConstGenParticlePtr p_des_ = p_des.genParticle(); // Get the underlying GenParticle of the descendant
                    printf("Descendant PID: %d, Status: %d\n", p_des.pid(), p_des_->status());

                    // Check if the PID of the descendant is the same as the PID of lep2
                    if(p_des.pid() == lep1.pid()){
                        printf("Found same lepton if :%d\n", isSame(p_des, lep1));
                        printf("Gen is same lepton if :%d\n", p_des.genParticle() == lep1.genParticle());
                        printf("Truth lepton 4-vector momentum: (E: %f, px: %f, py: %f, pz: %f)\n", p_des.mom().E(), p_des.mom().px(), p_des.mom().py(), p_des.mom().pz());

                        // Print the 4-vector momentum of lep1
                        printf("Rivet FS lepton 4-vector momentum: (E: %f, px: %f, py: %f, pz: %f)\n", lep1.mom().E(), lep1.mom().px(), lep1.mom().py(), lep1.mom().pz());

                    }
                }
            }
        } 
        //std::vector<int> ancestor_bis= GetAncestorPIDs(lep1);

/*         std::cout << "Contents of ancestor_bis: " << std::endl;
        for (int pid : ancestor_bis) {
            std::cout << pid << std::endl;
        } */


        
        // opposite charge and same flavour
/*         if (lep1.pid()+lep2.pid()!=0) vetoEvent; 
        _cutflows_merged.fillnext();
        _cutflows_resolved.fillnext(); */


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
        
        // Test the parents /anncestor pid of the particles in the signal jets and chek if come from W boson

/*         if (sjets_sig.size() > 0) {
            Jet sjets_signa = tag_jets[0];
            std::cout << "\n Tagging Jet pT: " << sjets_signa.mom().pT()/GeV << std::endl;
            std::cout << "\nJet nb tracks: " << CountChargedTracks(sjets_signa) << std::endl;
            int counter = 0;
            Particles sjets_signa_particles = sjets_signa.particles();
            std::sort(sjets_signa_particles.begin(), sjets_signa_particles.end(), [](Particle const &a, Particle const &b) {
                return a.pT() > b.pT(); // biggest pT will be first in array
                });
            printf("Tagging Jet 1 W candidate: %d",GetWCandidate(sjets_signa));
            for (const Particle& jet_p : sjets_signa_particles) {
                std::cout << "Jet particle pT: " << jet_p.mom().pT()/GeV << std::endl;
                std::cout << "Jet particle charge: " << jet_p.isCharged() << std::endl;
                printf("Jet particle nb %d PID: %d\n\n", counter, jet_p.pid());
                printf("Jet particle nb %d - Parent size PID: %d\n", counter, jet_p.parents().size());
                std::vector<int> ancestorPIDs = GetAncestorPIDs(jet_p);

                std::cout << "Ancestor PIDs: [";
                for (size_t i = 0; i < ancestorPIDs.size() - 1; ++i) {
                    if (ancestorPIDs[i] != ancestorPIDs[i + 1]) {
                        std::cout << ancestorPIDs[i];
                        if (i != ancestorPIDs.size() - 2) {
                            std::cout << ", ";
                        }
                    }
                }
                // Print the last PID if it's different from the second last
                if (ancestorPIDs.size() > 1 && ancestorPIDs[ancestorPIDs.size() - 1] != ancestorPIDs[ancestorPIDs.size() - 2]) {
                    std::cout << ", " << ancestorPIDs[ancestorPIDs.size() - 1];
                }
                std::cout << "]" << std::endl; 
                printf("\nJet particle nb %d - W candidate ancestor PID: %d\n\n\n", counter, GetWParent(jet_p).pid());
                counter++;
        }
        } */

        // Check if we are in the Merged signal region
        if (n_fjets > 0) {
            _cutflows_merged.fillnext();
            const FourMomentum fourvec_fjets = fjets[0].mom();
            if ((_docut==1 && (fourvec_fjets.mass() >= _jcuts["m_fjet_WZ"][0] && fourvec_fjets.mass() <= _jcuts["m_fjet_WZ"][1]))||_docut==0) {
                _cutflows_merged.fillnext();
                
                // Total cutflow of the merged region
                _cutflows_merged.fillnext();

                // Four vector of the QGC system in resolved region llJ
                const FourMomentum fourvec_fjets_ll = lep1.mom() + lep2.mom() + fjets[0].mom();

                // Four vector of the Full system in resolved region 
                const FourMomentum fourvec_fjets_full = lep1.mom() + lep2.mom() + fjets[0].mom() + tag1_jet + tag2_jet;

                // We are in the Merged signal region
                // Fill in the histogram of the merged region
                _h["merged_n_jets"]->fill(n_jets);
/*                 _h["merged_pt_tagjet1"]->fill(tag1_jet.pt());
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
                _h["merged_pt_WZ"]->fill(fourvec_fjets_ll.pt());                
                _h["merged_mass_Full"]->fill(fourvec_fjets_full.mass());                
                _h["merged_pt_Full"]->fill(fourvec_fjets_full.pt());    */             

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
            double signal_mjj = (signal_jet1 + signal_jet2).mass();
            //printf("signal_mjj: %f\n", signal_mjj);
            if ((_docut==1 && (signal_mjj < _jcuts["m_fjet_WZ"][0] || signal_mjj > _jcuts["m_fjet_WZ"][1]) )) {
                vetoEvent;
            }
            _cutflows_resolved.fillnext();

            const FourMomentum fourvec_signal_jets_ll = lep1.mom() + lep2.mom() + signal_jet1 + signal_jet2;

            // Four vector of the Full system in resolved region 
            const FourMomentum fourvec_signal_jets_full = lep1.mom() + lep2.mom() + signal_jet1 + signal_jet2 + tag1_jet + tag2_jet;
            // More than 2 signal jets candidates
            FourMomentum fourvec_signal_jjj;
/*             if (sjets_sig.size() > 2){
                const FourMomentum signal_jet3 = sjets_sig[2].mom();
                fourvec_signal_jjj = signal_jet1 + signal_jet2 + signal_jet3;
                double signal_mjjj = fourvec_signal_jjj.mass();
                // if ((_docut==1 && (signal_mjjj < _jcuts["m_jjj"] ))|| _docut==0) {vetoEvent;}
                _cutflows_resolved.fillnext();            
            } */

            printf("\n New event: \n");
            int WCount = 0;
            std::cout << "\n Tagging 1 Jet pT: " << tag_jets[0].mom().pT()/GeV << std::endl;
            std::cout << "Number of tracks in tagging jet 2: " << CountChargedTracks(tag_jets[1]) << std::endl;
            printf("Tagging Jet 1 W candidate: %d\n",GetWCandidate(tag_jets[0]));
            if (GetWCandidate(tag_jets[0])) {
                WCount++;
            }

            std::cout << "\n Tagging 2 Jet pT: " << tag_jets[1].mom().pT()/GeV << std::endl;
            std::cout << "Number of tracks in tagging jet 2: " << CountChargedTracks(tag_jets[1]) << std::endl;
            std::cout << "Tagging 2 W candidate: " << GetWCandidate(tag_jets[1]) << std::endl;
            if (GetWCandidate(tag_jets[1])) {
                WCount++;
            }

            for (int i = 0; i < sjets_sig.size(); i++) {
                std::cout << "\n Signal jet " << i+1 << " pT: " << sjets_sig[i].mom().pT()/GeV << std::endl;
                std::cout << "Number of tracks in Signal jet " << i+1 << ": " << CountChargedTracks(sjets_sig[i]) << std::endl;
                std::cout << "Signal jet " << i+1 << " W candidate: " << GetWCandidate(sjets_sig[i]) << std::endl;
                if (GetWCandidate(sjets_sig[i])) {
                    WCount++;
                }
            }
            std::cout << "Number of W candidates: " << WCount << std::endl;
            std::cout << "\n Invariant mass of the signal jets: " << signal_mjj << std::endl;



            _cutflows_resolved.fillnext();
            // Fill in the histograms for the resolved region
            //jet plots
            _h["resolved_n_jets"]->fill(n_jets);
            _h["resolved_pt_tagjet1"]->fill(tag1_jet.pt());
            _h["resolved_pt_tagjet2"]->fill(tag2_jet.pt());
/*             _h["resolved_eta_tagjets"]->fill(tag1_jet.eta()); _h["resolved_eta_tagjets"]->fill(tag2_jet.eta());
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
            _h["resolved_pt_Full"]->fill(fourvec_signal_jets_full.pt());   */  
        
            // More than 2 signal jets candidates
            if (sjets_sig.size() > 2){   
                _h["resolved_mjjj"]->fill(fourvec_signal_jjj.mass());
            }

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


    RIVET_DECLARE_PLUGIN(WpZ_llqq_pid);

}