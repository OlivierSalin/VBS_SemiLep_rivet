#include "Truth_matching.h"

namespace TruthMatching {



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

    int IsTruthTagJet(const std::vector<Particle>& all_particles, const Jet& tag_jet, double dR_cut = 0.2){
        std::vector<Particle> taggingQuarks = Tagging_quarks(all_particles);
        if(taggingQuarks.size() != 2) return -1; // Return -1 if the size of taggingQuarks is not 2

        for(const Particle& taggingQuark : taggingQuarks){
            const double dR = deltaR(taggingQuark, tag_jet);
            if(dR < dR_cut) return 1; // Return 1 if the jet matches a tagging quark
        }
        
        return 0; // Return 0 if the jet does not match any tagging quark
    }

    int NbTagJetMisID(const std::vector<Particle>& all_particles, const std::vector<Jet>& tag_jets, double dR_cut = 0.2){
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
     

    int IsTrueWboson(const std::vector<Particle>& all_particles, const Jet& jet, double dR_cut = 0.2, bool merged = false){
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

} // namespace TruthMatching