#ifndef TRUTH_MATCHING_H
#define TRUTH_MATCHING_H

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
#include <vector>
#include <limits>
#include <algorithm>

// Include other necessary headers here

namespace TruthMatching {

    const bool DescendantsWithParticle(const Particle& p, const Particle& q);
    const bool AncestorWithParticle(const Particle& p, const Particle& q);
    double MatchingJet(const Jet& jet, const Particle& particle);
    std::vector<Particle> ValidQuark_W(const Particle& W_boson);
    std::vector<Particle> Tagging_quarks(const Particles& all_particles);
    bool Check_VBS_event(const std::vector<Particle>& all_particles);
    const Particle GetWboson(const std::vector<Particle>& all_particles);
    std::vector<double> Truth_q_minDR_jets(const std::vector<Particle>& all_particles, const std::vector<Jet>& jets, bool merged = false);
    int IsTruthTagJet(const std::vector<Particle>& all_particles, const Jet& tag_jet, double dR_cut = 0.2);
    int NbTagJetMisID(const std::vector<Particle>& all_particles, const std::vector<Jet>& tag_jets, double dR_cut = 0.2);
    int IsTrueWboson(const std::vector<Particle>& all_particles, const Jet& jet, double dR_cut = 0.2, bool merged = false);

} // namespace TruthMatching

#endif // TRUTH_MATCHING_H