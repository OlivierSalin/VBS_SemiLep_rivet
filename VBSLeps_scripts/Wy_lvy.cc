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
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Tools/Cutflow.hh"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

namespace Rivet {


  /// @brief Add a short analysis description here
  class Wy_lvy : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(Wy_lvy);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      _cut_mode = getOption("cut");

      std::cout << "++++++received cut_mode" << _cut_mode << "will do them accordingly \n";

      std::ifstream json_file("Cuts/Wy_lvy_cuts.json");
      _jcuts = json::parse(json_file);
      std::cout << "++++++ to check json 1 var got photon pt min" << _jcuts["pt_photon"] << "\n";
      _electron_eta_cut = (Cuts::absetaIn(_jcuts["eta_lepton_electron"][0][0], _jcuts["eta_lepton_electron"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["eta_lepton_electron"][1][0], _jcuts["eta_lepton_electron"][1][1]));
      _photon_eta_cut = (Cuts::absetaIn(_jcuts["eta_photon"][0][0], _jcuts["eta_photon"][0][1])) || 
                                (Cuts::absetaIn(_jcuts["eta_photon"][1][0], _jcuts["eta_photon"][1][1]));
      _muon_eta_cut = Cuts::absetaIn(0.0, _jcuts["eta_lepton_muon"]);
      _lepton_pt_cut = Cuts::pT > dbl(_jcuts["pt_lepton"])*GeV; 

      // The basic final-state projection:
      // all final-state particles within
      const FinalState fs;

      // photons as separate particles for final state, not dressing
      DirectFinalState photons(Cuts::abspid == PID::PHOTON);
      declare(photons, "photons");

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

      // FS excluding the leading big pt photons, muons and neutrinos to calculate cone energy of bit pt photons
      // so basically electrons + jets left
      DirectFinalState bare_muon_for_ph_iso(Cuts::abspid == PID::MUON);
      DressedLeptons dressed_muons_for_ph_iso(photons_for_dressing, bare_muon_for_ph_iso, 0.1, _muon_eta_cut &&_lepton_pt_cut);
      declare(dressed_muons_for_ph_iso, "muons_for_ph_iso");
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(photons);
      vfs.addVetoOnThisFinalState(dressed_muons_for_ph_iso);
      vfs.addVetoOnThisFinalState(InvisibleFinalState());
      declare(vfs, "ElAndJetsForPhotonIsoCalc");

      declare(MissingMomentum(), "METFinder");


      // Book things and save names to normalize later
      
      // plots common with others
      std::ifstream jet_hist_file( "jet_hists.json");
      json jet_hist = json::parse(jet_hist_file);
      for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      std::ifstream y_hist_file( "photon_hists.json");
      json y_hist = json::parse(y_hist_file);
      for (json::iterator it = y_hist.begin(); it != y_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      std::ifstream lep_hist_file( "lepton_hists.json");
      json lep_hist = json::parse(lep_hist_file);
      for (json::iterator it = lep_hist.begin(); it != lep_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      // plots that are not in other ana
      std::ifstream ana_hist_file( "Wy_lvy_hists.json");
      json ana_hist = json::parse(ana_hist_file);
      for (json::iterator it = ana_hist.begin(); it != ana_hist.end(); ++it) {
        std::vector<double> h_bins = it.value();
        if (h_bins.size()==3) {
          book(_h[it.key()], it.key(), h_bins[0], h_bins[1], h_bins[2]);
        }
        else{
          book(_h[it.key()], it.key(), h_bins);
        }
        _hist_names.push_back(it.key());
      }

      //counter for efficiency
      book(_c["pos_w_initial"],"pos_w_initial");
      book(_c["neg_w_initial"],"neg_w_initial");
      for (std::string& i_clip : _clips){
        std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
        std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
        book(_c[i_pos_c_name], i_pos_c_name);
        book(_c[i_neg_c_name], i_neg_c_name);
      }

      // Cut-flows
      _cutflows.addCutflow("Wy_lvy_selections", {"one_lep","pt_MET","m_T_W","have_iso_photons_ok_pt_eta",
                          "m_ly_away_from_m_z","dR_lepton_photon","n_jets","no_b_jets","tagjets_pt_dR_j1j2", "dphi_MET_tagjets",
                          "m_tagjets","dy_tagjets","centrality_jjly", "n_gap_jets","two_HS_boson"});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // save weights before cuts
      double ev_nominal_weight =  event.weights()[0];
      if (ev_nominal_weight>=0){_c["pos_w_initial"]->fill();} // dont need anything in bracket as this will be weight on weight
      else {_c["neg_w_initial"]->fill();}

      _cutflows.fillinit();

      // Retrieve dressed leptons, sorted by pT
      Particles e_stable;
      Particles mu_stable;
      
      e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton_pt_cut);
      mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton_pt_cut);
      
      Particles leptons_stable = e_stable + mu_stable; 

      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=_jcuts["n_lepton_stable"])  vetoEvent; 
      _cutflows.fillnext();
      
      const FourMomentum lep = leptons_stable[0].mom();

      const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
      const double pt_MET = METfinder.missingPt()/GeV;
      if (pt_MET<_jcuts["pt_MET"]) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum fourvec_MET = METfinder.missingMomentum();
      const double d_phi_MET_lep = deltaPhi(fourvec_MET.phi(), lep.phi());
      const double m_T_W = sqrt( 2*lep.pT()*pt_MET*(1 - cos(d_phi_MET_lep)) );
      if (m_T_W<_jcuts["m_T_W"]) vetoEvent;
      _cutflows.fillnext();

      //photons
      Particles photons = apply<FinalState>(event, "photons").particlesByPt(_photon_eta_cut && Cuts::pT > 12*GeV);
      if (photons.empty())  vetoEvent;
      //photon cone calculation and photon OR with leptons 
      Particles isolated_photons;
      std::vector<double> cone_to_photon_fracs = {};
      Particles cone_sum_particles = apply<VetoedFinalState>(event, "ElAndJetsForPhotonIsoCalc").particles();
      for (const Particle &i_high_pt_photon : photons) {
        // check photon isolation
        double i_high_pt_photon_cone_E = 0.0;
        for (const Particle &i_p_cone : cone_sum_particles) {
          if (isSame(i_high_pt_photon,i_p_cone) || !isStable(i_p_cone)){
            continue;
          }
          if (deltaR(i_high_pt_photon, i_p_cone) < 0.4) { // etcone40
            i_high_pt_photon_cone_E += i_p_cone.Et();
          }
        }
        double i_cone_to_photon_frac = i_high_pt_photon_cone_E / i_high_pt_photon.pT(); 
        if (i_cone_to_photon_frac > dbl(_jcuts["cone_frac_photon"]))  continue;
        if (any(leptons_stable, deltaRLess(i_high_pt_photon, 0.4))) continue;
        isolated_photons += i_high_pt_photon;
        cone_to_photon_fracs.push_back(i_cone_to_photon_frac);

      }
      if (isolated_photons.size() < _jcuts["n_photons_iso"])  vetoEvent;
      const Particle& lead_iso_photon = isolated_photons[0];
      if (lead_iso_photon.pT() < dbl(_jcuts["pt_photon"])*GeV) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum fourvec_ly = lep + lead_iso_photon.mom();
      const double m_ly = fourvec_ly.mass();
      if (fabs(m_ly-91.18)<dbl(_jcuts["abs_diff_m_ly_m_z"])*GeV)  vetoEvent;
      _cutflows.fillnext();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 4.4 && Cuts::pT > 25*GeV); // then will do cut on two leading pt>50
      Jets btagging_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.5 && Cuts::pT > 20*GeV);
      // Remove all jets within certain dR of a dressed lepton
      double dR_cut = _jcuts["dR_particles"];
      idiscardIfAnyDeltaRLess(jets, leptons_stable, dR_cut);
      idiscardIfAnyDeltaRLess(jets, isolated_photons, dR_cut);
      //
      idiscardIfAnyDeltaRLess(btagging_jets, leptons_stable, dR_cut);
      idiscardIfAnyDeltaRLess(btagging_jets, isolated_photons, dR_cut);

      double dR_lepton_photon = deltaR(lep, lead_iso_photon); 
      if (dR_lepton_photon < dR_cut) vetoEvent; 
      _cutflows.fillnext();

      int n_jets = jets.size();
      if (n_jets < _jcuts["n_jets"])  vetoEvent;  
      _cutflows.fillnext();

      int n_b_jets = count(btagging_jets, hasBTag());
      if (n_b_jets>_jcuts["n_b_jets"]) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if ((tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"]) || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"]))) vetoEvent; 
      double dR_tagjets = deltaR(tag1_jet, tag2_jet); 
      if (dR_tagjets < dR_cut) vetoEvent; 
      _cutflows.fillnext();
      
      double dphi_MET_tagjet1 = fabs(deltaPhi(tag1_jet,fourvec_MET));
      double dphi_MET_tagjet2 = fabs(deltaPhi(tag2_jet,fourvec_MET));
      if ((dphi_MET_tagjet1<dbl(_jcuts["dphi_MET_tagjet"]) || dphi_MET_tagjet2<dbl(_jcuts["dphi_MET_tagjet"]))) vetoEvent; 
      _cutflows.fillnext();

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
      _cutflows.fillnext();

      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet));
      if (dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
      _cutflows.fillnext();

      const double centrality_jjly = fabs(0.5 * (fourvec_ly.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
      if (centrality_jjly > _jcuts["centrality_jjly"])  vetoEvent;
      _cutflows.fillnext();

      int n_gap_jets = 0;
      for (int i = 0; i < n_jets; i++) {
        const double i_jet_rap = jets[i].rap();
        if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
      }
      if (n_gap_jets > _jcuts["n_gap_jets"])  vetoEvent;
      _cutflows.fillnext();

      // do clipping - in case with 2+ bosons to avoid much work take two w with highest pt
      std::vector<FourMomentum> hs_bosons = {};
      const Particles all_particles = event.allParticles();
      for(const Particle& p_rivet : all_particles){
        ConstGenParticlePtr p_hepmc = p_rivet.genParticle();
        int status = p_hepmc->status();
        if (abs(status)==23 or abs(status)==22){
          int i_pid = p_hepmc->pid();
          FourMomentum i_mom = p_hepmc->momentum();
          if (abs(i_pid)==24 or abs(i_pid)==22){hs_bosons.push_back(i_mom);}
        }
      }
      std::sort(hs_bosons.begin(), hs_bosons.end(), [](FourMomentum const &a, FourMomentum const &b) {return a.pT() > b.pT(); });
      bool have_two_hs_bosons = false;
      double hs_diboson_mass = 0.0;
      int n_hs = hs_bosons.size();
      if (n_hs > 1){
        hs_diboson_mass = (hs_bosons[0] + hs_bosons[1]).mass() / GeV;
        have_two_hs_bosons = true;
      }
      if (!have_two_hs_bosons) vetoEvent; // just in case reject events where dont have wz somehow
      _cutflows.fillnext();

      //jet plots
      _h["n_jets"]->fill(n_jets);
      _h["pt_tagjet1"]->fill(tag1_jet.pt());
      _h["pt_tagjet2"]->fill(tag2_jet.pt());
      _h["m_tagjets"]->fill(m_tagjets);
      _h["dy_tagjets"]->fill(dy_tagjets);
      _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
      _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
      //lepton plots
      _h["n_lepton_stable"]->fill(nlep_stable);
      _h["pt_lepton"]->fill(lep.pT()); 
      _h["eta_lepton"]->fill(lep.eta());    
      //photon plots
      _h["n_photons_iso"]->fill(isolated_photons.size());
      _h["cone_frac_photon"]->fill(cone_to_photon_fracs[0]);
      _h["eta_photon"]->fill(lead_iso_photon.eta());
      _h["pt_photon"]->fill(lead_iso_photon.pT());
      // analysis-secific
      _h["pt_MET"]->fill(pt_MET);
      _h["m_T_W"]->fill(m_T_W);
      _h["m_ly"]->fill(m_ly);
      _h["centrality_jjly"]->fill(centrality_jjly);
      _h["dphi_MET_tagjet"]->fill(dphi_MET_tagjet1); _h["dphi_MET_tagjet"]->fill(dphi_MET_tagjet2);
      _h["dR_tagjets"]->fill(dR_tagjets);
      _h["dR_lepton_photon"]->fill(dR_lepton_photon);
      _h["n_gap_jets"]->fill(n_gap_jets);
      _h["n_b_jets"]->fill(n_b_jets);

      // clipping-related
      _h["m_Wy"]->fill(hs_diboson_mass);
      _h["pt_lepton_clip_inf"]->fill(lep.pT());
      if (ev_nominal_weight>=0){_c["pos_w_final_clip_inf"]->fill();}
      else {_c["neg_w_final_clip_inf"]->fill();}
      for (std::string& i_clip : _clips) {
        if (i_clip == "inf") continue; // as done above without cut
        int i_clip_num = std::stoi(i_clip);
        std::string i_pos_c_name = "pos_w_final_clip_" + i_clip;
        std::string i_neg_c_name = "neg_w_final_clip_" + i_clip;
        std::string i_hist_name = "pt_lepton_clip_" + i_clip;
        if (hs_diboson_mass < i_clip_num) {
          _h[i_hist_name]->fill(lep.pT());
          if (ev_nominal_weight >= 0) { _c[i_pos_c_name]->fill(); }
          else { _c[i_neg_c_name]->fill(); }
        }
      }

    }

    /// Normalise histograms etc., after the run
    void finalize() {
      std::cout << _cutflows.str();

      double pos_w_sum_initial = dbl(*_c["pos_w_initial"]); // from which also number of entries can be obtained
      double neg_w_sum_initial = dbl(*_c["neg_w_initial"]);
      double pos_w_sum_final = dbl(*_c["pos_w_final_clip_inf"]);
      double neg_w_sum_final = dbl(*_c["neg_w_final_clip_inf"]);
      MSG_INFO("\n pos weights initial final ratio " << pos_w_sum_initial <<" " << pos_w_sum_final <<" "<< pos_w_sum_final/pos_w_sum_initial << "\n" );
      MSG_INFO("\n neg weights initial final ratio " << neg_w_sum_initial <<" " << neg_w_sum_final <<" "<< neg_w_sum_final/neg_w_sum_initial << "\n" );

      // normalize all to 1 since in case of mostly negative weights not clear what it will do
      for (auto & i_name : _hist_names){ 
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
    std::string _cut_mode;
    Cut _electron_eta_cut;
    Cut _photon_eta_cut;
    Cut _muon_eta_cut;
    Cut _lepton_pt_cut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;
    std::vector<std::string> _clips {"inf", "3000", "2000", "1500", "1000", "700"};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(Wy_lvy);

}
