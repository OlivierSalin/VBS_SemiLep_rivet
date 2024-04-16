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
  class Zy_vvy : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(Zy_vvy);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      std::string txt_dir = "/exp/atlas/kurdysh/vbs_cross_terms_study/plotting/";
      
      std::string out_dir = getOption("OUTDIR");
      
      _docut = 0; // most cuts on number of particles are always applied to avoid segfault
      if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
      std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

      std::string jsonfilestr =  txt_dir + "Zy_vvy_cuts.json"; 
      std::cout << "++++++assume .json for this Zy_vvy" << "is " << jsonfilestr << "\n";
      std::ifstream json_file(jsonfilestr);
      
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
      std::ifstream jet_hist_file(txt_dir + "/jet_hists.json");      
      json jet_hist = json::parse(jet_hist_file);
      for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      std::ifstream y_hist_file(txt_dir + "/photon_hists.json");      
      json y_hist = json::parse(y_hist_file);
      for (json::iterator it = y_hist.begin(); it != y_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      // from lepton only take nleps to see 0
      std::ifstream lep_hist_file(txt_dir + "/lepton_hists.json"); 
      json lep_hist = json::parse(lep_hist_file);
      book(_h["n_lepton_stable"],  "n_lepton_stable", 
          lep_hist["n_lepton_stable"][0], lep_hist["n_lepton_stable"][1], lep_hist["n_lepton_stable"][2]);
      _hist_names.push_back("n_lepton_stable");
      // plots that are not in other ana
      std::ifstream ana_hist_file(txt_dir + "/Zy_vvy_hists.json");      
      json ana_hist = json::parse(ana_hist_file);
      for (json::iterator it = ana_hist.begin(); it != ana_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
      }
      
      //counter for efficiency
      book(_c["pos_w_initial"],"pos_w_initial");
      book(_c["pos_w_final"],"pos_w_final");
      book(_c["neg_w_initial"],"neg_w_initial");
      book(_c["neg_w_final"],"neg_w_final");

      // Cut-flows 
      _cutflows.addCutflow("Zy_vvy_selections", {"no_leptons", "have_iso_photons_ok_pt_eta",
                            "n_jets","pt_tagjet1_2","m_tagjets","centrality_jjy", 
                            "pt_MET", "dphi_MET_photon", "dphi_MET_tagjet"});


      // setup for  file used for drawing images
      if (_docut==1){
        std::vector<std::string> pic_particles = {"tagjet1", "tagjet2", "MET", "photon"};
        std::ofstream pic_csv (out_dir + "/info_for_image.csv", std::ofstream::out);
        for (auto & i_p : pic_particles){ 
          pic_csv << "eta_" + i_p +";";
          pic_csv << "phi_" + i_p +";";
          pic_csv << "pt_" + i_p +";";
        }
        pic_csv << "\n";
        pic_csv.close();
      }
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
      if (_docut==1){
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt(_electron_eta_cut && _lepton_pt_cut);
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt(_muon_eta_cut && _lepton_pt_cut);
      }
      else{
        e_stable = apply<FinalState>(event, "e_stable").particlesByPt();
        mu_stable = apply<FinalState>(event, "mu_stable").particlesByPt();
      }
      // will ahve either e or mu pair so don't need to sort combined e+mu array
      Particles leptons_stable = e_stable + mu_stable; 

      int nlep_stable = leptons_stable.size();
      if (nlep_stable!=_jcuts["n_lepton_stable"])  vetoEvent; 
      _cutflows.fillnext();

      //photons
      Particles photons = apply<FinalState>(event, "photons").particlesByPt(_photon_eta_cut);
      //photon cone calculation and photon OR with leptons 
      Particles isolated_photons;
      std::vector<double> cone_to_photon_fracs = {};
      Particles cone_sum_particles = apply<VetoedFinalState>(event, "ElAndJetsForPhotonIsoCalc").particles();
      for (const Particle &i_high_pt_photon : photons) {
        // check photon isolation
        double i_high_pt_photon_cone_E_40 = 0.0;
        double i_high_pt_photon_cone_pT_20 = 0.0;
        for (const Particle &i_p_cone : cone_sum_particles) {
          if (deltaR(i_high_pt_photon, i_p_cone) < 0.4) i_high_pt_photon_cone_E_40 += i_p_cone.Et();
          if (deltaR(i_high_pt_photon, i_p_cone) < 0.2) i_high_pt_photon_cone_pT_20 += i_p_cone.pT();
        }
        double i_cone_to_photon_frac_pT_20 = i_high_pt_photon_cone_pT_20 / i_high_pt_photon.pT(); 
        if (i_cone_to_photon_frac_pT_20 > _jcuts["cone_frac_photon"])  continue;
        if (i_high_pt_photon_cone_E_40 > 0.022*i_high_pt_photon.pT() + 2.45)  continue;
        isolated_photons += i_high_pt_photon;
        cone_to_photon_fracs.push_back(i_cone_to_photon_frac_pT_20);
      }
      if (isolated_photons.size()!=_jcuts["n_photons_iso"])  vetoEvent;
      const Particle& iso_photon = isolated_photons[0];
      if (_docut==1 && iso_photon.pT() < dbl(_jcuts["pt_photon"])*GeV) vetoEvent;
      _cutflows.fillnext();
      
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV); // then will do cut on two leading pt>50
      // Remove all jets within certain photon
      idiscardIfAnyDeltaRLess(jets, isolated_photons, 0.3);

      int n_jets = jets.size();
      if (n_jets < _jcuts["n_jets"])  vetoEvent;  
      _cutflows.fillnext();

      const FourMomentum tag1_jet = jets[0].mom();
      const FourMomentum tag2_jet = jets[1].mom();
      if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV || tag2_jet.pT()<dbl(_jcuts["pt_tagjet1"])*GeV)) vetoEvent; 
      _cutflows.fillnext();

      const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
      if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
      _cutflows.fillnext();
      const double dy_tagjets =  fabs(deltaRap(tag1_jet, tag2_jet)); // they dont do the cut

      const double centrality_jjy = fabs(0.5 * (iso_photon.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
      if (_docut==1 && centrality_jjy > _jcuts["centrality_jjy"])  vetoEvent;
      _cutflows.fillnext();

      // MET
      const MissingMomentum& METfinder = apply<MissingMomentum>(event, "METFinder");
      const double pt_MET = METfinder.missingPt()/GeV;
      if (_docut==1 && pt_MET<_jcuts["pt_MET"]) vetoEvent;
      _cutflows.fillnext();

      const FourMomentum fourvec_MET = METfinder.missingMomentum();
      double dphi_MET_photon = fabs(deltaPhi(iso_photon,fourvec_MET));
      if (_docut==1 && (dphi_MET_photon<_jcuts["dphi_MET_photon"])) vetoEvent;
      _cutflows.fillnext();

      double dphi_MET_tagjet1 = fabs(deltaPhi(tag1_jet,fourvec_MET));
      double dphi_MET_tagjet2 = fabs(deltaPhi(tag2_jet,fourvec_MET));
      if (_docut==1 && (dphi_MET_tagjet1<_jcuts["dphi_MET_tagjet"] || dphi_MET_tagjet2<_jcuts["dphi_MET_tagjet"])) vetoEvent;
      _cutflows.fillnext();

      //jet plots
      _h["n_jets"]->fill(n_jets);
      _h["pt_tagjet1"]->fill(tag1_jet.pt());
      _h["pt_tagjet2"]->fill(tag2_jet.pt());
      _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
      _h["phi_tagjets"]->fill(tag1_jet.phi()); _h["phi_tagjets"]->fill(tag2_jet.phi());
      _h["m_tagjets"]->fill(m_tagjets);
      _h["dy_tagjets"]->fill(dy_tagjets);
      _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
      //lepton plots - absence of 
      _h["n_lepton_stable"]->fill(nlep_stable);
      //photon plots
      _h["n_photons_iso"]->fill(isolated_photons.size());
      _h["cone_frac_photon"]->fill(cone_to_photon_fracs[0]);
      _h["eta_photon"]->fill(iso_photon.eta());
      _h["pt_photon"]->fill(iso_photon.pT());
      // analysis-secific
      _h["pt_MET"]->fill(pt_MET);
      _h["dphi_MET_photon"]->fill(dphi_MET_photon);
      _h["dphi_MET_tagjet"]->fill(dphi_MET_tagjet1); _h["dphi_MET_tagjet"]->fill(dphi_MET_tagjet2);
      _h["centrality_jjy"]->fill(centrality_jjy);
    
      // save weights after cuts
      if (ev_nominal_weight>=0){_c["pos_w_final"]->fill();}
      else {_c["neg_w_final"]->fill();}
    
      // file used for drawing images
      if (_docut==1){
        int ind_bigger_eta_tagjet = (tag1_jet.eta() >  tag2_jet.eta()) ? 0 : 1;
        int ind_smaller_eta_tagjet = static_cast<int>(!static_cast<bool>(ind_bigger_eta_tagjet));
        // pulling file into common with init() _fout didn't work so re-open
        std::ofstream pic_csv (getOption("OUTDIR") + "/info_for_image.csv", std::ofstream::app); 
        // tagjet1
        pic_csv << jets[ind_bigger_eta_tagjet].eta() << ";";
        pic_csv << jets[ind_bigger_eta_tagjet].phi() << ";";
        pic_csv << jets[ind_bigger_eta_tagjet].pt() << ";";
        //tagjet2
        pic_csv << jets[ind_smaller_eta_tagjet].eta() << ";";
        pic_csv << jets[ind_smaller_eta_tagjet].phi() << ";";
        pic_csv << jets[ind_smaller_eta_tagjet].pt() << ";";
        //MET        
        pic_csv << fourvec_MET.eta() << ";";
        pic_csv << fourvec_MET.phi() << ";";
        pic_csv << fourvec_MET.pt() << ";";
        //photon        
        pic_csv << iso_photon.eta() << ";";
        pic_csv << iso_photon.phi() << ";";
        pic_csv << iso_photon.pt() << ";";
        // terminate line
        pic_csv << "\n";
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      std::string cut_str = _cutflows.str();
      std::string cutflow_file = getOption("OUTDIR") + "/cutflow.txt";
      std::ofstream ofs (cutflow_file, std::ofstream::out); 
      ofs << cut_str;
      ofs.close();

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
    Cut _photon_eta_cut;
    Cut _muon_eta_cut;
    Cut _lepton_pt_cut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;

    /// @}


  };


  RIVET_DECLARE_PLUGIN(Zy_vvy);

}
