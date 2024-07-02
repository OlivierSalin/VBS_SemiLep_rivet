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
    class ZZ_llll : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZZ_llll);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
        std::string txt_dir = "/exp/atlas/salin/ATLAS/VBS_mc/plotting/";
        
        std::string out_dir = getOption("OUTDIR");
        
        _docut = 0; // most cuts on number of particles are always applied to avoid segfault
        if (out_dir.find("DOCUT_YES") != string::npos) _docut = 1;
        std::cout << "++++++received outidir" << out_dir << "meaning _docut is " << _docut << "\n";

        std::string jsonfilestr =  txt_dir + "ZZ_llll_cuts.json"; 
        std::cout << "++++++assume .json for this ZZ_llll" << "is " << jsonfilestr << "\n";
        std::ifstream json_file(jsonfilestr);
        
        _jcuts = json::parse(json_file);
        std::cout << "++++++ to check json 1 var got photon pt min" << _jcuts["m_tagjets"] << "\n";
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
        VetoedFinalState hadrons(FinalState(Cuts::absetaIn(0.0, _jcuts["eta_tagjets"])));
        hadrons.addVetoOnThisFinalState(dressed_e);
        hadrons.addVetoOnThisFinalState(dressed_mu);
        declare(hadrons, "hadrons");
        FastJets jetsfs(hadrons, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
        declare(jetsfs, "jets");

        declare(MissingMomentum(), "METFinder");

        // plots common with others
        std::ifstream jet_hist_file(txt_dir + "/jet_hists.json");      
        json jet_hist = json::parse(jet_hist_file);
        for (json::iterator it = jet_hist.begin(); it != jet_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        std::ifstream lep_hist_file(txt_dir + "/lepton_hists.json");      
        json lep_hist = json::parse(lep_hist_file);
        for (json::iterator it = lep_hist.begin(); it != lep_hist.end(); ++it) {
        book(_h[it.key()], it.key(), it.value()[0], it.value()[1], it.value()[2]);
        _hist_names.push_back(it.key());
        }
        // plots that are not in other ana
        std::ifstream ana_hist_file(txt_dir + "/ZZ_llll_hists.json");      
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
        _cutflows.addCutflow("Wmy_lvy_selections", {"have_four_lep","pt_lep1_2","dR_all_pairs","SFOC_2pairs_min",
                            "m_llll","n_jets","pt_tagjet1_2","m_tagjets","dy_tagjets","centrality_quadjj"});

        // setup for  file used for drawing images
        if (_docut==1){
        std::vector<std::string> pic_particles = {"tagjet1", "tagjet2", "lepton1", "lepton2", "lepton3", "lepton4"};
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

        int nlep = leptons.size();
        if (nlep<_jcuts["n_lepton_stable"])  vetoEvent; 
        _cutflows.fillnext();

        if (_docut==1 && (leptons[0].pT()<_jcuts["pt_lepton1"] || leptons[1].pT()<_jcuts["pt_lepton2"])) vetoEvent;
        _cutflows.fillnext();
        
        // deltaR between leptons for all pairs
        int num_pair_bad_dR = 0;
        for (int i = 0; i < nlep; i++) {
            const Particle& i_lep = leptons[i];
            for (int j = 0; j < nlep; j++) {
                // to avoid comparing to itself and double couting like (i,j),(j,i) as two separate pairs
                if (i==j || j>i) continue;
                const Particle& j_lep = leptons[j];
                double i_dR = deltaR(i_lep.rap(),i_lep.phi(), j_lep.rap(),j_lep.phi());
                if (i_dR<_jcuts["dR_all_pairs"]) num_pair_bad_dR+=1;
            }
        }
        if (_docut==1 && num_pair_bad_dR>0) vetoEvent;
        _cutflows.fillnext();

        // check SFOC and m_ll and save possible pairs
        std::vector<std::vector<int>> pairs_ind;
        for (int i = 0; i < nlep; i++) {
            const Particle& i_lep = leptons[i];
            for (int j = 0; j < nlep; j++) {
                // to avoid comparing to itself and double couting like (i,j),(j,i) as two separate pairs
                if (i==j || j>i) continue;
                const Particle& j_lep = leptons[j];
                int i_sum_pids = i_lep.pid() + j_lep.pid(); // to have SFOC will be 0
                double i_m_ll = (i_lep.mom() + j_lep.mom()).mass()/GeV;
                if (i_sum_pids==0 && i_m_ll>_jcuts["m_ll_all_pairs"]) pairs_ind.push_back({i,j});
            }
        }
        if (pairs_ind.size()<2) vetoEvent;
        _cutflows.fillnext();

        // order found pairs by how close they are m_z
        std::sort(pairs_ind.begin(), pairs_ind.end(), [this,leptons](std::vector<int> &ind_pair_1, std::vector<int> &ind_pair_2) {
            double dist_pair_1 = pair_m_dist_m_z(leptons[ind_pair_1[0]], leptons[ind_pair_1[1]]);
            double dist_pair_2 = pair_m_dist_m_z(leptons[ind_pair_2[0]], leptons[ind_pair_2[1]]);
            return dist_pair_1 < dist_pair_2; // closest to m_z will be first
            });
        // define two pairs as the ones with smallest dist
        const Particle& pair_1_lep_1 = leptons[pairs_ind[0][0]];
        const Particle& pair_1_lep_2 = leptons[pairs_ind[0][1]];
        const Particle& pair_2_lep_1 = leptons[pairs_ind[1][0]];
        const Particle& pair_2_lep_2 = leptons[pairs_ind[1][1]];
        int sum_abs_pids_quadruplet =  fabs(pair_1_lep_1.pid())+fabs(pair_1_lep_2.pid())+fabs(pair_2_lep_1.pid())+fabs(pair_2_lep_2.pid());
        double m_ll_pair_1 = (pair_1_lep_1.mom()+pair_1_lep_2.mom()).mass()/GeV;
        double m_ll_pair_2 = (pair_2_lep_1.mom()+pair_2_lep_2.mom()).mass()/GeV;

        const FourMomentum fourvec_llll = pair_1_lep_1.mom() + pair_1_lep_2.mom() + pair_2_lep_1.mom() + pair_2_lep_2.mom();
        double m_llll = fourvec_llll.mass()/GeV;
        if (_docut==1 && m_llll<_jcuts["m_llll"]) vetoEvent; 
        _cutflows.fillnext();

        // // Retrieve clustered jets, sorted by pT, with a minimum pT cut
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > dbl(_jcuts["pt_jet"])*GeV);
        idiscardIfAnyDeltaRLess(jets, leptons, 0.2);

        int n_jets = jets.size();
        if (n_jets < _jcuts["n_jets"])  vetoEvent;  
        _cutflows.fillnext();

        const FourMomentum tag1_jet = jets[0].mom();
        const FourMomentum tag2_jet = jets[1].mom();
        if (_docut==1 && (tag1_jet.pT()<dbl(_jcuts["pt_tagjet1"]) || tag2_jet.pT()<dbl(_jcuts["pt_tagjet2"]))) vetoEvent; 
        _cutflows.fillnext();

        const double m_tagjets = (tag1_jet + tag2_jet).mass()/GeV;
        if (_docut==1 && m_tagjets<_jcuts["m_tagjets"]) vetoEvent;
        _cutflows.fillnext();
        
        const FourMomentum fourvec_jj = tag1_jet+tag2_jet;
        const FourMomentum fourvec_lllljj = fourvec_llll + fourvec_jj;

        const double dy_tagjets = fabs(tag1_jet.rap() - tag2_jet.rap());
        if (_docut==1 && dy_tagjets<_jcuts["dy_tagjets"]) vetoEvent;
        if (_docut==1 && tag1_jet.eta()*tag2_jet.eta()>0) vetoEvent; 
        _cutflows.fillnext();

        const double centrality_quadjj = fabs(0.5 * (fourvec_llll.rap() - (tag1_jet.rap()+tag2_jet.rap())/2) / (tag1_jet.rap()-tag2_jet.rap()));
        if (_docut==1 && centrality_quadjj > _jcuts["centrality_quadjj"])  vetoEvent;
        _cutflows.fillnext();

        int n_gap_jets = 0;
        for (int i = 0; i < n_jets; i++) {
            const double i_jet_rap = jets[i].rap();
            if ((i_jet_rap < tag1_jet.rap() && i_jet_rap > tag2_jet.rap()) || (i_jet_rap < tag2_jet.rap() && i_jet_rap > tag1_jet.rap()))  ++n_gap_jets;
        }

        //jet plots
        _h["n_jets"]->fill(n_jets);
        _h["pt_tagjet1"]->fill(tag1_jet.pt());
        _h["pt_tagjet2"]->fill(tag2_jet.pt());
        _h["eta_tagjets"]->fill(tag1_jet.eta()); _h["eta_tagjets"]->fill(tag2_jet.eta());
        _h["phi_tagjets"]->fill(tag1_jet.phi()); _h["phi_tagjets"]->fill(tag2_jet.phi());
        _h["m_tagjets"]->fill(m_tagjets);
        _h["dy_tagjets"]->fill(dy_tagjets);
        _h["dphi_tagjets"]->fill(deltaPhi(tag1_jet,tag2_jet));
        //lepton plots
        _h["n_lepton_stable"]->fill(nlep);
        _h["pt_lepton"]->fill(pair_1_lep_1.pT()); _h["pt_lepton"]->fill(pair_1_lep_2.pT()); 
        _h["pt_lepton"]->fill(pair_2_lep_1.pT()); _h["pt_lepton"]->fill(pair_2_lep_2.pT()); 
        _h["eta_lepton"]->fill(pair_1_lep_1.eta()); _h["eta_lepton"]->fill(pair_1_lep_2.eta());    
        _h["eta_lepton"]->fill(pair_2_lep_1.eta()); _h["eta_lepton"]->fill(pair_2_lep_2.eta());    
        //ana-specific
        _h["sum_abs_pids_quadruplet"]->fill(sum_abs_pids_quadruplet);
        _h["m_ll_quadruplet"]->fill(m_ll_pair_1); _h["m_ll_quadruplet"]->fill(m_ll_pair_2);
        _h["m_llll"]->fill(m_llll);
        _h["n_gap_jets"]->fill(n_gap_jets);
        _h["pt_tagjets"]->fill(fourvec_jj.pT());
        _h["pt_llll"]->fill(fourvec_llll.pT());
        _h["centrality_quadjj"]->fill(centrality_quadjj);
        _h["pt_lllljj"]->fill(fourvec_lllljj.pT());

        // save weights after cuts
        if (ev_nominal_weight>=0){_c["pos_w_final"]->fill();}
        else {_c["neg_w_final"]->fill();}

        // file used for drawing images
        if (_docut==1){
            int ind_bigger_eta_tagjet = (tag1_jet.eta() >  tag2_jet.eta()) ? 0 : 1;
            int ind_smaller_eta_tagjet = static_cast<int>(!static_cast<bool>(ind_bigger_eta_tagjet));
            //
            int ind_of_ind_bigger_eta_pair_1 = (pair_1_lep_1.eta() >  pair_1_lep_2.eta()) ? 0 : 1;
            int ind_of_ind_smaller_eta_pair_1 = static_cast<int>(!static_cast<bool>(ind_of_ind_bigger_eta_pair_1));
            int ind_bigger_eta_pair_1 = pairs_ind[0][ind_of_ind_bigger_eta_pair_1]; 
            int ind_smaller_eta_pair_1 = pairs_ind[0][ind_of_ind_smaller_eta_pair_1]; 
            //
            int ind_of_ind_bigger_eta_pair_2 = (pair_2_lep_1.eta() >  pair_2_lep_2.eta()) ? 0 : 1;
            int ind_of_ind_smaller_eta_pair_2 = static_cast<int>(!static_cast<bool>(ind_of_ind_bigger_eta_pair_2));
            int ind_bigger_eta_pair_2 = pairs_ind[0][ind_of_ind_bigger_eta_pair_2]; 
            int ind_smaller_eta_pair_2 = pairs_ind[0][ind_of_ind_smaller_eta_pair_2]; 
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
            //
            //pair 1 lepton 1        
            pic_csv << leptons[ind_bigger_eta_pair_1].eta() << ";";
            pic_csv << leptons[ind_bigger_eta_pair_1].phi() << ";";
            pic_csv << leptons[ind_bigger_eta_pair_1].pt() << ";";
            //pair 1 lepton 2        
            pic_csv << leptons[ind_smaller_eta_pair_1].eta() << ";";
            pic_csv << leptons[ind_smaller_eta_pair_1].phi() << ";";
            pic_csv << leptons[ind_smaller_eta_pair_1].pt() << ";";
            //
            //pair 2 lepton 1        
            pic_csv << leptons[ind_bigger_eta_pair_2].eta() << ";";
            pic_csv << leptons[ind_bigger_eta_pair_2].phi() << ";";
            pic_csv << leptons[ind_bigger_eta_pair_2].pt() << ";";
            //pair 2 lepton 2        
            pic_csv << leptons[ind_smaller_eta_pair_2].eta() << ";";
            pic_csv << leptons[ind_smaller_eta_pair_2].phi() << ";";
            pic_csv << leptons[ind_smaller_eta_pair_2].pt() << ";";
            // terminate line
            pic_csv << "\n";
        }
    }

    double pair_m_dist_m_z(const Particle& lep_1, const Particle& lep_2){
        double i_m_ll = (lep_1.mom() + lep_2.mom()).mass()/GeV;
        double i_m_ll_dist_m_z = fabs(i_m_ll-91.18);
        return i_m_ll_dist_m_z;
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
    Cut _muon_eta_cut;
    Cut _electron_pt_cut;
    Cut _muon_pt_cut;
    json _jcuts;
    Cutflows _cutflows;
    std::vector<std::string> _hist_names;
    /// @}


    };


    RIVET_DECLARE_PLUGIN(ZZ_llll);

}
