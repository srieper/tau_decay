#include "tau_decay.h"
#include "deltaR.h"
//#include "Math/VectorUtil_Cint.h"

using std::vector;
using std::cout;
using std::endl;

tau_decay::tau_decay(MDSM_Mc_CL* mcs, int i, const MDSM_master_HCL* event) :n_reco_part(0), mother(0), reco_vector(0), deltaR_m(1000) {
  constr1(mcs, i);
  switch(decay_mode(1)) {
    case DECAY_MODE::MU:
      set_reco_vector(event->t_mu, event->t_mu_charge);
      break;
    case DECAY_MODE::E:
      set_reco_vector(event->t_el, event->t_el_charge);
      break;
    default:
      set_reco_vector(event->t_tau, event->t_tau_charge);
  }
  /*try {
    if(decay_mode()=="mu")
      set_reco_vector(event->t_mu, event->t_mu_charge);
    else if(decay_mode()=="e")
      set_reco_vector(event->t_el, event->t_el_charge);
    else {
      set_reco_vector(event->t_tau, event->t_tau_charge);
    }
  } catch(unknown_decay_mode) {
    set_reco_vector(event->t_tau, event->t_tau_charge);
  }*/
}

void tau_decay::set_reco_vector(const std::vector<ROOT::Math::XYZTVector>& t_lep, const vector<int>& t_lep_charge) {
  for(unsigned i=0; i<t_lep.size(); ++i)
    if((is_particle() && t_lep_charge[i]==-1) || (!is_particle() && t_lep_charge[i]==+1)) {
      double deltaR_temp = deltaR::DeltaR(t_lep[i], *mother_vec());
      if (deltaR_temp<deltaR_m) {
        deltaR_m=deltaR_temp;
        reco_vector=&t_lep[i];
      }
      ++n_reco_part;
    }
}

std::string tau_decay::decay_mode() const {
  switch(decay_mode(1)) {
    case DECAY_MODE::RHO:
      return "rho";
    case DECAY_MODE::A1:
      return "a1";
    case DECAY_MODE::MU:
      return "mu";
    case DECAY_MODE::E:
      return "e";
    case DECAY_MODE::PI:
      return "pi";
  }
  /*if(n_other())
    throw unknown_decay_mode();
  if(!n_mu()) {
    if(n_pi()==1 && !n_pi0())
      return "pi";
    else if(n_pi()==1 && n_pi0()==1)
      return "rho";
    else if((n_pi()+n_pi0())==3)
      return "a1";
  }
  if(is_hadronic())
    //throw unknown_decay_mode();
    return "had";
  else if(n_mu()==1 && n_e()==0)
    return "mu";
  else if(n_e()==1)
    return "e";*/

  throw unknown_decay_mode();
}

tau_decay::DECAY_MODE tau_decay::decay_mode(int) const {
  if(n_other())
    return DECAY_MODE::UNKNOWN;
  if(!n_mu()) {
    if(n_pi()==1 && !n_pi0())
      return DECAY_MODE::PI;
    else if(n_pi()==1 && n_pi0()==1)
      return DECAY_MODE::RHO;
    else if((n_pi()+n_pi0())==3)
      return DECAY_MODE::A1;
  }
  if(is_hadronic())
    return DECAY_MODE::HAD;
  else if(n_mu()==1 && n_e()==0)
    return DECAY_MODE::MU;
  else if(n_e()==1)
    return DECAY_MODE::E;

  return DECAY_MODE::UNKNOWN;
}

ROOT::Math::XYZTVector tau_decay::get_particle_vec(int pdgid) const {
  ROOT::Math::XYZTVector vec;
  switch (pdgid) {
    case 211: //pi+
      for(unsigned i=0; i<pip.size();i++) vec += *pip[i];
      return vec;
    case -211: //pi-
      for(unsigned i=0; i<pin.size();i++) vec += *pin[i];
      return vec;
    case 111: case -111: //pi0
      for(unsigned i=0; i<pi0.size(); i++) vec += *pi0[i];
      for(unsigned i=0; i<e.size(); i++) vec += *e[i];
    case 22: //gamma
      for(unsigned i=0; i<gamma.size(); i++) vec += *gamma[i];
      return vec;
    case 11: case -11: //e
      for(unsigned i=0; i<e.size(); i++) vec += *e[i];
      return vec;
    case 13: case -13:  //mu
      for(unsigned i=0; i<mu.size(); i++) vec += *mu[i];
      return vec;
    case 12: case -12: //nu
    case 14: case -14:
    case 16: case -16:
      for(unsigned i=0; i<nu.size(); i++) vec += *nu[i];
      return vec;
    default:
      //for(unsigned i=0; i<other_vis.size(); i++) vec += *other_vis[i];
      return vec;
  }
}

std::pair<double,double> tau_decay::get_Efrac_a1() const {
  if(!(n_pi()+n_pi0()==3)) throw wrong_particle();
  int charge = is_particle()?-1:1;
  double E1 = n_pi0()?get_particle_vec(charge*211).E():get_particle_vec(-charge*211).E();
  double E23 = n_pi0()?get_particle_vec(111).E():get_particle_vec(charge*211).E();
  return std::make_pair(E1, E23);
}

std::pair<double,double> tau_decay::get_Efrac_a1to3pi() const {
  if(n_pi()!=3 || n_pi0()!=0) throw wrong_particle();
  int charge = is_particle()?-1:1;
  double E1 =  get_particle_vec(-charge*211).E();
  double E23 =  get_particle_vec(charge*211).E();
  return std::make_pair(E1, E23);
}

void tau_decay::constr1(MDSM_Mc_CL* mcs, int i) {
  //reco_vector=0;
  //deltaR_m=0;
  particle_id = mcs->pdgid.akt_val;
  mother = mcs->get_vec_val(i);
  vector<int> daughters = mcs->getDaughters();
  for (unsigned k=0; k<daughters.size(); ++k) {
    mcs->set_akt_values(daughters[k]);
    if(mcs->status.akt_val!=1) {
      cout << "Status id: " << mcs->status.akt_val << ", pdgid: " << mcs->pdgid.akt_val << endl;
      continue;
    }
    switch (mcs->pdgid.akt_val) {
      case 111: case -111:
        pi0.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *pi0.back();
        break;
      case 22:
        gamma.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *gamma.back();
        break;
      case 211:
        pip.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *pip.back();
        break;
      case -211:
        pin.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *pin.back();
        break;
      case 13: case -13:
        mu.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *mu.back();
        break;
      case 11: case -11:
        e.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *e.back();
        break;
      case 310:
        K0s.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *K0s.back();
        break;
      case 130:
        K0l.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *K0l.back();
        break;
      case 321: case -321:
        Kc.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *Kc.back();
        break;
      case 12: case -12: //nu_e
      case 14: case -14: //nu_mu
      case 16: case -16: //nu_tau
        nu.push_back(mcs->get_vec_val(daughters[k]));
        break;
      default:
        other_vis.push_back(mcs->get_vec_val(daughters[k]));
        vis_vector += *other_vis.back();
        cout << "pdgid: " << mcs->pdgid.akt_val << endl;
        break;
    }
  }
  mcs->set_akt_values(i);
  deltaR_gen = mcs->mc_deltaR(vis_vector);
}
