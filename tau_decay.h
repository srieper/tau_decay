/*
tau_decay analyses tau decays on generator level. It has two constructors:

tau_decay(MDSM_Mc_CL* mcs, int i)

'i' is the number of the particle you want to analyse in the MDSM_Mc_CL object.
IMPORTANT:
This class works with pointers to elements of the MDSM_Mc_CL object handed over to the constructor;
so make sure the MDSM_Mc_CL object isn't destroyed before the object of this class.
The constructor won't check whether you called it with the number of a tau. (This class will also
work for other particles although most methods won't produce any usefull output.)


tau_decay(MDSM_Mc_CL*, int, const MDSM_master_HCL* event)

If you additionally hand over a MDSM_master_HCL object, the best matching reconstructed lepton
will be saved (based on charge, decay mode and smalest deltaR)
IMPORTANT:
This class works with pointers to elements of the MDSM_master_HCL object handed over to the constructor;
so make sure the MDSM_master_HCL object isn't destroyed before the object of this class.


std::string decay_mode() const

returns the decay mode; in case the decay mode cannot be identified it will throw an
unknown_decay_mode exception.
You can use

DECAY_MODE decay_mode(int) const

instead; this method won't thow an exception.


double get_deltaR() const

returns the deltaR between the generator and reco tau.


double get_deltaR_gen() const

returns the deltaR between the tau and its visible component.


Int_t n_reco_particles() const
returns the number of reco leptons matching generator tau by charge and decay mode.


ROOT::Math::XYZTVector get_particle_vec(int) const
If you hand over the pdgid of a particle this method should return the sum of all 

*/

#ifndef TAU_DECAY___
#define TAU_DECAY___

#include "MDSM_Mc_CL.h"
#include "Math/Vector4D.h"
#include <string>
#include <vector>
#include <utility>
#include "MDSM_master_HCL.h"

class tau_decay {
public:
  enum DECAY_MODE {UNKNOWN=0, RHO=1, E=2, MU=4, PI=8, A1=16, HAD=32};

  tau_decay(MDSM_Mc_CL* mcs, int i) :n_reco_part(0), mother(0), reco_vector(0), deltaR_m(1000) { constr1(mcs, i); }
  tau_decay(MDSM_Mc_CL*, int, const MDSM_master_HCL* event);
  //~tau_decay();
  bool is_particle() const { return particle_id>0; }
  bool is_hadronic() const { return (pin.size() || pip.size() || K0s.size() || Kc.size() || K0l.size() || other_vis.size()); }
  bool is_1prong() const { return ((pin.size()+pip.size())==1 && (other_vis.size()+mu.size())==0); }
  bool is_3prong() const { return ((pin.size()+pip.size())==3 && (other_vis.size()+mu.size())==0); }
  bool is_3pi() const { return ((pin.size()+pip.size()+pi0.size())==3 && (other_vis.size()+mu.size())==0); }
  const ROOT::Math::XYZTVector* mother_vec() const { return mother; }
  const ROOT::Math::XYZTVector* vis_vec() const { return &vis_vector; }
  const ROOT::Math::XYZTVector* reco_vec() const { return reco_vector; }
  int n_pi() const { return pip.size()+pin.size(); }
  int n_pip() const { return pip.size(); }
  int n_pin() const { return pin.size(); }
  int n_pi0() const { return (pi0.size()+(gamma.size()+int(e.size()/2))/2); }
  int n_mu() const { return mu.size(); }
  int n_e() const { return e.size(); }
  int n_other() const { return other_vis.size()+K0s.size()+K0l.size()+Kc.size(); }
  int n_K() const { return K0s.size()+K0l.size()+Kc.size(); }
  Int_t n_reco_particles() const { return n_reco_part; }
  std::pair<double,double> get_Efrac_a1() const;
  std::pair<double,double> get_Efrac_a1to3pi() const;
  double get_deltaR() const { return deltaR_m; }
  double get_deltaR_gen() const { return deltaR_gen; }
  ROOT::Math::XYZTVector get_particle_vec(int) const;
  std::string decay_mode() const;
  DECAY_MODE decay_mode(int) const;

  //enum DECAY_MODE {UNKNOWN=0, RHO=1, E=2, MU=4, PI=8, A1=16, HAD=32};
  struct wrong_particle {};
  struct unknown_decay_mode{};

private:
  int particle_id;
  Int_t n_reco_part;
  const ROOT::Math::XYZTVector* mother;
  const ROOT::Math::XYZTVector* reco_vector;
  std::vector<const ROOT::Math::XYZTVector*> pip, pin, pi0, gamma, mu, e, K0s, Kc, K0l, other_vis, nu;
  double deltaR_gen;
  double deltaR_m;
  ROOT::Math::XYZTVector vis_vector;

  void constr1(MDSM_Mc_CL*, int);
  void set_reco_vector(const std::vector<ROOT::Math::XYZTVector>&, const std::vector<int>&);
};

#endif
