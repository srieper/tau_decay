# tau_decay
For this class to work an interface has to be added to MDSM_Mc_CL:

const ROOT::Math::XYZTVector* MDSM_Mc_CL::get_vec_val(int i) {
  return &(*all_lvec_val[0]->pall_val)[i];
}
