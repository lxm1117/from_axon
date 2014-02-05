#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," Ca.mod");
    fprintf(stderr," KA_i0.mod");
    fprintf(stderr," KA_i1.mod");
    fprintf(stderr," KA_n1.mod");
    fprintf(stderr," KAa_i1.mod");
    fprintf(stderr," KDR_b0.mod");
    fprintf(stderr," KDR_i0.mod");
    fprintf(stderr," KDR_i1.mod");
    fprintf(stderr," KDR_mig.mod");
    fprintf(stderr," KDRa_i0.mod");
    fprintf(stderr," KDRa_i1.mod");
    fprintf(stderr," KDRb_i0.mod");
    fprintf(stderr," KDRb_i1.mod");
    fprintf(stderr," KDRc_i1.mod");
    fprintf(stderr," KIR_i0.mod");
    fprintf(stderr," Ktf_p1.mod");
    fprintf(stderr," Naf16_i1.mod");
    fprintf(stderr," Naf_b1.mod");
    fprintf(stderr," Naf_i0.mod");
    fprintf(stderr," Naf_i1.mod");
    fprintf(stderr," Naf_i_tst.mod");
    fprintf(stderr," Nafs_i1.mod");
    fprintf(stderr," Nap_i0.mod");
    fprintf(stderr," aabBK.mod");
    fprintf(stderr," ca1N.mod");
    fprintf(stderr," cdp.mod");
    fprintf(stderr," expsyn2b.mod");
    fprintf(stderr," expsyn2b_1.mod");
    fprintf(stderr," expsyn2c_1.mod");
    fprintf(stderr," fK_DR_n1.mod");
    fprintf(stderr," fNa_n1.mod");
    fprintf(stderr," h2_i0.mod");
    fprintf(stderr," h_2t.mod");
    fprintf(stderr," h_i0.mod");
    fprintf(stderr," h_in.mod");
    fprintf(stderr," h_n1.mod");
    fprintf(stderr," ipulse2.mod");
    fprintf(stderr," kext.mod");
    fprintf(stderr," kextna.mod");
    fprintf(stderr," morpho.mod");
    fprintf(stderr," pNa_n1.mod");
    fprintf(stderr," passive3.mod");
    fprintf(stderr, "\n");
  }
  _Ca_reg();
  _KA_i0_reg();
  _KA_i1_reg();
  _KA_n1_reg();
  _KAa_i1_reg();
  _KDR_b0_reg();
  _KDR_i0_reg();
  _KDR_i1_reg();
  _KDR_mig_reg();
  _KDRa_i0_reg();
  _KDRa_i1_reg();
  _KDRb_i0_reg();
  _KDRb_i1_reg();
  _KDRc_i1_reg();
  _KIR_i0_reg();
  _Ktf_p1_reg();
  _Naf16_i1_reg();
  _Naf_b1_reg();
  _Naf_i0_reg();
  _Naf_i1_reg();
  _Naf_i_tst_reg();
  _Nafs_i1_reg();
  _Nap_i0_reg();
  _aabBK_reg();
  _ca1N_reg();
  _cdp_reg();
  _expsyn2b_reg();
  _expsyn2b_1_reg();
  _expsyn2c_1_reg();
  _fK_DR_n1_reg();
  _fNa_n1_reg();
  _h2_i0_reg();
  _h_2t_reg();
  _h_i0_reg();
  _h_in_reg();
  _h_n1_reg();
  _ipulse2_reg();
  _kext_reg();
  _kextna_reg();
  _morpho_reg();
  _pNa_n1_reg();
  _passive3_reg();
}
