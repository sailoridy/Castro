#ifndef _Castro_error_F_H_
#define _Castro_error_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void laplac_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

  void denerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

  void velerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

  void temperror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

  void presserror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

  void nucerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

#ifdef RADIATION
  void raderror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);
#endif
  
  void enterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int tagval, const int clearval,
     amrex::Real* var, const int* var_lo, const int* var_hi,
     const int* lo, const int* hi,
     const int ncomp,
     const int level);

#ifdef AMREX_DIMENSION_AGNOSTIC
  void set_problem_tags
    (const int* lo, const int* hi,
     int* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* tagval, const int* clearval,
     const amrex::Real* dx, const amrex::Real* problo,
     const amrex::Real* time, const int* level);
#else
  void set_problem_tags
    (const int lo[], const int hi[],
     int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     BL_FORT_FAB_ARG(state),
     const int* tagval, const int* clearval,
     const amrex::Real dx[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);
#endif

#ifdef __cplusplus
}
#endif

#endif
