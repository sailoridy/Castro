
#include "Castro.H"
#include "Castro_error_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif
#include <Castro_error.H>
#include <Castro_prob_err_F.H>

using std::string;
using namespace amrex;

typedef StateDescriptor::BndryFunc BndryFunc;

void
Castro::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

    // This routine uses the special evaluation of the second derivative
    //   and can be called with any variable.  Note that two ghost cells are needed.
//  err_list.add("density",2,ErrorRec::Special,ca_laplac_error);
//  err_list.add("pressure",2,ErrorRec::Special,ca_laplac_error);

    err_list.add("density",1,ErrorRec::Special,ca_denerror);
    err_list.add("Temp",1,ErrorRec::Special,ca_temperror);
    err_list.add("pressure",1,ErrorRec::Special,ca_presserror);
    err_list.add("x_velocity",1,ErrorRec::Special,ca_velerror);
#if (BL_SPACEDIM >= 2)
    err_list.add("y_velocity",1,ErrorRec::Special,ca_velerror);
#endif
#if (BL_SPACEDIM == 3)
    err_list.add("z_velocity",1,ErrorRec::Special,ca_velerror);
#endif

//   err_list.add("entropy",1,ErrorRec::Special,ca_enterror);

#ifdef REACTIONS
    err_list.add("t_sound_t_enuc",0,ErrorRec::Special,ca_nucerror);
#endif

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
      err_list.add("rad",1,ErrorRec::Special,ca_raderror);
    }
#endif

    // Save the number of built-in functions; this will help us
    // distinguish between those, and the ones the user is about to add.

    num_err_list_default = err_list.size();

#include <Castro_prob_err_list.H>

}

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_laplac_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      laplac_error(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
                   *tagval, *clearval,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   *ncomp, *level);
  }

  void ca_denerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      denerror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
               *tagval, *clearval,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               *ncomp, *level);
  }

  void ca_velerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      velerror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
               *tagval, *clearval,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               *ncomp, *level);
  }

  void ca_temperror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      temperror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
                *tagval, *clearval,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
                AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                *ncomp, *level);
  }

  void ca_presserror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      presserror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
                 *tagval, *clearval,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 *ncomp, *level);
  }

  void ca_nucerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      nucerror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
               *tagval, *clearval,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               *ncomp, *level);
  }

#ifdef RADIATION
  void ca_raderror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      raderror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
               *tagval, *clearval,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               *ncomp, *level);
  }
#endif
  
  void ca_enterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     amrex::Real* data, const int* data_lo, const int* data_hi,
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level)
  {
#pragma gpu
      enterror(tag, AMREX_INT_ANYD(tag_lo), AMREX_INT_ANYD(tag_hi),
               *tagval, *clearval,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               *ncomp, *level);
  }

#ifdef __cplusplus
}
#endif
