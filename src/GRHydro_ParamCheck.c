 /*@@
   @file      GRHydro_ParamCheck.c
   @date      Sun Jun 26 14:07:35 EDT 2016
   @author
   @desc
   Parameter checking routine.
   @enddesc
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include <assert.h>

 /*@@
   @routine    GRHydro_CheckSpatialOrder
   @date       Sun Jun 26 14:08:21 EDT 2016
   @author     Roland Haas
   @desc
   Checks the parameters.
   @enddesc
   @calls
   @calledby
   @history

   @endhistory

@@*/

void CCTK_FNAME(GRHydro_CheckSpatialOrder)(void)
{
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_IsImplementationActive("ADMMacros")) {
    int type;
    const CCTK_INT *spatial_order =
      (CCTK_INT*)CCTK_ParameterGet("spatial_order", "ADMMacros", &type);
    assert(spatial_order != NULL);
    assert(type == PARAMETER_INT);
    if(sources_spatial_order != *spatial_order) {
      CCTK_VParamWarn(CCTK_THORNSTRING,
                      "ADMMacros::spatial_order (%d) differs from GRHydro::sources_spatial_order (%d), this is almost certainly incorrect since GRHydro no longer uses ADMMacros",
                      (int)*spatial_order, (int)sources_spatial_order);
    }
  }
}
