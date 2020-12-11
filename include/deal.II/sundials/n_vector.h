//-----------------------------------------------------------
//
//    Copyright (C) 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#ifndef dealii_sundials_n_vector_h
#define dealii_sundials_n_vector_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <sundials/sundials_nvector.h>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace internal
  {
    template <typename VectorType>
    VectorType *
    unwrap_nvector(N_Vector v);

    template <typename VectorType>
    N_Vector
    nvector_view(VectorType &vec);

    N_Vector_ID
    N_VGetVectorID_dealii(N_Vector v);

    N_Vector
    N_VCloneEmpty_dealii(N_Vector w);

    template <typename VectorType>
    N_Vector
    N_VClone_dealii(N_Vector w);


    template <typename VectorType>
    void
    N_VDestroy_dealii(N_Vector v);

    /**
     * Helper to create an empty vector with all operation set but no content.
     * @return a new N_Vector
     */
    template <typename VectorType>
    N_Vector
    N_VNewEmpty_dealii();


  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE
#endif
#endif
