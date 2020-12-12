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
    /**
     * Retrieve the underlying vector attached to N_Vector @p v.
     *
     * @note Users must ensure that they ask for the correct VectorType when
     *   calling this function and there are no type-safety checks in place.
     *
     * @tparam VectorType type of the vector that is stored in @p v
     * @param v vector to unwrap
     * @return the vector that is stored inside @p v
     */
    template <typename VectorType>
    VectorType *
    unwrap_nvector(N_Vector v);

    /**
     * Get a view of the vector @p v as SUNDIALS N_Vector.
     *
     * This function does not allocate any new memory or perform any
     * modifications on @p v. Instead it attaches the given vector @p v and its
     * vector operations to the N_Vector interface of SUNDIALS. The returned
     * vector can then be passed to SUNDIALS functions.
     *
     * The original vector @p v must be kept alive as long as the returned view
     * is in use. If the view is destroyed with N_VDestroy() (part of SUNDIALS)
     * the original vector @p v is not affected.
     *
     * @param v the vector to view as SUNDIALS N_Vector
     * @return the N_Vector view of the passed vector
     */
    template <typename VectorType>
    N_Vector
    nvector_view(VectorType &v);

  } // namespace internal
} // namespace SUNDIALS

DEAL_II_NAMESPACE_CLOSE
#endif
#endif
