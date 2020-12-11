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


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SUNDIALS
#  include <deal.II/base/exceptions.h>

#  include <deal.II/lac/vector_memory.h>

#  include <deal.II/sundials/n_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace SUNDIALS
{
  namespace internal
  {
    template <typename VectorType>
    class NVectorContent
    {
    public:
      /**
       * Create a non-owning content with an existing vector.
       * @param vector
       */
      NVectorContent(VectorType *vector);

      /**
       * Allocate a new vector wrapped in a new content object. The vector
       * will be deallocated automatically when this object is destroyed.
       */
      NVectorContent();

      VectorType *
      get();

    private:
      using PointerType =
        std::unique_ptr<VectorType, std::function<void(VectorType *)>>;
      GrowingVectorMemory<VectorType> mem;
      PointerType                     vector;
    };

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


    template <typename VectorType>
    void
    N_VLinearSum_dealii(realtype a,
                        N_Vector x,
                        realtype b,
                        N_Vector y,
                        N_Vector z);

    /**
     * Helper to create an empty vector with all operation set but no content.
     * @return a new N_Vector
     */
    template <typename VectorType>
    N_Vector
    N_VNewEmpty_dealii();

  } // namespace internal
} // namespace SUNDIALS



template <typename VectorType>
SUNDIALS::internal::NVectorContent<VectorType>::NVectorContent()
  : vector(typename VectorMemory<VectorType>::Pointer(mem))
{}



template <typename VectorType>
SUNDIALS::internal::NVectorContent<VectorType>::NVectorContent(
  VectorType *vector)
  : vector(vector, [](VectorType *) { /* not owning memory -> don't free*/ })
{}



template <typename VectorType>
VectorType *
SUNDIALS::internal::NVectorContent<VectorType>::get()
{
  return vector.get();
}



template <typename VectorType>
VectorType *
SUNDIALS::internal::unwrap_nvector(N_Vector v)
{
  Assert(v != nullptr, ExcInternalError());
  Assert(v->content != nullptr, ExcInternalError());
  auto *interface = reinterpret_cast<NVectorContent<VectorType> *>(v->content);
  return interface->get();
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::nvector_view(VectorType &vec)
{
  // Create an N_Vector with operators attached and empty content
  N_Vector v = N_VNewEmpty_dealii<VectorType>();
  if (v == nullptr)
    return (nullptr);

  // Create vector content using a pointer to the interface
  v->content = new NVectorContent<VectorType>(&vec);
  if (v->content == nullptr)
    {
      N_VDestroy(v);
      return nullptr;
    }

  return (v);
}



N_Vector_ID SUNDIALS::internal::N_VGetVectorID_dealii(N_Vector)
{
  return SUNDIALS_NVEC_CUSTOM;
}



N_Vector
SUNDIALS::internal::N_VCloneEmpty_dealii(N_Vector w)
{
  Assert(w != nullptr, ExcInternalError());
  N_Vector v = N_VNewEmpty();
  if (v == nullptr)
    return (nullptr);

  // Copy operations
  if (N_VCopyOps(w, v))
    {
      N_VDestroy(v);
      return (nullptr);
    }

  return (v);
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::N_VClone_dealii(N_Vector w)
{
  N_Vector v = N_VCloneEmpty_dealii(w);
  if (v == nullptr)
    return (nullptr);

  auto cloned = new NVectorContent<VectorType>();
  auto n      = unwrap_nvector<VectorType>(w)->size();

  cloned->get()->reinit(n);

  v->content = cloned;
  if (v->content == nullptr)
    {
      N_VDestroy(v);
      return nullptr;
    }

  return (v);
}



template <typename VectorType>
void
SUNDIALS::internal::N_VDestroy_dealii(N_Vector v)
{
  if (v == nullptr)
    return;

  if (v->content != nullptr)
    {
      auto *interface =
        reinterpret_cast<NVectorContent<VectorType> *>(v->content);
      // the NVectorInterface knows if it owns the memory and will free
      // correctly
      delete interface;
      v->content = nullptr;
    }

  /* free ops and vector */
  if (v->ops != nullptr)
    {
      free(v->ops);
      v->ops = nullptr;
    }
  free(v);
  v = nullptr;

  return;
}



template <typename VectorType>
void
SUNDIALS::internal::N_VLinearSum_dealii(realtype a,
                                        N_Vector x,
                                        realtype b,
                                        N_Vector y,
                                        N_Vector z)
{
  auto *x_dealii = unwrap_nvector<VectorType>(x);
  auto *y_dealii = unwrap_nvector<VectorType>(y);
  auto *z_dealii = unwrap_nvector<VectorType>(z);

  if (z_dealii == x_dealii)
    z_dealii->sadd(a, b, *y_dealii);
  else if (z_dealii == y_dealii)
    z_dealii->sadd(b, a, *x_dealii);
  else
    {
      *z_dealii = 0;
      z_dealii->add(a, *x_dealii, b, *y_dealii);
    }
}



/**
 * Helper to create an empty vector with all operation set but no content.
 * @return a new N_Vector
 */
template <typename VectorType>
N_Vector
SUNDIALS::internal::N_VNewEmpty_dealii()
{
  N_Vector v = N_VNewEmpty();
  if (v == nullptr)
    return (nullptr);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid = N_VGetVectorID_dealii;
  v->ops->nvclone       = N_VClone_dealii<VectorType>;
  v->ops->nvcloneempty  = N_VCloneEmpty_dealii;
  v->ops->nvdestroy     = N_VDestroy_dealii<VectorType>;

  /* standard vector operations */
  v->ops->nvlinearsum = N_VLinearSum_dealii<VectorType>;
  return (v);
}


template N_Vector
SUNDIALS::internal::nvector_view<Vector<double>>(Vector<double> &);

template Vector<double> *
  SUNDIALS::internal::unwrap_nvector<Vector<double>>(N_Vector);

DEAL_II_NAMESPACE_CLOSE



#endif
