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

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/la_parallel_vector.h>
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

    /**
     * Helper to create an empty vector with all operation set but no content.
     * @return a new N_Vector
     */
    template <typename VectorType>
    N_Vector
    create_empty_nvector();

    namespace NVectorOperations
    {
      N_Vector_ID
      get_vector_id(N_Vector v);

      N_Vector
      clone_empty(N_Vector w);

      template <typename VectorType>
      N_Vector
      clone(N_Vector w);


      template <typename VectorType>
      void
      destroy(N_Vector v);

      template <typename VectorType>
      sunindextype
      get_global_length(N_Vector v);

      template <typename VectorType>
      void
      linear_sum(realtype a, N_Vector x, realtype b, N_Vector y, N_Vector z);

      template <typename VectorType>
      void
      set_constant(realtype c, N_Vector v);

    } // namespace NVectorOperations

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
  N_Vector v = create_empty_nvector<VectorType>();
  Assert(v != nullptr, ExcInternalError());

  // Create a non-owning vector content using a pointer and attach to N_Vector
  v->content = new NVectorContent<VectorType>(&vec);
  Assert(v->content != nullptr, ExcInternalError());
  return (v);
}



N_Vector_ID SUNDIALS::internal::NVectorOperations::get_vector_id(N_Vector)
{
  return SUNDIALS_NVEC_CUSTOM;
}



N_Vector
SUNDIALS::internal::NVectorOperations::clone_empty(N_Vector w)
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
SUNDIALS::internal::NVectorOperations::clone(N_Vector w)
{
  N_Vector v = clone_empty(w);

  // the corresponding delete is called in N_VDestroy_dealii
  auto  cloned   = new NVectorContent<VectorType>();
  auto *w_dealii = unwrap_nvector<VectorType>(w);
  cloned->get()->reinit(*w_dealii);

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
SUNDIALS::internal::NVectorOperations::destroy(N_Vector v)
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
}



template <typename VectorType>
sunindextype
SUNDIALS::internal::NVectorOperations::get_global_length(N_Vector v)
{
  return unwrap_nvector<VectorType>(v)->size();
}



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::linear_sum(realtype a,
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



template <typename VectorType>
void
SUNDIALS::internal::NVectorOperations::set_constant(realtype c, N_Vector v)
{
  auto *v_dealii = unwrap_nvector<VectorType>(v);
  *v_dealii      = c;
}



template <typename VectorType>
N_Vector
SUNDIALS::internal::create_empty_nvector()
{
  N_Vector v = N_VNewEmpty();
  if (v == nullptr)
    return (nullptr);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid = NVectorOperations::get_vector_id;
  v->ops->nvclone       = NVectorOperations::clone<VectorType>;
  v->ops->nvcloneempty  = NVectorOperations::clone_empty;
  v->ops->nvdestroy     = NVectorOperations::destroy<VectorType>;
  //  v->ops->nvspace           = undef;
  //  v->ops->nvgetcommunicator = undef;
  v->ops->nvgetlength = NVectorOperations::get_global_length<VectorType>;

  /* standard vector operations */
  v->ops->nvlinearsum = NVectorOperations::linear_sum<VectorType>;
  v->ops->nvconst     = NVectorOperations::set_constant<VectorType>;
  //  v->ops->nvprod         = undef;
  //  v->ops->nvdiv          = undef;
  //  v->ops->nvscale        = undef;
  //  v->ops->nvabs          = undef;
  //  v->ops->nvinv          = undef;
  //  v->ops->nvaddconst     = undef;
  //  v->ops->nvdotprod      = undef;
  //  v->ops->nvmaxnorm      = undef;
  //  v->ops->nvwrmsnorm     = undef;
  //  v->ops->nvwrmsnormmask = undef;
  //  v->ops->nvmin          = undef;
  //  v->ops->nvwl2norm      = undef;
  //  v->ops->nvl1norm       = undef;
  //  v->ops->nvcompare      = undef;
  //  v->ops->nvinvtest      = undef;
  //  v->ops->nvconstrmask   = undef;
  //  v->ops->nvminquotient  = undef;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  //  v->ops->nvdotprodlocal     = undef;
  //  v->ops->nvmaxnormlocal     = undef;
  //  v->ops->nvminlocal         = undef;
  //  v->ops->nvl1normlocal      = undef;
  //  v->ops->nvinvtestlocal     = undef;
  //  v->ops->nvconstrmasklocal  = undef;
  //  v->ops->nvminquotientlocal = undef;
  //  v->ops->nvwsqrsumlocal     = undef;
  //  v->ops->nvwsqrsummasklocal = undef;
  return (v);
}


template N_Vector
SUNDIALS::internal::nvector_view<Vector<double>>(Vector<double> &);

template Vector<double> *
  SUNDIALS::internal::unwrap_nvector<Vector<double>>(N_Vector);

template N_Vector
SUNDIALS::internal::nvector_view<BlockVector<double>>(BlockVector<double> &);

template BlockVector<double> *
  SUNDIALS::internal::unwrap_nvector<BlockVector<double>>(N_Vector);

template LinearAlgebra::distributed::Vector<double> *
  SUNDIALS::internal::unwrap_nvector<
    LinearAlgebra::distributed::Vector<double>>(N_Vector);

template N_Vector
SUNDIALS::internal::nvector_view<LinearAlgebra::distributed::Vector<double>>(
  LinearAlgebra::distributed::Vector<double> &);

DEAL_II_NAMESPACE_CLOSE



#endif
