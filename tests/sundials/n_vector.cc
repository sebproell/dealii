//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
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

#include <deal.II/sundials/n_vector.h>

#include "../tests.h"

using namespace SUNDIALS::internal;

DeclExceptionMsg(NVectorTestError,
                 "The internal N_Vector implementation didn't pass a test");



template <typename VectorType>
void
test_nvector_view_unwrap()
{
  VectorType vector(2);
  auto       n_vector = nvector_view<VectorType>(vector);

  Assert(n_vector != nullptr, NVectorTestError());
  Assert(n_vector->content != nullptr, NVectorTestError());

  auto *vector_unwrapped = unwrap_nvector<VectorType>(n_vector);
  Assert(vector_unwrapped == &vector, NVectorTestError())

      deallog
    << "test_nvector_view_unwrap OK" << std::endl;
}



template <typename VectorType>
void
test_get_vector_id()
{
  VectorType vector(2);
  auto       n_vector = nvector_view<VectorType>(vector);
  auto       id       = N_VGetVectorID(n_vector);
  Assert(id == SUNDIALS_NVEC_CUSTOM, NVectorTestError());
  deallog << "test_get_vector_id OK" << std::endl;
}



template <typename VectorType>
void
test_clone()
{
  VectorType vector(3);
  auto       n_vector = SUNDIALS::internal::nvector_view<VectorType>(vector);
  auto       cloned   = N_VClone(n_vector);

  Assert(cloned != nullptr, NVectorTestError());
  AssertDimension(unwrap_nvector<VectorType>(cloned)->size(), 3);

  N_VDestroy(cloned);
  deallog << "test_clone OK" << std::endl;
}



template <typename VectorType>
void
test_destroy()
{
  GrowingVectorMemory<VectorType>                   mem;
  typename GrowingVectorMemory<VectorType>::Pointer vector(mem);
  auto n_vector = SUNDIALS::internal::nvector_view<VectorType>(*vector);

  Assert(n_vector != nullptr, NVectorTestError());
  auto cloned = N_VClone(n_vector);
  Assert(cloned != nullptr, NVectorTestError());

  // destroy memory-owning vector
  N_VDestroy(cloned);

  // clone empty(=without an actual vector attached) and destroy
  cloned = N_VCloneEmpty(n_vector);
  N_VDestroy(cloned);

  // destroy non-owning vector
  N_VDestroy(n_vector);

  // check that nullptr destroy is accepted
  n_vector = nullptr;
  N_VDestroy(n_vector);

  deallog << "test_destroy OK" << std::endl;
}



template <typename VectorType>
void
test_linear_sum()
{
  VectorType va({1.0, 2.0});
  VectorType vb({-1.0, -1.0});


  auto n_va = SUNDIALS::internal::nvector_view<VectorType>(va);
  auto n_vb = SUNDIALS::internal::nvector_view<VectorType>(vb);

  auto  n_vc = N_VClone(n_va);
  auto *p_vc = unwrap_nvector<VectorType>(n_vc);

  // test sum into third vector
  N_VLinearSum(1.0, n_va, 2.0, n_vb, n_vc);
  Assert(*p_vc == Vector<double>({-1.0, 0.0}), NVectorTestError());
  // repeat to test that sum overwrites initial content
  N_VLinearSum(1.0, n_va, 2.0, n_vb, n_vc);
  Assert(*p_vc == Vector<double>({-1.0, 0.0}), NVectorTestError());

  // test store sum into one of the summands
  N_VLinearSum(1.0, n_va, 2.0, n_vb, n_va);
  Assert(va == Vector<double>({-1.0, 0.0}), NVectorTestError());
  va = VectorType({1.0, 2.0});

  N_VLinearSum(1.0, n_va, 2.0, n_vb, n_vb);
  Assert(vb == Vector<double>({-1.0, 0.0}), NVectorTestError());

  deallog << "test_linear_sum OK" << std::endl;
}



int
main(int argc, char **argv)
{
  initlog();

  using VectorType = Vector<double>;

  // test conversion between vectors
  test_nvector_view_unwrap<VectorType>();
  test_get_vector_id<VectorType>();

  // test vector operations
  test_clone<VectorType>();
  test_destroy<VectorType>();
  test_linear_sum<VectorType>();
}