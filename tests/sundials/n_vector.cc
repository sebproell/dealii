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

#include <deal.II/base/logstream.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/sundials/n_vector.h>

#include "../tests.h"

#include "../test_grids.h"

using namespace SUNDIALS::internal;

// anonymous namespace groups helper functions for testing
namespace
{
  DeclExceptionMsg(NVectorTestError,
                   "The internal N_Vector implementation didn't pass a test.");



  template <typename VectorType>
  VectorType
  create_test_vector();



  template <>
  Vector<double>
  create_test_vector()
  {
    return Vector<double>(3 /*size*/);
  }



  template <>
  BlockVector<double>
  create_test_vector()
  {
    const int num_blocks = 2;
    const int size_block = 3;
    return BlockVector<double>(num_blocks, size_block);
  }

  template <>
  LinearAlgebra::distributed::Vector<double>
  create_test_vector()
  {
    return LinearAlgebra::distributed::Vector<double>(3 /*size*/);
  }

  bool
  operator==(const LinearAlgebra::distributed::Vector<double> &a,
             const LinearAlgebra::distributed::Vector<double> &b)
  {
    return std::equal(a.begin(), a.end(), b.begin());
  }
} // namespace


template <typename VectorType>
void
test_nvector_view_unwrap()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = nvector_view<VectorType>(vector);

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
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = nvector_view<VectorType>(vector);
  auto id       = N_VGetVectorID(n_vector);
  Assert(id == SUNDIALS_NVEC_CUSTOM, NVectorTestError());
  deallog << "test_get_vector_id OK" << std::endl;
}



template <typename VectorType>
void
test_clone()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = nvector_view<VectorType>(vector);
  auto cloned   = N_VClone(n_vector);

  Assert(cloned != nullptr, NVectorTestError());
  AssertDimension(unwrap_nvector<VectorType>(cloned)->size(), vector.size());

  N_VDestroy(cloned);
  deallog << "test_clone OK" << std::endl;
}



template <typename VectorType>
void
test_destroy()
{
  GrowingVectorMemory<VectorType>                   mem;
  typename GrowingVectorMemory<VectorType>::Pointer vector(mem);
  auto n_vector = nvector_view<VectorType>(*vector);

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
test_length()
{
  auto vector   = create_test_vector<VectorType>();
  auto n_vector = nvector_view<VectorType>(vector);
  Assert(N_VGetLength(n_vector) == vector.size(), NVectorTestError());

  deallog << "test_length OK" << std::endl;
}



template <typename VectorType>
void
test_linear_sum()
{
  auto va       = create_test_vector<VectorType>();
  auto vb       = create_test_vector<VectorType>();
  auto vc       = create_test_vector<VectorType>();
  auto expected = create_test_vector<VectorType>();

  auto n_va = SUNDIALS::internal::nvector_view<VectorType>(va);
  auto n_vb = SUNDIALS::internal::nvector_view<VectorType>(vb);
  auto n_vc = SUNDIALS::internal::nvector_view<VectorType>(vc);

  va = 1.0;
  vb = 2.0;

  expected = -3.0;

  // test sum into third vector
  N_VLinearSum(1.0, n_va, -2.0, n_vb, n_vc);
  Assert(vc == expected, NVectorTestError());
  // repeat to test that sum overwrites initial content
  N_VLinearSum(1.0, n_va, -2.0, n_vb, n_vc);
  Assert(vc == expected, NVectorTestError());

  // test store sum into one of the summands
  N_VLinearSum(1.0, n_va, -2.0, n_vb, n_va);
  Assert(va == expected, NVectorTestError());
  va = 1.0;

  N_VLinearSum(1.0, n_va, -2.0, n_vb, n_vb);
  Assert(vb == expected, NVectorTestError());

  deallog << "test_linear_sum OK" << std::endl;
}



template <typename VectorType>
void
test_set_constant()
{
  auto vector   = create_test_vector<VectorType>();
  auto expected = create_test_vector<VectorType>();
  expected      = 1.0;

  auto *n_vector = nvector_view(vector);

  N_VConst(1.0, n_vector);
  Assert(vector == expected, NVectorTestError());

  N_VConst(1.0, n_vector);
  Assert(vector == expected, NVectorTestError());

  deallog << "test_set_constant OK" << std::endl;
}



template <typename VectorType>
void
run_all_tests(const std::string &prefix)
{
  LogStream::Prefix p(prefix);
  // test conversion between vectors
  test_nvector_view_unwrap<VectorType>();
  test_get_vector_id<VectorType>();

  // test vector operations
  test_clone<VectorType>();
  test_destroy<VectorType>();
  test_length<VectorType>();
  test_linear_sum<VectorType>();
  test_set_constant<VectorType>();
}

int
main(int argc, char **argv)
{
  initlog();

  using VectorType = Vector<double>;

  run_all_tests<Vector<double>>("Vector<double>");
  run_all_tests<BlockVector<double>>("BlockVector<double>");
  run_all_tests<LinearAlgebra::distributed::Vector<double>>(
    "LinearAlgebra::distributed::Vector<double>");
}