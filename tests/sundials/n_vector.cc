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

void
test_nvector_view(N_Vector v)
{
  Assert(v != nullptr, NVectorTestError());
  Assert(v->content != nullptr, NVectorTestError());
  deallog << "test_nvector_view OK" << std::endl;
}

template <typename VectorType>
void
test_unwrap_vector(N_Vector v)
{
  auto *vector = unwrap_nvector<VectorType>(v);
  Assert(vector != nullptr, NVectorTestError());
  deallog << "test_unwrap_vector OK" << std::endl;
}

void
test_get_vector_id(N_Vector v)
{
  auto id = N_VGetVectorID(v);
  Assert(id == SUNDIALS_NVEC_CUSTOM, NVectorTestError());
  deallog << "test_get_vector_id OK" << std::endl;
}

void
test_clone(N_Vector v)
{
  auto cloned = N_VClone(v);
  Assert(cloned != nullptr, NVectorTestError());
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


int
main(int argc, char **argv)
{
  initlog();

  using VectorType = Vector<double>;
  GrowingVectorMemory<VectorType>          mem;
  GrowingVectorMemory<VectorType>::Pointer vector(mem);
  auto n_vector = SUNDIALS::internal::nvector_view<VectorType>(*vector);

  test_nvector_view(n_vector);
  test_unwrap_vector<VectorType>(n_vector);
  test_get_vector_id(n_vector);
  test_clone(n_vector);


  test_destroy<VectorType>();
}