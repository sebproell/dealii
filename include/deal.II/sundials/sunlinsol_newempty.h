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

/**
 * This file contains a verbatim copy of code snippets distributed within the
 * SUNDIALS package, see the license below.
 */

/* -----------------------------------------------------------------
 * Programmer(s): Daniel Reynolds @ SMU
 *                David J. Gardner, Carol S. Woodward, and
 *                Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2020, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.

  BSD 3-Clause License

  Copyright (c) 2002-2019, Lawrence Livermore National Security and Southern
  Methodist University. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------*/

#ifndef dealii_sundials_sunlinsol_newempty_h
#define dealii_sundials_sunlinsol_newempty_h

#include <deal.II/base/config.h>
#ifdef DEAL_II_WITH_SUNDIALS

DEAL_II_NAMESPACE_OPEN
namespace SUNDIALS
{
  namespace internal
  {
    /* -----------------------------------------------------------------
     * Create a new empty SUNLinearSolver object
     * ----------------------------------------------------------------- */
    SUNLinearSolver
    SUNLinSolNewEmpty()
    {
      SUNLinearSolver     LS;
      SUNLinearSolver_Ops ops;

      /* create linear solver object */
      LS = NULL;
      LS = (SUNLinearSolver)malloc(sizeof *LS);
      if (LS == NULL)
        return (NULL);

      /* create linear solver ops structure */
      ops = NULL;
      ops = (SUNLinearSolver_Ops)malloc(sizeof *ops);
      if (ops == NULL)
        {
          free(LS);
          return (NULL);
        }

      /* initialize operations to NULL */
      ops->gettype           = NULL;
      ops->setatimes         = NULL;
      ops->setpreconditioner = NULL;
      ops->setscalingvectors = NULL;
      ops->initialize        = NULL;
      ops->setup             = NULL;
      ops->solve             = NULL;
      ops->numiters          = NULL;
      ops->resnorm           = NULL;
      ops->resid             = NULL;
      ops->lastflag          = NULL;
      ops->space             = NULL;
      ops->free              = NULL;

      /* attach ops and initialize content to NULL */
      LS->ops     = ops;
      LS->content = NULL;

      return (LS);
    }

    /* -----------------------------------------------------------------
     * Free a generic SUNLinearSolver (assumes content is already empty)
     * ----------------------------------------------------------------- */
    void
    SUNLinSolFreeEmpty(SUNLinearSolver S)
    {
      if (S == NULL)
        return;

      /* free non-NULL ops structure */
      if (S->ops)
        free(S->ops);
      S->ops = NULL;

      /* free overall N_Vector object and return */
      free(S);
      return;
    }

  } // namespace internal
} // namespace SUNDIALS
DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SUNDIALS
#endif // dealii_sundials_sunlinsol_newempty_h
