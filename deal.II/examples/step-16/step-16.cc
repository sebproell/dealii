/* $Id$ */
/* Author: Guido Kanschat, University of Heidelberg, 2003 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2003, 2004 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

				 // Again, the first few include files
				 // are already known, so we won't
				 // comment on them:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
				 // These are the new include files
				 // required for multi-level methods.
				 // First, the file defining the
				 // multigrid method itself.
#include <multigrid/multigrid.h>
				 // The DoFHandler is replaced by an
				 // MGDoFHandler which is defined
				 // here.
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>

				 // Then, we need some pre-made
				 // transfer routines between grids.
#include <multigrid/mg_transfer.h>
                                 // And a file in which equivalents to the
                                 // DoFTools class are declared:
#include <multigrid/mg_dof_tools.h>

				 // This is C++ ... see step 5 for
				 // further comments.
#include <fstream>
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


				 // This class is based on the same
				 // class in step 5. Remark that we
				 // replaced the DoFHandler by
				 // MGDoFHandler. since this inherits
				 // from DoFHandler, the new object
				 // incorporates the old functionality
				 // plus the new functions for degrees
				 // of freedom on different
				 // levels. Furthermore, we added
				 // MGLevelObjects for sparsity
				 // patterns and matrices.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem ();
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
				     // We add this function for
				     // assembling the multilevel
				     // matrices.
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    MGDoFHandler<dim>      mg_dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

				     // Here are the new objects for
				     // handling level matrices.
    MGLevelObject<SparsityPattern> mg_sparsity;
				     // We use number type float to
				     // save memory. It's only a
				     // preconditioner!
    MGLevelObject<SparseMatrix<float> > mg_matrices;
    
    Vector<double>       solution;
    Vector<double>       system_rhs;
};


				 // This function is as before.
template <int dim>
LaplaceProblem<dim>::LaplaceProblem () :
                fe (1),
		mg_dof_handler (triangulation)
{}



				 // This is the function of step 5
				 // augmented by the setup of the
				 // multi-level matrices in the end.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  mg_dof_handler.distribute_dofs (fe);

  std::cout << "   Number of degrees of freedom: "
	    << mg_dof_handler.n_dofs()
	    << std::endl;

  sparsity_pattern.reinit (mg_dof_handler.n_dofs(),
			   mg_dof_handler.n_dofs(),
			   mg_dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (mg_dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

				   // The multi-level objects are
				   // resized to hold matrices for
				   // every level. The coarse level is
				   // zero (this is mandatory right
				   // now but may change in a future
				   // revision). Remark, that the
				   // finest level is nlevels-1.
  const unsigned int nlevels = triangulation.n_levels();
  mg_sparsity.resize(0, nlevels-1);
  mg_matrices.resize(0, nlevels-1);

				   // Now, we have to build a matrix
				   // on each level. Technically, we
				   // could use the matrix initialized
				   // above on the finest
				   // level. Beware that this is not
				   // true anymore with local
				   // refinement!
  for (unsigned int level=0;level<nlevels;++level)
    {
      mg_sparsity[level].reinit (mg_dof_handler.n_dofs(level),
				 mg_dof_handler.n_dofs(level),
				 mg_dof_handler.max_couplings_between_dofs());
      MGTools::make_sparsity_pattern (mg_dof_handler, mg_sparsity[level], level);
      mg_sparsity[level].compress();
      mg_matrices[level].reinit(mg_sparsity[level]);
    }
}

				 // This is the standard assemble
				 // function you have seen a lot of
				 // times before.
				 //
				 // A small difference, though: we
				 // assemble the matrix for Helmholtz'
				 // equation so we can solve the
				 // Neumann boundary value problem.
template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{  
  QGauss2<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_values    |
				       update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active(),
						 endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix.clear ();
      cell_rhs.clear ();

				       // As before, we want the
				       // FEValues object to compute
				       // the quantities which we told
				       // him to compute in the
				       // constructor using the update
				       // flags.
      fe_values.reinit (cell);
				       // It should be noted that the
				       // creation of the
				       // coefficient_values object is
				       // done outside the loop over
				       // all cells to avoid memory
				       // allocation each time we
				       // visit a new cell.
      
				       // With all this, the loops
				       // then look like this (the
				       // parentheses around the
				       // product of the two gradients
				       // are needed to indicate the
				       // dot product; we have to
				       // overrule associativity of
				       // the operator* here, since
				       // the compiler would otherwise
				       // complain about an undefined
				       // product of double*gradient
				       // since it parses
				       // left-to-right):
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
				   * fe_values.shape_grad(j,q_point)
				   * fe_values.JxW(q_point));

					     // For the right hand
					     // side, a constant value
					     // is used again:
	    cell_rhs(i) += (fe_values.shape_value(i,q_point)
			    * 1.0 * fe_values.JxW(q_point));
	  };


      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    system_matrix.add (local_dof_indices[i],
			       local_dof_indices[j],
			       cell_matrix(i,j));
	  
	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
	};
    };

				   // Again use zero boundary values:
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}


				 // Here is another assemble
				 // function. The integration core is
				 // the same as above. Only the loop
				 // goes over all existing cells now
				 // and the results must be entered
				 // into the right matrix.

				 // Since we only do multi-level
				 // preconditioning, no right-hand
				 // side is assembled here.
template <int dim>
void LaplaceProblem<dim>::assemble_multigrid () 
{  
  QGauss2<dim>  quadrature_formula;

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   UpdateFlags(update_gradients |
				       update_q_points  |
				       update_JxW_values));

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.n_quadrature_points;

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

				   // 
  typename DoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					  endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
				       // Remember the level of the
				       // current cell.
      const unsigned int level = cell->level();
      cell_matrix.clear ();

				       // Compute the values specified
				       // by update flags above.
      fe_values.reinit (cell);

				       // This is exactly the
				       // integration loop of the cell
				       // matrix above.
      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point)
				   * fe_values.shape_grad(j,q_point)
				   * fe_values.JxW(q_point));
	  };


				       // Oops! This is a tiny
				       // difference easily
				       // forgotten. The indices we
				       // want here are the ones for
				       // that special level, not for
				       // the global
				       // matrix. Therefore, a little
				       // 'mg' entered into the
				       // function call.
      cell->get_mg_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
					     // And now add everything
					     // to the matrix on the
					     // right level.
	    mg_matrices[level].add (local_dof_indices[i],
                                    local_dof_indices[j],
                                    cell_matrix(i,j));
	};
    };

				   // Again use zero boundary values:
  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}



				 // The solution process again looks
				 // mostly like in the previous
				 // examples. However, we will now use
				 // a preconditioned conjugate
				 // gradient algorithm. It is not very
				 // difficult to make this change:
template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverCG<>              cg (solver_control, vector_memory);

				   // The only thing we have to alter
				   // is that we need an object which
				   // will act as a preconditioner. We
				   // will use SSOR (symmetric
				   // successive overrelaxation), with
				   // a relaxation factor of 1.2. For
				   // this purpose, the SparseMatrix
				   // class has a function which does
				   // one SSOR step, and we need to
				   // package the address of this
				   // function together with the
				   // matrix on which it should act
				   // (which is the matrix to be
				   // inverted) and the relaxation
				   // factor into one object. This can
				   // be done like this:
  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);
				   // (Note that we did not have to
				   // explicitely pass the address of
				   // the SSOR function of the matrix
				   // to this objects, rather it is
				   // hardcoded into the object, thus
				   // the name.)
				   //
				   // The default template parameters
				   // of the ``PreconditionRelaxation''
				   // class is the matrix type, which
				   // defaults to the types used in
				   // this program.

				   // Calling the solver now looks
				   // mostly like in the example
				   // before, but where there was an
				   // object of type
				   // PreconditionIdentity before,
				   // there now is the newly generated
				   // preconditioner object.
  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  std::cout << "   " << solver_control.last_step()
	    << " CG iterations needed to obtain convergence."
	    << std::endl;
}



				 // Here is the data output, which is
				 // a simplified version of step 5. We
				 // do a standard gnuplot output for
				 // each grid produced in the
				 // refinement process.
template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
				   // Construct and initialize a DataOut object
  DataOut<dim> data_out;

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

				   // The following block generates
				   // the file name and opens the
				   // file. This looks awkward because
				   // of a change in the handling of
				   // string streams (See step 5 for explanation).
  
#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream filename;
#else
  std::ostrstream filename;
#endif
  filename << "solution-"
	   << cycle
	   << ".gnuplot";

#ifndef HAVE_STD_STRINGSTREAM
  filename << std::ends;
#endif

#ifdef HAVE_STD_STRINGSTREAM
  std::ofstream output (filename.str().c_str());
#else
  std::ofstream output (filename.str());
#endif

  data_out.write_gnuplot (output);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<6; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
					   // Generate a simple hypercube grid.
	  GridGenerator::hyper_cube(triangulation);
	  
	}
				       // If this is not the first
				       // cycle, then simply refine
				       // the grid once globally.
      else
	triangulation.refine_global (1);

				       // Write some output and do all
				       // the things that we have
				       // already seen in the previous
				       // examples.
      std::cout << "   Number of active cells: "
		<< triangulation.n_active_cells()
		<< std::endl
		<< "   Total number of cells: "
		<< triangulation.n_cells()
		<< std::endl;

      setup_system ();
      assemble_system ();
      solve ();
      output_results (cycle);
    };
}

    

				 // The main function looks mostly
				 // like the one in the previous
				 // example, so we won't comment on it
				 // further.
int main () 
{
  deallog.depth_console (0);

  LaplaceProblem<2> laplace_problem_2d;
  laplace_problem_2d.run ();

				   // Finally, we have promised to
				   // trigger an exception in the
				   // Coefficient class. For this, we
				   // have to call its ``value_list''
				   // function with two arrays of
				   // different size (the number in
				   // parentheses behind the name of
				   // the object). We have commented
				   // out these lines in order to
				   // allow the program to exit
				   // gracefully in normal situations
				   // (we use the program in
				   // day-to-day testing of changes to
				   // the library as well), so you
				   // will only get the exception by
				   // un-commenting the following
				   // lines.
/*  
  Coefficient<2>    coefficient;
  std::vector<Point<2> > points (2);
  std::vector<double>    coefficient_values (1);
  coefficient.value_list (points, coefficient_values);
*/
  
  return 0;
}
