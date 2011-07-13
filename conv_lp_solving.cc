/**** written by Thomas Schoenemann as an employee of Lund University, August 2010 ****/

#include "conv_lp_solving.hh"
#include "vector.hh"
#include "projection.hh"

//#define USE_EXPLICIT_GRADIENT

#ifdef USE_CUDA
#include "cuda_conv_lp.cuh"
#endif

void eq_constrained_lp_solving_auglagrange(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
					   const SparseMatrixDescription<double> matrix_descr, const double* rhs,  
					   double* solution) {

  Math1D::NamedVector<double> lagrange_multiplier(nConstraints,0.0,MAKENAME(lagrange_multiplier));

  double penalty = 100.0;

  double alpha = 1e-7 / penalty;

  Math1D::Vector<double> ax(nConstraints,0.0);

  for (uint outer_iter = 1; outer_iter <= 15; outer_iter++) {

    if (outer_iter != 1) {
      penalty *= 10.0;
      alpha *= 0.01;
    }

    std::cerr << "############# penalty " << penalty << std::endl;

    double last_energy = 1e50;


    for (uint iter = 1; iter <= 200; iter++) {

      /*** 1. calculate current energy ****/
      double energy = 0.0;
      for (uint v=0; v < nVars; v++)
	energy += solution[v] * cost[v];

      ax.set_constant(0.0);
      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += entry*solution[col]; 
      }      
      for (uint c=0; c < nConstraints; c++) {
	double temp = ax[c] - rhs[c];

	energy += lagrange_multiplier[c] *  temp;
	energy += 0.5*penalty * temp*temp;
      }

      std::cerr << "energy: " << energy << std::endl;

      if (energy > last_energy) 
	alpha *= 0.1;

      last_energy = energy;

      /*** 2. calculate gradient ***/

      Math1D::Vector<double> grad(nVars);
      for (uint v=0; v < nVars; v++)
	grad[v] = cost[v];

      ax.set_constant(0.0);

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += solution[col] * entry;
      }

      for (uint c=0; c < nConstraints; c++)
	ax[c] -= rhs[c];

      ax *= penalty;
      ax += lagrange_multiplier;

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	grad[col] += entry*ax[row];
      }

      /*** 3. go in the direction of the negative gradient  ***/
      for (uint v=0; v < nVars; v++) {

	solution[v] -= alpha*grad[v];
      }

      /*** 4. reproject to the convex set defined by the variable bounds ***/
      for (uint v=0; v < nVars; v++) {

	if (solution[v] < var_lb[v])
	  solution[v] = var_lb[v];
	else if (solution[v] > var_ub[v])
	  solution[v] = var_ub[v];
      }

    } //end of inner iterations


    double lp_energy = 0.0;
    for (uint v=0; v < nVars; v++)
      lp_energy += solution[v]*cost[v];

    std::cerr << "lp-energy: " << lp_energy << std::endl;


    /*** update lagrange multipliers ****/
#if 1
    ax.set_constant(0.0);
    for (uint k=0; k < matrix_descr.nEntries(); k++) {
      
      const uint row = matrix_descr.row_indices()[k];
      const uint col = matrix_descr.col_indices()[k];
      const double entry = matrix_descr.value()[k];
      
      ax[row] += entry*solution[col]; 
    }      
    for (uint c=0; c < nConstraints; c++) {

      double temp = ax[c] - rhs[c];
      
      lagrange_multiplier[c] += penalty * temp;
    }
#endif
    
  }
}

void eq_and_simplex_constrained_lp_solve_auglagr(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
						 const SparseMatrixDescription<double> matrix_descr, const double* rhs,  
						 uint nSimplices, const uint* simplex_starts, double* solution) {


  Math1D::NamedVector<double> lagrange_multiplier(nConstraints,0.0,MAKENAME(lagrange_multiplier));

  double penalty = 100.0;

  //double alpha = 1e-1 / penalty;

  Math1D::Vector<double> ax(nConstraints,0.0);

  for (uint outer_iter = 1; outer_iter <= 15; outer_iter++) {

    if (outer_iter != 1) {
      penalty *= 1.25; //1.5;
      //alpha *= 0.01;
    }

    double alpha = 1e-1 / penalty;
  

    std::cerr << "############# penalty " << penalty << std::endl;

    double last_energy = 1e50;

    uint iter_since_restart = 0;

    double save_energy = 1e50;

    for (uint iter = 1; iter <= 1000; iter++) {

      /*** 1. calculate current energy ****/
      double energy = 0.0;
      for (uint v=0; v < nVars; v++)
	energy += solution[v] * cost[v];

      ax.set_constant(0.0);
      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += entry*solution[col]; 
      }      
      for (uint c=0; c < nConstraints; c++) {
	double temp = ax[c] - rhs[c];

	energy += lagrange_multiplier[c] *  temp;
	energy += 0.5*penalty * temp*temp;
      }

      std::cerr << "energy: " << energy << std::endl;

      if ((iter_since_restart % 15)) {

	if (iter_since_restart != 0 && fabs(save_energy-energy) < 1e-6) {

	  std::cerr << "iter converged -> CUTOFF" << std::endl;
	  break;
	}

	save_energy = energy;
      }

      if (iter_since_restart > 7 && energy > last_energy) {
	alpha *= 0.5;

	iter_since_restart = 0;

	std::cerr << "RESTART" << std::endl;
      }
      else 
	iter_since_restart++;

      last_energy = energy;

      /*** 2. calculate gradient ***/

      Math1D::Vector<double> grad(nVars);
      for (uint v=0; v < nVars; v++)
	grad[v] = cost[v];

      ax.set_constant(0.0);

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += solution[col] * entry;
      }

      for (uint c=0; c < nConstraints; c++) 
	ax[c] -= rhs[c];

      ax *= penalty;
      ax += lagrange_multiplier;

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	grad[col] += entry*ax[row];
      }

      /*** 3. go in the direction of the negative gradient  ***/
      for (uint v=0; v < nVars; v++) {

	solution[v] -= alpha*grad[v];
      }

      /*** 4. reproject to the convex set defined by the simplices and variable bounds ***/
      for (uint s=0; s < nSimplices; s++) {
	
	uint start = simplex_starts[s];
	uint nPoints = simplex_starts[s+1] - start;

	projection_on_simplex(solution+start, nPoints);
      }

      for (uint v=0; v < nVars; v++) {

	if (solution[v] < var_lb[v])
	  solution[v] = var_lb[v];
	else if (solution[v] > var_ub[v])
	  solution[v] = var_ub[v];
      }

    } //end of inner iterations


    /*** update lagrange multipliers ****/
#if 1
    ax.set_constant(0.0);
    for (uint k=0; k < matrix_descr.nEntries(); k++) {
      
      const uint row = matrix_descr.row_indices()[k];
      const uint col = matrix_descr.col_indices()[k];
      const double entry = matrix_descr.value()[k];
      
      ax[row] += entry*solution[col]; 
    }      
    for (uint c=0; c < nConstraints; c++) {

      double temp = ax[c] - rhs[c];
      
      lagrange_multiplier[c] += penalty * temp;
    }
#endif
    
  }
}



