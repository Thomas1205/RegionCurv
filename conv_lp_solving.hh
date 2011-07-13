/**** written by Thomas Schoenemann as an employee of Lund University, August 2010 ****/

#ifndef CONV_LP_SOLVING_HH
#define CONV_LP_SOLVING_HH

#include "sparse_matrix_description.hh"
#include "projection.hh"

void eq_constrained_lp_solving_auglagrange(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
					   const SparseMatrixDescription<double> matrix_descr, const double* rhs,  
					   double* solution);


template<typename T>
void eq_constrained_lp_solving_auglagrange_nesterov(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
						    const SparseMatrixDescription<T> matrix_descr, const double* rhs,  
						    double* solution, double start_penalty = 100.0, double stepsize_coeff = 1.0,
						    uint inner_iter = 1000, uint outer_iter = 15, double penalty_factor = 1.25, int nUsableEntries = -1);

//NOTE: simplices MAY NOT overlap
void eq_and_simplex_constrained_lp_solve_auglagr(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
						 const SparseMatrixDescription<double> matrix_descr, const double* rhs,  
						 uint nSimplices, const uint* simplex_starts, double* solution);

//NOTE: simplices MAY NOT overlap
template<typename T>
void eq_and_simplex_constrained_lp_solve_auglagr_nesterov(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
							  const SparseMatrixDescription<T> matrix_descr, const double* rhs,  
							  uint nSimplices, const uint* simplex_starts, double* solution,
							  double start_penalty = 100.0, double stepsize_coeff = 1.0, 
							  uint inner_iter = 1000, uint outer_iter = 15, double penalty_factor = 1.25);




/************************************************* implementation *************************************/

template<typename T>
void eq_constrained_lp_solving_auglagrange_nesterov(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
						    const SparseMatrixDescription<T> matrix_descr, const double* rhs,  
						    double* solution, double start_penalty, double stepsize_coeff, 
						    uint nInnerIter, uint nOuterIter, double penalty_factor, int nUsableEntries) {

#ifdef USE_CUDA
  //#if 0
  //NOTE: here we are not checking if the matrix is really sorted
  Math1D::Vector<uint> row_start(nConstraints+1,MAX_UINT);
  for (uint e=0; e < matrix_descr.nEntries(); e++) {

    uint r = matrix_descr.row_indices()[e];
    row_start[r] = std::min(row_start[r],e);
  }
  row_start[nConstraints] = matrix_descr.nEntries();
  
  cuda_eq_constrained_lp_solving_auglagrange_nesterov(nVars, nConstraints, cost, var_lb, var_ub, matrix_descr.value(), matrix_descr.col_indices(), 
						      row_start.direct_access(), rhs, solution, start_penalty, nInnerIter, nOuterIter, stepsize_coeff,
						      penalty_factor);

  return;
#endif

  Math1D::NamedVector<double> aux_solution(nVars,MAKENAME(aux_solution));

#ifdef USE_EXPLICIT_GRADIENT
  Math1D::Vector<double> grad(nVars);
#endif

  Math1D::NamedVector<double> lagrange_multiplier(nConstraints,0.0,MAKENAME(lagrange_multiplier));

  double penalty = start_penalty;
  //double alpha = 1e-1 / penalty;

  Math1D::Vector<double> ax(nConstraints,0.0);

  uint nMatrixEntries = matrix_descr.nEntries();
  if (nUsableEntries > 0)
    nMatrixEntries = (uint) nUsableEntries;

  for (uint outer_iter = 1; outer_iter <= nOuterIter; outer_iter++) {

    if (outer_iter != 1) {
      penalty *= penalty_factor;
    }

    double alpha = 1e-1 * stepsize_coeff / penalty;

    std::cerr << "############# iteration " << outer_iter << ", penalty " << penalty << std::endl;

    double last_energy = 1e50;

    double prev_t = 1.0;

    for (uint v=0; v < nVars; v++)      
      aux_solution[v] = solution[v];

    uint iter_since_restart = 0;

    double saved_energy = 1e50;

    for (uint iter = 1; iter <= nInnerIter; iter++) {

      /*** 1. calculate current energy ****/
      double energy = 0.0;
      for (uint v=0; v < nVars; v++)
	energy += solution[v] * cost[v];

      ax.set_constant(0.0);
      for (uint k=0; k < nMatrixEntries; k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += entry*solution[col]; 
      }      
      //double infeas = 0.0;
      for (uint c=0; c < nConstraints; c++) {
	const double temp = ax[c] - rhs[c];

	energy += lagrange_multiplier[c] *  temp;
	energy += 0.5*penalty * temp*temp;
	//infeas += temp*temp;
      }

      std::cerr << "iter " << iter << ",energy: ";
      std::cerr.precision(12);
      std::cerr << energy << std::endl;
      //std::cerr << "infeas: " << infeas << std::endl;

      if ((iter_since_restart % 15) == 0) {

	
	if (fabs(saved_energy - energy) < 1e-6) {

	  std::cerr << "inner iteration converged. -> CUTOFF" << std::endl;
	  break;
	}	  

	saved_energy = energy;
      }

      //std::cerr << "last:" << last_energy << std::endl;
#if 1
      if (energy > pow(1.25,sign(last_energy))*last_energy || (energy > last_energy && iter_since_restart > 5)) {
 	alpha *= 0.5;
	iter_since_restart = 0;

	prev_t = 1.0;

	std::cerr << "RESTART" << std::endl;
	
	for (uint v=0; v < nVars; v++)
	  aux_solution[v] = solution[v];
      }
      else
	iter_since_restart++;
#endif

      last_energy = energy;

      /*** 2. calculate gradient ***/

#ifdef USE_EXPLICIT_GRADIENT
      for (uint v=0; v < nVars; v++)
	grad[v] = cost[v];
#endif

      ax.set_constant(0.0);

      for (uint k=0; k < nMatrixEntries; k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += aux_solution[col] * entry;
      }

      for (uint c=0; c < nConstraints; c++)
	ax[c] -= rhs[c];
      ax *= penalty;

      ax += lagrange_multiplier;

#ifndef USE_EXPLICIT_GRADIENT
      for (uint v=0; v < nVars; v++)
	aux_solution[v] -= alpha*cost[v];
#endif

      for (uint k=0; k < nMatrixEntries; k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

#ifdef USE_EXPLICIT_GRADIENT
	grad[col] += entry*ax[row];
#else
	aux_solution[col] -= alpha*entry*ax[row];
#endif
      }

      /*** 3. go in the direction of the negative gradient  ***/
#ifdef USE_EXPLICIT_GRADIENT
      for (uint v=0; v < nVars; v++) {

	aux_solution[v] -= alpha*grad[v];
      }
#endif
      
      /*** 4. reproject to the convex set defined by the variable bounds ***/
      for (uint v=0; v < nVars; v++) {

	if (aux_solution[v] < var_lb[v])
	  aux_solution[v] = var_lb[v];
	else if (aux_solution[v] > var_ub[v])
	  aux_solution[v] = var_ub[v];
      }

      /*** 5. update variables according to Nesterov scheme ***/
      const double new_t = 0.5 * (1 + sqrt(1+4*prev_t*prev_t));
      const double nesterov_fac = (prev_t - 1) / new_t;
      
      for (uint i=0; i < nVars; i++) {
	
	const double old_aux = aux_solution.direct_access(i);
	aux_solution.direct_access(i) = old_aux + nesterov_fac*(old_aux - solution[i]) ;
	solution[i] = old_aux;
      }
      
      prev_t = new_t;

    } //end of inner iterations


    double lp_energy = 0.0;
    for (uint v=0; v < nVars; v++)
      lp_energy += solution[v]*cost[v];

    std::cerr.precision(12);
    std::cerr << "lp-energy: " << lp_energy << std::endl;

    double penalty_term = 0.0;

    /*** update lagrange multipliers ****/
#if 1
    ax.set_constant(0.0);
    for (uint k=0; k < nMatrixEntries; k++) {
      
      const uint row = matrix_descr.row_indices()[k];
      const uint col = matrix_descr.col_indices()[k];
      const double entry = matrix_descr.value()[k];
      
      ax[row] += entry*solution[col]; 
    }      
    for (uint c=0; c < nConstraints; c++) {

      const double temp = ax[c] - rhs[c];
      
      lagrange_multiplier[c] += penalty * temp;

      penalty_term += 0.5*penalty*temp*temp;
    }
#endif

    std::cerr << "penalty term: " << penalty_term << std::endl;
  }
}



template<typename T>
void eq_and_simplex_constrained_lp_solve_auglagr_nesterov(uint nVars, uint nConstraints, const double* cost, const double* var_lb, const double* var_ub,
							  const SparseMatrixDescription<T> matrix_descr, const double* rhs,  
							  uint nSimplices, const uint* simplex_starts, double* solution,
							  double start_penalty, double stepsize_coeff, uint nInnerIter, uint nOuterIter,
							  double penalty_factor) {

  //TODO: check if the initial <code> solution </code> satisfies the simplex constraints. If not, do a projection


#ifdef USE_CUDA
  //NOTE: here we are not checking if the matrix is really sorted
  //      and that the simplices are all of same width
  Math1D::Vector<uint> row_start(nConstraints+1,MAX_UINT);
  for (uint e=0; e < matrix_descr.nEntries(); e++) {

    uint r = matrix_descr.row_indices()[e];
    row_start[r] = std::min(row_start[r],e);
  }
  row_start[nConstraints] = matrix_descr.nEntries();
  
  cuda_eq_and_simplex_constr_lp_solving_auglag_nesterov(nVars, nConstraints, cost, nSimplices, simplex_starts[1] - simplex_starts[0],
							var_lb, var_ub, matrix_descr.value(), matrix_descr.col_indices(), 
							row_start.direct_access(), rhs, solution, start_penalty, nInnerIter, nOuterIter, 
							stepsize_coeff, penalty_factor);

  return;
#endif

  Math1D::NamedVector<double> aux_solution(nVars,MAKENAME(aux_solution));

  Math1D::NamedVector<double> lagrange_multiplier(nConstraints,0.0,MAKENAME(lagrange_multiplier));

  double penalty = start_penalty;

  //double alpha = 1e-1 / penalty;

  Math1D::Vector<double> ax(nConstraints,0.0);

  for (uint outer_iter = 1; outer_iter <= nOuterIter; outer_iter++) {

    double best_energy = 1e300;

    if (outer_iter != 1) {
      penalty *= penalty_factor;
    }

    double alpha = 1e-1 * stepsize_coeff / penalty;

    std::cerr << "############# outer iteration " << outer_iter << ", penalty " << penalty << std::endl;

    double last_energy = 1e50;

    double prev_t = 1.0;

    for (uint v=0; v < nVars; v++)      
      aux_solution[v] = solution[v];

    uint iter_since_restart = 0;

    double save_energy = 1e50;

    double energy_landmark = 1e50;

    const double cutoff = (outer_iter == nOuterIter) ? 0.5e-7: 1e-5;

    for (uint iter = 1; iter <= nInnerIter; iter++) {

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

	const double cur_rhs = (rhs != 0) ? rhs[c] : 0.0;
	const double temp = ax[c] - cur_rhs;

	energy += lagrange_multiplier[c] * temp;
	energy += 0.5*penalty * temp*temp;
      }

      std::cerr << "iter: " << iter << ", energy: ";
      std::cerr.precision(12);
      std::cerr << energy << std::endl;

      if ((iter % 10) == 1) {
	
	if (fabs(energy_landmark-energy) < cutoff) {
	  std::cerr << "apparently near convergence -> CUTOFF" << std::endl;
	  break;
	}
	energy_landmark = energy;
      }


      if ((iter_since_restart % 15) == 0) {

	if (iter_since_restart >= 15 && fabs(energy - save_energy) < 1e-6) {
	  
	  std::cerr << "iter converged -> CUTOFF" << std::endl;
	  break;
	}

	save_energy = energy;
      }

#if 1
      if ((energy > 1.25*best_energy && iter_since_restart > 1) || (energy > last_energy && iter_since_restart > 5)) {
 	alpha *= 0.5;
	iter_since_restart = 0;

	prev_t = 1.0;

	std::cerr << "RESTART" << std::endl;
	
	for (uint v=0; v < nVars; v++)
	  aux_solution[v] = solution[v];
      }
      else
	iter_since_restart++;
#endif

      best_energy = std::min(best_energy,energy);

      last_energy = energy;

      /*** 2. calculate gradient ***/

#ifdef USE_EXPLICIT_GRADIENT
      //std::cerr << "setting up the gradient" << std::endl;
      Math1D::Vector<double> grad(nVars);
      for (uint v=0; v < nVars; v++)
	grad[v] = cost[v];
#endif

      ax.set_constant(0.0);

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	ax[row] += aux_solution[col] * entry;
      }

      if (rhs != 0) {
	for (uint c=0; c < nConstraints; c++)
	  ax[c] -= rhs[c];
      }

      ax *= penalty;

      ax += lagrange_multiplier;

#ifndef USE_EXPLICIT_GRADIENT
      for (uint v=0; v < nVars; v++)
	aux_solution[v] -= alpha * cost[v];
#endif

      for (uint k=0; k < matrix_descr.nEntries(); k++) {

	const uint row = matrix_descr.row_indices()[k];
	const uint col = matrix_descr.col_indices()[k];
	const double entry = matrix_descr.value()[k];

	//NOTE: the multiplication with <code> penalty </code> was already included in <code> ax </code>
#ifdef USE_EXPLICIT_GRADIENT
	grad[col] += entry*ax[row];
#else
	aux_solution[col] -= alpha * entry*ax[row];
#endif
      }

      /*** 3. go in the direction of the negative gradient  ***/
#ifdef USE_EXPLICIT_GRADIENT
      for (uint v=0; v < nVars; v++) {

	aux_solution[v] -= alpha*grad[v];
      }
#endif

      /*** 4. reproject to the convex set defined by the simplices and the variable bounds ***/
      for (uint s=0; s < nSimplices; s++) {
	
	const uint start = simplex_starts[s];
	const uint nPoints = simplex_starts[s+1] - start;

	projection_on_simplex(aux_solution.direct_access()+start, nPoints);
      }

      for (uint v=0; v < nVars; v++) {

	double lb = (var_lb != 0) ? var_lb[v] : 0.0;
	double ub = (var_ub != 0) ? var_ub[v] : 1.0;

	if (aux_solution[v] < lb)
	  aux_solution[v] = lb;
	else if (aux_solution[v] > ub)
	  aux_solution[v] = ub;
      }

      /*** 5. update variables according to Nesterov scheme ***/
      const double new_t = 0.5 * (1 + sqrt(1+4*prev_t*prev_t));
      const double nesterov_fac = (prev_t - 1) / new_t;
      
      for (uint i=0; i < nVars; i++) {
	
	const double old_aux = aux_solution.direct_access(i);
	aux_solution.direct_access(i) = old_aux + nesterov_fac*(old_aux - solution[i]) ;
	solution[i] = old_aux;
      }
      
      prev_t = new_t;

    } //end of inner iterations

    double lp_energy = 0.0;
    for (uint v=0; v < nVars; v++)
      lp_energy += solution[v]*cost[v];

    std::cerr.precision(12);
    std::cerr << "lp-energy: " << lp_energy << std::endl;

    double penalty_term = 0.0;

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

      double cur_rhs = (rhs != 0) ? rhs[c] : 0.0; 

      const double temp = ax[c] - cur_rhs;
      
      lagrange_multiplier[c] += penalty * temp;

      penalty_term += 0.5*penalty*temp*temp;
    }
#endif
    
    std::cerr << "penalty term: " << penalty_term << std::endl;
  }

}







#endif
