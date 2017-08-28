//solveIncrement() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"
#include "../../include/varTypeEnums.h"
#include "../../include/matrixFreeImplicitOperator.h"


//solve each time increment
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::solveIncrement(){

    //log time
    computing_timer.enter_section("matrixFreePDE: solveIncrements");
    Timer time;
    char buffer[200];

    //compute residual vectors
    computeRHS();

    //solve for each field
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
        currentFieldIndex = fieldIndex; // Used in computeLHS()

        // Add Neumann BC terms to the residual vector for the current field, if appropriate
        // Currently commented out because it isn't working yet
        //applyNeumannBCs();

        // --------------------------------------------------------------------
        // Explicit parabolic (first order derivatives in time) fields
        // --------------------------------------------------------------------
        if (fields[fieldIndex].pdetype==PARABOLIC){

            // Explicit-time step each DOF
            // Takes advantage of knowledge that the length of solutionSet and residualSet is an integer multiple of the length of invM for vector variables
            unsigned int invM_size = invM.local_size();
            for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
                solutionSet[fieldIndex]->local_element(dof)=			\
                invM.local_element(dof%invM_size)*residualSet[fieldIndex]->local_element(dof);
            }

            // Set the Dirichelet values (hanging node constraints don't need to be distributed every time step, only at output)
            constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
            solutionSet[fieldIndex]->update_ghost_values();

            // Print update to screen
            if (currentIncrement%userInputs.skip_print_steps==0){
                sprintf(buffer, "field '%2s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
                fields[fieldIndex].name.c_str(),				\
                solutionSet[fieldIndex]->l2_norm(),			\
                residualSet[fieldIndex]->l2_norm());
                pcout<<buffer;
            }
        }

        // --------------------------------------------------------------------
        // Implicit parabolic (first order derivatives in time) fields
        // --------------------------------------------------------------------
        else if (fields[fieldIndex].pdetype==IMPLICIT_PARABOLIC){

            // Initially, I'm not going to worry about Dirichlet BCs, just trying to get the right matrix to be solved

            // Allen-Cahn
            // some hard-coded constants for testing


            ACOperator<dim,degree,double> system_matrix(userInputs,matrixFreeObject);

            system_matrix.invM = invM;

            vectorType X, Mu;

            Mu.reinit (invM);
            X.reinit  (Mu);

            system_matrix.X = &X;

            //Begin solve
            //compute fn(n0)+lambda*M^(-1)*K*n0

            system_matrix.invMK(X,*solutionSet[fieldIndex]); //M^(-1)*K*n0
            double n0;
            for (unsigned int k=0; k<solutionSet[fieldIndex]->local_size(); ++k){
              n0=solutionSet[fieldIndex]->local_element(k);
              double fnV = 4.0*n0*(n0-1.0)*(n0-0.5);

              X.local_element(k) = n0 - system_matrix.dt * system_matrix.mobility * fnV; // n0 - dt*L*fn(n0)

              //X.local_element(k)=fnV+ system_matrix.lambda*X.local_element(k); //fn(n0)+lambda*M^(-1)*K*n0
            }

            //(1 + mobility*dt*lambda*M^(-1)*K*M^(-1)*K)Mu=f(c0)+lambda*M^(-1)*K*c0

            //solver controls
            double tol_value;
            if (userInputs.abs_tol == true){
                tol_value = userInputs.solver_tolerance;
            }
            else {
                tol_value = userInputs.solver_tolerance*residualSet[fieldIndex]->l2_norm();
            }
            SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

            // Perform the actual conjugate gradient solve
            SolverCG<vectorType>              solver (solver_control);
            solver.solve(system_matrix,*solutionSet[fieldIndex],X,PreconditionIdentity());

            sprintf(buffer, "field '%2s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e, solution: %12.6e\n", \
            fields[fieldIndex].name.c_str(),			\
            residualSet[fieldIndex]->l2_norm(),			\
            solver_control.last_value(),				\
            solver_control.last_step(), solver_control.tolerance(), solutionSet[fieldIndex]->l2_norm());
            pcout<<buffer;

            // Cahn-Hilliard, modified from Shiva's code from the first CHiMaD Hackathon
            // some hard-coded constants for testing
            /*
            CHOperator<dim,degree,double> system_matrix(userInputs,matrixFreeObject);

            system_matrix.invM = invM;

            vectorType X, Mu;

            Mu.reinit (invM);
            X.reinit  (Mu);

            system_matrix.X = &X;

            //Begin solve
            //compute fc(c0)+lambda*M^(-1)*K*c0

            system_matrix.invMK(X,*solutionSet[fieldIndex]); //M^(-1)*K*c0
            double c0;
            for (unsigned int k=0; k<solutionSet[fieldIndex]->local_size(); ++k){
              c0=solutionSet[fieldIndex]->local_element(k);
              double fcV = 4.0*c0*(c0-1.0)*(c0-0.5);
              X.local_element(k)=fcV+ system_matrix.lambda*X.local_element(k); //f(c0)+lambda*M^(-1)*K*c0
            }

            //(1 + mobility*dt*lambda*M^(-1)*K*M^(-1)*K)Mu=f(c0)+lambda*M^(-1)*K*c0

            //solver controls
            double tol_value;
            if (userInputs.abs_tol == true){
                tol_value = userInputs.solver_tolerance;
            }
            else {
                tol_value = userInputs.solver_tolerance*residualSet[fieldIndex]->l2_norm();
            }
            SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

            // Perform the actual conjugate gradient solve
            SolverCG<vectorType>              solver (solver_control);
            solver.solve(system_matrix,Mu,X,PreconditionIdentity());


            //c=c0-mobility*dt*M^(-1)*K*mu
            system_matrix.invMK(X,Mu);
            X*=(system_matrix.mobility*system_matrix.dt);
            *solutionSet[fieldIndex]-=X;

            */

        }

        // --------------------------------------------------------------------
        //Elliptic (time-independent) fields
        // --------------------------------------------------------------------
        else if (fields[fieldIndex].pdetype==ELLIPTIC){

            //implicit solve
            //apply Dirichlet BC's
            // Loops through all DoF to which ones have Dirichlet BCs applied, replace the ones that do with the Dirichlet value
            // Is this needed? Why are we applying BCs to the residualSet?
            // This clears the residual where we want to apply Dirichlet BCs, otherwise the solver sees a positive residual
            for (std::map<types::global_dof_index, double>::const_iterator it=valuesDirichletSet[fieldIndex]->begin(); it!=valuesDirichletSet[fieldIndex]->end(); ++it){
                if (residualSet[fieldIndex]->in_local_range(it->first)){
                    (*residualSet[fieldIndex])(it->first) = 0.0; //it->second; //*jacobianDiagonal(it->first);
                }
            }

            //solver controls
            double tol_value;
            if (userInputs.abs_tol == true){
                tol_value = userInputs.solver_tolerance;
            }
            else {
                tol_value = userInputs.solver_tolerance*residualSet[fieldIndex]->l2_norm();
            }

            SolverControl solver_control(userInputs.max_solver_iterations, tol_value);

            // Currently the only allowed solver is SolverCG, the SolverType input variable is a dummy
            SolverCG<vectorType> solver(solver_control);

            //solve
            try{
                if (fields[fieldIndex].type == SCALAR){
                    dU_scalar=0.0;
                    solver.solve(*this, dU_scalar, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
                }
                else {
                    dU_vector=0.0;
                    solver.solve(*this, dU_vector, *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
                }
            }
            catch (...) {
                pcout << "\nWarning: implicit solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing solverTolerance.\n";
            }
            if (fields[fieldIndex].type == SCALAR){
                *solutionSet[fieldIndex]+=dU_scalar;
            }
            else {
                *solutionSet[fieldIndex]+=dU_vector;
            }

            if (currentIncrement%userInputs.skip_print_steps==0){
                double dU_norm;
                if (fields[fieldIndex].type == SCALAR){
                    dU_norm = dU_scalar.l2_norm();
                }
                else {
                    dU_norm = dU_vector.l2_norm();
                }
                sprintf(buffer, "field '%2s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e, solution: %12.6e, dU: %12.6e\n", \
                fields[fieldIndex].name.c_str(),			\
                residualSet[fieldIndex]->l2_norm(),			\
                solver_control.last_value(),				\
                solver_control.last_step(), solver_control.tolerance(), solutionSet[fieldIndex]->l2_norm(), dU_norm);
                pcout<<buffer;
            }

        }

        //Hyperbolic (second order derivatives in time) fields and general
        //non-linear PDE types not yet implemented
        else{
            pcout << "PRISMS-PF Error: solveIncrement.cc: unknown equation type, given pdetype is " << fields[fieldIndex].pdetype << "\n";
            exit(-1);
        }

        //check if solution is nan
        if (!numbers::is_finite(solutionSet[fieldIndex]->l2_norm())){
            sprintf(buffer, "PRISMS-PF ERROR: field '%s' solution is NAN. exiting.\n\n",
            fields[fieldIndex].name.c_str());
            pcout<<buffer;
            exit(-1);
        }
    }
    if (currentIncrement%userInputs.skip_print_steps==0){
        pcout << "wall time: " << time.wall_time() << "s\n";
    }
    //log time
    computing_timer.exit_section("matrixFreePDE: solveIncrements");

}

#include "../../include/matrixFreePDE_template_instantiations.h"
