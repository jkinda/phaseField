#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
    // Constructor
    customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};

    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

private:
	#include "../../include/typeDefs.h"

	const userInputParameters<dim> userInputs;

	// Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set the LHS of the governing equations (in equations.h)
	void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

	// Function to set postprocessing expressions (in postprocess.h)
	#ifdef POSTPROCESS_FILE_EXISTS
	void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
					variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
					const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
	#endif

	// Function to set the nucleation probability (in nucleation.h)
	#ifdef NUCLEATION_FILE_EXISTS
	double getNucleationProbability(variableValueContainer variable_value, double dV) const;
	#endif

	// ================================================================
	// Methods specific to this subclass
	// ================================================================

    void applyNeumannBCs();

	// ================================================================
	// Model constants specific to this subclass
	// ================================================================

	double McV = userInputs.get_model_constant_double("McV");
	double KcV = userInputs.get_model_constant_double("KcV");
    double rho = userInputs.get_model_constant_double("rho");
    double c_alpha = userInputs.get_model_constant_double("c_alpha");
    double c_beta= userInputs.get_model_constant_double("c_beta");
    double k = userInputs.get_model_constant_double("k");
    double epsilon = userInputs.get_model_constant_double("epsilon");

	// ================================================================

};

// =================================================================================
// Methods to apply non-zero Neumann BCs
// =================================================================================
template <int dim, int degree>
void customPDE<dim,degree>::applyNeumannBCs(){
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// NOTE: Currently this function doesn't work and it's call is commented out in solveIncrement.
	// The result is off by almost exactly a factor of 100,000. I don't know what the issue is.
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	// Check to the BC for the current field
	unsigned int starting_BC_list_index = 0;
	for (unsigned int i=0; i<this->currentFieldIndex; i++){
		if (userInputs.var_type[i] == SCALAR){
			starting_BC_list_index++;
		}
		else {
			starting_BC_list_index+=dim;
		}
	}

	if (userInputs.var_type[this->currentFieldIndex] == SCALAR){
		for (unsigned int direction = 0; direction < 2*dim; direction++){
			if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] == NEUMANN){

				typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet[0]->begin_active(), endc = this->dofHandlersSet[0]->end();
				FESystem<dim>* fe= this->FESet[this->currentFieldIndex];
				QGaussLobatto<dim-1>  face_quadrature_formula(degree+1);
				FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula, update_values | update_JxW_values | update_quadrature_points);
				const unsigned int n_face_q_points  = face_quadrature_formula.size(), dofs_per_cell=fe->dofs_per_cell;
				Vector<double> cell_rhs(dofs_per_cell);
				std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

				// Loop over each face on a boundary
				for (;cell!=endc; ++cell){
					for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
						if (cell->face(f)->at_boundary()){
							if (cell->face(f)->boundary_id() == direction){
								fe_face_values.reinit (cell, f);
								cell_rhs=0.0;
								for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point){
									double neumann_value = userInputs.BC_list[starting_BC_list_index].var_BC_val[direction];
									for (unsigned int i=0; i<dofs_per_cell; ++i){

                                        dealii::Point<dim, double> p = fe_face_values.quadrature_point(q_point);

                                        double A = 0.0002;
                                        double B = -0.01;
                                        double C = 0.02;

                                        if (direction == 0){
                                            neumann_value = (A*p(1)+B);
                                        }
                                        else if (direction = 1){
                                            neumann_value = -(A*p(1)+B);
                                        }
                                        else if (direction = 2){
                                            neumann_value = (A*p(0)+C);
                                        }
                                        else{
                                            neumann_value = -(A*p(0)+C);
                                        }



										cell_rhs(i) += (neumann_value * fe_face_values.shape_value(i,q_point) * fe_face_values.JxW(q_point));
									}
								}
								cell->get_dof_indices (local_dof_indices);
								//assemble
								for (unsigned int i=0; i<dofs_per_cell; ++i){
									(*(this->residualSet[this->currentFieldIndex]))[local_dof_indices[i]] += cell_rhs(i);
								}
							}
						}
					}
				}
			}

		}
	}

}
