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


	// ================================================================
	// Model constants specific to this subclass
	// ================================================================
    
// Matrial Properties constant
    //double ml = userInputs.get_model_constant_double("ml");
    double c0 = userInputs.get_model_constant_double("c0");
    double r0 = userInputs.get_model_constant_double("r0");
    double omega = userInputs.get_model_constant_double("omega");
    double p0 = userInputs.get_model_constant_double("p0");
    //double gamma = userInputs.get_model_constant_double("gamma");
    //double d0 = userInputs.get_model_constant_double("d0");
    double Dl = userInputs.get_model_constant_double("Dl");
    double Ds = userInputs.get_model_constant_double("Ds");
    
//Dimensionless parameters
	
	//double Tl = userInputs.get_model_constant_double("Tl");
	//double G = userInputs.get_model_constant_double("G");
	//double R = userInputs.get_model_constant_double("R");
	//double Vp = userInputs.get_model_constant_double("Vp");
    
// New input
	//double W0 = userInputs.get_model_constant_double("W0");
    //double tau0 = userInputs.get_model_constant_double("tau0");
    
    double epsilon = userInputs.get_model_constant_double("epsilon");
    double k = userInputs.get_model_constant_double("k");
    double lamda = userInputs.get_model_constant_double("lamda");
    
    double Dtilt = userInputs.get_model_constant_double("Dtilt");
    double Vtilt = userInputs.get_model_constant_double("Vtilt");
    double ltilt = userInputs.get_model_constant_double("ltilt");
    

	// ================================================================

};
