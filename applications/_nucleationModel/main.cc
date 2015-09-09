//Coupled Cahn-Hilliard and Allen-Cahn implementation with nucleation
//general headers
#include "../../include/dealIIheaders.h"

//coupled Cahn-Hilliard and Allen-Cahn problem headers
#include "parameters.h"
#include "../../src/models/diffusion/coupledCHAC.h"

//structure representing each nucleus 
struct nucleus{
  unsigned int index;
  dealii::Point<problemDIM> center;
  double radius;
  double seededTime, seedingTime;
};
//vector of all nucleus seeded in the problem
std::vector<nucleus> nuclei;

//initial condition function for concentration
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
    //return the value of the initial concentration field at point p 
    return 0.03 + 1.0e-3*(2*(0.5 - (double)(std::rand() % 100 )/100.0));
  }
};

//apply initial conditions
template <int dim>
void CoupledCHACProblem<dim>::applyInitialConditions()
{
  unsigned int fieldIndex;
  //call initial condition function for c
  fieldIndex=this->getFieldIndex("c");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    InitialConditionC<dim>(),			\
			    *this->solutionSet[fieldIndex]);
  //call initial condition function for n
  fieldIndex=this->getFieldIndex("n");
  VectorTools::interpolate (*this->dofHandlersSet[fieldIndex],		\
			    ZeroFunction<dim>(1),			\
			    *this->solutionSet[fieldIndex]);
}

//nucleation model implementation
template <int dim>
void CoupledCHACProblem<dim>::modifySolutionFields()
{
  //current time
  double t=this->currentTime;
  unsigned int inc=this->currentIncrement;
  double dx=spanX/std::pow(2.0,refineFactor);

  //nucleation parameters
  double minDistBetwenNuclei=spanX/2.0;
  unsigned int maxNumberNuclei=5;

  //get the list of node points in the domain
  std::map<dealii::types::global_dof_index, dealii::Point<dim> > support_points;
  dealii::DoFTools::map_dofs_to_support_points (dealii::MappingQ1<dim>(), *this->dofHandlersSet[0], support_points);   
  //fields
  vectorType* n=this->solutionSet[this->getFieldIndex("n")];
  vectorType* c=this->solutionSet[this->getFieldIndex("c")];

  //populate nuclei vector
  unsigned int index=-1;
  if (inc==1){
    nucleus* temp;
    //first nucleus
    /*
    temp = new nucleus;
    temp->index=++index; 
    temp->center=dealii::Point<dim>(spanX/2.0,spanY/2.0);
    temp->radius=spanX/32.0;
    temp->seededTime=t;
    temp->seedingTime=this->finalTime/10.0;
    nuclei.push_back(*temp);
    */
    //add nuclei based on concentration field values
    //loop over all points in the domain
    for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
      unsigned int dof=it->first;
      //set only local owned values of the parallel vector
      if (n->locally_owned_elements().is_element(dof)){
	dealii::Point<dim> nodePoint=it->second;
	double nValue=(*n)(dof);
	double cValue=(*c)(dof);
	if ((cValue>0.0305) and (nuclei.size()<maxNumberNuclei)){
	  //loop over all existing nuclei to check if they are in the vicinity
	  bool isClose=false;
	  for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
	    if (thisNuclei->center.distance(nodePoint)<minDistBetwenNuclei){
	      isClose=true;
	    }
	  }
	  if (!isClose){
	    temp = new nucleus;
	    temp->index=++index; 
	    temp->center=nodePoint;
	    temp->radius=spanX/32.0;
	    temp->seededTime=t;
	    temp->seedingTime=this->finalTime;
	    nuclei.push_back(*temp);
	  }
	}
      }
    }
    this->pcout << "number of nuclei seeded: " << nuclei.size() << std::endl;
  }

  //seed nuclei
  unsigned int fieldIndex=this->getFieldIndex("n");
  for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
    //
    dealii::Point<dim> center=thisNuclei->center;
    double radius=thisNuclei->radius;
    double seededTime=thisNuclei->seededTime;
    double seedingTime=thisNuclei->seedingTime;
    //loop over all points in the domain
    for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
      unsigned int dof=it->first;
      //set only local owned values of the parallel vector
      if (n->locally_owned_elements().is_element(dof)){
	dealii::Point<dim> nodePoint=it->second;
	//check conditions and seed nuclei
	double r=nodePoint.distance(center);
	if (r<=(radius+8*dx)){
	  if ((t>seededTime) && (t<(seededTime+seedingTime))){
	    (*n)(dof)=0.5*(1.0-std::tanh((r-spanX/16.0)/(3*dx)));
	  }
	}
      }
    }
  }
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      CoupledCHACProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "n"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
      problem.init (); 
      problem.solve();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  return 0;
}
