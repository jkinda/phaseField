#ifndef MATRIXFREEIMPLICITOPERATOR_H
#define MATRIXFREEIMPLICITOPERATOR_H



using namespace dealii;

//Matrix-free implementation
template <int dim, int degree, typename number>
class CHOperator : public Subscriptor
{
    #include "typeDefs.h"


    public:
        // Temporary hard-coded values for testing purposes
        double dt = 1.0e-1;
        double lambda = 2.0;
        double mobility = 2.0/(0.9);



    CHOperator (userInputParameters<dim>, MatrixFree<dim,number>);
    userInputParameters<dim> userInputs;

    void clear();
    void reinit (const DoFHandler<dim>  &dof_handler, const ConstraintMatrix  &constraints, const unsigned int level = numbers::invalid_unsigned_int);
    unsigned int m () const;
    unsigned int n () const;
    void invMK (vectorType &dst, const vectorType &src) const;
    void vmult (vectorType &dst, const vectorType &src) const;
    number el (const unsigned int row, const unsigned int col) const;
    void set_diagonal (const vectorType &diagonal);
    std::size_t memory_consumption () const;
    vectorType *X;
    vectorType invM;
    private:
    void local_apply (const MatrixFree<dim,number> &data, vectorType &dst,
    	      const vectorType &src,
    	      const std::pair<unsigned int,unsigned int> &cell_range) const;
    MatrixFree<dim,number>      data;



};

//Constructor
template <int dim, int degree, typename number>
CHOperator<dim,degree,number>::CHOperator (userInputParameters<dim> _userInputs, MatrixFree<dim,number> _data):Subscriptor(), userInputs(_userInputs), data(_data)
{}

//Matrix free data structure size
template <int dim, int degree, typename number>
unsigned int CHOperator<dim,degree,number>::m () const
{
    return data.get_vector_partitioner()->size();
}
template <int dim, int degree, typename number>
unsigned int CHOperator<dim,degree,number>::n () const
{
    return data.get_vector_partitioner()->size();
}
//Reset Matrix free data structure
template <int dim, int degree, typename number>
void CHOperator<dim,degree,number>::clear ()
{
    data.clear();
}
//Initialize Matrix Free data structure
template <int dim, int degree, typename number>
void CHOperator<dim,degree,number>::reinit (const DoFHandler<dim>  &dof_handler, const ConstraintMatrix  &constraint, const unsigned int level)
{
    typename MatrixFree<dim,number>::AdditionalData additional_data;
    additional_data.mpi_communicator = MPI_COMM_WORLD;
    additional_data.tasks_parallel_scheme = MatrixFree<dim,number>::AdditionalData::partition_partition;
    additional_data.mapping_update_flags = (update_gradients | update_JxW_values);
    QGaussLobatto<1> quadrature (degree+1);
    data.reinit (dof_handler, constraint, quadrature, additional_data);

    //Compute  invM
    data.initialize_dof_vector (invM);
    VectorizedArray<double> one = make_vectorized_array (1.0);
    //Select gauss lobatto quad points which are suboptimal but give diogonal M
    FEEvaluation<dim,degree> fe_eval(data);
    const unsigned int            n_q_points = fe_eval.n_q_points;
    for (unsigned int cell=0; cell<data.n_macro_cells(); ++cell)
    {
        fe_eval.reinit(cell);
        for (unsigned int q=0; q<n_q_points; ++q)
        fe_eval.submit_value(one,q);
        fe_eval.integrate (true,false);
        fe_eval.distribute_local_to_global (invM);
    }
    invM.compress(VectorOperation::add);
    //
    for (unsigned int k=0; k<invM.local_size(); ++k)
    if (std::abs(invM.local_element(k))>1e-15)
    invM.local_element(k) = 1./invM.local_element(k);
    else
    invM.local_element(k) = 0;
}

//Implement finite element operator application
template <int dim, int degree, typename number>
void CHOperator<dim,degree,number>::
local_apply (const MatrixFree<dim,number>               &data,
           vectorType                                 &dst,
           const vectorType                           &src,
           const std::pair<unsigned int,unsigned int> &cell_range) const
{
    FEEvaluation<dim,degree> mat(data);
    //loop over all "cells"
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
    {
        mat.reinit (cell);
        mat.read_dof_values(src);
        mat.evaluate (false,true,false);
        for (unsigned int q=0; q<mat.n_q_points; ++q){
            mat.submit_gradient(mat.get_gradient(q),q);
        }
        mat.integrate (false,true);
        mat.distribute_local_to_global (dst);
    }
}

// Matrix free data structure vmult operations.
template <int dim, int degree, typename number>
void CHOperator<dim,degree,number>::vmult (vectorType &dst, const vectorType &src) const
{
    invMK(*X,src); //X=M^(-1)*K*src
    invMK(dst,*X); //dst=M^(-1)*K*X
    dst*=(mobility*dt*lambda); //dst*=(M*dt*Lambda)
    dst+=src;
}

template <int dim, int degree, typename number>
void CHOperator<dim,degree,number>::invMK (vectorType &dst, const vectorType &src) const
{
    dst=0.0;
    //perform w=K*v
    data.cell_loop (&CHOperator::local_apply, this, dst, src);
    //perform q=M^-1*w...to give q=M^-1*K*v
    for (unsigned int k=0; k<invM.local_size(); ++k)
      dst.local_element(k)*=invM.local_element(k);
}

//datastructure memory consumption estimation
template <int dim, int degree, typename number>
std::size_t CHOperator<dim,degree,number>::memory_consumption () const
{
    return (data.memory_consumption ());
}



#endif /* MATRIXFREEIMPLICITOPERATOR_H */
