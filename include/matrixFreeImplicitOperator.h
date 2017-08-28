#ifndef MATRIXFREEIMPLICITOPERATOR_H
#define MATRIXFREEIMPLICITOPERATOR_H



using namespace dealii;

//Matrix-free implementation for Allen-Cahn
template <int dim, int degree, typename number>
class ACOperator : public Subscriptor
{
    #include "typeDefs.h"


    public:
        // Temporary hard-coded values for testing purposes
        double dt = 1.0e-1;
        double lambda = 2.0;
        double mobility = 1.0;



    ACOperator (userInputParameters<dim>, MatrixFree<dim,number>);
    userInputParameters<dim> userInputs;

    void clear();
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
ACOperator<dim,degree,number>::ACOperator (userInputParameters<dim> _userInputs, MatrixFree<dim,number> _data):Subscriptor(), userInputs(_userInputs), data(_data)
{}

//Matrix free data structure size
template <int dim, int degree, typename number>
unsigned int ACOperator<dim,degree,number>::m () const
{
    return data.get_vector_partitioner()->size();
}
template <int dim, int degree, typename number>
unsigned int ACOperator<dim,degree,number>::n () const
{
    return data.get_vector_partitioner()->size();
}
//Reset Matrix free data structure
template <int dim, int degree, typename number>
void ACOperator<dim,degree,number>::clear ()
{
    data.clear();
}

//Implement finite element operator application
template <int dim, int degree, typename number>
void ACOperator<dim,degree,number>::
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
void ACOperator<dim,degree,number>::vmult (vectorType &dst, const vectorType &src) const
{
    // (1 + mobility*dt*lambda * M^(-1)*K) n

    invMK(dst,src); //dst = M^(-1)*K*src
    dst*=(mobility*dt*lambda); //dst* = (M*dt*Lambda)
    dst+=src;
}

template <int dim, int degree, typename number>
void ACOperator<dim,degree,number>::invMK (vectorType &dst, const vectorType &src) const
{
    dst=0.0;
    //perform w=K*v
    data.cell_loop (&ACOperator::local_apply, this, dst, src);
    //perform q=M^-1*w...to give q=M^-1*K*v
    for (unsigned int k=0; k<invM.local_size(); ++k)
      dst.local_element(k)*=invM.local_element(k);
}

//datastructure memory consumption estimation
template <int dim, int degree, typename number>
std::size_t ACOperator<dim,degree,number>::memory_consumption () const
{
    return (data.memory_consumption ());
}

// -----------------------------------------------------------------------------------------------
//Matrix-free implementation for Cahn-Hilliard
template <int dim, int degree, typename number>
class CHOperator : public Subscriptor
{
    #include "typeDefs.h"


    public:
        // Temporary hard-coded values for testing purposes
        double dt = 1.0e0;
        double lambda = 1.5;
        double mobility = 1.0;



    CHOperator (userInputParameters<dim>, MatrixFree<dim,number>);
    userInputParameters<dim> userInputs;

    void clear();
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
