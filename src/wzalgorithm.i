%module wzalgorithm

%{
    #include "REvol.h"
    #include "Result.h"
    #include "ParticleSwarmOptimization.h"

    using namespace std;
    using namespace wzalgorithm;
%}


%feature("flatnested");

%include <std_vector.i>

namespace std {
  %template(VecDouble) vector<double>;
}

%include "config.h"
%include "Result.h"

%rename(REvolIndividual) wzalgorithm::REvol::Individual;
%ignore *::operator=;
%ignore wzalgorithm::REvol::Individual::Individual(wzalgorithm::REvol::Individual&&);
%include "REvol.h"

%rename(ParticleSwarmOptimizationResult) wzalgorithm::ParticleSwarmOptimization::Result;
%rename(ParticleSwarmOptimizationParticle) wzalgorithm::ParticleSwarmOptimization::Particle;
%include "ParticleSwarmOptimization.h"


%inline %{
#include <Python.h>
#include <iostream>

namespace wzalgorithm {
    struct REvolSuccessPredicate
    {
        PyObject* evaluator_;


        REvolSuccessPredicate(PyObject* evaluatorObject) :
                evaluator_(evaluatorObject)
        {
        }


        ~REvolSuccessPredicate()
        {
            Py_DECREF(evaluator_);
        }


        bool operator()(REvol::Individual& i)
        {
            auto individual = SWIG_NewPointerObj(
                    &i,
                    SWIGTYPE_p_wzalgorithm__REvol__Individual,
                    0 /* do not own */);
            auto tupleArgument = PyTuple_Pack(1, individual);
            auto result = PyObject_CallObject(evaluator_, tupleArgument);
            Py_DECREF(tupleArgument);
            return (PyBool_Check(result) && PyLong_AsLong(result) == 1);
        }
    };
}
%}


%feature("director") wzalgorithm::REvolSuccessPredicate;
%extend wzalgorithm::REvol {
   %template(runPredicated) run<REvolSuccessPredicate>;
}


