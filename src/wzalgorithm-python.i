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
%ignore wzalgorithm::REvol::Individual::Individual(wzalgorithm::REvol::Individual &&);
%include "REvol.h"

%rename(ParticleSwarmOptimizationResult) wzalgorithm::ParticleSwarmOptimization::Result;
%rename(ParticleSwarmOptimizationParticle) wzalgorithm::ParticleSwarmOptimization::Particle;
%include "ParticleSwarmOptimization.h"


%inline %{
#include <Python.h>
#include <stdexcept>
#include <iostream>

namespace wzalgorithm {
    struct REvolSuccessPredicate
    {
        REvolSuccessPredicate(PyObject* evaluatorObject) :
                evaluator_(evaluatorObject)
        {
            assert(nullptr != evaluator_);
            Py_INCREF(evaluator_);
        }

        ~REvolSuccessPredicate()
        {
        }

        bool operator()(REvol::Individual& i)
        {
            auto individual = SWIG_NewPointerObj(
                    &i,
                    SWIGTYPE_p_wzalgorithm__REvol__Individual,
                    0 /* No Python GC */ | 0 /* Not an opaque pointer */);
            auto tupleArgument = PyTuple_Pack(1, individual);
            auto result = PyObject_CallObject(evaluator_, tupleArgument);
            Py_DECREF(tupleArgument);

            if (!result) {
                throw std::invalid_argument("Python-code evaluator object "
                        "did not return a result!");
                return false;
            }
            return (PyBool_Check(result) && PyLong_AsLong(result) == 1);
        }


    private:
        PyObject* evaluator_;
    };
}
%}


%feature("director") wzalgorithm::REvolSuccessPredicate;
%extend wzalgorithm::REvol {
    %template(run) run<REvolSuccessPredicate>;
}


%feature("shadow") wzalgorithm::REvol::run(REvol::Individual const&, REvolSuccessPredicate&) %{
    def run(self, individual, evaluator):
        print("Hello, World!")
        $action
%}
