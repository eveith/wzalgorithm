#ifndef WZALGORITHM_RESULT_H
#define WZALGORITHM_RESULT_H


#include "config.h"


namespace wzalgorithm {


    /*!
     * \brief A general result structure for optimization algorithms
     *
     * Although the wzalgorithm namespace contains more than one algorithm,
     * the overall goal remains the same: to minimize an error value and
     * to return the solution that calculates to the corresponding error
     * value.
     *
     * The Result structure serves to return exactly that:
     *
     * 1. the final error value
     * 2. the parameters to obtain this error value.
     */
    struct Result
    {
        double error;
        vector_t parameters;
    };


} // namespace wzalgorithm

#endif // WZALGORITHM_RESULT_H
