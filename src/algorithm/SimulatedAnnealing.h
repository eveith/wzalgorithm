#ifndef WINZENT_ALGORITHM_SIMULATEDANNEALING_H
#define WINZENT_ALGORITHM_SIMULATEDANNEALING_H


#include <cstddef>
#include <functional>

#include <boost/random.hpp>

#include "agent_global.h"


namespace Winzent {
    namespace Algorithm {

        class AGENTSHARED_EXPORT SimulatedAnnealing
        {
        public:


            typedef std::function<qreal(QVector<qreal>)> Evaluator;


            //! Constructs a new instance of the Simulated Annealing optimizer
            SimulatedAnnealing();


            /*!
             * \brief Retrieves the maximum number of iterations the algorithm
             *  will run for
             *
             * \return The maximum number of iterations
             */
            size_t maxIterations() const;


            /*!
             * \brief Sets the maximum number of iterations the algorithm will
             *  run for.
             *
             * \param[in] iterations The maximum number of iterations
             *
             * \return `*this`
             */
            SimulatedAnnealing &maxIterations(const size_t &iterations);


            /*!
             * \brief Computes the temperature at the current state
             *
             * \param[in] timeFraction Fraction of runtime that has been
             *  expended so far, e.g., `k / kmax`
             *
             * \return The temperature
             */
            qreal temperature(const qreal &timeFraction) const;


            /*!
             * \brief Generates a neighbor of the given state at the given
             *  temperature
             *
             * \param[in] state The state for which a neighbor should be
             *  generated
             *
             * \param[in] temperature The current temperature
             *
             * \return A neighbor of `state`
             */
            QVector<qreal> neighbor(
                    const QVector<qreal> &state,
                    const qreal &temperature) const;


            /*!
             * \brief Applies the Simulated Annealing algorithm to a problem
             *  represented by the evaluator function
             *
             * \param[in] origin The initial state
             *
             * \param[in] evaluator The evaluation ("fitness") function that
             *  calculates the fitness of a given state.
             *
             * \return The optimized state
             */
            QVector<qreal> run(
                    const QVector<qreal> &origin,
                    Evaluator &evaluator);


        private:


            size_t m_maxIterations;


            //! Random Number Generator
            boost::random::mt11213b m_randomNumberGenerator;


            //! An uniform distribution `[0; 1)` from which we draw
            boost::random::uniform_01<qreal> m_uniformDistribution;
        };

    } // namespace Algorithm
} // namespace Winzent

#endif // WINZENT_ALGORITHM_SIMULATEDANNEALING_H
