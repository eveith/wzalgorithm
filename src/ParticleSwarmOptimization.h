#ifndef WINZENT_ALGORITHM_PARTICLESWARMOPTIMIZATION_H
#define WINZENT_ALGORITHM_PARTICLESWARMOPTIMIZATION_H


#include <tuple>
#include <cstddef>
#include <ostream>
#include <functional>

#include <boost/random.hpp>

#include "config.h"


namespace Winzent {
    namespace Algorithm {
        namespace detail {


            /*!
             * \brief A particle in the particle swarm
             *
             * The particle is the basic maintenance structure of the
             * Particle Swarm Optimization algorithm. It carries the particles
             * current position, the fitness at the position, its best
             * position and fitness, as well as its velocity.
             *
             * The Particle structure is passed to the user-supplied
             * evaluation function. Although this function can access the
             * data structure as a whole, the typical modus operandi would be
             * to use the particle's current posistion as input parameters to
             * the objective function, and place the function's result in the
             * #currentFitness variable.
             */
            struct Particle
            {
                //! The vector type used here.
                typedef Winzent::Algorithm::vector_t Vector;


                //! The best fitness this particle ever obtained
                double bestFitness;


                /*!
                 * \brief The fitness at the current position, i.e.,
                 *  the result of applying the objective function to the
                 *  current position
                 */
                double currentFitness;


                //! The best position, i.e., the source of the #bestFitness
                Vector bestPosition;


                /*!
                 * \brief The current position, i.e., the parameters to the
                 *  objective function
                 */
                Vector currentPosition;


                //! The particle's velocity
                Vector velocity;


                /*!
                 * \brief The smaller-than operator, applied to the particle's
                 *  fitness.
                 */
                bool operator <(Particle const& rhs) const {
                    return bestFitness < rhs.bestFitness;
                }
            };


            struct ParticleSwarmOptimizationResult
            {
                //! The best particle
                Particle bestParticle;

                //! Number of iterations the algorithm ran
                size_t iterationsUsed;
            };
        }


        class ParticleSwarmOptimization
        {
        public:


            //! Pre-defined default swarm size (cp. SPSO 2007)
            static const size_t DEFAULT_SWARM_SIZE = 40;


            //! The contant C, ~1.193
            static const double C;


            //! The contant W, ~0.721
            static const double W;


            /*!
             * \brief The Evaluator designates the user-supplied fitness
             *  function
             *
             * This functional object represents the objective function that
             * evaluates each particle's fitness. The goal of the PSO is to
             * optimize the parameters in such a way that its result is as low
             * as possible.
             *
             * \param[in] particle The particle that is currently being
             *  evaluated
             *
             * \return `true` if the particle matches the user-determined
             *  success criterion, `false` if the PSO algorithm should
             *  continue
             *
             * \sa detail::Particle
             *
             * \sa #run()
             */
            typedef std::function<bool(detail::Particle&)> Evaluator;


            //! The Swarm data type
            typedef std::vector<detail::Particle> Swarm;


            //! Three neighbors to a particle
            typedef std::tuple<size_t, size_t, size_t> NeighborIndices;


            /*!
             * \brief A particle neighborhood
             */
            typedef std::vector<detail::Particle const*> Neighborhood;


            ParticleSwarmOptimization();


            /*!
             * \brief Gets the swarm size
             *
             * \return The current swarm size
             *
             * \sa DEFAULT_SWARM_SIZE
             */
            size_t swarmSize() const;


            /*!
             * \brief Sets the size of the swarm
             *
             * \param[in] size The new size
             *
             * \return `*this`
             */
            ParticleSwarmOptimization &swarmSize(size_t size);


            /*!
             * \brief Gets the maximum number of iterations the algorithm will
             *  run for.
             *
             * \return The maximum number of iterations
             */
            size_t maxIterations() const;


            /*!
             * \brief Sets the maximum number of iterations the algorithm will
             *  try to optimize the evaluation function's parameters
             *
             * \param[in] iterations The number of iterations
             *
             * \return `*this`
             */
            ParticleSwarmOptimization &maxIterations(size_t iterations);


            /*!
             * \brief Gets the lower boundary of the search space
             *
             * \return The lower boundary
             */
            double lowerBoundary() const;


            /*!
             * \brief Sets the lower boundary
             *
             * \param[in] boundary The new lower boundary
             *
             * \return `*this`
             */
            ParticleSwarmOptimization& lowerBoundary(double boundary);


            /*!
             * \brief Gets the upper boundary of the search space
             *
             * \return The upper boundary
             */
            double upperBoundary() const;


            /*!
             * \brief Sets the upper boundary
             *
             * \param[in] boundary The new upper boundary
             *
             * \return `*this`
             */
            ParticleSwarmOptimization& upperBoundary(double boundary);


            /*!
             * \brief Calculates the neighbors of a particle
             *
             * Uses the ring topology methodology to calculate the neighbors
             * of a given patricle in the swarm.
             *
             * \param[in] particleIndex The particle index to find the
             *  neighbors for
             *
             * \return The indices of the particle's neighbors
             */
            NeighborIndices neighbors(size_t particleIndex, size_t swarmSize)
                    const;


            /*!
             * \brief bestPreviousBestPosition retrieves the argmin of all
             *  previous best positions in a neighboorhood
             *
             * \param[in] neighborhood The particles that form a neighborhood
             *
             * \return The position vector
             *
             * \sa #neighbors()
             */
            vector_t bestPreviousBestPosition(Neighborhood const& neighborhood)
                    const;


            /*!
             * \brief Creates the initial swarm
             *
             * \param[in] dimension Dimension of the search space, i.e., the
             *  number of parameters to the evaluation function
             *
             * \param[in] evaluator The evaluation function that initially
             *  evaluates each particle's fitness
             *
             * \return The new swarm
             *
             * \sa Evaluator
             */
            Swarm createSwarm(size_t dimension, Evaluator const& evaluator);


            /*!
             * \brief Applies the algorithm with the configured parameters
             *
             * \param[in] dimensions The dimensions of the search space.
             *
             * \param[in] evalutor The objective evaluation function, i.e.,
             *  the fitness function
             *
             * \return The best set of parameters found
             *
             * \sa Evaluator
             */
            detail::ParticleSwarmOptimizationResult run(
                    size_t dimension,
                    Evaluator const& evaluator);


        private:


            //! The size of the swarm
            size_t m_swarmSize;


            //! Maximum number of iterations we search for
            size_t m_maxIterations;


            //! The lower boundary of the search space
            double m_lowerBoundary;


            //! The upper boundary of the search space
            double m_upperBoundary;


            //! Our RNG
            boost::random::mt11213b m_randomNumberGenerator;


            //! A uniform distribution for the RNG, [0;1)
            boost::random::uniform_01<double> m_uniformDistribution;
        };
    } // namespace Algorithm
} // namespace Winzent


namespace std {
    static ostream& operator <<(
            ostream& os,
            Winzent::Algorithm::detail::Particle const& particle)
    {
        os
                << "Particle("
                << "velocity = " << particle.velocity << ", "
                << "bestFitness = " << particle.bestFitness << ", "
                << "bestPosition = " << particle.bestPosition << ", "
                << "currentFitness = " << particle.currentFitness << ", "
                << "currentPosition = " << particle.currentPosition
                << ")";

        return os;
    }


    static ostream& operator <<(
            ostream& os,
            Winzent::Algorithm::ParticleSwarmOptimization::Swarm const& swarm)
    {
        os << "Swarm(";

        for (auto const& p: swarm) {
            os << p;
            if (&p != &(swarm.back())) {
                os << ", ";
            }
        }

        os << ")";
        return os;
    }
}

#endif // WINZENT_ALGORITHM_PARTICLESWARMOPTIMIZATION_H
