#include <cstddef>
#include <functional>

#include <QVector>

#include <boost/random.hpp>

#include "ParticleSwarmOptimization.h"


namespace Winzent {
    namespace Algorithm {

        ParticleSwarmOptimization::ParticleSwarmOptimization():
                m_swarmSize(DEFAULT_SWARM_SIZE),
                m_maxIterations(0)
        {
        }


        size_t ParticleSwarmOptimization::swarmSize() const
        {
            return m_swarmSize;
        }


        ParticleSwarmOptimization &ParticleSwarmOptimization::swarmSize(
                const size_t &size)
        {
            m_swarmSize = size;
            return *this;
        }


        size_t ParticleSwarmOptimization::maxIterations() const
        {
            return m_maxIterations;
        }


        ParticleSwarmOptimization &ParticleSwarmOptimization::maxIterations(
                const size_t &iterations)
        {
            m_maxIterations = iterations;
            return *this;
        }


        qreal ParticleSwarmOptimization::lowerBoundary() const
        {
            return m_lowerBoundary;
        }


        ParticleSwarmOptimization &ParticleSwarmOptimization::lowerBoundary(
                const qreal &boundary)
        {
            m_lowerBoundary = boundary;
            return *this;
        }


        Swarm ParticleSwarmOptimization::createSwarm(
                const size_t &dimension,
                const Evaluator &evaluator)
                const
        {
            Swarm swarm;

            for (size_t i = 0; i != swarmSize(); ++i) {
                detail::Particle particle;

                for (size_t d = 0; d != dimension; ++d) {
                    qreal param =
                            lowerBoundary()
                            + m_uniformDistribution(m_randomNumberGenerator)
                                * (upperBoundary() - lowerBoundary());

                    particle.currentPosition.push_back(param);

                    // U[min_d-x_{i,d}; max_d - x_{i,d}]
                    particle.velocity = (lowerBoundary() - param)
                            + m_uniformDistribution(m_randomNumberGenerator)
                                * ((lowerBoundary() - param)
                                   - (upperBoundary() - param));
                }

                particle.currentFitness = evaluator(
                        particle.currentParameters);
                particle.bestFitness = particle.currentFitness;
                particle.bestParameters = particle.currentParameters;
            }

            return swarm;
        }


        QVector<qreal> ParticleSwarmOptimization::run(
                const size_t &dimension,
                const Evaluator &evaluator)
        {
            Swarm swarm = createSwarm(dimension, evaluator);
            detail::Particle &best;

            std::sort(swarm.begin(), swarm.end());
            best = swarm.front();

            for (size_t i = 0; i != maxIterations(); ++i) {

            }

            return best.bestParameters;
        }
    } // namespace Algorithm
} // namespace Winzent

