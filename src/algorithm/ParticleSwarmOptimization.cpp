#include <cmath>
#include <tuple>
#include <cstddef>
#include <functional>

#include <QVector>

#include <boost/random.hpp>

#include "ParticleSwarmOptimization.h"


namespace Winzent {
    namespace Algorithm {


        const qreal ParticleSwarmOptimization::C = 0.5 + std::log(2);
        const qreal ParticleSwarmOptimization::W = 1.0 / (2.0 * std::log(2));


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


        qreal ParticleSwarmOptimization::upperBoundary() const
        {
            return m_upperBoundary;
        }


        ParticleSwarmOptimization &ParticleSwarmOptimization::upperBoundary(
                const qreal &boundary)
        {
            m_upperBoundary = boundary;
            return *this;
        }



        ParticleSwarmOptimization::NeighborIndices
        ParticleSwarmOptimization::neighbors(
                const size_t &particleIndex,
                const size_t &swarmSize)
                const
        {
            return std::make_tuple(
                    (particleIndex - 1) % swarmSize,
                    particleIndex,
                    (particleIndex + 1) % swarmSize);
        }


        ParticleSwarmOptimization::Swarm
        ParticleSwarmOptimization::createSwarm(
                const size_t &dimension,
                const Evaluator &evaluator)
        {
            Swarm swarm;

            for (size_t i = 0; i != swarmSize(); ++i) {
                detail::Particle particle;

                for (size_t d = 0; d != dimension; ++d) {
                    qreal p = lowerBoundary()
                            + m_uniformDistribution(m_randomNumberGenerator)
                                * (upperBoundary() - lowerBoundary());
                    qreal v = (lowerBoundary() - p)
                            + m_uniformDistribution(m_randomNumberGenerator)
                                * ((lowerBoundary() - p)
                                   - (upperBoundary() - p));

                    particle.currentPosition.push_back(p);
                    particle.velocity.push_back(v);
                }

                particle.currentFitness = evaluator(
                        particle.currentPosition);
                particle.bestFitness = particle.currentFitness;
                particle.bestPosition = particle.currentPosition;

                swarm.push_back(particle);
            }

            // Now we need to compute the neighborhood information. We need to
            // iterate again over the swarm, sadly.

            for (size_t i = 0; i != swarmSize(); ++i) {
                NeighborIndices ni = neighbors(i, swarmSize());
                qreal fitness = swarm.at(std::get<0>(ni)).currentFitness;
                size_t index = 0;

                if (swarm.at(std::get<1>(ni)).currentFitness < fitness) {
                    index = 1;
                    fitness = swarm.at(std::get<1>(ni)).currentFitness;
                }
                if (swarm.at(std::get<2>(ni)).currentFitness < fitness) {
                    index = 2;
                    fitness = swarm.at(std::get<2>(ni)).currentFitness;
                }

                swarm[i].bestPreviousBestPosition =
                        swarm[index].currentPosition;
            }

            return swarm;
        }


        QVector<qreal> ParticleSwarmOptimization::run(
                const size_t &dimension,
                const Evaluator &evaluator)
        {
            Swarm swarm = createSwarm(dimension, evaluator);
            detail::Particle *best;

            std::sort(swarm.begin(), swarm.end());
            best = &(swarm.front());

            for (size_t i = 0; i != maxIterations(); ++i) {
                for (auto &particle: swarm) {
                    for (size_t j = 0; j != particle.currentPosition.size();
                            ++j) {
                        qreal v = particle.velocity[j];
                        qreal x = particle.currentPosition[j];
                        qreal p = particle.bestPosition[j];
                        qreal g = x + C * ((p - x)/2.0);
                        qreal r = g + ((g-x)-g) * m_uniformDistribution(
                                m_randomNumberGenerator);

                        particle.velocity[j] = W*v + r - x;
                        particle.currentPosition[j] =
                                W * particle.velocity[j] + r;
                    }

                    particle.currentFitness = evaluator(
                            particle.currentPosition);
                    if (particle.currentFitness > particle.bestFitness) {
                        particle.bestFitness = particle.currentFitness;
                        particle.bestPosition = particle.currentPosition;
                    }

                    if (particle.bestFitness > best->bestFitness) {
                        best = &particle;
                    }
                }
            }

            return best->bestPosition;
        }
    } // namespace Algorithm
} // namespace Winzent

