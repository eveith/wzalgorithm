#include <cmath>
#include <tuple>
#include <limits>
#include <cstddef>
#include <functional>

#include <QList>
#include <QVector>

#include <boost/random.hpp>

#include <log4cxx/logger.h>
#include <log4cxx/logmanager.h>

#include "ParticleSwarmOptimization.h"


namespace Winzent {
    namespace Algorithm {


        const qreal ParticleSwarmOptimization::C = 0.5 + std::log(2);
        const qreal ParticleSwarmOptimization::W = 1.0 / (2.0 * std::log(2));


        ParticleSwarmOptimization::ParticleSwarmOptimization():
                logger(log4cxx::LogManager::getLogger(
                    "Winzent.Algorithm.ParticleSwarmOptimization")),
                m_swarmSize(DEFAULT_SWARM_SIZE),
                m_maxIterations(std::numeric_limits<size_t>::max()),
                m_randomNumberGenerator(0xBEEFu)
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


        QVector<qreal> ParticleSwarmOptimization::bestPreviousBestPosition(
                const ParticleSwarmOptimization::Neighborhood &neighborhood)
                const
        {
            auto *best = neighborhood.front();

            for (auto &particle: neighborhood) {
                if (particle->bestFitness < best->bestFitness) {
                    best = particle;
                }
            }

            return best->bestPosition;
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

                evaluator(particle);
                particle.bestFitness = particle.currentFitness;
                particle.bestPosition = particle.currentPosition;

                swarm.push_back(particle);
            }

            LOG4CXX_DEBUG(logger, "Created swarm: " << swarm);

            return swarm;
        }


        detail::ParticleSwarmOptimizationResult
        ParticleSwarmOptimization::run(
                const size_t &dimension,
                const Evaluator &evaluator)
        {
            Swarm swarm = createSwarm(dimension, evaluator);
            detail::Particle *best;

            std::sort(swarm.begin(), swarm.end());
            best = &(swarm.front());
            bool success = false;
            size_t i;

            for (i = 0; i < maxIterations() && !success; ++i) {
                for (int k = 0; k != swarm.size(); ++k) {
                    auto &particle = swarm[k];
                    auto neighborIndices = neighbors(k, swarm.size());
                    auto bestNeighborhoodPosition = bestPreviousBestPosition({
                            &(swarm[std::get<0>(neighborIndices)]),
                            &(swarm[std::get<1>(neighborIndices)]),
                            &(swarm[std::get<2>(neighborIndices)]) });

                    for (detail::Particle::Vector::size_type j = 0;
                            j != particle.currentPosition.size(); ++j) {
                        bool isBestInNeighborhood = (particle.bestPosition
                                == bestNeighborhoodPosition);
                        qreal v = particle.velocity[j];
                        qreal x = particle.currentPosition[j];
                        qreal p = particle.bestPosition[j];
                        qreal l = bestNeighborhoodPosition[j];
                        qreal g = (isBestInNeighborhood
                                ? x + C * ((p - x) / 2.0)
                                : x + C * ((p + l - 2.0 * x) / 3.0));
                        qreal r = g + ((g-x)-g) * m_uniformDistribution(
                                m_randomNumberGenerator);

                        qreal newV = W*v + r - x;
                        qreal newX = W * particle.velocity[j] + r;

                        if (newX < lowerBoundary()) {
                            newX = lowerBoundary();
                            newV *= -0.5;
                        }
                        if (newX > upperBoundary()) {
                            newX = upperBoundary();
                            newV *= -0.5;
                        }

                        particle.velocity[j] = newV;
                        particle.currentPosition[j] = newX;
                    }

                    success = evaluator(particle);

                    if (particle.currentFitness < particle.bestFitness) {
                        particle.bestFitness = particle.currentFitness;
                        particle.bestPosition = particle.currentPosition;
                    }

                    if (particle.bestFitness < best->bestFitness) {
                        best = &particle;
                    }
                }
            }

            return { *best, i };
        }
    } // namespace Algorithm
} // namespace Winzent

