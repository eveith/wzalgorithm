#include <cmath>
#include <tuple>
#include <limits>
#include <cstddef>
#include <functional>

#include <boost/random.hpp>

#include "ParticleSwarmOptimization.h"


using std::begin;
using std::end;


namespace wzalgorithm {
    const double ParticleSwarmOptimization::C = 0.5 + std::log(2);
    const double ParticleSwarmOptimization::W = 1.0 / (2.0 * std::log(2));
    const ParticleSwarmOptimization::Swarm::size_type
            ParticleSwarmOptimization::DEFAULT_SWARM_SIZE;


    ParticleSwarmOptimization::ParticleSwarmOptimization():
            m_swarmSize(DEFAULT_SWARM_SIZE),
            m_maxIterations(std::numeric_limits<size_t>::max()),
            m_randomNumberGenerator(0xBEEFu)
    {
    }


    ParticleSwarmOptimization::Swarm::size_type
    ParticleSwarmOptimization::swarmSize() const
    {
        return m_swarmSize;
    }


    ParticleSwarmOptimization& ParticleSwarmOptimization::swarmSize(
            ParticleSwarmOptimization::Swarm::size_type size)
    {
        m_swarmSize = size;
        return *this;
    }


    ParticleSwarmOptimization::epoch_t
    ParticleSwarmOptimization::maxIterations() const
    {
        return m_maxIterations;
    }


    ParticleSwarmOptimization& ParticleSwarmOptimization::maxIterations(
            ParticleSwarmOptimization::epoch_t iterations)
    {
        m_maxIterations = iterations;
        return *this;
    }


    double ParticleSwarmOptimization::lowerBoundary() const
    {
        return m_lowerBoundary;
    }


    ParticleSwarmOptimization& ParticleSwarmOptimization::lowerBoundary(
            double boundary)
    {
        m_lowerBoundary = boundary;
        return *this;
    }


    double ParticleSwarmOptimization::upperBoundary() const
    {
        return m_upperBoundary;
    }


    ParticleSwarmOptimization& ParticleSwarmOptimization::upperBoundary(
            double boundary)
    {
        m_upperBoundary = boundary;
        return *this;
    }



    ParticleSwarmOptimization::NeighborIndices
    ParticleSwarmOptimization::neighbors(
            ParticleSwarmOptimization::Swarm::size_type particleIndex,
            ParticleSwarmOptimization::Swarm::size_type swarmSize)
            const
    {
        return std::make_tuple(
                (particleIndex - 1) % swarmSize,
                particleIndex,
                (particleIndex + 1) % swarmSize);
    }


    vector_t ParticleSwarmOptimization::bestPreviousBestPosition(
            ParticleSwarmOptimization::Neighborhood const& neighborhood)
            const
    {
        auto *best = neighborhood.front();

        for (auto const& particle: neighborhood) {
            if (particle->bestFitness < best->bestFitness) {
                best = particle;
            }
        }

        return best->bestPosition;
    }


    ParticleSwarmOptimization::Swarm
    ParticleSwarmOptimization::createSwarm(
            vector_t::size_type dimension,
            Evaluator const& evaluator)
    {
        Swarm swarm;

        for (size_t i = 0; i != swarmSize(); ++i) {
            ParticleSwarmOptimization::Particle particle;

            for (vector_t::size_type d = 0; d != dimension; ++d) {
                double p = lowerBoundary()
                        + m_uniformDistribution(m_randomNumberGenerator)
                            * (upperBoundary() - lowerBoundary());
                double v = (lowerBoundary() - p)
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

        return swarm;
    }


    ParticleSwarmOptimization::Result
    ParticleSwarmOptimization::run(
            vector_t::size_type dimension,
            Evaluator const& evaluator)
    {
        Swarm swarm = createSwarm(dimension, evaluator);
        std::sort(begin(swarm), end(swarm));

        auto* best      = &(swarm.front());
        bool success    = false;
        epoch_t i       = 0;

        while (i < maxIterations() && ! success) {
            for (Swarm::size_type k = 0; k != swarm.size(); ++k) {
                auto &particle = swarm[k];
                auto neighborIndices = neighbors(k, swarm.size());
                auto bestNeighborhoodPosition = bestPreviousBestPosition({
                        &(swarm[std::get<0>(neighborIndices)]),
                        &(swarm[std::get<1>(neighborIndices)]),
                        &(swarm[std::get<2>(neighborIndices)]) });

                for (vector_t::size_type j = 0;
                        j != particle.currentPosition.size();
                        ++j) {
                    bool isBestInNeighborhood = (particle.bestPosition
                            == bestNeighborhoodPosition);
                    const double v = particle.velocity[j];
                    const double x = particle.currentPosition[j];
                    const double p = particle.bestPosition[j];
                    const double l = bestNeighborhoodPosition[j];
                    const double g = (isBestInNeighborhood
                            ? x + C * ((p - x) / 2.0)
                            : x + C * ((p + l - 2.0 * x) / 3.0));
                    const double r = g + ((g-x)-g) * m_uniformDistribution(
                            m_randomNumberGenerator);

                    double newV = W*v + r - x;
                    double newX = W * particle.velocity[j] + r;

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

                success |= evaluator(particle);

                if (particle.currentFitness < particle.bestFitness) {
                    particle.bestFitness = particle.currentFitness;
                    particle.bestPosition = particle.currentPosition;
                }

                if (particle.bestFitness < best->bestFitness) {
                    best = &particle;
                }
            }

            i += 1;
        }

        return { *best, i };
    }
} // namespace wzalgorithm

