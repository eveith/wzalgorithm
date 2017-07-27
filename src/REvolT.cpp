#include <cfenv>
#include <cmath>
#include <limits>
#include <future>
#include <cassert>
#include <cstdlib>
#include <ostream>
#include <algorithm>
#include <functional>

#include <CTPL/ctpl.h>

#include <boost/random.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "REvol.h"
#include "REvolT.h"


using std::exp;
using std::fabs;
using std::ceil;
using std::size_t;
using std::distance;
using std::ptrdiff_t;
using std::make_pair;
using std::numeric_limits;


namespace wzalgorithm {
    void REvolT::agePopulation(REvolT::Population& population)
    {
        for (auto i : population) {
            i->age();
        }
    }


    static void sort(
            ptrdiff_t l,
            ptrdiff_t r,
            REvolT::Population& population)

    {
        auto i = l;
        auto j = r;
        auto x = population[vector_t::size_type(l + r) / 2];

        do {
            while (population[vector_t::size_type(i)]->isBetterThan(*x)) {
                ++i;
            }

            while (x->isBetterThan(*(population[vector_t::size_type(j)]))) {
                --j;
            }

            if (i <= j) {
                std::swap(
                        population[vector_t::size_type(i)],
                        population[vector_t::size_type(j)]);
                ++i;
                --j;
            }
        } while (i <= j);

        if (l < j) {
            sort(l, j, population);
        }

        if (i < r) {
            sort(i, r, population);
        }
    }


    void REvolT::sort(REvolT::Population& population)
    {
        wzalgorithm::sort(
                0,
                static_cast<ptrdiff_t>(population.size()) - 1,
                population);
    }


    REvolT::REvolT():
            m_maxNoSuccessEpochs(std::numeric_limits<size_t>::max()),
            m_populationSize(0),
            m_eliteSize(0),
            m_gradientWeight(1.0),
            m_successWeight(1.0),
            m_eamin(std::numeric_limits<double>::min()),
            m_ebmin(1e-12),
            m_ebmax(1e-1),
            m_startTTL(0),
            m_measurementEpochs(5000),
            m_targetSuccess(0.25),
            m_randomNumberGenerator(0xCAFEu)
    {
    }


    REvolT::~REvolT()
    {
    }


    double REvolT::frandom()
    {
        return frandom(m_randomNumberGenerator);
    }


    double REvolT::frandom(REvolT::rng_t& rng)
    {
        return m_uniformDistribution(rng);
    }



    size_t REvolT::maxEpochs() const
    {
        return m_maxEpochs;
    }


    REvolT& REvolT::maxEpochs(std::size_t epochs)
    {
        m_maxEpochs = epochs;
        return *this;
    }

    size_t REvolT::maxNoSuccessEpochs() const
    {
        return m_maxNoSuccessEpochs;
    }


    REvolT& REvolT::maxNoSuccessEpochs(std::size_t epochs)
    {
        m_maxNoSuccessEpochs = epochs;
        return *this;
    }


    size_t REvolT::populationSize() const
    {
        return m_populationSize;
    }


    REvolT& REvolT::populationSize(size_t size)
    {
        m_populationSize = size;

        if (0 == eliteSize()) {
            eliteSize(static_cast<size_t>(std::ceil(size * 0.1)));
        }

        if (0 == startTTL()) {
            startTTL(static_cast<ptrdiff_t>(size) * 3);
        }

        return *this;
    }


    size_t REvolT::eliteSize() const
    {
        return m_eliteSize;
    }


    REvolT& REvolT::eliteSize(size_t size)
    {
        m_eliteSize = size;
        return *this;
    }


    double REvolT::gradientWeight() const
    {
        return m_gradientWeight;
    }


    REvolT& REvolT::gradientWeight(double weight)
    {
        m_gradientWeight = weight;
        return *this;
    }


    double REvolT::successWeight() const
    {
        return m_successWeight;
    }


    REvolT& REvolT::successWeight(double weight)
    {
        m_successWeight = weight;
        return *this;
    }


    double REvolT::targetSuccess() const
    {
        return m_targetSuccess;
    }


    REvolT& REvolT::targetSuccess(double targetSuccess)
    {
        m_targetSuccess = targetSuccess;
        return *this;
    }


    double REvolT::eamin() const
    {
        return m_eamin;
    }


    REvolT& REvolT::eamin(double eamin)
    {
        m_eamin = eamin;
        return *this;
    }


    double REvolT::ebmin() const
    {
        return m_ebmin;
    }


    REvolT& REvolT::ebmin(double ebmin)
    {
        m_ebmin = ebmin;
        return *this;
    }


    double REvolT::ebmax() const
    {
        return m_ebmax;
    }


    REvolT& REvolT::ebmax(double ebmax)
    {
        m_ebmax = ebmax;
        return *this;
    }


    std::ptrdiff_t REvolT::startTTL() const
    {
        return m_startTTL;
    }


    REvolT& REvolT::startTTL(std::ptrdiff_t ttl)
    {
        m_startTTL = ttl;
        return *this;
    }


    size_t REvolT::measurementEpochs() const
    {
        return m_measurementEpochs;
    }


    REvolT& REvolT::measurementEpochs(size_t epochs)
    {
        m_measurementEpochs = epochs;
        return *this;
    }


    double REvolT::clamp(double dx, double parameter) const
    {
        if (std::fetestexcept(FE_UNDERFLOW)) {
            dx = eamin();
        }

        if (std::fetestexcept(FE_OVERFLOW)) {
            dx = ebmax() * fabs(parameter);
        }

        if (dx < ebmin() * fabs(parameter)) {
            dx = ebmin() * fabs(parameter);
        }

        if (dx > ebmax() * fabs(parameter)) {
            dx = ebmax() * fabs(parameter);
        }

        if (dx < eamin()) {
            dx = eamin();
        }

        std::feclearexcept(FE_ALL_EXCEPT);
        return dx;
    }


    bool REvolT::hasSensibleTrainingParameters() const
    {
        bool ok = true;

        if (0 >= populationSize()) {
            ok = false;
        }

        if (0 >= eliteSize()) {
            ok = false;
        }

        if (eliteSize() >= populationSize()) {
            ok = false;
        }

        if (startTTL() <= 0) {
            ok = false;
        }

        if (0 >= measurementEpochs()) {
            ok = false;
        }

        return ok;
    }


    REvolT::Population REvolT::generateInitialPopulation(
            detail::Individual const& origin)
    {
        assert(origin.parameters.size() == origin.scatter.size());

        Population population;
        population.reserve(populationSize() + 1);
        auto numParameters = origin.parameters.size();


        // Handle base (origin) individual:

        auto* baseIndividual = new detail::Individual(origin);
        baseIndividual->timeToLive = startTTL();
        baseIndividual->restrictions.push_back(
                std::numeric_limits<double>::max());
        population.push_back(baseIndividual);

        thread_local rng_t rng;
        FuturePopulation futurePopulation;
        futurePopulation.reserve(populationSize());
        ctpl::thread_pool pool(static_cast<int>(
                std::thread::hardware_concurrency()));

        for (size_t i = 0; i < populationSize(); ++i) {
                futurePopulation.push_back(pool.push(
                        [this, &numParameters, &origin](int) {
                auto individual = new detail::Individual();
                individual->timeToLive = startTTL();
                individual->parameters.reserve(numParameters);
                individual->scatter.reserve(numParameters);
                individual->restrictions.push_back(
                        numeric_limits<double>::max());

                for (size_t j = 0; j != numParameters; ++j) {
                    const double r = origin.scatter.at(j)
                            * exp(0.4 * (0.5 - frandom(rng)));
                    individual->scatter.push_back(r);
                    individual->parameters.push_back(
                            origin.parameters.at(j) + r * (
                                    frandom(rng)
                                    - frandom(rng)
                                    + frandom(rng)
                                    - frandom(rng)));
                }

                return individual;
            }));
        }

        for (auto& f: futurePopulation) {
            population.push_back(f.get());
        }

        assert(population.size() == populationSize() + 1);
        return population;
    }


    void REvolT::modifyIndividual(
            detail::Individual& individual,
            PopulationRange populationRange,
            double currentSuccess)
    {
        thread_local rng_t rng;

        // Select proper individuals:

        const auto nIndividuals = distance(
                populationRange.first,
                populationRange.second);
        assert(nIndividuals >= 2);
        auto eliteIdx = abs(static_cast<ptrdiff_t>(
                m_rnDistribution(rng) % eliteSize()
                    - (m_rnDistribution(rng) % eliteSize())));
        auto otherIdx = m_rnDistribution(rng) % nIndividuals;

        if (**(populationRange.first + otherIdx)
                < **(populationRange.first + eliteIdx)) {
            std::swap(eliteIdx, otherIdx);
        }

        auto eliteIndividual = *(populationRange.first + eliteIdx);
        auto const otherIndividual = *(populationRange.first + otherIdx);


        // Determine influence of current success rate and gradient:

        double xlp = 0.0;
        const double successRate = currentSuccess / targetSuccess() - 1.0;
        const int gradientSwitch(m_rnDistribution(rng) % 3);
        const double expvar = exp(frandom(rng) - frandom(rng));

        if (2 == gradientSwitch) {
            xlp = (frandom(rng) + frandom(rng) + frandom(rng) + frandom(rng)
                        + frandom(rng) + frandom(rng) + frandom(rng)
                        + frandom(rng) + frandom(rng) + frandom(rng)
                        - frandom(rng) - frandom(rng) - frandom(rng)
                        - frandom(rng) - frandom(rng) - frandom(rng))
                * gradientWeight();

            if (xlp > 0.0) {
                xlp *= 0.5;
            }

            xlp *= exp(gradientWeight() * successRate);
        }

        // Now modify the new individual:

        const auto numParameters = eliteIndividual->parameters.size();

        assert(numParameters == eliteIndividual->parameters.size());
        assert(numParameters == otherIndividual->parameters.size());

        for (vector_t::size_type i = 0; i != numParameters; ++i) {
            std::feclearexcept(FE_ALL_EXCEPT);

            double dx = eliteIndividual->scatter.at(i) * exp(
                    successWeight() * successRate);
            dx = clamp(dx, eliteIndividual->parameters.at(i));

            // Mutate scatter:

            eliteIndividual->scatter[i] = dx;

            if (frandom(rng) < 0.5) {
                dx = eliteIndividual->scatter.at(i);
            } else {
                dx = 0.5 * (eliteIndividual->scatter.at(i)
                        + otherIndividual->scatter.at(i));
            }

            dx *= expvar;
            dx = clamp(dx, eliteIndividual->parameters.at(i));

            // Generate new scatter:

            individual.scatter[i] = dx;

            dx *= (frandom(rng) + frandom(rng) + frandom(rng) + frandom(rng)
                    + frandom(rng) - frandom(rng) - frandom(rng)
                    - frandom(rng) - frandom(rng) - frandom(rng));

            if (0 == gradientSwitch) { // Everything from the elite, p=2/3
                if (m_rnDistribution(rng) % 3 < 2) {
                    dx += eliteIndividual->parameters.at(i);
                } else {
                    dx += otherIndividual->parameters.at(i);
                }
            } else if (1 == gradientSwitch) { // use eliteIndividual
                dx += eliteIndividual->parameters.at(i);
            } else if (2 == gradientSwitch) { // use elite & gradient
                dx += eliteIndividual->parameters.at(i);
                dx += xlp * (eliteIndividual->parameters.at(i)
                        - otherIndividual->parameters.at(i));
            }

            individual.parameters[i] = dx;
        }

        individual.timeToLive = startTTL();
        individual.restrictions[0] = numeric_limits<double>::max();
    }


    detail::REvolResult REvolT::run(
            const detail::Individual& origin,
            REvol::Evaluator const& succeeds)
    {
        if (0 == startTTL()) {
            startTTL(static_cast<ptrdiff_t>(populationSize()) * 3);
        }

        if (0 == measurementEpochs()) {
            measurementEpochs(static_cast<size_t>(ceil(maxEpochs() / 200.0)));
        }

        if (!hasSensibleTrainingParameters()) {
            throw "Training parameters have no sensible values, "
                        "won't train.";
            return { origin, 0 };
        }

        size_t lastSuccess = 0;
        size_t epoch       = 0;
        auto currentSuccess= targetSuccess();

        detail::Individual* bestIndividual;
        Population population = generateInitialPopulation(origin);
        ctpl::thread_pool pool(static_cast<int>(
                std::thread::hardware_concurrency()));
        std::vector<std::future<bool>> evaluations;
        evaluations.reserve(population.size());

        // Do first-time evaluation of the population:

        for (auto individual: population) {
            evaluations.push_back(pool.push([&succeeds, &individual](int) {
                return succeeds(*individual);
            }));
        }

        for (size_t i = 0; i != evaluations.size(); ++i) {
            if (evaluations[i].get()) {
                bestIndividual = population.at(i);
                goto out;
            }
        }

        sort(population);
        bestIndividual = population.front();


        do {
            evaluations.clear();

            // Modify the worst individual(s):
            for (int i = 1; i <= pool.size(); ++i) {
                evaluations.push_back(pool.push([
                        &i,
                        this,
                        &pool,
                        &succeeds,
                        &population,
                        &currentSuccess,
                        &bestIndividual](int) -> bool {
                    auto individual = *(population.end() - i);
                    modifyIndividual(
                            *individual,
                            make_pair(
                                population.begin(),
                                population.end() - pool.size()),
                            currentSuccess);

                    return succeeds(*individual);
                }));
            }

            // Check for addition of a new individual, or global success:

            auto worstIndividual = population.at(vector_t::size_type(
                    population.size() - pool.size() - 1));

            for (auto it = evaluations.begin(); it != evaluations.end();
                    it++) {
                auto newIndividual = *(population.end() - 1
                            - (it - evaluations.begin()));

                assert(it->valid());
                if (it->get()) {
                    bestIndividual = newIndividual;
                    break;
                }

                // Check for global or, at least, local improvement:

                if (! worstIndividual->isBetterThan(*newIndividual)) {
                    if (worstIndividual->timeToLive >= 0) {
                        currentSuccess = REvol::pt1(
                                currentSuccess,
                                1.0,
                                measurementEpochs());
                    } else {
                        currentSuccess = REvol::pt1(
                                currentSuccess,
                                -1.0,
                                measurementEpochs());
                    }
                }

                if (newIndividual->isBetterThan(*bestIndividual)) {
                    lastSuccess = epoch;
                    bestIndividual = newIndividual;
                    bestIndividual->timeToLive = ptrdiff_t(epoch);
                }
            }

            // Sort the list and do a bit of caretaking:

            agePopulation(population);
            sort(population);
            currentSuccess = REvol::pt1(
                    currentSuccess,
                    0.0,
                    measurementEpochs());
            epoch += 1;
        } while (epoch < maxEpochs()
                && epoch - lastSuccess <= maxNoSuccessEpochs());

out:
        detail::REvolResult result = { *bestIndividual, epoch };
        for (auto i : population) {
            delete i;
        }
        return result;
    }
} // namespace wzalgorithm


namespace std {
    ostream &operator<<(
            ostream &os,
            const wzalgorithm::REvolT& algorithm)
    {
        os
                << "REvolT("
                << "maxNoSuccessEpochs = " << algorithm.maxNoSuccessEpochs()
                << ", populationSize = " << algorithm.populationSize()
                << ", eliteSize = " << algorithm.eliteSize()
                << ", startTTL = " << algorithm.startTTL()
                << ", eamin = " << algorithm.eamin()
                << ", ebmin = " << algorithm.ebmin()
                << ", ebmax = " << algorithm.ebmax();
        os << ")";
        return os;
    }
} // namespace std

