#include <cfenv>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <functional>

#include <boost/random.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "REvol.h"


using std::exp;
using std::fabs;
using std::ceil;
using std::size_t;
using std::ptrdiff_t;
using std::numeric_limits;


namespace wzalgorithm {
    detail::Individual::Individual(): timeToLive(0)
    {
    }


    detail::Individual::Individual(detail::Individual const& other):
            parameters(other.parameters),
            scatter(other.scatter),
            timeToLive(other.timeToLive),
            restrictions(other.restrictions)
    {
    }


    detail::Individual::Individual(detail::Individual&& other):
            parameters(std::move(other.parameters)),
            scatter(std::move(other.scatter)),
            timeToLive(std::move(other.timeToLive)),
            restrictions(std::move(other.restrictions))
    {
    }


    detail::Individual& detail::Individual::age()
    {
        timeToLive -= 1;
        return *this;
    }


    bool detail::Individual::isBetterThan(
            detail::Individual const& other)
            const
    {
        if (this->timeToLive < 0 && other.timeToLive >= 0) {
            return false;
        }

        if (this->timeToLive >= 0 && other.timeToLive < 0) {
            return true;
        }

        if (0 == restrictions.size() || 0 == other.restrictions.size()) {
            return false;
        }

        auto size = (
                restrictions.size() > other.restrictions.size()
                ? other.restrictions.size()
                : this->restrictions.size());
        for (vector_t::size_type i = 1; i < size; ++i) {
            if (this->restrictions.at(i) < other.restrictions.at(i)) {
                return true;
            } else {
                if (this->restrictions.at(i) > other.restrictions.at(i)) {
                    return false;
                }
            }
        }

        if (this->restrictions.at(0) == other.restrictions.at(0)) {
            if (this->timeToLive > other.timeToLive) {
                return true;
            } else {
                return false;
            }
        } else {
            if (this->restrictions.at(0) < other.restrictions.at(0)) {
                return true;
            } else {
                return false;
            }
        }

        return false;
    }


    bool detail::Individual::isIndividual1Better(
            detail::Individual const& i1,
            detail::Individual const& i2)
    {
        return i1.isBetterThan(i2);
    }


    bool detail::Individual::operator ==(detail::Individual const& other)
            const
    {
            return (other.timeToLive == this->timeToLive
                    && other.restrictions == this->restrictions
                    && other.scatter == this->scatter
                    && other.parameters == this->parameters);
    }


    bool detail::Individual::operator !=(detail::Individual const& other)
            const
    {
        return !(*this == other);
    }


    detail::Individual& detail::Individual::operator =(
            detail::Individual const& other)
    {
        if (this == &other) {
            return *this;
        }

        this->timeToLive    = other.timeToLive;
        this->parameters    = other.parameters;
        this->scatter       = other.scatter;
        this->restrictions  = other.restrictions;

        return *this;
    }


    bool detail::Individual::operator <(detail::Individual const& other)
            const
    {
        return this->isBetterThan(other);
    }


    double REvol::pt1(double y, double u, double t)
    {
        return (t + 1.0 != 1.0) ? (y + ((u - y) / t)) : u;
    }


    void REvol::agePopulation(Population& population)
    {
        for (auto& i: population) {
            i.age();
        }
    }


    double REvol::frandom()
    {
        return m_uniformDistribution(m_randomNumberGenerator);
    }


    REvol::REvol():
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


    size_t REvol::maxEpochs() const
    {
        return m_maxEpochs;
    }


    REvol& REvol::maxEpochs(std::size_t epochs)
    {
        m_maxEpochs = epochs;
        return *this;
    }

    size_t REvol::maxNoSuccessEpochs() const
    {
        return m_maxNoSuccessEpochs;
    }


    REvol& REvol::maxNoSuccessEpochs(std::size_t epochs)
    {
        m_maxNoSuccessEpochs = epochs;
        return *this;
    }


    size_t REvol::populationSize() const
    {
        return m_populationSize;
    }


    REvol& REvol::populationSize(size_t size)
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


    size_t REvol::eliteSize() const
    {
        return m_eliteSize;
    }


    REvol& REvol::eliteSize(size_t size)
    {
        m_eliteSize = size;
        return *this;
    }


    double REvol::gradientWeight() const
    {
        return m_gradientWeight;
    }


    REvol& REvol::gradientWeight(double weight)
    {
        m_gradientWeight = weight;
        return *this;
    }


    double REvol::successWeight() const
    {
        return m_successWeight;
    }


    REvol& REvol::successWeight(double weight)
    {
        m_successWeight = weight;
        return *this;
    }


    double REvol::targetSuccess() const
    {
        return m_targetSuccess;
    }


    REvol& REvol::targetSuccess(double targetSuccess)
    {
        m_targetSuccess = targetSuccess;
        return *this;
    }


    double REvol::eamin() const
    {
        return m_eamin;
    }


    REvol& REvol::eamin(double eamin)
    {
        m_eamin = eamin;
        return *this;
    }


    double REvol::ebmin() const
    {
        return m_ebmin;
    }


    REvol& REvol::ebmin(double ebmin)
    {
        m_ebmin = ebmin;
        return *this;
    }


    double REvol::ebmax() const
    {
        return m_ebmax;
    }


    REvol& REvol::ebmax(double ebmax)
    {
        m_ebmax = ebmax;
        return *this;
    }


    std::ptrdiff_t REvol::startTTL() const
    {
        return m_startTTL;
    }


    REvol& REvol::startTTL(std::ptrdiff_t ttl)
    {
        m_startTTL = ttl;
        return *this;
    }


    size_t REvol::measurementEpochs() const
    {
        return m_measurementEpochs;
    }


    REvol& REvol::measurementEpochs(size_t epochs)
    {
        m_measurementEpochs = epochs;
        return *this;
    }


    double REvol::clamp(double dx, double parameter) const
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


    bool REvol::hasSensibleTrainingParameters() const
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


    REvol::Population REvol::generateInitialPopulation(
            const detail::Individual &origin)
    {
        assert(origin.parameters.size() == origin.scatter.size());

        Population population;
        population.reserve(populationSize() + 1);
        auto numParameters = origin.parameters.size();


        // Handle base (origin) individual:

        auto baseIndividual = new detail::Individual(origin);
        baseIndividual->timeToLive = startTTL();
        baseIndividual->restrictions.push_back(
                std::numeric_limits<double>::max());
        population.push_back(baseIndividual);


        for (size_t i = 0; i < populationSize(); ++i) {
            auto individual = new detail::Individual();
            individual->timeToLive = startTTL();
            individual->parameters.reserve(numParameters);
            individual->scatter.reserve(numParameters);
            individual->restrictions.push_back(
                    numeric_limits<double>::max());

            for (size_t j = 0; j != numParameters; ++j) {
                const double r = origin.scatter.at(j)
                        * exp(0.4 * (0.5 - frandom()));
                individual->scatter.push_back(r);
                individual->parameters.push_back(
                        origin.parameters.at(j)
                            + r * (frandom()-frandom()+frandom()-frandom()));
            }

            population.push_back(individual);
        }


        assert(population.size() == populationSize() + 1);
        return population;
    }


    void REvol::modifyWorstIndividual(
            Population& population,
            double currentSuccess)
    {
        assert(population.size() >= 3);

        // Select proper individuals:

        Population::size_type eliteIdx(
                static_cast<Population::size_type>(abs(static_cast<long long>(
                    m_rnDistribution(m_randomNumberGenerator)
                        % eliteSize())
                - static_cast<long long>(
                    m_rnDistribution(m_randomNumberGenerator)
                        % eliteSize()))));
        Population::size_type otherIdx(
                m_rnDistribution(m_randomNumberGenerator)
                    % (population.size() - 1));

        if (population.at(otherIdx) < population.at(eliteIdx)) {
            std::swap(eliteIdx, otherIdx);
        }

        auto &individual = population.back();
        auto &eliteIndividual = population.at(eliteIdx);
        auto &otherIndividual = population.at(otherIdx);


        // Determine influence of current success rate and gradient:

        double xlp = 0.0;
        const double successRate = currentSuccess / targetSuccess() - 1.0;
        const int gradientSwitch(m_rnDistribution(m_randomNumberGenerator)%3);
        const double expvar = exp(frandom() - frandom());

        if (2 == gradientSwitch) {
            xlp = (frandom() + frandom() + frandom() + frandom()
                        + frandom() + frandom() + frandom() + frandom()
                        + frandom() + frandom() - frandom() - frandom()
                        - frandom() - frandom() - frandom() - frandom())
                * gradientWeight();

            if (xlp > 0.0) {
                xlp *= 0.5;
            }

            xlp *= exp(gradientWeight() * successRate);
        }

        // Now modify the new individual:

        const auto numParameters = eliteIndividual.parameters.size();

        assert(numParameters == eliteIndividual.parameters.size());
        assert(numParameters == otherIndividual.parameters.size());

        for (vector_t::size_type i = 0; i != numParameters; ++i) {
            std::feclearexcept(FE_ALL_EXCEPT);

            double dx = eliteIndividual.scatter.at(i) * exp(
                    successWeight() * successRate);
            dx = clamp(dx, eliteIndividual.parameters.at(i));

            // Mutate scatter:

            eliteIndividual.scatter[i] = dx;

            if (frandom() < 0.5) {
                dx = eliteIndividual.scatter.at(i);
            } else {
                dx = 0.5 * (eliteIndividual.scatter.at(i)
                        + otherIndividual.scatter.at(i));
            }

            dx *= expvar;
            dx = clamp(dx, eliteIndividual.parameters.at(i));

            // Generate new scatter:

            individual.scatter[i] = dx;

            dx *= (frandom() + frandom() + frandom() + frandom()
                    + frandom() - frandom() - frandom() - frandom()
                    - frandom() - frandom());

            if (0 == gradientSwitch) { // Everything from the elite, p=2/3
                if (m_rnDistribution(m_randomNumberGenerator) % 3 < 2) {
                    dx += eliteIndividual.parameters.at(i);
                } else {
                    dx += otherIndividual.parameters.at(i);
                }
            } else if (1 == gradientSwitch) { // use eliteIndividual
                dx += eliteIndividual.parameters.at(i);
            } else if (2 == gradientSwitch) { // use elite & gradient
                dx += eliteIndividual.parameters.at(i);
                dx += xlp * (eliteIndividual.parameters.at(i)
                        - otherIndividual.parameters.at(i));
            }

            individual.parameters[i] = dx;
        }

        individual.timeToLive = startTTL();
        individual.restrictions[0] = numeric_limits<double>::max();
    }


    detail::REvolResult REvol::run(
            const detail::Individual &origin,
            const Evaluator &succeeds)
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
        Population population = generateInitialPopulation(origin);
        detail::Individual* bestIndividual = &(population.front());

        do {
            // Modify the worst individual:

            if (0 < epoch) {
                modifyWorstIndividual(population, currentSuccess);

                if (succeeds(population.back())) {
                    bestIndividual = &(population.back());
                    break;
                }
            } else {
                for (auto &individual: population) {
                    if (succeeds(individual)) {
                        bestIndividual = &individual;
                        break;
                    }
                }

                population.sort();
                bestIndividual = &(population.front());
            }

            // Check for addition of a new individual:

            auto &newIndividual = population.back();
            auto &worstIndividual = population.at(population.size() - 2);

            // Check for global or, at least, local improvement:

            if (! worstIndividual.isBetterThan(newIndividual)) {
                if (worstIndividual.timeToLive >= 0) {
                    currentSuccess = pt1(
                            currentSuccess,
                            1.0,
                            measurementEpochs());
                } else {
                    currentSuccess = pt1(
                            currentSuccess,
                            -1.0,
                            measurementEpochs());
                }
            }

            if (newIndividual.isBetterThan(*bestIndividual)) {
                lastSuccess = epoch;
                bestIndividual = &newIndividual;
                bestIndividual->timeToLive = static_cast<ptrdiff_t>(epoch);
            }

            // Sort the list and do a bit of caretaking:

            agePopulation(population);
            population.sort();
            currentSuccess = pt1(
                    currentSuccess,
                    0.0,
                    measurementEpochs());
            epoch += 1;
        } while (epoch < maxEpochs()
                && epoch - lastSuccess <= maxNoSuccessEpochs());

        detail::REvolResult result = { *bestIndividual, epoch };
        return result;
    }
} // namespace wzalgorithm


namespace std {
    ostream &operator<<(
            ostream &os,
            const wzalgorithm::detail::Individual &individual)
    {
        os << "Individual(";

        os << "TTL = " << individual.timeToLive << ", ";

        os << "Parameters = (";
        for (wzalgorithm::vector_t::size_type i = 0;
                i < individual.parameters.size(); ++i) {
            os << individual.parameters.at(i);
            if (i < individual.parameters.size() - 1) {
                os << ", ";
            }
        }

        os << "), Scatter = (";
        for (wzalgorithm::vector_t::size_type i = 0;
                i < individual.scatter.size(); ++i) {
            os << individual.scatter.at(i);
            if (i < individual.scatter.size() - 1) {
                os << ", ";
            }
        }

        os << "), Restrictions = (";
        for (wzalgorithm::vector_t::size_type i = 0;
                i < individual.restrictions.size(); ++i) {
            os << individual.restrictions.at(i);
            if (i < individual.restrictions.size() - 1) {
                os << ", ";
            }
        }

        os << "))";

        return os;
    }


    ostream &operator<<(
            ostream &os,
            const wzalgorithm::REvol &algorithm)
    {
        os
                << "REvol("
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


    ostream &operator<<(
            ostream &os,
            const wzalgorithm::REvol::Population &v)
    {
        os << "(";
        for (auto const& i: v) {
            os << i;
            if (&i != &(v.back())) {
                os << ", ";
            }
        }
        os << ")";

        return os;
    }
} // namespace std

