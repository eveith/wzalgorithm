#include <limits>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <functional>

#include <QPair>
#include <QList>
#include <QVector>

#include <log4cxx/logger.h>
#include <log4cxx/logmanager.h>

#include <boost/random.hpp>

#include "REvol.h"


using std::exp;
using std::fabs;
using std::numeric_limits;


namespace Winzent {
    namespace Algorithm {


        detail::Individual::Individual(): timeToLive(0)
        {
        }


        detail::Individual::Individual(const detail::Individual &other):
                parameters(other.parameters),
                scatter(other.scatter),
                timeToLive(other.timeToLive),
                restrictions(other.restrictions)
        {
        }


        detail::Individual &detail::Individual::age()
        {
            timeToLive -= 1;
            return *this;
        }


        bool detail::Individual::isBetterThan(
                const detail::Individual &other)
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
            for (auto i = 1; i < size; ++i) {
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


        bool detail::Individual::operator ==(const detail::Individual &other)
                const
        {
#ifdef      QT_DEBUG
                bool ok = true;
                ok &= timeToLive == other.timeToLive;
                ok &= other.restrictions == restrictions;
                ok &= other.scatter == scatter;
                ok &= other.parameters == parameters;
                return ok;
#else
                return (other.timeToLive == timeToLive
                        && other.restrictions == restrictions
                        && other.scatter == scatter
                        && other.parameters == parameters);
#endif
        }


        detail::Individual &detail::Individual::operator =(
                const detail::Individual &other)
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


        bool detail::Individual::operator <(const detail::Individual &other)
                const
        {
            return this->isBetterThan(other);
        }


        qreal REvol::pt1(const qreal &y, const qreal &u, const qreal &t)
        {
            if (t + 1.0 != 1.0) {
                return y + ((u - y) / t);
            } else {
                return u;
            }
        }


        void REvol::agePopulation(Population &population)
        {
            for (auto &i: population) {
                i->age();
            }
        }


        qreal REvol::frandom()
        {
            return m_uniformDistribution(m_randomNumberGenerator);
        }


        REvol::REvol():
                logger(log4cxx::LogManager::getLogger(
                    "Winzent.Algorithm.REvol")),
                m_maxNoSuccessEpochs(std::numeric_limits<size_t>::max()),
                m_populationSize(0),
                m_eliteSize(0),
                m_gradientWeight(1.0),
                m_successWeight(1.0),
                m_eamin(1e-30),
                m_ebmin(1e-7),
                m_ebmax(1e-1),
                m_startTTL(0),
                m_measurementEpochs(5000),
                m_success(0.25),
                m_targetSuccess(0.25),
                m_randomNumberGenerator(0xCAFEu)
        {
        }


        size_t REvol::maxEpochs() const
        {
            return m_maxEpochs;
        }


        REvol &REvol::maxEpochs(const std::size_t &epochs)
        {
            m_maxEpochs = epochs;
            return *this;
        }

        size_t REvol::maxNoSuccessEpochs() const
        {
            return m_maxNoSuccessEpochs;
        }


        REvol &REvol::maxNoSuccessEpochs(const std::size_t &epochs)
        {
            m_maxNoSuccessEpochs = epochs;
            return *this;
        }


        size_t REvol::populationSize() const
        {
            return m_populationSize;
        }


        REvol &REvol::populationSize(const std::size_t &size)
        {
            m_populationSize = size;

            if (0 == eliteSize()) {
                eliteSize(std::ceil(size * 0.1));
            }

            if (0 == startTTL()) {
                startTTL(size * 3);
            }

            return *this;
        }


        size_t REvol::eliteSize() const
        {
            return m_eliteSize;
        }


        REvol &REvol::eliteSize(const std::size_t &size)
        {
            m_eliteSize = size;
            return *this;
        }


        qreal REvol::gradientWeight() const
        {
            return m_gradientWeight;
        }


        REvol &REvol::gradientWeight(const qreal &weight)
        {
            m_gradientWeight = weight;
            return *this;
        }


        qreal REvol::successWeight() const
        {
            return m_successWeight;
        }


        REvol &REvol::successWeight(const qreal &weight)
        {
            m_successWeight = weight;
            return *this;
        }


        qreal REvol::targetSuccess() const
        {
            return m_targetSuccess;
        }


        REvol &REvol::targetSuccess(const qreal &targetSuccess)
        {
            m_targetSuccess = targetSuccess;
            return *this;
        }


        qreal REvol::eamin() const
        {
            return m_eamin;
        }


        REvol &REvol::eamin(const qreal &eamin)
        {
            m_eamin = eamin;
            return *this;
        }


        qreal REvol::ebmin() const
        {
            return m_ebmin;
        }


        REvol &REvol::ebmin(const qreal &ebmin)
        {
            m_ebmin = ebmin;
            return *this;
        }


        qreal REvol::ebmax() const
        {
            return m_ebmax;
        }


        REvol &REvol::ebmax(const qreal &ebmax)
        {
            m_ebmax = ebmax;
            return *this;
        }


        std::ptrdiff_t REvol::startTTL() const
        {
            return m_startTTL;
        }


        REvol &REvol::startTTL(const std::ptrdiff_t &ttl)
        {
            m_startTTL = ttl;
            return *this;
        }


        size_t REvol::measurementEpochs() const
        {
            return m_measurementEpochs;
        }


        REvol &REvol::measurementEpochs(const std::size_t &epochs)
        {
            m_measurementEpochs = epochs;
            return *this;
        }


        qreal REvol::clamp(const qreal &dx, const qreal &parameter) const
        {
            qreal cdx = dx;

            if (std::fetestexcept(FE_UNDERFLOW)) {
                LOG4CXX_DEBUG(logger, "Underflow detected");
                cdx = eamin();
            }

            if (std::fetestexcept(FE_OVERFLOW)) {
                LOG4CXX_DEBUG(logger, "Overflow detected");
                cdx = ebmax() * fabs(parameter);
            }

            if (cdx < ebmin() * fabs(parameter)) {
                cdx = ebmin() * fabs(parameter);
            }

            if (cdx > ebmax() * fabs(parameter)) {
                cdx = ebmax() * fabs(parameter);
            }

            if (cdx < eamin()) {
                cdx = eamin();
            }

            std::feclearexcept(FE_ALL_EXCEPT);
            return cdx;
        }


        bool REvol::hasSensibleTrainingParameters() const
        {
            bool ok = true;

            if (0 >= populationSize()) {
                LOG4CXX_ERROR(logger, "Population size is 0");
                ok = false;
            }

            if (0 >= eliteSize()) {
                LOG4CXX_ERROR(logger, "Elite size is 0");
                ok = false;
            }

            if (eliteSize() >= populationSize()) {
                LOG4CXX_ERROR(
                        logger,
                        "Elite is bigger or equal to population");
                ok = false;
            }

            if (startTTL() <= 0) {
                LOG4CXX_ERROR(logger, "No sensible start TTL (<= 0)");
                ok = false;
            }

            if (0 >= measurementEpochs()) {
                LOG4CXX_ERROR(
                        logger,
                        "Invalid number of epochs for measurement,"
                            "must be > 0");
                ok = false;
            }

            return ok;
        }


        REvol::Population REvol::generateInitialPopulation(
                const detail::Individual &origin)
        {
            Q_ASSERT(origin.parameters.size() == origin.scatter.size());

            auto baseIndividual = new detail::Individual(origin);
            baseIndividual->timeToLive = startTTL();
            baseIndividual->restrictions.push_front(
                    std::numeric_limits<qreal>::infinity());

            Population population;
            population.push_back(baseIndividual);
            auto numParameters = baseIndividual->parameters.size();

            for (Population::size_type i = 1; i < populationSize() + 1; ++i) {
                auto individual = new detail::Individual();
                individual->timeToLive = startTTL();
                individual->parameters.reserve(numParameters);
                individual->scatter.reserve(numParameters);
                individual->restrictions.push_back(
                        numeric_limits<qreal>::infinity());

                for (auto j = 0; j != numParameters; ++j) {
                    qreal r = baseIndividual->scatter.at(j) * exp(
                            0.4 * (0.5 - frandom()));
                    individual->scatter.push_back(r);
                    individual->parameters.push_back(
                            baseIndividual->parameters.at(j)
                                +r*(frandom()-frandom()+frandom()-frandom()));
                }

                population.push_back(individual);
            }

            Q_ASSERT(population.size() == populationSize() + 1);
            return population;
        }


        void REvol::modifyWorstIndividual(Population &population)
        {
            Q_ASSERT(population.size() >= 3);

            auto *individual = population.back();
            auto *eliteIndividual = population.at(abs(
                    (m_rnDistribution(m_randomNumberGenerator) % eliteSize())
                        - (m_rnDistribution(m_randomNumberGenerator)
                            % eliteSize())));
            auto *otherIndividual = population.at(
                    m_rnDistribution(m_randomNumberGenerator)
                        % (population.size()-1));

            if (otherIndividual->isBetterThan(*eliteIndividual)) {
                auto *tmp = otherIndividual;
                otherIndividual = eliteIndividual;
                eliteIndividual = tmp;
            }

            qreal xlp = 0.0;
            qreal successRate = m_success / m_targetSuccess - 1.0;
            int gradientSwitch = m_rnDistribution(m_randomNumberGenerator)%3;
            qreal expvar = exp(frandom() - frandom());

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

            auto numParameters = eliteIndividual->parameters.size();

            Q_ASSERT(numParameters == eliteIndividual->parameters.size());
            Q_ASSERT(numParameters == otherIndividual->parameters.size());

            for (auto i = 0; i != numParameters; ++i) {
                std::feclearexcept(FE_ALL_EXCEPT);

                qreal dx = eliteIndividual->scatter.at(i) * exp(
                        successWeight() * successRate);

                dx = clamp(dx, eliteIndividual->parameters.at(i));

                // Mutate scatter:

                eliteIndividual->scatter[i] = dx;

                if (frandom() < 0.5) {
                    dx = eliteIndividual->scatter.at(i);
                } else {
                    dx = 0.5 * (eliteIndividual->scatter.at(i)
                            + otherIndividual->scatter.at(i));
                }

                dx *= expvar;
                dx = clamp(dx, eliteIndividual->parameters.at(i));

                // Generate new scatter:

                individual->scatter[i] = dx;

                dx *= (frandom() + frandom() + frandom() + frandom()
                        + frandom() - frandom() - frandom() - frandom()
                        - frandom() - frandom());

                if (0 == gradientSwitch) { // Everything from the elite, p=2/3
                    if (m_rnDistribution(m_randomNumberGenerator) % 3 < 2) {
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

                individual->parameters[i] = dx;
            }

            individual->timeToLive = startTTL();
            individual->restrictions[0] = numeric_limits<qreal>::infinity();

            LOG4CXX_DEBUG(logger, "Modified " << individual);
        }


        detail::REvolResult REvol::run(
                const detail::Individual &origin,
                const Evaluator &succeeds)
        {
            if (0 == startTTL()) {
                startTTL(populationSize() * 3);
            }

            if (0 == measurementEpochs()) {
                measurementEpochs(std::ceil(maxEpochs() / 200.0));
            }

            if (!hasSensibleTrainingParameters()) {
                LOG4CXX_ERROR(
                        logger,
                        "Training parameters have no sensible values, "
                            "won't train.");
                return { origin, 0 };
            }

            size_t lastSuccess = 0;
            size_t epoch       = 0;
            Population population = generateInitialPopulation(origin);
            detail::Individual *bestIndividual = population.front();

            do {
                // Modify the worst individual:

                if (0 < epoch) {
                    modifyWorstIndividual(population);

                    if (succeeds(*(population.back()))) {
                        bestIndividual = population.back();
                        break;
                    }
                } else {
                    for (auto &individual: population) {
                        if (succeeds(*individual)) {
                            bestIndividual = individual;
                            break;
                        }
                    }

                    std::sort(
                            population.begin(),
                            population.end(),
                            [](detail::Individual *a, detail::Individual *b) {
                        return *a < *b;
                    });
                    bestIndividual = population.front();
                }

                // Check for addition of a new individual:

                auto &newIndividual = population.back();
                auto &worstIndividual = population.at(population.size() - 2);

                // Check for global or, at least, local improvement:

                if (! worstIndividual->isBetterThan(*newIndividual)) {
                    if (worstIndividual->timeToLive >= 0) {
                        m_success = pt1(
                                m_success,
                                1.0,
                                measurementEpochs());
                    } else {
                        m_success = pt1(
                                m_success,
                                -1.0,
                                measurementEpochs());
                    }
                }

                if (newIndividual->isBetterThan(*bestIndividual)) {
                    lastSuccess = epoch;
                    bestIndividual = newIndividual;
                    bestIndividual->timeToLive = epoch;
                }

                // Sort the list and do a bit of caretaking:

                agePopulation(population);
                std::sort(
                        population.begin(),
                        population.end(),
                        [](detail::Individual *a, detail::Individual *b) {
                    return *a < *b;
                });

                Q_ASSERT(detail::Individual::isIndividual1Better(
                        *population.front(),
                        *population.back()));

                m_success = pt1(m_success, 0.0, measurementEpochs());

                LOG4CXX_DEBUG(
                        logger,
                        "Iteration(epoch = " << epoch
                            << ", maxEpochs = " << maxEpochs()
                            << ", success = " << m_success
                            << ", targetSuccess = " << m_targetSuccess
                            << ", lastSuccess = " << lastSuccess
                            << ", maxNoSuccessEpochs = "
                                << maxNoSuccessEpochs()
                            << ", population = " << population
                            << ")");

                epoch += 1;
            } while (epoch < maxEpochs()
                    && epoch - lastSuccess <= maxNoSuccessEpochs());

            LOG4CXX_DEBUG(
                    logger,
                    "Training ended after "
                        << epoch
                        << " epochs; "
                        << "winner: "
                        << *bestIndividual);

            detail::REvolResult result = { *bestIndividual, epoch };
            for (auto &individual: population) {
                delete individual;
            }
            return result;
        }
    } // namespace ANN
} // namespace Winzent


namespace std {
    ostream &operator<<(
            ostream &os,
            const Winzent::Algorithm::detail::Individual &individual)
    {
        os << "Individual(";

        os << "TTL = " << individual.timeToLive << ", ";

        os << "Parameters = (";
        for (int i = 0; i < individual.parameters.size(); ++i) {
            os << individual.parameters.at(i);
            if (i < individual.parameters.size() - 1) {
                os << ", ";
            }
        }

        os << "), Scatter = (";
        for (int i = 0; i < individual.scatter.size(); ++i) {
            os << individual.scatter.at(i);
            if (i < individual.scatter.size() - 1) {
                os << ", ";
            }
        }

        os << "), Restrictions = (";
        for (int i = 0; i < individual.restrictions.size(); ++i) {
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
            const Winzent::Algorithm::REvol &algorithm)
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
            const Winzent::Algorithm::REvol::Population &v)
    {
        os << "(";
        for (const auto &i: v) {
            os << *i;
            if (i != v.back()) {
                os << ", ";
            }
        }
        os << ")";

        return os;
    }
}

