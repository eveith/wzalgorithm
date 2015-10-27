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
#include <boost/ptr_container/ptr_vector.hpp>

#include "REvol.h"


#pragma STDC FENV_ACCESS ON


using std::exp;
using std::fabs;
using std::numeric_limits;


namespace Winzent {
    namespace Algorithm {


        log4cxx::LoggerPtr REvol::logger =
                log4cxx::LogManager::getLogger("Winzent.Algorithm.REvol");


        Individual::Individual(): timeToLive(0)
        {
        }


        Individual::Individual(const Individual &other):
            parameters(other.parameters),
            scatter(other.scatter),
            timeToLive(other.timeToLive),
            restrictions(other.restrictions)
        {
        }


        Individual &Individual::age()
        {
            timeToLive -= 1;
            return *this;
        }


        int Individual::compare(const Individual &other) const
        {
            if (this->timeToLive < 0 && other.timeToLive >= 0) {
                return -1;
            }

            if (this->timeToLive >= 0 && other.timeToLive < 0) {
                return 1;
            }

            if (this->restrictions.at(0) < other.restrictions.at(0)) {
                return 1;
            } else if (restrictions.at(0) > other.restrictions.at(0)) {
                return -1;
            } else if (restrictions.at(0) == other.restrictions.at(0)) {
                if (timeToLive > other.timeToLive) {
                    return 1;
                }
            } else {
                auto size = (
                        restrictions.size() > other.restrictions.size()
                        ? other.restrictions.size()
                        : this->restrictions.size());
                for (auto i = 1; i < size; ++i) {
                    if (restrictions.at(i) < other.restrictions.at(i)) {
                        return 1;
                    } else if (other.restrictions.at(i)
                            < this->restrictions.at(i)) {
                        return -1;
                    }
                }
            }

            return 0;
        }


        bool Individual::isBetterThan(const Individual &other) const
        {
            return (1 == compare(other));
        }


        bool Individual::isIndividual1Better(
                Individual const& i1,
                Individual const& i2)
        {
            return i1.isBetterThan(i2);
        }


        bool Individual::operator ==(const Individual &other) const
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
                        && other.parameters() == parameters());
#endif
        }


        Individual &Individual::operator =(const Individual &other)
        {
            if (this == &other) {
                return *this;
            }

            timeToLive    = other.timeToLive;
            parameters    = other.parameters;
            scatter       = other.scatter;
            restrictions  = other.restrictions;

            return *this;
        }


        qreal REvol::pt1(const qreal &y, const qreal &u, const qreal &t)
        {
            qreal r = 0.0;

            if (t + 1.0 != 1.0) {
                r = y + ((u - y) / t);
            } else {
                r = u;
            }

            return r;
        }


        void REvol::sortPopulation(
                Population& population)
        {
            population.sort(&Individual::isIndividual1Better);
        }


        void REvol::agePopulation(Population &population)
        {
            for (Individual &i: population) {
                i.age();
            }
        }


        qreal REvol::frandom()
        {
            return m_uniformDistribution(m_randomNumberGenerator);
        }


        REvol::REvol():
                m_maxNoSuccessEpochs(0),
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
                m_normalDistributionZero(0.0, 0.5),
                m_normalDistributionMinusTwo(-2.0, 0.5)
        {
        }


        size_t REvol::maxEpochs() const
        {
            return m_maxEpochs;
        }


        REvol &REvol::maxEpochs(const size_t &epochs)
        {
            m_maxEpochs = epochs;
            if (0 == maxNoSuccessEpochs()) {
                maxNoSuccessEpochs(epochs);
            }
            return *this;
        }

        size_t REvol::maxNoSuccessEpochs() const
        {
            return m_maxNoSuccessEpochs;
        }


        REvol &REvol::maxNoSuccessEpochs(const size_t &epochs)
        {
            m_maxNoSuccessEpochs = epochs;
            return *this;
        }


        size_t REvol::populationSize() const
        {
            return m_populationSize;
        }


        REvol &REvol::populationSize(const size_t &size)
        {
            m_populationSize = size;

            if (0 == eliteSize()) {
                eliteSize(std::ceil(size * 0.1));
            }

            return *this;
        }


        size_t REvol::eliteSize() const
        {
            return m_eliteSize;
        }


        REvol &REvol::eliteSize(const size_t &size)
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


        int REvol::startTTL() const
        {
            return m_startTTL;
        }


        REvol &REvol::startTTL(const int &ttl)
        {
            m_startTTL = ttl;
            return *this;
        }


        size_t REvol::measurementEpochs() const
        {
            return m_measurementEpochs;
        }


        REvol &REvol::measurementEpochs(const size_t &epochs)
        {
            m_measurementEpochs = epochs;
            return *this;
        }


        qreal REvol::applyDxBounds(
                const qreal &dx,
                const qreal &parameter)
                const
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
                const Individual &origin)
        {
            Q_ASSERT(origin.parameters.size() == origin.scatter.size());

            Individual *baseIndividual = new Individual(origin);
            baseIndividual->timeToLive = startTTL();
            baseIndividual->restrictions.push_front(
                    std::numeric_limits<qreal>::infinity());

            Population population;
            population.push_back(baseIndividual);
            auto numParameters = baseIndividual->parameters.size();

            for (auto i = 1; i < populationSize() + 1; ++i) {
                Individual *individual = new Individual();
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
                                + r * m_normalDistributionZero(
                                    m_randomNumberGenerator));
                }

                population.push_back(individual);
            }

            Q_ASSERT(population.size() == populationSize() + 1);
            return population;
        }


        QPair<Individual &, Individual &> REvol::modifyIndividual(
                Individual &individual,
                Population &population)
        {
            Individual &eliteIndividual = population.at(abs(
                    (m_rnDistribution(m_randomNumberGenerator) % eliteSize())
                        - (m_rnDistribution(m_randomNumberGenerator)
                            % eliteSize())));
            Individual &otherIndividual = population.at(
                    m_rnDistribution(m_randomNumberGenerator)
                        % population.size());

            if (otherIndividual.isBetterThan(eliteIndividual)) {
                Individual &tmp = eliteIndividual;
                eliteIndividual = otherIndividual;
                otherIndividual = tmp;
            }

            qreal xlp = 0.0;
            qreal successRate = m_success / m_targetSuccess - 1.0;
            int gradientSwitch = m_rnDistribution(m_randomNumberGenerator) % 3;
            qreal expvar = exp(frandom() - frandom());

            if (2 == gradientSwitch) {
                xlp = m_normalDistributionMinusTwo(m_randomNumberGenerator)
                    * gradientWeight();

                if (xlp > 0.0) {
                    xlp *= 0.5;
                }

                xlp *= exp(gradientWeight() * successRate);
            }

            // Now modify the new individual:

            auto numParameters = eliteIndividual.parameters.size();
            QVector<qreal> newParameters;
            newParameters.reserve(numParameters);

            Q_ASSERT(numParameters == eliteIndividual.parameters.size());
            Q_ASSERT(numParameters == otherIndividual.parameters.size());

            for (auto i = 0; i != numParameters; ++i) {
                std::feclearexcept(FE_ALL_EXCEPT);

                qreal dx = eliteIndividual.scatter.at(i) * exp(
                        successWeight() * successRate);

                dx = applyDxBounds(dx, eliteIndividual.parameters.at(i));

                // Mutate scatter:

                eliteIndividual.scatter[i] = dx;

                if (frandom() < 0.5) {
                    dx = eliteIndividual.scatter.at(i);
                } else {
                    dx = 0.5 * (eliteIndividual.scatter.at(i)
                            + otherIndividual.scatter.at(i));
                }

                dx *= expvar;
                dx = applyDxBounds(dx, eliteIndividual.parameters.at(i));

                // Generate new scatter:

                individual.scatter[i] = dx;

                dx = dx * (frandom() + frandom()
                        + frandom() + frandom() + frandom() - frandom()
                        - frandom() - frandom() - frandom() - frandom());

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

                newParameters.push_back(dx);
            }

            Q_ASSERT(newParameters.size() == individual.parameters.size());
            individual.parameters = newParameters;
            individual.timeToLive = startTTL();
            individual.restrictions[0] = numeric_limits<qreal>::infinity();

#ifdef      QT_DEBUG
                for (auto i = 0; i != newParameters.size(); ++i) {
                    Q_ASSERT(newParameters.at(i)
                             == individual.parameters.at(i));
                }
#endif

            LOG4CXX_DEBUG(logger, "Created " << individual);

            return QPair<Individual &, Individual &>(
                    eliteIndividual,
                    otherIndividual);
        }


        Individual REvol::run(const Individual &origin, Evaluator evaluator)
        {
            if (0 == startTTL()) {
                startTTL(std::ceil(maxEpochs() * 0.1));
            }

            if (0 == measurementEpochs()) {
                measurementEpochs(std::ceil(maxEpochs() / 200.0));
            }

            if (!hasSensibleTrainingParameters()) {
                LOG4CXX_ERROR(
                        logger,
                        "Training parameters have no sensible values, "
                            "won't train.");
                return origin;
            }

            size_t lastSuccess = 0;
            size_t epoch       = 0;
            bool successful    = false;
            Population population = generateInitialPopulation(origin);

            do {
                // Modify the worst individual:

                if (0 != epoch) {
                    auto individual = population.pop_back();
                    auto srcIndividuals = modifyIndividual(
                            *individual,
                            population);

                    successful |= evaluator(*individual);
                    successful |= evaluator(srcIndividuals.first);
                    successful |= evaluator(srcIndividuals.second);

                    population.push_back(individual.release());
                } else {
                    for (auto &individual: population) {
                        successful |= evaluator(individual);
                    }

                    sortPopulation(population);
                }

                // Check for addition of a new individual:

                Individual &newIndividual = population.back();
                Individual &bestIndividual = population.front();
                Individual &worstIndividual = population.at(
                        population.size() - 2);

                // Check for global or, at least, local improvement:

                if (! worstIndividual.isBetterThan(newIndividual)) {
                    if (worstIndividual.timeToLive >= 0) {
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

                if (newIndividual.isBetterThan(bestIndividual)) {
                    lastSuccess = epoch;
                    //newIndividual.timeToLive = 0;
                    newIndividual.timeToLive = epoch;
                }

                // Sort the list and do a bit of caretaking:

                sortPopulation(population);
                agePopulation(population);
                m_success = pt1(m_success, 0.0, measurementEpochs());

                LOG4CXX_DEBUG(
                        logger,
                        "Iteration(epoch = " << epoch
                            << ", success = " << m_success
                            << ", targetSuccess = " << m_targetSuccess
                            << ", lastSuccess = " << lastSuccess
                            << ", population = " << population
                            << ")");

                epoch++;
            } while (! successful
                    && epoch < maxEpochs()
                    && epoch - lastSuccess <= maxNoSuccessEpochs());

            LOG4CXX_DEBUG(
                    logger,
                    "Training ended after " << epoch << " epochs; "
                        << population.front());

            return population.front();
        }
    } // namespace ANN
} // namespace Winzent


namespace std {
    ostream &operator<<(
            ostream &os,
            const Winzent::Algorithm::Individual &individual)
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
                << ", eamin = " << algorithm.eamin()
                << ", ebmin = " << algorithm.ebmin()
                << ", ebmax = " << algorithm.ebmax();
        os << ")";
        return os;
    }


    ostream &operator<<(
            ostream &os,
            const boost::ptr_vector<Winzent::Algorithm::Individual> &v)
    {
        os << "(";
        for (const Winzent::Algorithm::Individual &i: v) {
            os << i;
            if (&i != &(v.back())) {
                os << ", ";
            }
        }
        os << ")";

        return os;
    }
}

