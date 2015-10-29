#include <QPair>
#include <QVector>

#include "SimulatedAnnealing.h"

namespace Winzent {
    namespace Algorithm {
        typedef QPair<qreal, QVector<qreal>> State;

        SimulatedAnnealing::SimulatedAnnealing()
        {
        }


        size_t SimulatedAnnealing::maxIterations() const
        {
            return m_maxIterations;
        }


        void SimulatedAnnealing::maxIterations(const size_t &iterations)
        {
            m_maxIterations = iterations;
            return *this;
        }


        qreal SimulatedAnnealing::temperature(const qreal &timeFraction) const
        {
            return timeFraction;
        }


        QVector<qreal> SimulatedAnnealing::neighbor(
                const QVector<qreal> &state,
                const qreal &temperature)
                const
        {
            return state;
        }


        QVector<qreal> SimulatedAnnealing::run(
                const QVector<qreal> &origin,
                Evaluator &evaluator)
        {
            State bestState = qMakePair(evaluator(origin), origin);
            qreal currentTemperature;

            for (size_t k = 0; k != m_maxIterations; ++k) {
                currentTemperature = temperature(
                        k / static_cast<qreal>(m_maxIterations));
                auto sNew = neighbor(bestState.second, currentTemperature);
                auto eNew = evaluator(sNew);
                qreal de = eNew - bestState.first;

                if (eNew < bestState.first) {
                    bestState = qMakePair(eNew, sNew);
                }
                if (std::exp(- de / currentTemperature)
                        > m_uniformDistribution(m_randomNumberGenerator)) {
                    bestState = qMakePair(eNew, sNew);
                }
            }

            return bestState.second;
        }
    } // namespace Algorithm
} // namespace Winzent

