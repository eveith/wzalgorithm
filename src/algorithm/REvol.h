#ifndef WINZENT_ALGORITHM_REVOL_H_
#define WINZENT_ALGORITHM_REVOL_H_


#include <vector>
#include <cstddef>
#include <functional>

#include <QPair>
#include <QVector>

#include <log4cxx/logger.h>

#include <boost/random.hpp>

#include "algorithm_global.h"


namespace Winzent {
    namespace Algorithm {
        namespace detail {


            /*!
             * \brief The Individual class represents an individual during the
             *  training phase of the evolutionary algorithm
             */
            struct ALGORTHMSHARED_EXPORT Individual
            {
                //! The parameters vector of this individual
                QVector<qreal> parameters;


                //! The scatter vector for the parameters
                QVector<qreal> scatter;


                //! The individual's time to live: How long may it exist?
                std::ptrdiff_t timeToLive;


                //! Restrictions specifing the fitness of the individual
                QVector<qreal> restrictions;


                /*!
                 * \brief Constructs an empty Individual
                 *
                 * Initializes the Individual's TTL to 0.
                 */
                Individual();


                /*!
                 * \brief Copy constructor
                 *
                 * \param[in] The other individual to copy
                 */
                Individual(const Individual &other);


                /*!
                 * \brief Ages the individuum, i.e., reduces its time to live.
                 *
                 * \return `*this`
                 */
                Individual &age();


                /*!
                 * \brief Compares one individual to another.
                 *
                 * \param other The other individual
                 *
                 * \return -1 if the other individual is better than this one,
                 *  0 if they are equal, or 1 if this one is better than
                 *  the other individual.
                 */
                int compare(const Individual &other) const;


                /*!
                 * \brief Checks whether this object is better than
                 *  another one
                 *
                 * \param[in] other The object to check againts
                 *
                 * \return true if the object is better, false otherwise.
                 */
                bool isBetterThan(const Individual &other) const;


                /*!
                 * \brief Checks whether an individual is better than another
                 *
                 * \param[in] i1 The first individual
                 *
                 * \param[in] i2 The second individual
                 *
                 * \return True, iff `i1` is better than `i2`.
                 */
                static bool isIndividual1Better(
                        Individual const& i1,
                        Individual const& i2);


                /*!
                 * \brief Compares two Individuals for equality in
                 *  all vectors.
                 *
                 * \param[in] other The Individual to compare the current
                 *  one to
                 *
                 * \return `true` iff equal, `false` otherwise.
                 */
                bool operator==(const Individual &other) const;


                /*!
                 * \brief Deep copy operator
                 *
                 * \param[in] other The other Individual to copy
                 *
                 * \return A deep copy
                 */
                Individual &operator=(const Individual &other);


                /*!
                 * \brief Compare two Individuals, applying the semantics of
                 *  #isBetterThan
                 *
                 * \param[in] other The other Individual to compare against
                 *
                 * \return `true` if this Individual is better than the other
                 *  one, `false` otherwise
                 */
                bool operator<(const Individual &other) const;
            };


            //! The result of a run of the REvol algorithm
            struct ALGORTHMSHARED_EXPORT REvolResult
            {
                //! The best individual
                Individual bestIndividual;

                //! Number of iterations the algorithm took
                std::size_t iterationsUsed;
            };
        }


        /*!
         * \brief Implements an optimization algorithm called "REvol" that
         *  implements a multi-part evolutionary algorithm with implicit
         *  gradient information and dynamic reproduction probabilty spread
         *  of the individuals.
         */
        class ALGORTHMSHARED_EXPORT REvol
        {
        public:


            //! Auto-deleting vector for the population
            typedef std::vector<detail::Individual> Population;


            /*!
             * \brief User-supplied function for evaluation of individuals
             *
             * This function is supplied by the caller and used to evaluate
             * individuals. It has write-access to the individual it shall
             * evaluate.
             *
             * Its return type signifies the individual's success or failure
             * to solve the problem presented to it by the user. As long as
             * this function returns `false`, the algorithm continues. Once
             * the function returns `true`, REvol terminates, since the
             * user-supplied restrictions have all been fulfilled.
             *
             * \param[inout] individual The individual that is being evaluated
             *
             * \return `true` if the individual's evaluation satisfies the
             *  user-determined success criterion, `false` otherwise
             */
            typedef std::function<bool(detail::Individual &)> Evaluator;


            //! A time-discrete LTI system of first order
            static qreal pt1(const qreal &y, const qreal &u, const qreal &t);


            /*!
             * \brief Reduces each individual's time to live.
             *
             * \param population The population
             */
            static void agePopulation(Population &population);


            //! Creates a new, un-initialized instance of this algorithm
            REvol();


            /*!
             * \brief Creates a random number in the interval [0.0, 1.0)
             *
             * \return A random number between 0.0 and 1.0 (exclusive)
             */
            qreal frandom();


            /*!
             * \brief Returns the maximum number of iterations this algorithm
             *  will try to optimize
             *
             * \return The maximum number of epochs
             */
            std::size_t maxEpochs() const;


            /*!
             * \brief Sets the maximum number of iterations (epochs) this
             *  algorithm will run
             *
             * \param[in] epochs The maximum number of epochs
             *
             * \return `*this`
             */
            REvol &maxEpochs(const std::size_t &epochs);


            /*!
             * \brief Returns the maximum number of epochs that may pass
             *  without a global improvement before the algorithm stops
             *
             * \return The number of epochs
             */
            std::size_t maxNoSuccessEpochs() const;


            /*!
             * \brief Sets the maximum number of epochs that may pass without
             *  a global improvement
             *
             * \param[in] epochs The new maximum number of epochs
             *
             * \return `*this`
             */
            REvol &maxNoSuccessEpochs(const std::size_t &epochs);


            /*!
             * \brief Returns the size of the population
             *
             * \return The population's size
             */
            std::size_t populationSize() const;


            /*!
             * \brief Sets the size of the population
             *
             * \param[in] size The new population size
             *
             * \return `*this`
             */
            REvol &populationSize(const std::size_t &size);


            /*!
             * \brief Returns the size of the elite
             *
             * \return The number of elite individuals within the population
             */
            std::size_t eliteSize() const;


            /*!
             * \brief Sets the number of elite objects within the population
             *
             * \param[in] size The new size of the elite
             *
             * \return `*this`
             */
            REvol &eliteSize(const std::size_t &size);


            /*!
             * \brief Returns the currently set gradient weight
             *
             * \return The gradient weight
             */
            qreal gradientWeight() const;


            /*!
             * \brief Sets the new gradient weight
             *
             * This weight factor sets the influence of the
             * implicit gradient information on the training process.
             * Setting it to 0.0 completeley disables this feature.
             * Values between [1.0, 3.0] typically yield the best
             * results. Values > 5.0 are probably not useful.
             *
             * \param weight The new weight
             *
             * \return `*this`
             */
            REvol &gradientWeight(const qreal &weight);


            /*!
             * \brief Retrieves the weight factor applied to the error metric
             *
             * \return The error weight
             */
            qreal successWeight() const;


            /*!
             * \brief Sets the new error weight
             *
             * This value influences how much a given reproduction success
             * is weighted in when creating a new offspring and its parameters.
             *
             * \param weight The new weight
             *
             * \return `*this`
             */
            REvol &successWeight(const qreal &weight);


            /*!
             * \brief The smallest absolute change applied during object
             *  creation
             *
             * \return The smallest absolute delta
             */
            qreal eamin() const;


            /*!
             * \brief Sets the smallest absolute delta for parameters
             *  or scatter
             *
             * This is the smallest absolute change we apply; smaller values
             * are not accepted. Typically, this is the smallest floating
             * point number the architecture accepts and acts as a safety net.
             *
             * \param eamin The new smallest absolute delta
             *
             * \return `*this`
             */
            REvol &eamin(const qreal &eamin);


            /*!
             * \brief Retrieves the current smallest relative change
             *
             * \return The smallest relative change for scatter/parameters
             */
            qreal ebmin() const;


            /*!
             * \brief Sets the new smallest relative delta applied to scatter
             *  and parameters
             *
             * We default to an `ebmin` so that `|(1.0 + ebmin) > 1.0|`.
             *
             * \param ebmin The new relative minimum
             *
             * \return `*this`
             */
            REvol &ebmin(const qreal &ebmin);


            /*!
             * \brief The relative maximum delta for scatter and parameters
             *
             * \return The relative maximum
             */
            qreal ebmax() const;


            /*!
             * \brief Sets the new relative maximum delta for scatter and
             *  parameter changes
             *
             * We default to an ebmax that does not lead to a too big spread
             * of individuals through reproduction: `ebmax < 10.0`.
             *
             * \param ebmax The new maximum delta
             *
             * \return `*this`
             */
            REvol &ebmax(const qreal &ebmax);


            /*!
             * \brief The default start TTL for new Individuals
             *
             * This variable is used when generating new objects. It should
             * be set before starting the training. A fail-safe default is
             * applied if not set: `ceil(0.1 * maxEpochs)`.
             *
             * \return The initial Time-To-Live
             */
            std::ptrdiff_t startTTL() const;


            /*!
             * \brief Sets the new initial TTL for new Individuals
             *
             * \param[in] ttl The new TTL
             *
             * \return `*this`
             */
            REvol &startTTL(const std::ptrdiff_t &ttl);


            /*!
             * \brief Mean number of epochs for application in the pt1 method.
             *
             * \return Number of epochs
             */
            size_t measurementEpochs() const;


            /*!
             * \brief Number of epochs applied to the dc1 method
             *
             * \param[in] epochs The number of epochs
             *
             * \return `*this`
             */
            REvol &measurementEpochs(const size_t &epochs);


            /*!
             * \brief Generates the initial population from the supplied base
             *  network
             *
             * \param[in] origin The origin individual the user supplied for
             *  training
             *
             * \return The population, including the elite
             */
            Population generateInitialPopulation(
                    const detail::Individual &origin);


            /*!
             * \brief Creates a new individual as part of the combination
             *  and crossover process
             *
             * This method generates a new offspring by crossover and
             * mutation. It selects an of the elite and one normal object.
             *
             * In the process, it also modifies existing objects in order to
             * re-train them.
             *
             * \param population The current population; must be sorted
             *
             * \return A pair containing the two individuals out of the
             *  population that were used to modifiy the given target
             *  individual. The frist reference denotes the better, the second
             *  the other individual that were chosen.
             */
            void modifyWorstIndividual(Population &population);


            /*!
             * \brief Executes the multi-part evoluationary algorithm
             *
             * \param[in] origin The origin individual
             *
             * \param[in] succeeds The evaluating predicate that returns
             *  `true` if the given individual succeeds at the user-defined
             *  goal, or `false` if the search must go on.
             *
             * \return The best individual
             */
            detail::REvolResult run(
                    const detail::Individual &origin,
                    const Evaluator &succeeds);


        private:


            //! Our logger
            log4cxx::LoggerPtr logger;


            //! The maximum number of iterations this alorithm runs
            std::size_t m_maxEpochs;


            /*!
             * \brief Maximum number of epochs that may pass without a global
             *  improvement
             */
            std::size_t m_maxNoSuccessEpochs;


            //! \brief Overall size of the population
            std::size_t m_populationSize;


            //! \brief Size of the elite, contained in the population
            std::size_t m_eliteSize;


            //! \brief Weight of implicit gradient information.
            qreal m_gradientWeight;


            //! \brief Weight of the reproduction success rate
            qreal m_successWeight;


            /*!
             * \brief Smallest absolute delta; typically the smallest number
             *  we can store
             */
            qreal m_eamin;


            /*!
             * \brief Smallest relative delta
             */
            qreal m_ebmin;


            /*!
             * \brief The biggest relative change
             */
            qreal m_ebmax;


            /*!
             * \brief Initial Time-To-Live for new individuals
             */
            std::ptrdiff_t m_startTTL;


            /*!
             * \brief Number of epochs to apply to the dc1 method
             */
            std::size_t m_measurementEpochs;


            /*!
             * \brief Success of reproduction
             */
            qreal m_success;


            /*!
             * \brief Target value on which the population has reached
             *  equilibrium
             */
            qreal m_targetSuccess;


            /*!
             * \brief Our Random Number Generator
             */
            boost::random::mt11213b m_randomNumberGenerator;


            /*!
             * \brief A uniform distribution `[0; 1)` from which we draw
             *  random numbers
             */
            boost::random::uniform_01<qreal> m_uniformDistribution;


            //! The uniform integer distribution used to select individuals
            boost::uniform_int<> m_rnDistribution;


            /*!
             * \brief Applies the bounds defined in ebmin, eamin and eamax
             *  given another object's parameter
             *
             * \param[in] dx The delta X that shall be checked and corrected
             *
             * \param[in] parameter Another object's parameter
             *
             * \return The corrected delta X
             */
            qreal applyDxBounds(const qreal &dx, const qreal &parameter)
                    const;


            /*!
             * \brief Checks that all parameters are within safe bounds
             *
             * \return true iff all parameters are in sensible bounds
             */
            bool hasSensibleTrainingParameters() const;
        };
    }
}


namespace std {
    ostream &operator<<(
            ostream &os,
            const Winzent::Algorithm::detail::Individual &individual);
    ostream &operator<<(
            ostream &os,
            const Winzent::Algorithm::REvol::Population &v);
    ostream &operator<<(
            ostream &os,
            const Winzent::Algorithm::REvol &algorithm);
}

#endif // WINZENT_ALGORITHM_REVOL_H_
