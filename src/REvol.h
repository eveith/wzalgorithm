#ifndef WINZENT_ALGORITHM_REVOL_H_
#define WINZENT_ALGORITHM_REVOL_H_


#include <cstddef>
#include <functional>

#include <boost/random.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "config.h"


namespace wzalgorithm {
    namespace detail {


        /*!
         * \brief The Individual class represents an individual during the
         *  training phase of the evolutionary algorithm
         */
        struct Individual
        {
            //! The parameters vector of this individual
            vector_t parameters;


            //! The scatter vector for the parameters
            vector_t scatter;


            //! The individual's time to live: How long may it exist?
            std::ptrdiff_t timeToLive;


            //! Restrictions specifing the fitness of the individual
            vector_t restrictions;


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
            Individual(Individual const& other);


            //! \brief Move constructor
            Individual(Individual&& other);


            /*!
             * \brief Ages the individuum, i.e., reduces its time to live.
             *
             * \return `*this`
             */
            Individual& age();


            /*!
             * \brief Checks whether this object is better than
             *  another one
             *
             * \param[in] other The object to check againts
             *
             * \return true if the object is better, false otherwise.
             */
            bool isBetterThan(Individual const& other) const;


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
            bool operator ==(Individual const& other) const;


            /*!
             * \brief Compares two individuals for unequality
             *
             * \param[in] other The other individual to compare the
             *  current one to
             *
             * \return `! (*this == other)`
             *
             * \sa operator==()
             */
            bool operator !=(Individual const& other) const;


            /*!
             * \brief Deep copy operator
             *
             * \param[in] other The other Individual to copy
             *
             * \return A deep copy
             */
            Individual& operator=(Individual const& other);


            /*!
             * \brief Compare two Individuals, applying the semantics of
             *  #isBetterThan
             *
             * \param[in] other The other Individual to compare against
             *
             * \return `true` if this Individual is better than the other
             *  one, `false` otherwise
             */
            bool operator <(Individual const& other) const;
        };


        //! The result of a run of the REvol algorithm
        struct REvolResult
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
    class REvol
    {
    public:


        //! Auto-deleting vector for the population
        typedef boost::ptr_vector<detail::Individual> Population;


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
        typedef std::function<bool (detail::Individual&)> Evaluator;


        //! A time-discrete LTI system of first order
        static double pt1(double y, double u, double t);


        /*!
         * \brief Reduces each individual's time to live.
         *
         * \param population The population
         */
        static void agePopulation(Population& population);


        //! Creates a new, un-initialized instance of this algorithm
        REvol();


        /*!
         * \brief Creates a random number in the interval [0.0, 1.0)
         *
         * \return A random number between 0.0 and 1.0 (exclusive)
         */
        double frandom();


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
        REvol& maxEpochs(size_t epochs);


        /*!
         * \brief Returns the maximum number of epochs that may pass
         *  without a global improvement before the algorithm stops
         *
         * \return The number of epochs
         */
        size_t maxNoSuccessEpochs() const;


        /*!
         * \brief Sets the maximum number of epochs that may pass without
         *  a global improvement
         *
         * \param[in] epochs The new maximum number of epochs
         *
         * \return `*this`
         */
        REvol& maxNoSuccessEpochs(size_t epochs);


        /*!
         * \brief Returns the size of the population
         *
         * \return The population's size
         */
        size_t populationSize() const;


        /*!
         * \brief Sets the size of the population
         *
         * \param[in] size The new population size
         *
         * \return `*this`
         */
        REvol& populationSize(size_t size);


        /*!
         * \brief Returns the size of the elite
         *
         * \return The number of elite individuals within the population
         */
        size_t eliteSize() const;


        /*!
         * \brief Sets the number of elite objects within the population
         *
         * \param[in] size The new size of the elite
         *
         * \return `*this`
         */
        REvol& eliteSize(size_t size);


        /*!
         * \brief Returns the currently set gradient weight
         *
         * \return The gradient weight
         */
        double gradientWeight() const;


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
        REvol& gradientWeight(double weight);


        /*!
         * \brief Retrieves the weight factor applied to the error metric
         *
         * \return The error weight
         */
        double successWeight() const;


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
        REvol& successWeight(double weight);


        /*!
         * \brief Accesses the target success rate.
         *
         * \return The target sucess rate
         */
        double targetSuccess() const;


        /*!
         * \brief Sets the targeted success rate
         *
         * \param[in] targetSuccess The targeted success rate
         *
         * \return `*this`
         *
         * \sa REvol::measurementEpochs()
         */
        REvol &targetSuccess(double targetSuccess);


        /*!
         * \brief The smallest absolute change applied during object
         *  creation
         *
         * \return The smallest absolute delta
         */
        double eamin() const;


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
        REvol& eamin(double eamin);


        /*!
         * \brief Retrieves the current smallest relative change
         *
         * \return The smallest relative change for scatter/parameters
         */
        double ebmin() const;


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
        REvol& ebmin(double ebmin);


        /*!
         * \brief The relative maximum delta for scatter and parameters
         *
         * \return The relative maximum
         */
        double ebmax() const;


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
        REvol& ebmax(double ebmax);


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
        REvol& startTTL(std::ptrdiff_t ttl);


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
        REvol& measurementEpochs(size_t epochs);


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
                detail::Individual const& origin);


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
         * \param[in] currentSuccess The momentary success rate
         *
         * \return A pair containing the two individuals out of the
         *  population that were used to modifiy the given target
         *  individual. The frist reference denotes the better, the second
         *  the other individual that were chosen.
         */
        void modifyWorstIndividual(
                Population& population,
                double currentSuccess);


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
                detail::Individual const& origin,
                Evaluator const& succeeds);


    private:


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
        double m_gradientWeight;


        //! \brief Weight of the reproduction success rate
        double m_successWeight;


        /*!
         * \brief Smallest absolute delta; typically the smallest number
         *  we can store
         */
        double m_eamin;


        //! \brief Smallest relative delta
        double m_ebmin;


        //! \brief The biggest relative change
        double m_ebmax;


        //! \brief Initial Time-To-Live for new individuals
        std::ptrdiff_t m_startTTL;


        //! \brief Number of epochs to apply to the dc1 method
        std::size_t m_measurementEpochs;


        /*!
         * \brief Target value on which the population has reached
         *  equilibrium
         */
        double m_targetSuccess;


        /*!
         * \brief Our Random Number Generator
         */
        boost::random::mt11213b m_randomNumberGenerator;


        /*!
         * \brief A uniform distribution `[0; 1)` from which we draw
         *  random numbers
         */
        boost::random::uniform_01<double> m_uniformDistribution;


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
        double clamp(double dx, double parameter) const;


        /*!
         * \brief Checks that all parameters are within safe bounds
         *
         * \return true iff all parameters are in sensible bounds
         */
        bool hasSensibleTrainingParameters() const;
    };
} // namespace wzalgorithm


namespace std {
    ostream& operator<<(
            ostream& os,
            wzalgorithm::detail::Individual const& individual);
    ostream& operator<<(
            ostream &os,
            wzalgorithm::REvol::Population const& v);
    ostream& operator<<(
            ostream& os,
            wzalgorithm::REvol const& algorithm);
} // namespace std

#endif // WINZENT_ALGORITHM_REVOL_H_
