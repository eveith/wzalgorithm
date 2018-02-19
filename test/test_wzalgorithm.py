import math
import unittest
import wzalgorithm



class CrossinTrayEvaluator:
    def crossin_tray(x, y):
        return -0.0001 * math.pow( \
                (abs( \
                    math.sin(x) \
                    * math.sin(y) \
                    * math.exp( \
                        abs( \
                            100 - math.sqrt(math.pow(x, 2) * math.pow(y, 2))/ math.pi))) \
                + 1), \
                0.1)

    def __init__(self):
        print("Welcome!")

    def __call__(self, individual):
        individual.restrictions.resize(1)
        individual.restrictions[0] = CrossinTrayEvaluator.crossin_tray(individual.parameters[0],
                                                                       individual.parameters[1])
        return individual.restrictions[0] < -2.07

    def __del__(self):
        print("Good-bye says the evaluator.")


class TestREvol(unittest.TestCase):
    def test_crossin_tray(self):
        revol = wzalgorithm.REvol()
        origin = revol.generateOrigin(2)

        revol.maxEpochs(5000)
        result = revol.runPredicated(origin,
                                     wzalgorithm.REvolSuccessPredicate(CrossinTrayEvaluator()))
        print("Result is: {}".format(result.error))
