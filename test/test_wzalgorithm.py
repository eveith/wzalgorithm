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

    def __call__(self, individual):
        individual.restrictions.resize(1)
        individual.restrictions[0] = CrossinTrayEvaluator.crossin_tray(individual.parameters[0],
                                                                       individual.parameters[1])
        return individual.restrictions[0] < -2.07


class EggholderEvaluator:
    def eggholder(x1, x2):
        return -(x2 + 47) * math.sin(math.sqrt(abs(x2 + 0.5 * x1 +47))) \
            - x1 * math.sin(math.sqrt(abs(x1 - (x2 + 47))))

    def __call__(self, individual):
        individual.restrictions.resize(1)
        individual.restrictions[0] = EggholderEvaluator.eggholder(individual.parameters[0],
                                                                  individual.parameters[1])
        return individual.restrictions[0] < -959


class TestREvol(unittest.TestCase):
    def test_crossin_tray(self):
        revol = wzalgorithm.REvol()
        origin = revol.generateOrigin(2)

        revol.maxEpochs(10000)
        result = revol.run(origin,
                           wzalgorithm.REvolSuccessPredicate(CrossinTrayEvaluator()))
        self.assertTrue(result.error < -2.07)

    def test_eggholder(self):
        revol = wzalgorithm.REvol()
        origin = revol.generateOrigin(2)

        revol.maxEpochs(10000)
        result = revol.run(origin,
                           wzalgorithm.REvolSuccessPredicate(EggholderEvaluator()))
        self.assertTrue(result.error < -959)



if __name__ == '__main__':
    unittest.main()
