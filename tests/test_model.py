import unittest
from goodwin_model.model import goodwin_model_simplified
from goodwin_model.integration import integrate_model

class TestGoodwinModel(unittest.TestCase):
    
    def test_goodwin_model_simplified(self):
        y0 = [1, 0, 0, 0, 0, 0, 0]
        t = [0, 1]
        parameters = {
            'nu1': 0.7, 'nu2': 0.5, 'nu3': 0.40, 'nu4': 0.3, 'nu5': 0.70,
            'nu6': 0.35, 'nu7': 0.3, 'nu8': 0.2, 'nu9': 0.1, 'nu10': 0.2,
            'nu11': 0.03, 'nu12': 0.05, 'nu13': 0.3, 'nu14': 0.2, 'K1': 1.00,
            'K2': 1.00, 'K3': 1.00, 'K4': 1.00, 'K5': 0.80, 'K6': 1.00,
            'K7': 1.00, 'K8': 1.00, 'K9': 1.00, 'K10': 1.00, 'hill': 3,
            'hill_S': 1.00, 'hill_W': 1.00, 'b': 1.0, 'c': 1.0, 'd': 0.3
        }

        result = goodwin_model_simplified(y0, t, parameters)
        self.assertEqual(len(result), 7)
    
    def test_integrate_model(self):
        y0 = [1, 0, 0, 0, 0, 0, 0]
        t = [0, 1]
        parameters = {
            'nu1': 0.7, 'nu2': 0.5, 'nu3': 0.40, 'nu4': 0.3, 'nu5': 0.70,
            'nu6': 0.35, 'nu7': 0.3, 'nu8': 0.2, 'nu9': 0.1, 'nu10': 0.2,
            'nu11': 0.03, 'nu12': 0.05, 'nu13': 0.3, 'nu14': 0.2, 'K1': 1.00,
            'K2': 1.00, 'K3': 1.00, 'K4': 1.00, 'K5': 0.80, 'K6': 1.00,
            'K7': 1.00, 'K8': 1.00, 'K9': 1.00, 'K10': 1.00, 'hill': 3,
            'hill_S': 1.00, 'hill_W': 1.00, 'b': 1.0, 'c': 1.0, 'd': 0.3
        }

        result = integrate_model(y0, t, parameters)
        self.assertEqual(result.shape, (2, 7))

if __name__ == '__main__':
    unittest.main()
