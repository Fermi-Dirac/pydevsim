import unittest



class LightsourceTestCase(unittest.TestCase):
    def test_am0(self):
        from devsim.light_sources import AM0
        src = AM0()
        self.assertEqual(src.lambda_min, 280)
        self.assertEqual(src.lambda_max, 4000)
        self.assertEqual(src.irradiance(280), 4000)
        self.assertEqual(src.irradiance(4000), 4000)