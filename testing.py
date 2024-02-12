import unittest
import numpy as np
from optimization import PolygonSpringSimulation
from setup import DiskEnergyCalculator

class TestPolygonSpringSimulation(unittest.TestCase):

    def setUp(self):
        # Define test parameters
        vertices = [[-5, -5], [3, -6], [3, 2], [-2, 3]]
        ideal_distances = [10, 4, 6, 4]
        ideal_angles = [np.deg2rad(90), np.deg2rad(90), np.deg2rad(120), np.deg2rad(60)]
        k_edges = 1
        k_angle = 1
        max_steps = 100
        dt = 0.05
        optimizer = "gd"
        use_adam_optimizer = False 
        log_energy = "log_energy.txt"
        log_geom = "log_geom.txt"
        n_geom = 10
        # Initialize the simulation object
        self.sim = PolygonSpringSimulation(vertices, ideal_distances, ideal_angles, k_angle, k_edges, max_steps, dt,
            optimizer, use_adam_optimizer, log_energy, log_geom, n_geom)

    def test_calculate_angle(self):
        # Test calculate_angle method
        angle = self.sim.calculate_angle(np.array([0,0]), np.array([1,0]), np.array([1,1]))
        self.assertAlmostEqual(angle, np.pi/2)  # Expected angle for these three points

    def test_calculate_energy(self):
        # Test calculate_energy method
        energy = self.sim.calculate_energy(self.sim.vertices)
        # You may need to define an expected value based on your specific simulation
        # For this example, I'll just assert that the energy is a float
        self.assertIsInstance(energy, float)

    def test_bfgs_optimizer(self):
        # Test BFGS optimizer
        optimized_vertices = self.sim.bfgs_optimizer()
        self.assertEqual(optimized_vertices.shape, (self.sim.num_vertices, 2))

    def test_update_position(self):
        # Test update_position method
        initial_vertices = np.copy(self.sim.vertices)
        self.sim.update_position()
        updated_vertices = self.sim.vertices
        # Make sure the vertices are updated after calling update_position
        self.assertFalse(np.allclose(initial_vertices, updated_vertices))

    def test_should_stop(self):
        # Test should_stop method
        self.assertFalse(self.sim.should_stop())  # As it should not stop initially

class TestDiskEnergyCalculator(unittest.TestCase):

    def test_calc_comb_energy(self):
        # Test calc_comb_energy method
        energy1 = DiskEnergyCalculator.calc_comb_energy([1, 1, 2], 4)
        energy2 = DiskEnergyCalculator.calc_comb_energy([1, 2, 1], 4)
        self.assertAlmostEqual(energy1, energy2, places=5)

    def test_modify_one_length(self):
        # Test modify_one_length method
        combination = DiskEnergyCalculator.ProteinConfiguration([1, 2, 6, 12], 17)
        length = combination.map_combination_to_lengths()
        modified_length = DiskEnergyCalculator.ProteinAnalyzer.modify_one_length(length, disk_radius=7, L=2)
        expected_modified_length = [14, 4, 56, 4, 84, 4, 84, 4]
        self.assertEqual(modified_length, expected_modified_length)

    def test_is_valid_inequality(self):
        # Test is_valid_inequality method
        lengths_true = [14, 4, 56, 4, 84, 4, 84, 4]
        self.assertTrue(DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(lengths_true, L=2))

        lengths_false = [14, 4, 14, 4, 210, 4]
        self.assertFalse(DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(lengths_false, L=2))
    import unittest

    def test_cyclic_filter(self):
        # Define some sample data
        combinations = [[1, 3, 5, 6], [2, 4, 6, 8], [3, 6, 9, 12]]
        n_disks = 4

        # Call the function under test
        filtered_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(combinations, n_disks)

        # Define the expected output after filtering
        expected_filtered_combinations = [[1, 3, 5, 6], [2, 4, 6, 8]]

        # Assert that the filtered combinations match the expected output
        self.assertEqual(filtered_combinations, expected_filtered_combinations)

    def test_cyclic_filter_empty_input(self):
        # Test with empty input
        combinations = []
        n_disks = 4

        # Call the function under test
        filtered_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(combinations, n_disks)

        # The result should be an empty list
        self.assertEqual(filtered_combinations, [])

    def test_cyclic_filter_single_combination(self):
        # Test with a single combination
        combinations = [[1, 2, 3, 4]]
        n_disks = 4

        # Call the function under test
        filtered_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(combinations, n_disks)

        # The result should be the same as the input since there are no cyclic duplicates
        self.assertEqual(filtered_combinations, combinations)

    def test_cyclic_filter_duplicates(self):
        # Test with duplicate combinations
        combinations = [[1, 2, 3, 4], [1, 3, 2, 4], [4, 3, 2, 1]]
        n_disks = 4

        # Call the function under test
        filtered_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(combinations, n_disks)

        # The result should contain only one unique combination
        self.assertEqual(len(filtered_combinations), 1)
        self.assertEqual(filtered_combinations[0], [1, 2, 3, 4])


if __name__ == '__main__':
    unittest.main()
