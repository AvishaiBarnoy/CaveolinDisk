import numpy as np
import itertools as it
import math
import sys
import matplotlib.pyplot as plt

class DiskEnergyCalculator:
    class ProteinConfiguration:
        def __init__(self, combination: list, n_disks: int):
            self.combination = combination
            self.n_disks = n_disks 

        def map_combination_to_lengths(self) -> list:
            combination = np.array(self.combination)
            lengths = np.diff(combination).tolist()
            last_length = self.n_disks - combination[-1] + combination[0]
    
            if combination[-1] == self.n_disks:
                lengths.append(1)
            else:
                lengths.append(last_length)
    
            if combination[0] != 1:
                lengths[-1] += combination[0] - 1
    
            return lengths
        # old implementation:
            lengths = [self.combination[i + 1] - self.combination[i] for i in range(len(self.combination) - 1)]
            lengths.append(self.n_disks - self.combination[-1] + self.combination[0])
            if self.combination[-1] == self.n_disks:
                lengths[-1] = 1
            if self.combination[0] != 1:
                lengths[-1] += self.combination[0] - 1
            return lengths

    class ProteinAnalyzer:
        @staticmethod
        def generate_combinations(n_disks):
            disks = list(range(1, n_disks + 1))
            all_combinations = [list(it.combinations(disks, r)) for r in range(3, len(disks) + 1)]
            valid_combinations = [combination for sublist in all_combinations for combination in sublist]
            return valid_combinations
        # old implementation:
            disks = list(range(1, n_disks + 1))
            all_combinations = []
            for r in range(1, len(disks) + 1):
                combinations_object = it.combinations(disks, r)
                combinations_list = list(combinations_object)
                all_combinations += combinations_list
            valid_combinations = [combination for combination in all_combinations if len(combination) >= 3]
            return valid_combinations
            # return  DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(valid_combinations, n_disks)

        @staticmethod
        def map_many_combinations(combinations: list, n_disks: int) -> list:
            lengths = [DiskEnergyCalculator.ProteinConfiguration(combination, n_disks).map_combination_to_lengths() for combination in combinations]
            return lengths
        # old implementation:
            lengths = []
            for i in combinations:
                configuration = DiskEnergyCalculator.ProteinConfiguration(i, n_disks)
                lengths.append(configuration.map_combination_to_lengths())
            return lengths
        
        @staticmethod
        def modify_one_length(length_list, disk_radius=7, L=2) -> list:
            modified_length = []
            L *= 2
            disk_radius *= 2
            for l in length_list:
                modified_length.extend([l * disk_radius, L])
                #modified_length.append(l * disk_radius * 2)
                #modified_length.append(L * 2)
            return modified_length

        @staticmethod
        def modify_many_lengths(lengths: list, disk_radius=7, L=2) -> list:
            modified_lengths = []
            for length in lengths:
                modified_length = []
                L *= 2
                disk_radius *= 2
                for l in length:
                    modified_length.extend([l * disk_radius, L])
                    modified_lengths.append(modified_length)
            return modified_lengths
            # new implementation
            modified_lengths = [DiskEnergyCalculator.ProteinAnalyzer.modify_one_length(length) for length in lengths]
            return modified_lengths
            # old implementation:
            for length in lengths:
                modified_length = []
                for l in length:
                    modified_length.append(l * disk_radius * 2)
                    if l != length[-1]:
                        modified_length.append(2 * L)
                modified_lengths.append(modified_length)
            return modified_lengths

        @staticmethod
        def is_valid_inequality(length_list: list, L=2) -> bool:
            total_length = sum(length_list)
            for i in length_list:
                sum_lengths = 0.5 * (total_length - i)
                if i > sum_lengths:
                    return False
            return True
            # old implementation:
            total_length = sum(length_list)
            for i in length_list:
                sum_lengths = 0.5 * (total_length - i)
                if i <= sum_lengths:
                    pass
                elif i > sum_lengths:
                    return False
            return True

        @staticmethod
        def rotate_list(combination):
            """
            Warning: Depracted 
            Rotate a list until it returns to its original configuration.
            """
            rotations = [combination]
            for _ in range(1, len(combination)):
                combination = combination[1:] + [combination[0]]
                rotations.append(combination)
            return rotations

        @staticmethod
        def cyclic_filter(combinations, n_disks):
            """
            Takes combinations, converts to lengths, and removes cyclic degenerate combinations
            """
            lengths = DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(combinations, n_disks)
            filtered = []
            already = set()  # Use a set for faster membership checks

            for i, length in enumerate(lengths):
                # Generate all rotations of the current combination
                rotations = {tuple(length[j:] + length[:j]) for j in range(len(length))}
                mirror_rotations = {tuple(length[::-1][j:] + length[::-1][:j]) for j in range(len(length))}  # Mirror rotations
                all_rotations = rotations.union(mirror_rotations)  # Union of original and mirror rotations
                # Check if any rotation has been seen before
                if not any(rot in already for rot in all_rotations):
                    already.update(all_rotations)  # Update the set with new rotations
                    filtered.append(combinations[i])

            return filtered
            lengths = DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(combinations, n_disks)
            filtered = []
            already = set()  # Use a set for faster membership checks

            for i, length in enumerate(lengths):
                # Generate all rotations of the current combination
                rotations = {tuple(length[j:] + length[:j]) for j in range(len(length))}
                # Check if any rotation has been seen before
                if not any(rot in already for rot in rotations):
                    already.update(rotations)  # Update the set with new rotations
                    filtered.append(combinations[i])

            return filtered

        @staticmethod
        def filter_valid_combinations(combinations: list, length_list: list, n_disks: int, L=2) -> list:
            valid_combinations = []
            combinations = [] 
            for combination in combinations:
                if combination[0] == 1:
                    if DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(length_list, L):
                        valid_combinations.append(combination)

            valid_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(valid_combinations, n_disks)

            return valid_combinations

    @staticmethod
    def F_el(L=2, R=7, l=2):
        l = 2  # decay length of tilt
        L = L  # half distance between proteins
        C = 0.4
        F = -C / (1/np.tanh(R/l) + 0.5 * 1/np.tanh(L/l))
        print('elastic', F)
        return F

    @staticmethod
    def F_dev(L, phi, phi_s):
        return 2 * (0.5 * 10 * 1/np.tanh(L) * (phi - phi_s) ** 2)

    @staticmethod
    def F_vdw(L, R=7, l=2):
        L *= 2  # half distance between proteins
        A = 1 # Hamaker constant
        h = 2 # monolayer height
        B = round(1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l), 3)
        F = -B / (L / l) ** (3 / 2)
        print('vdw', F)
        return F

    @staticmethod
    def F_tot(L, R, l=2):
        return DiskEnergyCalculator.F_el(L, R, l) + DiskEnergyCalculator.F_vdw(L, R, l)

    @staticmethod
    def calc_comb_energy(length_list: list, L=2, eta=0.09) -> float:
        comb_energy = 0.0
        F_disk_14  = DiskEnergyCalculator.F_tot(L=L, R=7, l=2)
        F_disk_eta = DiskEnergyCalculator.F_tot(L=eta, R=7, l=2)
        for i in length_list:
            if i >= 14:
                print(f"len segments {i//14 -1}")
                comb_energy += F_disk_14 + F_disk_eta * (i//14 - 1) 
                print(f"comb energy {comb_energy}") 
        return comb_energy

        # old implementatin: 
        for i in length_list:
            if i == 14:
                F_disk = DiskEnergyCalculator.F_tot(L=L, R=7, l=2)
                combination_total_energy += F_disk
            elif i > 14:
                F_disks = DiskEnergyCalculator.F_tot(L=L, R=7, l=2) + DiskEnergyCalculator.F_tot(L=eta, R=7, l=2) * (i - 1)
                combination_total_energy += F_disks
        return combination_total_energy

    @staticmethod
    def calc_many_combinations(lengths: list, combinations: list, eta) -> float:
        combinations_energies = []
        for i, j in enumerate(lengths):
            energy_i = DiskEnergyCalculator.calc_comb_energy(j, eta)
            combinations_energies.append((combinations[i], energy_i))
        return combinations_energies

    @staticmethod
    def check_circle(lengths: list[float], radius: float) -> float:
        sum_of_angles = 0
        for length in lengths:
            angle = 2 * np.arcsin(length / (2 * radius))
            sum_of_angles += angle
        return sum_of_angles

    @staticmethod
    def calculate_circle_points(lengths: list[float]) -> list[tuple[float, float]]:
        min_radius = max(lengths) / 2
        max_radius = sum(lengths) / 2

        if min_radius > max_radius:
            return []
            raise ValueError("Segments are too long to fit in a circle")

        x_0 = (min_radius + max_radius) / 2
        d_0 = ((max_radius - min_radius) / 2) * 0.9

        while d_0 > 1e-9:
            current_angle = DiskEnergyCalculator.check_circle(lengths, x_0)
            if current_angle < 2 * math.pi:
                x_0 -= d_0
            else:
                x_0 += d_0
            d_0 /= 1.8
        print(f"The smallest radius that can fit all segments is {x_0}")
        angle = 0.0
        points = []
        for length in lengths:
            x = x_0 * np.cos(angle)
            y = x_0 * np.sin(angle)
            points.append((x, y))
            angle += 2 * np.arcsin(length / (2 * x_0))

        return points

if __name__ == "__main__":
    combinations = [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 2, 3, 4]]
    lengths = DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(combinations, 4)
    phi = 0.1206
    ''' 
    assert round(DiskEnergyCalculator.calc_comb_energy([1, 1, 2]),5) == round(DiskEnergyCalculator.calc_comb_energy([1, 2, 1]),5), "Fail for degenerate combinations"
    
    lengths = [[1, 1, 2], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
    modified_lengths = DiskEnergyCalculator.ProteinAnalyzer.modify_many_lengths(lengths, disk_radius=7, L=2)
    # print(modified_lengths)

    combination = DiskEnergyCalculator.ProteinConfiguration([1, 2, 6, 12], 17)
    length = combination.map_combination_to_lengths()
    assert DiskEnergyCalculator.ProteinAnalyzer.modify_one_length(length, disk_radius=7, L=2) == [14, 4, 56, 4, 84, 4, 84, 4]

    lengths = [14, 4, 56, 4, 84, 4, 84, 4]
    mod_lengths = DiskEnergyCalculator.ProteinAnalyzer.modify_one_length(length, disk_radius=7, L=2)
    assert DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(mod_lengths, L=2) == True
    lengths = [14, 4, 14, 4, 210, 4]
    assert DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(lengths, L=2) == False
    '''

    lengths = [14]
    energy = DiskEnergyCalculator.calc_comb_energy(lengths, 1)
    print(energy)
    assert round(energy,2) == -0.24

    import timeit
    new = """
def is_valid_inequality(length_list: list, L=2):
    total_length = sum(length_list)
    for i in length_list:
        sum_lengths = 0.5 * (total_length - i)
        if i > sum_lengths:
            return False
    return True
lengths = [14, 4, 56, 4, 84, 4, 84, 4]
is_valid_inequality(lengths)
lengths = [14, 4, 14, 4, 210, 4]
is_valid_inequality(lengths)
"""

    old = """
def is_valid_inequality(length_list: list, L=2):
    for i in length_list:
        sum_lengths = 0.5 * (sum(length_list) - i)
        if i <= sum_lengths - i:
            pass
        elif i > sum_lengths:
            return False
    return True
lengths = [14, 4, 56, 4, 84, 4, 84, 4]
is_valid_inequality(lengths)
lengths = [14, 4, 14, 4, 210, 4]
is_valid_inequality(lengths)
"""

    # print("new:", timeit.timeit(new, number = 100000))
    # print("old:", timeit.timeit(old, number = 100000))
