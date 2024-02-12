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
                if DiskEnergyCalculator.ProteinAnalyzer.is_valid_inequality(length_list, L) and combination[0] == 1:
                    valid_combinations.append(combination)

            valid_combinations = DiskEnergyCalculator.ProteinAnalyzer.cyclic_filter(valid_combinations, n_disks)

            return valid_combinations

    @staticmethod
    def F_dev(L, phi, phi_s):
        return 2 * (0.5 * 10 * 1/np.tanh(L) * (phi - phi_s) ** 2)

    @staticmethod
    def F_tot(L, R=7, l=2):
        # VdW energy
        A = 1  # Hamaker constant   [kT]
        h = 2  # monolayer height   [nm]
        B = 1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l)
        F_vdw = - B / (L*2 / l) ** (3/2)
        
        # Elastic energy TODO: explicit C
        C = 0.4
        F_el = - C / (1/np.tanh(R/l) + 0.5 / np.tanh(L/l))
        return F_vdw + F_el

    @staticmethod
    def calc_comb_energy(length_list: list, L=2, eta=0.09) -> float:
        comb_tot_energy = 0.0

        length_list = np.array(length_list) // 14
        non_zero = np.count_nonzero(length_list)
        comb_tot_energy += non_zero * DiskEnergyCalculator.F_tot(L=2)
        
        n_close_inter = length_list[length_list != 0] - 1
        comb_tot_energy += n_close_inter.sum() * DiskEnergyCalculator.F_tot(L=eta/2)
        return comb_tot_energy

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
    e_tot = DiskEnergyCalculator.F_tot
    assert round(e_tot(L=2),3) == -0.242
    
    lengths = [14, 4, 14, 4, 14, 4]
    e_lengths = DiskEnergyCalculator.calc_comb_energy(lengths)
    assert round(e_lengths / 3, 3) == -0.242

    combinations = [[1, 2, 3], [1, 2, 4], [2, 3, 4], [1, 2, 3, 4]]
    lengths = DiskEnergyCalculator.ProteinAnalyzer.map_many_combinations(combinations, 4)
    phi = 0.1206

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