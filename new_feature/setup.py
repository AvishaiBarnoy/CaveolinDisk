import numpy as np
import itertools as it
import matplotlib.pyplot as plt

class Combinatorics:
    def __init__(self, n_disks, raw_combinations=[], disk_radius=7, L=2):
        self.n_disks = n_disks
        self.raw_combinations = []
    
    def generate_combinations(self) -> list:
        disks = list(range(1, self.n_disks + 1))
        all_combinations = [list(it.combinations(disks, r)) for r in range(3, len(disks) + 1)]
        self.raw_combinations = [combination for sublist in all_combinations for combination in sublist]
        return self.raw_combinations

class Combination:
    def __init__(self, combination: list, n_disks: int, lengths=[],
                 mod_length=[], energy=0, points=[]):
        self.combination = combination
        self.n_disks = n_disks 
        self.lengths = lengths
        self.energy = energy
        self.mod_length = mod_length
        self.points = points

    def map_combination_to_lengths(self) -> list:
        """ 
        takes combination [1, 3, 5, 6] and n_disks=6
        return length list unmodified [2, 2, 1, 1]
        """
        # TODO: unite with modify_one_length
        combination = np.array(self.combination)
        self.lengths = np.diff(combination).tolist()
        last_length = self.n_disks - combination[-1] + combination[0]
    
        if combination[-1] == self.n_disks:
            self.lengths.append(1)
        else:
            self.lengths.append(last_length)
    
        if combination[0] != 1:
            self.lengths[-1] += combination[0] - 1
        # print(self.combination, self.lengths)
        # print(f"lengths {combination}\t{self.lengths}")
        assert sum(self.lengths) == self.n_disks, f"{self.lengths}, {self.n_disks}, something went wrong with conversion to length list"
        return self.lengths
    
    def modify_one_length(self, disk_radius=7, L=2) -> list:
        L *= 2
        disk_radius *= 2
        self.mod_length = []
        for l in self.lengths:
            self.mod_length.extend([l * disk_radius, L])
        return self.mod_length

    def F_tot(self, L, R=7, l=2):
        A = 1  
        h = 2  
        B = 1 / 270 * A * h / l ** 2 * 1 / np.sqrt(R / l)
        F_vdw = - B / (L*2 / l) ** (3/2)
        
        C = 0.4
        F_el = - C / (1/np.tanh(R/l) + 0.5 / np.tanh(L/l))
        return F_vdw + F_el

    def calc_comb_energy(self) -> float:
        L = 2
        eta = 0.09
        length_list = np.array(self.mod_length) // 14
        non_zero = np.count_nonzero(length_list)
        self.energy += non_zero * self.F_tot(L=2)
        
        n_close_inter = length_list[length_list != 0] - 1
        self.energy += n_close_inter.sum() * self.F_tot(L=eta)
        return self.energy

    def is_valid_inequality(self, L=2) -> bool:
        total_length = sum(self.mod_length)
        for i in self.mod_length:
            sum_lengths = 0.5 * (total_length - i)
            if i > sum_lengths:
                return False
        return True
    
    def check_circle(self, radius: float) -> float:
        sum_of_angles = 0
        for length in self.mod_length:
            angle = 2 * np.arcsin(length / (2 * radius))
            sum_of_angles += angle
        return sum_of_angles

    def calculate_circle_points(self) -> list[tuple[float, float]]:
        min_radius = max(self.mod_length) / 2
        max_radius = sum(self.mod_length) / 2

        if min_radius > max_radius:
            return []
            raise ValueError("Segments are too long to fit in a circle")

        x_0 = (min_radius + max_radius) / 2
        d_0 = ((max_radius - min_radius) / 2) * 0.9

        while d_0 > 1e-9:
            current_angle = self.check_circle(x_0)
            if current_angle < 2 * np.pi:
                x_0 -= d_0
            else:
                x_0 += d_0
            d_0 /= 1.8
        # print(f"The smallest radius that can fit all segments is {x_0}")
        angle = 0.0
        self.points = []
        for length in self.mod_length:
            x = x_0 * np.cos(angle)
            y = x_0 * np.sin(angle)
            self.points.append((x, y))
            angle += 2 * np.arcsin(length / (2 * x_0))

        return self.points

class group_operations:
    @staticmethod
    def map_many_combinations(combinations: list, n_disks: int) -> list:
        lengths = [Combination(combination, n_disks).map_combination_to_lengths() for combination in combinations]
        return lengths
        
    @staticmethod
    def cyclic_filter(combinations, n_disks):
        lengths = group_operations.map_many_combinations(combinations, n_disks)
        filtered = []
        already = set()  

        for i, length in enumerate(lengths):
            rotations = {tuple(length[j:] + length[:j]) for j in range(len(length))}
            mirror_rotations = {tuple(length[::-1][j:] + length[::-1][:j]) for j in range(len(length))}
            all_rotations = rotations.union(mirror_rotations)  
            if not any(rot in already for rot in all_rotations):
                already.update(all_rotations)  
                filtered.append(combinations[i])

        return filtered

    @staticmethod
    def filter_valid_combinations(combinations, length_list: list, n_disks: int, L=2) -> list:
        valid_combinations = []
        combinations = [] 
        for combination in combinations:
            if Combination.is_valid_inequality(length_list, L) and combination[0] == 1:
                valid_combinations.append(combination)

        valid_combinations = group_operations.cyclic_filter(valid_combinations, n_disks)

        return valid_combinations

    @staticmethod
    def calc_many_combinations(lengths: list, combinations: list, eta) -> float:
        combinations_energies = []
        for i, j in enumerate(lengths):
            energy_i = Combination.calc_comb_energy(j, eta)
            combinations_energies.append((combinations[i], energy_i))
        return combinations_energies

if __name__ == "__main__":
    # check energy of 1 protein
    comb = Combination([1], n_disks=1)
    lengths = comb.map_combination_to_lengths()
    mod_length = comb.modify_one_length()
    energy = comb.calc_comb_energy()
    assert round(comb.energy,3) == -0.242
   
    # check energy of 3 proteins
    comb = Combination([1,2,3], n_disks=3)
    lengths = comb.map_combination_to_lengths()
    mod_length = comb.modify_one_length()
    energy_lengths = comb.calc_comb_energy()
    assert round(energy_lengths / 3, 3) == -0.242
    
    comb = Combination(list(range(1,17)), n_disks=17)
    comb.map_combination_to_lengths()
    comb.modify_one_length()
    comb.calc_comb_energy()
    # print(comb.mod_length)
    # print(comb.energy)

    # check conversion of combination and inequality
    comb = Combination([1, 2, 6, 12], 17)
    length = comb.map_combination_to_lengths()
    mod_length = comb.modify_one_length()
    assert length == [1, 4, 6, 6]
    assert mod_length == [14, 4, 56, 4, 84, 4, 84, 4]
    assert comb.is_valid_inequality(L=2) == True
    
    # check a False inequality
    comb = Combination([1, 2, 3], 17)
    length = comb.map_combination_to_lengths()
    mod_length = comb.modify_one_length()
    assert mod_length == [14, 4, 14, 4, 210, 4]
    assert comb.is_valid_inequality(L=2) == False

    combinations = [[1, 2, 3], [1, 2, 4], [1, 2, 3, 4]]
    lengths = group_operations.map_many_combinations(combinations, 4)
    assert lengths == [[1, 1, 2], [1, 2, 1], [1, 1, 1, 1]]
    cyclic = group_operations.cyclic_filter(combinations, 4)
    assert cyclic == [[1,2,3], [1,2,3,4]]
    # assert round(group_operations.calc_many_combinations(lengths, combinations)[0][1], 5) == round(group_operations.calc_many_combinations(lengths, combinations)[1][1], 5), "Fail for degenerate combinations"
   

    print("Partition:")
    n = 17 
    partitions = Combinatorics(n_disks=n)
    partitions.generate_combinations()
    combinations = partitions.raw_combinations
    filtered_combinations = [combination for combination in combinations if combination[0] == 1]
 
    lengths = group_operations.map_many_combinations(filtered_combinations, n_disks=n)
    
    cyclic = group_operations.cyclic_filter(filtered_combinations, n_disks=n)
    cyclic_lengths = group_operations.map_many_combinations(cyclic, n_disks=n)
    # print(cyclic)
    # print(cyclic_lengths)
    # print(f"len cyclics {len(cyclic_lengths)}")
    
    cyclic_mod = []
    for i in cyclic:
        comb = Combination(i, n)
        length = comb.map_combination_to_lengths()
        mod_length = comb.modify_one_length()
        cyclic_mod.append(mod_length)
    # print(cyclic_mod)
    dn4 = []
    for i in cyclic_lengths:
        if len(i) == 13:
            dn4.append(i)
    for i in dn4:
        # print(i)
        pass
    # print(len(dn4))

    
    l = [28,4,28,4,14,4,14,4,14,4,14,4,28,4,28,4,14,4,14,4,14,4,14,4,14,4]
    comb = [1,3,5,6,7,8,9,11,13,14,15,16,17]
    print('l', len(l), len(l)/2, (34-len(l))/2)
    comb = Combination(comb, n_disks=17)
    comb.map_combination_to_lengths()
    comb.modify_one_length()
    print(comb.is_valid_inequality())
    print(comb.lengths)
    print(comb.mod_length)
    print(comb.calculate_circle_points())
