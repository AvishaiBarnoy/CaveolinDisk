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
    
    def partitions(n):
        """
        generates all unique integer partitions
        """
        if n == 0:
            yield ()
        else:
            for i in range(1, n + 1):
                for p in Combinatorics.partitions(n - i):
                    yield(i,) + p

class Combination:
    def __init__(self, combination=[], n_disks=0, lengths=[],
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
        assert sum(self.lengths) == self.n_disks, f"{self.lengths}, {self.n_disks}, something went wrong with conversion to length list"
        return self.lengths
    
    def modify_one_length(self, disk_radius=7, L=2) -> list:
        L *= 2
        disk_radius *= 2
        self.mod_length = []
        for l in self.lengths:
            self.mod_length.extend([l * disk_radius, L])
        return self.mod_length

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
    # TODO: possibly restrutcutre, though Combination class is for a specific combinatio
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
    def cyclic_partitions(lsts):
        # TODO: should change to work on lengths and not on combinations
        unique = set()

        # print(lsts, type(lsts))
        lsts = [list(lst) for lst in lsts]
        for sublist in lsts:
            original_lst = tuple(sublist)
        
            rotations = []
            sublist = list(sublist)
            while True:
                sublist = [sublist[-1]] + sublist[:-1]  # Rotate the list to the right
                rotations.append(tuple(sublist))
            
                if tuple(sublist) in unique: # or tuple(sublist) == original_lst: # rotations:
                    break

                if tuple(sublist) == original_lst and len(rotations) > 1:
                    unique.add(original_lst)
                    break

        return unique

    @staticmethod
    def filter_partitions(n):
        """
        filters partitions
        """
        partition_generator = Combinatorics.partitions(n) 
        partition_count = 0 
        partLen3 = []
        while True:
            try:
                partition = next(partition_generator)
                if partition == [2,2,2]:
                    print("found [2,2,2]")
                partition_count += 1
                if len(partition) >= 3: # toss-out partitions shorter than 3 vertex 
                    partLen3.append(partition)
            except StopIteration:
                break
    
        unique = []
        for i in partLen3:
            permutations = set(it.permutations(i))
            sub_uniques = group_operations.cyclic_partitions(permutations)
            unique.extend(sub_uniques)
        return unique 

if __name__ == "__main__":
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

    print("Partition:")
    n = 17 
    partitions = Combinatorics(n_disks=n)
    partitions.generate_combinations()
    combinations = partitions.raw_combinations
    filtered_combinations = [combination for combination in combinations if combination[0] == 1]
 
    lengths = group_operations.map_many_combinations(filtered_combinations, n_disks=n)
    
    cyclic = group_operations.cyclic_filter(filtered_combinations, n_disks=n)
    cyclic_lengths = group_operations.map_many_combinations(cyclic, n_disks=n)
    
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
        pass
    
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
