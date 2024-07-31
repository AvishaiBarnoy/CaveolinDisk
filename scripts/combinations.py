import numpy as np
import itertools as it

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
    def __init__(self, combination=[], n_disks=0, lengths=[], mod_length=[], energy=0, points=[]):
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

    def modify_one_length(self, L, disk_radius=7) -> list:
        """
        converts from combination to actual lengths
        """
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
    def cyclic_partitions(lsts):
        unique = set()

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
        partLen3 = []
        while True:
            try:
                partition = next(partition_generator)
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

def calc_ideal_angle(L, R=7, k=0.4e-19, kt=20e-3, h=2, a=0.7):
    '''
    f_param = Delta_epsilon / sqrt(k*k_t) * h / a
    D_epsilon - energy diff of lipids on-top of protein and in membrane ~ 0.3
    h - monolayer thickness ~ 2 nm, a - lipid length ~ 0.7 nm
    k - monolayer rigidiy ~ 1e-9 J, k_t - tilt modulus ~ 30 mJ/m
    '''
    xi = np.sqrt(k/kt) * 1e9 # J / nm
    depsilon = 4 # kT/nm
    f_param = h/a * depsilon * 1/np.sqrt(k*kt/1e18) * 4.11e-21
    return np.pi - f_param * 1 / (2/np.tanh(R/xi) + 1/np.tanh(L/xi))

def caveolin_radius(L, R=7, k=0.4e-19, kt=20e-3):
    phi = calc_ideal_angle(L=L, R=R, k=k, kt=kt)
    R_c = (R+L*np.cos(np.pi-phi))/np.sin(np.pi-phi)
    return R_c

def calc_n_disks(L, R):
    Rc = caveolin_radius(L=L, R=R)
    circumference = 2*np.pi*Rc
    n_disks = circumference/(2*R + 2*L)
    return round(n_disks)

if __name__ == "__main__":
    R = 7
    L = 2

    phi_s = calc_ideal_angle(L=L, R=R)
    print("phi*:", phi_s)
    Rc = caveolin_radius(L=L, R=R)
    print("Rc:", Rc)
    n_disks = calc_n_disks(L=L, R=R)
    print(n_disks, round(n_disks))
