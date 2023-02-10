
# Import general usage libraries
import os
from typing import Callable

# Import general usage scientific libraries
from numpy import linspace, ndarray

# Import local modules
from base_integrals import (L1, L1_dA, L1_dB, L2, L3, L3_dA, L3_dB,
                            M1, M1_dA, M1_dB, M2, M3, M3_dA, M3_dB)


class Bounds:
    def __init__(self) -> None:
        self.a_max = 0.0
        self.a_min = 0.0
        self.b_max = 0.0
        self.b_min = 0.0
        self.h_max = 0.0
        self.h_min = 0.0
        self.num_a = 10
        self.num_b = 10
        self.num_h = 10


def create_reference_database(folder_path: str)->None:
    # Create L1 reference data
    bnd_l1 = Bounds()
    bnd_l1.a_max = 1.0
    bnd_l1.a_min = 0.0
    bnd_l1.b_max = 1.0
    bnd_l1.b_min = 0.0
    bnd_l1.h_max = 1.0
    bnd_l1.h_min = 1e-16
    write_reference_file_3d(folder_path, "L1", L1, bnd_l1)

    # Create L1_dA reference data
    bnd_l1_dA = Bounds()
    bnd_l1_dA.a_max = 1.0
    bnd_l1_dA.a_min = 0.0
    bnd_l1_dA.b_max = 1.0
    bnd_l1_dA.b_min = 0.0
    bnd_l1_dA.h_max = 1.0
    bnd_l1_dA.h_min = 1e-16
    write_reference_file_3d(folder_path, "L1_dA", L1_dA, bnd_l1_dA)

    # Create L1_dB reference data
    bnd_l1_dB = Bounds()
    bnd_l1_dB.a_max = 1.0
    bnd_l1_dB.a_min = 0.0
    bnd_l1_dB.b_max = 1.0
    bnd_l1_dB.b_min = 0.0
    bnd_l1_dB.h_max = 1.0
    bnd_l1_dB.h_min = 1e-16
    write_reference_file_3d(folder_path, "L1_dB", L1_dB, bnd_l1_dB)

    # Create L2 reference data
    bnd_l2 = Bounds()
    bnd_l2.h_max = 1.0
    bnd_l2.h_min = 1e-16
    write_reference_file_1d(folder_path, "L2", L2, bnd_l2)

    # Create L3 reference data
    bnd_l3 = Bounds()
    bnd_l3.a_max = 1.0
    bnd_l3.a_min = 0.0
    bnd_l3.b_max = 1.0
    bnd_l3.b_min = 0.0
    bnd_l3.h_max = 1000.0
    bnd_l3.h_min = 1.0
    write_reference_file_3d(folder_path, "L3", L3, bnd_l3)

    # Create L3_dA reference data
    bnd_l3_dA = Bounds()
    bnd_l3_dA.a_max = 1.0
    bnd_l3_dA.a_min = 0.0
    bnd_l3_dA.b_max = 1.0
    bnd_l3_dA.b_min = 0.0
    bnd_l3_dA.h_max = 1000.0
    bnd_l3_dA.h_min = 1.0
    write_reference_file_3d(folder_path, "L3_dA", L3_dA, bnd_l3_dA)

    # Create L3_dB reference data
    bnd_l3_dB = Bounds()
    bnd_l3_dB.a_max = 1.0
    bnd_l3_dB.a_min = 0.0
    bnd_l3_dB.b_max = 1.0
    bnd_l3_dB.b_min = 0.0
    bnd_l3_dB.h_max = 1000.0
    bnd_l3_dB.h_min = 1.0
    write_reference_file_3d(folder_path, "L3_dB", L3_dB, bnd_l3_dB)

    # Create M1 reference data
    bnd_m1 = Bounds()
    bnd_m1.a_max = 1.0
    bnd_m1.a_min = 0.0
    bnd_m1.b_max = 2.0
    bnd_m1.b_min = 1.0
    bnd_m1.h_max = 1.0
    bnd_m1.h_min = 1e-16
    write_reference_file_3d(folder_path, "M1", M1, bnd_m1)

    # Create M1_dA reference data
    bnd_m1_dA = Bounds()
    bnd_m1_dA.a_max = 1.0
    bnd_m1_dA.a_min = 0.0
    bnd_m1_dA.b_max = 2.0
    bnd_m1_dA.b_min = 1.0
    bnd_m1_dA.h_max = 1.0
    bnd_m1_dA.h_min = 1e-16
    write_reference_file_3d(folder_path, "M1_dA", M1_dA, bnd_m1_dA)

    # Create M1_dB reference data
    bnd_m1_dB = Bounds()
    bnd_m1_dB.a_max = 1.0
    bnd_m1_dB.a_min = 0.0
    bnd_m1_dB.b_max = 2.0
    bnd_m1_dB.b_min = 1.0
    bnd_m1_dB.h_max = 1.0
    bnd_m1_dB.h_min = 1e-16
    write_reference_file_3d(folder_path, "M1_dB", M1_dB, bnd_m1_dB)

    # Create M2 reference data
    bnd_m2 = Bounds()
    bnd_m2.h_max = 1.0
    bnd_m2.h_min = 1e-16
    write_reference_file_1d(folder_path, "M2", M2, bnd_m2)

    # Create M3 reference data
    bnd_m3 = Bounds()
    bnd_m3.a_max = 1.0
    bnd_m3.a_min = 0.0
    bnd_m3.b_max = 2.0
    bnd_m3.b_min = 1.0
    bnd_m3.h_max = 1000.0
    bnd_m3.h_min = 1.0
    write_reference_file_3d(folder_path, "M3", M3, bnd_m3)

    # Create M3_dA reference data
    bnd_m3_dA = Bounds()
    bnd_m3_dA.a_max = 1.0
    bnd_m3_dA.a_min = 0.0
    bnd_m3_dA.b_max = 2.0
    bnd_m3_dA.b_min = 1.0
    bnd_m3_dA.h_max = 1000.0
    bnd_m3_dA.h_min = 1.0
    write_reference_file_3d(folder_path, "M3_dA", M3_dA, bnd_m3_dA)

    # Create M3_dB reference data
    bnd_m3_dB = Bounds()
    bnd_m3_dB.a_max = 1.0
    bnd_m3_dB.a_min = 0.0
    bnd_m3_dB.b_max = 2.0
    bnd_m3_dB.b_min = 1.0
    bnd_m3_dB.h_max = 1000.0
    bnd_m3_dB.h_min = 1.0
    write_reference_file_3d(folder_path, "M3_dB", M3_dB, bnd_m3_dB)


def write_reference_file_1d(folder_path: str, file_name: str, f_def: Callable, bounds: Bounds)->None:
    file_path = os.path.join(folder_path, file_name+".dat")
    with open(file_path, "w") as fid:
        # Write number of test points
        fid.writelines(f"NumPoints {bounds.num_h:d}\n")
        
        # Write header
        header_labels = ["H", "G_int"]
        fid.writelines(" ".join([f"{si:^22s}" for si in header_labels]) + "\n")

        # Loop over bounds to write the reference data
        h_values = linspace(bounds.h_min, bounds.h_max, bounds.num_h)
        for k in range(bounds.num_h):
            vtw = [
                    h_values[k],
                    f_def(h_values[k])[0].real
                    ]
            fid.writelines("      ".join(f"{fi:^0.16f}" for fi in vtw) + "\n")  


def write_reference_file_3d(folder_path: str, file_name: str, f_def: Callable, bounds: Bounds)->None:
    file_path = os.path.join(folder_path, file_name+".dat")
    with open(file_path, "w") as fid:
        # Write number of test points
        num_points = bounds.num_a*bounds.num_b*bounds.num_h
        fid.writelines(f"NumPoints {num_points:d}\n")

        # Write header
        header_labels = ["A", "B", "H", "G_int"]
        fid.writelines(" ".join([f"{si:^22s}" for si in header_labels]) + "\n")

        # Loop over bounds to write the reference data
        a_values = linspace(bounds.a_min, bounds.a_max, bounds.num_a)
        b_values = linspace(bounds.b_min, bounds.b_max, bounds.num_b)
        h_values = linspace(bounds.h_min, bounds.h_max, bounds.num_h)
        for i in range(bounds.num_a):
            for j in range(bounds.num_b):
                for k in range(bounds.num_h):
                    vtw = [
                            a_values[i],
                            b_values[j],
                            h_values[k],
                            f_def(a_values[i], b_values[j], h_values[k])[0].real
                            ]
                    fid.writelines("      ".join(f"{fi:^0.16f}" for fi in vtw) + "\n")

if __name__ == "__main__":
    this_path = os.path.dirname(os.path.abspath(__file__))
    target_folder = os.path.join(os.path.dirname(this_path), "tests", "tests_data", "green_fin_depth_tables")
    create_reference_database(target_folder)