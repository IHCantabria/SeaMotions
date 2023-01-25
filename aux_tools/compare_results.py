
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from numpy import array, fromstring, pi, sqrt, zeros

# Import local modules
from hydlib.waves.airy import w2k
from base_integrals import G_integral


class RefData(object):

    def __init__(self) -> None:
        self.A = array([])
        self.grav_acc = 0.0
        self.green = array([])
        self.H = array([])
        self.num_A = 0
        self.num_H = 0
        self.num_z = 0
        self.water_depth = 0.0
        self.z = array([])
        self.zeta = array([])

    def load_data(self, file_path: str)->None:
        # Read data from disk
        with open(file_path, "r") as fid:
            # Read gravitational acceleration
            self.grav_acc = float(fid.readline().strip().split()[1])

            # Read water depth
            self.water_depth = float(fid.readline().strip().split()[1])

            # Read A parameter
            self.num_A = int(fid.readline().strip().split()[1])
            self.A = fromstring(fid.readline().strip(), dtype=float, sep=" ")

            # Read H parameter
            self.num_H = int(fid.readline().strip().split()[1])
            self.H = fromstring(fid.readline().strip(), dtype=float, sep=" ")

            # Read Z parameter
            self.num_z = int(fid.readline().strip().split()[1])
            self.z = fromstring(fid.readline().strip(), dtype=float, sep=" ")
            self.zeta = fromstring(fid.readline().strip(), dtype=float, sep=" ")
            
            # Read Green function data header
            fid.readline()
            
            # Read Greeen function data
            lines = fid.readlines()
            self.green = zeros((self.num_A*self.num_H*self.num_z, ), dtype=complex)
            for i, line_i in enumerate(lines):
                a = fromstring(line_i.strip(), dtype=float, sep=" ")
                self.green[i] = complex(a[0], a[1])

    def print_out(self)->None:
        print(f"Grav.Acc: {self.grav_acc}")
        print(f"Water Depth: {self.water_depth}")
        print(f"Num.A: {self.num_A}")
        print("A: ", self.A)
        print(f"Num.H: {self.num_H}")
        print("H: ", self.H)
        print(f"Num.z: {self.num_z}")
        print("z: ", self.z)
        print("zeta: ", self.zeta)
        print("G_integral")
        for i in range(self.green.shape[0]):
            print(f"{self.green[i].real:0.16f} {self.green[i].imag:0.16f}i")


def compare(file_path: str, f_def: Callable, f_name: str)->None:
    tol = 1e-7

    # Read reference data
    ref_data = RefData()
    ref_data.load_data(file_path)

    # Loop over reference data to check the value of 
    # the integrals performed
    for i in range(ref_data.num_A):
        R = ref_data.A[i]*ref_data.water_depth

        for j in range(ref_data.num_H):
            T = 2*pi*sqrt(ref_data.water_depth/ref_data.H[j]/ref_data.grav_acc)

            for k in range(ref_data.num_z):
                # Calculate integral values
                gc = G_integral(
                                R,
                                ref_data.z[k],
                                ref_data.zeta[k],
                                T,
                                ref_data.water_depth
                                )

                # Get reference data
                index_rd = (
                            + i*(ref_data.num_H*ref_data.num_z)
                            + j*ref_data.num_z
                            + k
                            )
                gr = ref_data.green[index_rd]

                # Compare results
                abs_diff = gc-gr
                rel_diff = (gc-gr)/abs(abs_diff)
                if (abs(abs_diff)>tol) and (abs(rel_diff)>tol):
                    k0 = w2k(2*pi/T, ref_data.water_depth, method="bisection")
                    B2 = (ref_data.z[k]+ref_data.zeta[k]+2*ref_data.water_depth)/ref_data.water_depth
                    raise ValueError(
                                    "The difference between the actual Green function \n"
                                    "calculation and the reference one exceed the predefined \n"
                                    f"tolerance for the function: {f_name}.\n"
                                    + "Input Parameters:\n"
                                    + f"-> A: {ref_data.A[i]}\n"
                                    + f"-> B2: {B2}\n"
                                    + f"-> H: {ref_data.H[j]}\n"
                                    + f"-> R: {R}\n"
                                    + f"-> Z: {ref_data.z[k]}\n"
                                    + f"-> Zeta: {ref_data.zeta[k]}\n"
                                    + f"-> Wave Period: {T}\n" 
                                    + f"-> Water Depth: {ref_data.water_depth}\n"
                                    + f"-> Grav.Acc: {ref_data.grav_acc}\n"
                                    + f"-> nu: {(2*pi/T)**2.0/ref_data.grav_acc}\n"
                                    + f"-> k0: {k0:0.16f}\n"
                                    + "Numerical Error Details:\n"
                                    + f"---> Gc: {gc.real:0.16f} {gc.imag:0.16f}i\n"
                                    + f"---> Gr: {gr.real:0.16f} {gr.imag:0.16f}i\n"
                                    + f"---> Abs.Diff: {abs_diff.real} {abs_diff.imag}i - |Abs.Diff|: {abs(abs_diff)}\n"
                                    + f"---> Rel.Diff: {rel_diff.real} {rel_diff.imag}i - |Rel.Diff|: {abs(rel_diff)}\n"
                                    )


if __name__ == "__main__":
    file_path = r"E:\sergio\developments\SeaMotions\tests\tests_data\green_fin_depth_tables\fin_depth_integral.dat"
    compare(file_path, G_integral, "G")