
# Import general usage libraries
from typing import Callable

# Import general usage scientific libraries
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from numpy import exp, flip, linspace, log, log10, meshgrid, pi, sign, zeros
from numpy import abs as np_abs
from scipy.special import jv, roots_chebyt, yv

# Import local modules
from base_integrals import fxy
from fit_cheby import eval_chebyshev, eval_chebyshev_filter, fit_chebyshev


class FitProperties:

    def __init__(self) -> None:
        self.x_max = 0
        self.x_map = None
        self.x_min = 0.0
        self.y_max = 0.0
        self.y_map = None
        self.y_min = 0.0
        self.cheby_order = 0
        self.cheby_tol = 0.0


def fit_residual_region_14()->None:
    # Set up domain extension and fit properties
    fit_props = FitProperties()
    fit_props.x_max = 3.0
    fit_props.x_min = 0.3
    fit_props.y_max = 4.0
    fit_props.y_min = 0.3
    fit_props.cheby_order = 20
    fit_props.cheby_tol = 1e-7

    # Launch fit
    fit_residual(resiudal_region_1, fit_props)


def resiudal_region_1(X: float, Y: float)->float:
    c0 = (
        -2*jv(0, X)*log(X+Y)
        -(pi*yv(0, X)-2*jv(0, X)*log(X))
        )*exp(-Y)
    
    return (fxy(X, Y)-c0)*exp(Y)


def fit_residual(f_residual: Callable,
                fit_props: FitProperties
                )->None:
    # Define parametric space limits
    x_max = fit_props.x_max
    x_min = fit_props.x_min
    y_max = fit_props.y_max
    y_min = fit_props.y_min
    cheby_order = fit_props.cheby_order
    cheby_tol = fit_props.cheby_tol

    # Define parametric space for function fit
    num_x_fit = 100
    num_y_fit = 100
    x_fit_poly = roots_chebyt(num_x_fit)[0]
    y_fit_poly = roots_chebyt(num_y_fit)[0]
    x_fit = (x_fit_poly+1)*(x_max-x_min)/2 + x_min
    y_fit = (x_fit_poly+1)*(y_max-y_min)/2 + y_min
    Xfp,Yfp = meshgrid(x_fit_poly, y_fit_poly)
    Xf,Yf = meshgrid(x_fit, y_fit)

    # Calculate residual function over the fit points
    Ff = zeros((num_x_fit, num_y_fit))
    for i in range(num_x_fit):
        for j in range(num_y_fit):
            Ff[i, j] = f_residual(x_fit[i], y_fit[j])

    # Calculate polynomial fit coefficients
    C = fit_chebyshev(Xfp.ravel(), Yfp.ravel(), Ff.ravel(), cheby_order)
    
    # Filter coefficients to the precision required
    n_cheby = linspace(0, cheby_order-1, cheby_order, dtype=int)
    NCX,NCY = meshgrid(n_cheby, n_cheby)
    pos = np_abs(C) > cheby_tol
    C_filter = C[pos]
    NCX_filter = (NCX.ravel())[pos]
    NCY_filter = (NCY.ravel())[pos]

    # Define parametric space
    num_x = 150
    num_y = 100
    xe = linspace(-1.0, 1.0, num_x)
    ye = linspace(-1.0, 1.0, num_y)
    x = (xe+1)*(x_max-x_min)/2.0 + x_min
    y = (ye+1)*(y_max-y_min)/2.0 + y_min
    Xe,Ye = meshgrid(xe, ye)
    X,Y = meshgrid(x, y)
    X0 = x_min
    Y0 = y_min

    # Calculate residual function over the evaluation points
    Fr = zeros((num_x*num_y, ))
    for i in range(num_y):
        for j in range(num_x):
            index = i*num_x+j
            Fr[index] = f_residual(x[j], y[i])
    Fr = Fr.reshape(num_y, num_x)
    
    Ffe = eval_chebyshev_filter(Xe.ravel(), Ye.ravel(), C_filter, NCX_filter, NCY_filter)
    Ffe = Ffe.reshape(num_y, num_x)
    

    FY0 = zeros((num_x, ))
    for i in range(num_x):
        FY0[i] = f_residual(x[i], Y0)

    FX0 = zeros((num_y, ))
    for i in range(num_y):
        FX0[i] = f_residual(X0, y[i])

    # Plot residual function
    fig = plt.figure()
    ax0 = fig.add_subplot(231)
    ax1 = fig.add_subplot(232, projection='3d')
    ax4 = fig.add_subplot(233)
    ax2 = fig.add_subplot(234)
    ax3 = fig.add_subplot(235)
    ax5 = fig.add_subplot(236)

    cfx = ax0.contourf(X, Y, Fr)
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    plt.colorbar(cfx, ax=ax0)

    psf = ax1.plot_surface(X, Y, Ffe, cmap="jet")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    plt.colorbar(psf, ax=ax1)

    ax2.set_title(f"Y = {y_min:0.1E}")
    ax2.plot(x, FY0)
    ax2.set_xlabel("X")

    ax3.set_title(f"X = {x_min:0.1E}")
    ax3.plot(y, FX0)
    ax3.set_xlabel("Y")

    ax4.set_title("log10(|Ffe-Fr|)")
    cntf = ax4.contourf(X, Y, log10(np_abs(Ffe-Fr)), cmap="jet")
    # cntf = ax4.plot_surface(X, Y, Ffe, cmap="jet")
    # ax4.scatter(Xf.ravel(), Yf.ravel(), Ff.ravel())
    ax4.set_xlabel("X")
    ax4.set_ylabel("Y")
    plt.colorbar(cntf, ax=ax4)

    ax5.set_title("log10(|Ffe-Fr|/|Fr|)")
    cntf = ax5.contourf(X, Y, log10(np_abs(Ffe-Fr)/np_abs(Fr)))
    ax5.set_xlabel("X")
    ax5.set_ylabel("Y")
    plt.colorbar(cntf, ax=ax5)

    plt.show()


if __name__ == "__main__":
    fit_residual_region_14()