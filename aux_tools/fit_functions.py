
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


class FcnParameters:

    def __init__( self ) -> None:
        self.x0      = 0
        self.x1      = 10
        self.x2      = 20
        self.x3      = 30
        self.f0      = 0
        self.f1      = 2
        self.df0dx   = 0
        self.df1dx   = 0


def calculate_step_function_coeffs( fcn_params: FcnParameters  ) -> np.ndarray:
    # Get local copies
    x1      = fcn_params.x1
    x2      = fcn_params.x2
    f0      = fcn_params.f0
    f1      = fcn_params.f1
    df0dx   = fcn_params.df0dx
    df1dx   = fcn_params.df1dx

    # Define system matrix
    A       = np.array([
                            [ x1**3.0, x1**2.0, x1, 1.0 ],
                            [ x2**3.0, x2**2.0, x2, 1.0 ],
                            [ 3.0*x1**2.0, 2.0*x1, 1.0, 0.0 ],
                            [ 3.0*x2**2.0, 2.0*x2, 1.0, 0.0 ]
                        ])
    B       = np.array([ f0, f1, df0dx, df1dx ])

    # Calculate coefficients
    C       = np.linalg.solve( A, B )

    return C


def besselj0( x: float ) -> float:
    return sp.special.jv(0, x)


def exp_squared( x: float, y: float ) -> float:
    xc = 15.0
    yc = 15.0

    return np.exp( -( ( x - xc )**2.0 + ( y - yc )**2.0 ) )


def sin_wave( x: float ) -> float:
    w = 2.0*np.pi/0.5
    return np.sin( w*x )


def sin_wave_exp( x: float, y:float ) -> float:
    w = 2.0*np.pi/0.5
    return np.sin( w*x ) * np.exp( -np.abs(y) )


def step_function( fcn_params: FcnParameters, c: np.ndarray, x: float ) -> np.ndarray:
    f = 0.0
    if x < fcn_params.x1:
        f = fcn_params.f0

    elif x > fcn_params.x2:
        f = fcn_params.f1
    
    else:
        f = c[0]*x**3.0 + c[1]*x**2.0 + c[2]*x + c[3]
    
    return f


def step_function_1d( fcn_params: FcnParameters, x: float ) -> float:
    # Calculate smoothing coefficients
    coeffs      = calculate_step_function_coeffs( fcn_params )

    # Calculate step function
    step_fcn    = step_function( fcn_params, coeffs, x )

    return step_fcn


def step_function_2d( X, Y ) -> np.ndarray:
    a = 0


if __name__ == "__main__":
    # X, Y = np.meshgrid( x, x, indexing="ij" )

    # plt.plot(  )

    fcn_params  = FcnParameters( )
    x           = np.linspace( fcn_params.x0, fcn_params.x3, 1000 )
    y           = np.linspace( fcn_params.x0, fcn_params.x3, 1000 )
    X, Y        = np.meshgrid( x, y, indexing="ij" )
    F           = np.zeros( X.shape )

    for i in range( x.shape[0] ):
        for j in range( y.shape[0] ):
            # F[i, j] = step_function_1d( fcn_params, x[i] )
            # F[i, j] = exp_squared( x[i], y[j] )
            # F[i, j] = sin_wave( x[i], y[j] )
            # F[i, j] = sin_wave_exp( x[i], y[j] )
            F[i, j] = besselj0( x[i] )

    plt.contourf( X, Y, F )
    plt.colorbar( )
    plt.show( )