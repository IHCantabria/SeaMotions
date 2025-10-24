
# Import general usage libraries
import argparse
import h5py
import multiprocessing
import os
import time

# Import general usage scientific libraries
import numpy as np
import scipy as sp

# Import general usage scientific plotting libraries
import matplotlib.pyplot as plt


def besselj0( u: np.ndarray ) -> np.ndarray:
    # Add low u values
    f   =   (
                0.999999 
                - 
                2.24999239 * ( u / 3.0 )**2.0
                +
                1.26553572 * ( u / 3.0 )**4.0
                -
                0.31602189 * ( u / 3.0 )**6.0
                +
                0.04374224 * ( u / 3.0 )**8.0
                -
                0.00331563 * ( u / 3.0 )**10.0
            )

    # Add asymptotic values
    if type( u ) is np.ndarray:
        pos         = u > 3.0
        f0          =   ( 
                            0.79788454
                            -
                            0.00553897 * ( 3.0 / u[ pos ] )**2.0
                            +
                            0.00099336 * ( 3.0 / u[ pos ] )**4.0
                            -
                            0.00044346 * ( 3.0 / u[ pos ] )**6.0
                            +
                            0.00020445 * ( 3.0 / u[ pos ] )**8.0
                            -
                            0.00004959 * ( 3.0 / u[ pos ] )**10.0
                        )
        theta0      =   (
                            u[ pos ]
                            -
                            np.pi / 4.0
                            -
                            0.04166592 * ( 3.0 / u[ pos ] )
                            +
                            0.00239399 * ( 3.0 / u[ pos ] )**3.0
                            -
                            0.00073984 * ( 3.0 / u[ pos ] )**5.0
                            +
                            0.00031099 * ( 3.0 / u[ pos ] )**7.0
                            -
                            0.00007605 * ( 3.0 / u[ pos ] )**9.0
                        )
        f[ pos ]    =   f0 * np.cos( theta0 ) / np.sqrt( u[ pos ] )


        # f[ pos ]    =  (
                            # np.sqrt( 2.0 / np.pi ) * np.cos( np.pi/4.0 - u[ pos ] ) / np.sqrt( u[ pos ] )
                            # -
                            # np.sqrt( 2.0 / np.pi ) * np.sin( np.pi/4.0 - u[ pos ] ) / 8.0 / u[ pos ]**1.5
                        # )
    else:
        if u > 3.0:
            f0          =   ( 
                                0.79788454
                                -
                                0.00553897 * ( 3.0 / u )**2.0
                                +
                                0.00099336 * ( 3.0 / u )**4.0
                                -
                                0.00044346 * ( 3.0 / u )**6.0
                                +
                                0.00020445 * ( 3.0 / u )**8.0
                                -
                                0.00004959 * ( 3.0 / u )**10.0
                            )
            theta0      =   (
                                u
                                -
                                np.pi / 4.0
                                -
                                0.04166592 * ( 3.0 / u )
                                +
                                0.00239399 * ( 3.0 / u )**3.0
                                -
                                0.00073984 * ( 3.0 / u )**5.0
                                +
                                0.00031099 * ( 3.0 / u )**7.0
                                -
                                0.00007605 * ( 3.0 / u )**9.0
                            )
            f           =   f0 * np.cos( theta0 ) / np.sqrt( u )
            # f           =   (
                                # np.sqrt( 2.0 / np.pi ) * np.cos( np.pi/4.0 - u ) / np.sqrt( u )
                                # -
                                # np.sqrt( 2.0 / np.pi ) * np.sin( np.pi/4.0 - u ) / 8.0 / u**1.5
                            # )
            
    return f


def besselj0_fit( u: np.ndarray ) -> np.ndarray:
    # Add low u values
    f   =   (
                0.999999 
                - 
                2.24999239 * ( u / 3.0 )**2.0
                +
                1.26553572 * ( u / 3.0 )**4.0
                -
                0.31602189 * ( u / 3.0 )**6.0
                +
                0.04374224 * ( u / 3.0 )**8.0
                -
                0.00331563 * ( u / 3.0 )**10.0
            )

    # Add asymptotic values
    if type( u ) is np.ndarray:
        pos         = u > 3.0
        f0          =   ( 
                            0.79788454
                            -
                            0.00553897 * ( 3.0 / u[ pos ] )**2.0
                            +
                            0.00099336 * ( 3.0 / u[ pos ] )**4.0
                            -
                            0.00044346 * ( 3.0 / u[ pos ] )**6.0
                            +
                            0.00020445 * ( 3.0 / u[ pos ] )**8.0
                            -
                            0.00004959 * ( 3.0 / u[ pos ] )**10.0
                        )
        theta0      =   u[ pos ] - np.pi / 4.0
        theta1      =   (
                            -0.04166592 * ( 3.0 / u[ pos ] )
                            +
                            0.00239399 * ( 3.0 / u[ pos ] )**3.0
                            -
                            0.00073984 * ( 3.0 / u[ pos ] )**5.0
                            +
                            0.00031099 * ( 3.0 / u[ pos ] )**7.0
                            -
                            0.00007605 * ( 3.0 / u[ pos ] )**9.0
                        )
        f[ pos ]    =   f0 * ( np.cos( theta0 ) - np.sin( theta0 ) * theta1 ) / np.sqrt( u[ pos ] )

    else:
        if u > 3.0:
            f0          =   ( 
                                0.79788454
                                -
                                0.00553897 * ( 3.0 / u )**2.0
                                +
                                0.00099336 * ( 3.0 / u )**4.0
                                -
                                0.00044346 * ( 3.0 / u )**6.0
                                +
                                0.00020445 * ( 3.0 / u )**8.0
                                -
                                0.00004959 * ( 3.0 / u )**10.0
                            )
            theta0      =   u - np.pi / 4.0
            theta1      =   (
                                -0.04166592 * ( 3.0 / u )
                                +
                                0.00239399 * ( 3.0 / u )**3.0
                                -
                                0.00073984 * ( 3.0 / u )**5.0
                                +
                                0.00031099 * ( 3.0 / u )**7.0
                                -
                                0.00007605 * ( 3.0 / u )**9.0
                            )
            f           =   f0 * ( np.cos( theta0 ) - np.sin( theta0 ) * theta1 ) / np.sqrt( u )
            
    return f



def cvar_fcn_0( u: np.ndarray, b1: float, b2: float, A: float ) -> np.ndarray:
    return (
                1 / ( 1 + np.exp( -2.0 * u ) )
                *
                (
                    np.exp( -u * ( 2 - b1 ) )
                    +
                    np.exp( -u * ( 2 + b1 ) )
                    +
                    np.exp( -u * ( 2 - b2 ) )
                    +
                    np.exp( -u * ( 2 + b2 ) )
                )
                *
                sp.special.jv( 0, u * A )
            )


def cvar_fcn_1( u: np.ndarray, T: float, b1: float, b2: float, A: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = -u * ( 2 - b1 )
    arg_2 = -u * ( 2 + b1 )
    arg_3 = -u * ( 2 - b2 )
    arg_4 = -u * ( 2 + b2 )

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0
    arg_2_nozeros = arg_2 < 300.0
    arg_3_nozeros = arg_3 < 300.0
    arg_4_nozeros = arg_4 < 300.0

    total = 0.0
    if type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        term_2  = np.zeros( ( u.shape[0], ) )
        term_3  = np.zeros( ( u.shape[0], ) )
        term_4  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - b1 ) )
        term_2[ arg_2_nozeros ] = np.exp( -u[ arg_2_nozeros ] * ( 2 + b1 ) )
        term_3[ arg_3_nozeros ] = np.exp( -u[ arg_3_nozeros ] * ( 2 - b2 ) )
        term_4[ arg_4_nozeros ] = np.exp( -u[ arg_4_nozeros ] * ( 2 + b2 ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    ( 1 / ( 1 - term_0[ u_zero_pos ] ) )
                                    *
                                    np.sqrt( g * u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) ) )
                                    *
                                    (
                                        term_1[ u_zero_pos ]
                                        +
                                        term_2[ u_zero_pos ]
                                        +
                                        term_3[ u_zero_pos ]
                                        +
                                        term_4[ u_zero_pos ]
                                    )
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )
    else:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )

        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - b1 ) )
        
        term_2 = 0.0
        if arg_2_nozeros:
            term_2 = np.exp( -u * ( 2 + b1 ) )

        term_3 = 0.0
        if arg_3_nozeros:
            term_3 = np.exp( -u * ( 2 - b2 ) )

        term_4 = 0.0
        if arg_4_nozeros:
            term_4 = np.exp( -u * ( 2 + b2 ) )

        # Calculate total function value
        print( "---" )
        print( u, arg_0, term_0, arg_0_nozeros )
        print( u, arg_1, term_1, arg_1_nozeros )
        print( u, arg_2, term_2, arg_2_nozeros )
        print( u, arg_3, term_3, arg_3_nozeros )
        print( u, arg_4, term_4, arg_4_nozeros )
        total  = 0.0
        if u > 1e-16:
            total   = (
                            ( 1 / ( 1 - term_0 ) )
                            *
                            np.sqrt( g * u * tanhs( u ) )
                            *
                            np.sin( T * np.sqrt( u * tanhs( u) ) )
                            *
                            (
                                term_1
                                +
                                term_2
                                +
                                term_3
                                +
                                term_4
                            )
                            *
                            sp.special.jv( 0, u * A )
                        )
        print( u, ( 1 / ( 1 - term_0 ) ), np.sqrt( g * u * tanhs( u ) ), sp.special.jv( 0, u * A ), total )
        print( "---" )

    return total


def G_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * ( 2 - B )
    arg_2 = - u * ( 2 + B )

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0
    arg_2_nozeros = arg_2 < 300.0

    if type( u ) is float:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )

        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )
        
        term_2 = 0.0
        if arg_2_nozeros:
            term_2 = np.exp( -u * ( 2 + B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            ( 1 / ( 1 - term_0 ) )
                            *
                            np.sqrt( u * tanhs( u) )
                            *
                            np.sin( T * np.sqrt( u * tanhs( u) ) )
                            *
                            (
                                term_1
                                +
                                term_2
                            )
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        term_2  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )
        term_2[ arg_2_nozeros ] = np.exp( -u[ arg_2_nozeros ] * ( 2 + B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    ( 1 / ( 1 - term_0[ u_zero_pos ] ) )
                                    *
                                    np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) ) )
                                    *
                                    (
                                        term_1[ u_zero_pos ]
                                        +
                                        term_2[ u_zero_pos ]
                                    )
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_fcn1( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * ( 2 - B )
    arg_2 = - u * ( 2 + B )

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0
    arg_2_nozeros = arg_2 < 300.0

    if type( u ) is float:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )

        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )
        
        term_2 = 0.0
        if arg_2_nozeros:
            term_2 = np.exp( -u * ( 2 + B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            ( 1 / ( 1 - term_0 ) )
                            *
                            np.sqrt( u * tanhs( u) )
                            *
                            np.sin( T * np.sqrt( u * tanhs( u) ) )
                            *
                            (
                                term_1
                            )
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        term_2  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )
        term_2[ arg_2_nozeros ] = np.exp( -u[ arg_2_nozeros ] * ( 2 + B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    ( 1 / ( 1 - term_0[ u_zero_pos ] ) )
                                    *
                                    np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) ) )
                                    *
                                    (
                                        term_1[ u_zero_pos ]
                                    )
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_fcn_v6( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * ( 2 - B )
    arg_2 = - u * ( 2 + B )

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0
    arg_2_nozeros = arg_2 < 300.0

    if type( u ) is float:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )

        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )
        
        term_2 = 0.0
        if arg_2_nozeros:
            term_2 = np.exp( -u * ( 2 + B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            ( 1 / ( 1 - term_0 ) )
                            *
                            np.sqrt( u * tanhs( u) )
                            *
                            np.sin( T * np.sqrt( u * tanhs( u) ) )
                            *
                            (
                                term_1
                            )
                            *
                            sp.special.jv( 6, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        term_2  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )
        term_2[ arg_2_nozeros ] = np.exp( -u[ arg_2_nozeros ] * ( 2 + B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    ( 1 / ( 1 - term_0[ u_zero_pos ] ) )
                                    *
                                    np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] * tanhs( u[ u_zero_pos ] ) ) )
                                    *
                                    (
                                        term_1[ u_zero_pos ]
                                    )
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_fcn_clement( u: np.ndarray, mu: float, beta: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * mu

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = np.abs( arg_1 ) < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( arg_1 )
        
        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            np.sqrt( u )
                            *
                            np.sin( beta * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * np.sqrt( 1.0 - mu**2.0 ) )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( arg_1[ arg_1_nozeros ] )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    np.sqrt( u[ u_zero_pos ] )
                                    *
                                    np.sin( beta * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * np.sqrt( 1 - mu**2.0 ) )
                                )

    return total


def G_fcn_clement_dbeta( u: np.ndarray, mu: float, beta: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * mu

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = np.abs( arg_1 ) < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( arg_1 )
        
        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            u
                            *
                            np.cos( beta * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * np.sqrt( 1.0 - mu**2.0 ) )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( arg_1[ arg_1_nozeros ] )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = (
                                    u[ u_zero_pos ]
                                    *
                                    np.cos( beta * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * np.sqrt( 1 - mu**2.0 ) )
                                )

    return total


def G_fcn_clement_dbeta2( u: np.ndarray, mu: float, beta: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * mu

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = np.abs( arg_1 ) < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( arg_1 )
        
        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = -(
                            u**(3.0/2.0)
                            *
                            np.sin( beta * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * np.sqrt( 1.0 - mu**2.0 ) )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( arg_1[ arg_1_nozeros ] )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = -(
                                    u[ u_zero_pos ]**(3.0/2.0)
                                    *
                                    np.sin( beta * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * np.sqrt( 1 - mu**2.0 ) )
                                )

    return total


def G_fcn_clement_dbeta3( u: np.ndarray, mu: float, beta: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * mu

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = np.abs( arg_1 ) < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( arg_1 )
        
        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = -(
                            u**2.0
                            *
                            np.cos( beta * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * np.sqrt( 1.0 - mu**2.0 ) )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( arg_1[ arg_1_nozeros ] )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   = -(
                                    u[ u_zero_pos ]**2.0
                                    *
                                    np.cos( beta * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * np.sqrt( 1 - mu**2.0 ) )
                                )

    return total


def G_G_inf_diff( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    return G_fcn1( u, A, B, T, g ) - G_inf_asymp_fcn( u, A, B, T, g )


def G_inf_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )
        
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            ( 1 / ( 1 - term_0 ) )
                            *
                            np.sqrt( u )
                            *
                            np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   =  (
                                    ( 1.0 / ( 1 - term_0[ u_zero_pos ] ) )
                                    *
                                    np.sqrt( u[ u_zero_pos ] )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_inf_asymp_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_0 = 0.0
        if arg_0_nozeros:
            term_0 = np.exp( - 4 * u )
        
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate sqrt term
        sq_term = np.sqrt( u )
        if u < 5.0:
            # sq_term = -1.30752811e-01 * u**2.0 + 9.68179363e-01 * u -4.37969484e-16
            # sq_term =   (
            #                 -4.18837940e-05 * u**5.0 
            #                 -2.54911636e-03 * u**4.0  
            #                 +4.63769389e-02 * u**3.0 
            #                 -2.99080237e-01 * u**2.0 
            #                 +1.12798792e+00 * u
            #                 -2.61408726e-15
            #             )
            # sq_term = np.sqrt( u * u * 0.331 )
            sq_term = np.sqrt( u * np.tanh( u ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            # total   = (
            #                 np.sqrt( 2.0 / np.pi / A )
            #                 *
            #                 np.sin( T * np.sqrt( u ) )
            #                 *
            #                 term_1
            #                 *
            #                 np.cos( u * A - np.pi/4 )
            #             )
            total   = (
                            1 / ( 1 - term_0 )
                            *
                            sq_term
                            # np.sqrt( u )
                            *
                            np.sin( T * sq_term )
                            # np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                            *
                            # sp.special.jv( 0, u * A )
                            besselj0_fit( u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_0[ arg_0_nozeros ] = np.exp( - 4 * u[ arg_0_nozeros ] )
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Calculate sqrt term
        sq_term = np.sqrt( u )
        pos_3   = u < 5.0
        # sq_term[ pos_3 ] = sq_term = -1.30752811e-01 * u[ pos_3 ]**2.0 + 9.68179363e-01 * u[ pos_3 ] -4.37969484e-16
        # sq_term[ pos_3 ]    =   (
        #                             -4.18837940e-05 * u[ pos_3 ]**5.0 
        #                             -2.54911636e-03 * u[ pos_3 ]**4.0  
        #                             +4.63769389e-02 * u[ pos_3 ]**3.0 
        #                             -2.99080237e-01 * u[ pos_3 ]**2.0 
        #                             +1.12798792e+00 * u[ pos_3 ]
        #                             -2.61408726e-15
        #                         )
        # sq_term[ pos_3 ] = np.sqrt( u[ pos_3 ] * u[ pos_3] * 0.331 )
        sq_term[ pos_3 ] = np.sqrt( u[ pos_3 ] * np.tanh( u[ pos_3] ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        # total[u_zero_pos]   =  (
        #                             np.sqrt( 2.0 / np.pi / A )
        #                             *
        #                             np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
        #                             *
        #                             term_1[ u_zero_pos ]
        #                             *
        #                             np.cos( u * A - np.pi/4 )
        #                         )
        total[u_zero_pos]   =  (
                                    1 / ( 1 - term_0[ u_zero_pos ] )
                                    *
                                    sq_term[ u_zero_pos ]
                                    # np.sqrt( u[ u_zero_pos ] )
                                    *
                                    np.sin( T * sq_term[ u_zero_pos ] )
                                    # np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    # sp.special.jv( 0, u * A )
                                    besselj0_fit( u_zero_pos * A )
                                )

    return total


def G_inf_asymp2_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_0 = - 4 * u
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_0_nozeros = arg_0 < 300.0
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            # total   = (
            #                 np.sqrt( 2.0 / np.pi / A )
            #                 *
            #                 np.sin( T * np.sqrt( u ) )
            #                 *
            #                 term_1
            #                 *
            #                 np.cos( u * A - np.pi/4 )
            #             )
            total   = (
                            np.sqrt( u )
                            *
                            np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                            *
                            # sp.special.jv( 0, u * A )
                            besselj0( u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_0  = np.zeros( ( u.shape[0], ) )
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        # total[u_zero_pos]   =  (
        #                             np.sqrt( 2.0 / np.pi / A )
        #                             *
        #                             np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
        #                             *
        #                             term_1[ u_zero_pos ]
        #                             *
        #                             np.cos( u * A - np.pi/4 )
        #                         )
        total[u_zero_pos]   =  (
                                    np.sqrt( u_zero_pos )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    # sp.special.jv( 0, u * A )
                                    besselj0( u_zero_pos * A )
                                )

    return total


def G_inf_simple_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    B       = 2
    arg_1   = - u * ( 2 - B )

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = np.abs( arg_1 ) < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( arg_1 )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            np.sqrt( g * u )
                            *
                            np.sin( T * np.sqrt( u ) )
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * B )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   =  (
                                    np.sqrt( g * u )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_inf_star_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            np.sqrt( g * u )
                            *
                            np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   =  (
                                    np.sqrt( g * u[ u_zero_pos ] )
                                    *
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


def G_inf_star2_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   =  (
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                )

    return total


def G_inf_star3_fcn( u: np.ndarray, A: float, B: float, T: float, g: float ) -> np.ndarray:
    # Calculate exponential function arguments
    arg_1 = - u * B

    # Define non-zero ranges for the exponential terms
    arg_1_nozeros = arg_1 < 300.0

    if type( u ) is float:
        term_1 = 0.0
        if arg_1_nozeros:
            term_1 = np.exp( -u * ( 2 - B ) )

        # Calculate total function value
        total  = 0.0
        if u > 1e-16:
            total   = (
                            np.sin( T * np.sqrt( u ) )
                            *
                            term_1
                            *
                            sp.special.jv( 0, u * A )
                        )

    elif type( u ) is np.ndarray:
        # Allocate space for the terms
        term_1  = np.zeros( ( u.shape[0], ) )
        total   = np.zeros( ( u.shape[0], ) )

        # Calculate non-zero function value
        term_1[ arg_1_nozeros ] = np.exp( -u[ arg_1_nozeros ] * ( 2 - B ) )

        # Check for singularity due to u=0
        u_zero_pos          = u > 1e-16
        total[u_zero_pos]   =  (
                                    np.sin( T * np.sqrt( u[ u_zero_pos ] ) )
                                    *
                                    term_1[ u_zero_pos ]
                                    *
                                    sp.special.jv( 0, u[ u_zero_pos ] * A )
                                )

    return total


# def integrate_fcn( fcn, h: float ) -> float:
#     I1AA = sp.integrate.quad(
#                                 fcn,
#                                 0 * h,
#                                 0.1 * h                              
#                             )
#     I1AB = sp.integrate.quad(
#                                 fcn,
#                                 0.1 * h,
#                                 0.2 * h                              
#                             )
#     I1AC = sp.integrate.quad(
#                                 fcn,
#                                 0.2 * h,
#                                 0.3 * h                              
#                             )
#     I1AD = sp.integrate.quad(
#                                 fcn,
#                                 0.3 * h,
#                                 0.4 * h                              
#                             )
#     I1AE = sp.integrate.quad(
#                                 fcn,
#                                 0.4 * h,
#                                 0.5 * h                              
#                             )
#     I1AF = sp.integrate.quad(
#                                 fcn,
#                                 0.5 * h,
#                                 0.6 * h                              
#                             )
#     I1AG = sp.integrate.quad(
#                                 fcn,
#                                 0.6 * h,
#                                 0.7 * h                              
#                             )
#     I1AH = sp.integrate.quad(
#                                 fcn,
#                                 0.7 * h,
#                                 0.8 * h                              
#                             )
#     I1AI = sp.integrate.quad(
#                                 fcn,
#                                 0.8 * h,
#                                 0.9 * h                              
#                             )
#     I1AJ = sp.integrate.quad(
#                                 fcn,
#                                 0.9 * h,
#                                 1.0 * h                              
#                             )
    
#     I1B = sp.integrate.quad(
#                                 fcn,
#                                 1 * h,
#                                 2 * h                              
#                             )
#     I1C = sp.integrate.quad(
#                                 fcn,
#                                 2 * h,
#                                 3 * h                              
#                             )
#     I1D = sp.integrate.quad(
#                                 fcn,
#                                 3 * h,
#                                 4 * h                              
#                             )
#     I1E = sp.integrate.quad(
#                                 fcn,
#                                 4 * h,
#                                 5 * h                              
#                             )
    
#     return ( I1AA[0] + I1AB[0] + I1AC[0] + I1AD[0] + I1AE[0] + I1AF[0] + I1AG[0] + I1AH[0] + I1AI[0] + I1AJ[0] + I1B[0] + I1C[0] + I1D[0] + I1E[0] )


def integrate_adapt_fcn( fcn, a: float, b:float, eps=1e-6, ref_step=10 ) -> float:
    # Define first interval
    # lt  = ( b - a ) / 1000

    # if lt > ref_step:
        # lt = ref_step
    lt = ref_step

    # Loop over refinement levels
    is_conv = False
    for _ in range( 5 ):
        # Integrate function
        values = integrate_fcn( fcn, a, b, lt, eps=eps )

        # Check convergence
        if np.abs( values[1] ) < eps:
            is_conv = True
            break
    
        # Define next interval
        lt /= 10.0

    if not is_conv:
        raise ValueError( "Convergence not found!" )
    
    return values


def integrate_fcn( fcn, a: float, b:float, du: float, eps=1e-6 ) -> float:
    # Define intervals
    interval_np = int( np.ceil( ( b - a ) / du ) )
    intervals   = np.linspace( a, b, interval_np )
    integrals   = np.zeros_like( intervals )

    # Loop over integral to find the overall function integral value
    int_total   = 0.0
    total_err   = 0.0
    hist_np     = 100
    for i in range( intervals.shape[0] - 1 ):
        int_value   = sp.integrate.quad(
                                            fcn,
                                            intervals[i],
                                            intervals[i+1]
                                        )
        if np.abs( int_value[1] ) > eps:
            break
        
        integrals[i]    = int_value[0]
        int_total       += int_value[0]
        total_err       += np.abs( int_value[1] )

        # Check for convergence
        if i > hist_np:
            if np.abs( integrals[i-hist_np:i] ).max( ) < eps:
                break

    return int_total, total_err


def integrate_full_matrix( procs_np: int, input_fipath: str ) -> None:
    # Define integration matrix
    # A = np.array( [ 0.0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 100.0 ] )
    # B = np.array( [ 0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.95, 1.99 ] )
    # T = np.arange( 0, 100, 0.01 )

    # fcn = G_fcn1
    fcn = G_fcn_v6
    # fcn = G_G_inf_diff
    A = np.array( [ 10.0 ] )
    # B = np.array( [ 0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.95, 1.99 ] )
    B = np.array( [ 1.7, 1.8, 1.9, 1.95, 1.99, 1.9999 ] )
    T = np.arange( 0, 100, 0.1 )
    # fcn, A, B, T = read_input_file( input_fipath )

    # Define input arguments for the pool
    input_args  = [ ]
    for ai in A:
        for bi in B:
            input_args.append( ( fcn, ai, bi, T ) )

    # Define pool for computing
    jobs_np = len( A ) * len( B )
    data    = [ ]
    with multiprocessing.Pool( procs_np ) as p:
        for i, d in enumerate( p.imap_unordered( integrate_over_time, input_args ) ):
            t_exec = np.array( d[2] ).sum( )
            print( "Job: ", i, " - out of: ", jobs_np, " - A: ", input_args[i][1], " - B: ", input_args[i][2], " - Exec Time: ", t_exec, flush=True )
            data.append( d )

    # Allocate space to allocate data
    fcn_val     = np.zeros( ( A.shape[0], B.shape[0], T.shape[0] ) )
    fcn_err     = np.zeros( ( A.shape[0], B.shape[0], T.shape[0] ) )
    fcn_time    = np.zeros( ( A.shape[0], B.shape[0], T.shape[0] ) )

    # Unpack data into separate matrixes
    for i in range( A.shape[0] ):
        for j in range( B.shape[0] ):
            index               = i * B.shape[0] + j
            fcn_val[i, j, :]    = np.array( data[ index ][0] )
            fcn_err[i, j, :]    = np.array( data[ index ][1] )
            fcn_time[i, j, :]   = np.array( data[ index ][2] )

    # Storage data into disk
    fopath = os.path.dirname( input_fipath )
    fipath = os.path.join( fopath, "database.h5" )
    with h5py.File( fipath, "w" ) as fid:
        fid.create_dataset( "A", data=A )
        fid.create_dataset( "B", data=B )
        fid.create_dataset( "T", data=T )
        fid.create_dataset( "fcn_val", data=fcn_val )
        fid.create_dataset( "fcn_err", data=fcn_err )
        fid.create_dataset( "fcn_time", data=fcn_time )


def integrate_full_matrix_clement( procs_np: int, input_fipath: str ) -> None:
    # fcn = G_fcn1
    # mu = np.array( [ 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0 ] )
    mu      = np.linspace( 1e-3, 9e-3, 10 )
    beta    = np.arange( 0, 20.0, 0.01 )
    # fcn, A, B, T = read_input_file( input_fipath )

    # Storage data into disk
    fopath = os.path.dirname( input_fipath )
    fipath = os.path.join( fopath, "database.h5" )
    with h5py.File( fipath, "w" ) as fid:
        fid.create_dataset( "A", data=mu )
        fid.create_dataset( "B", data=beta )

    # Define input arguments for the pool
    fcn = G_fcn_clement
    input_args  = [ ]
    for i, ai in enumerate( mu ):
        input_args.append( ( fcn, ai, beta, i ) )

    with h5py.File( fipath, "a" ) as fid:
        integrate_full_matrix_clement_f( procs_np, mu, beta, input_args, fid, "f" )

    # Define input arguments for the pool
    fcn = G_fcn_clement_dbeta
    input_args  = [ ]
    for i, ai in enumerate( mu ):
        input_args.append( ( fcn, ai, beta, i ) )

    with h5py.File( fipath, "a" ) as fid:
        integrate_full_matrix_clement_f( procs_np, mu, beta, input_args, fid, "f_dbeta" )

    # Define input arguments for the pool
    fcn = G_fcn_clement_dbeta2
    input_args  = [ ]
    for i, ai in enumerate( mu ):
        input_args.append( ( fcn, ai, beta, i ) )

    with h5py.File( fipath, "a" ) as fid:
        integrate_full_matrix_clement_f( procs_np, mu, beta, input_args, fid, "f_dbeta2" )

    # Define input arguments for the pool
    fcn = G_fcn_clement_dbeta3
    input_args  = [ ]
    for i, ai in enumerate( mu ):
        input_args.append( ( fcn, ai, beta, i ) )

    with h5py.File( fipath, "a" ) as fid:
        integrate_full_matrix_clement_f( procs_np, mu, beta, input_args, fid, "f_dbeta3" )

    

def integrate_full_matrix_clement_f( procs_np: int, mu: np.ndarray, beta: np.ndarray, input_args: list, fid: h5py.File, fcn_name: str ) -> None:
    # Define pool for computing
    jobs_np = len( mu )
    data    = [ ]
    with multiprocessing.Pool( procs_np ) as p:
        for i, d in enumerate( p.imap_unordered( integrate_over_time_clement, input_args ) ):
            t_exec = np.array( d[2] ).sum( )
            print( "Job: ", i, " - out of: ", jobs_np, " - A: ", input_args[i][1], " - Exec Time: ", t_exec, flush=True )
            data.append( d )

    # Allocate space to allocate data
    fcn_val     = np.zeros( ( mu.shape[0], beta.shape[0] ) )
    fcn_err     = np.zeros( ( mu.shape[0], beta.shape[0] ) )
    fcn_time    = np.zeros( ( mu.shape[0], beta.shape[0] ) )

    # Unpack data into separate matrixes
    for index in range( mu.shape[0] ):
        i               = int( data[ index ][3] )
        fcn_val[i, :]   = np.array( data[ index ][0] )
        fcn_err[i, :]   = np.array( data[ index ][1] )
        fcn_time[i, :]  = np.array( data[ index ][2] )

    fid.create_dataset( f"{fcn_name:s}_val", data=fcn_val )
    fid.create_dataset( f"{fcn_name:s}_err", data=fcn_err )
    fid.create_dataset( f"{fcn_name:s}_time", data=fcn_time )


def integrate_over_time( *args ) -> list:
    # Unpack input arguments
    fcn = args[0][0]
    A   = args[0][1]
    B   = args[0][2]
    T   = args[0][3]

    # Loop over time to integrate function
    fcn_storage         = np.zeros( ( T.shape[0], ) )
    fcn_err_storage     = np.zeros( ( T.shape[0], ) )
    fcn_time_storage    = np.zeros( ( T.shape[0], ) )
    hist_np             = 100
    for i in range( T.shape[0] ):
        t0          = time.perf_counter( )
        G_int       = integrate_adapt_fcn( lambda ul: fcn( ul, A, B, T[i], 9.81 ), 0.0, 40.0, ref_step=10 )
        t1          = time.perf_counter( )

        # Storage results
        fcn_storage[i]      = G_int[0]
        fcn_err_storage[i]  = G_int[1]
        fcn_time_storage[i] = t1 - t0

        # Check if latest results are zero so the function is damped out
        if i > hist_np:
            if np.abs( fcn_storage[i-hist_np:i] ).max( ) < 1e-6:
                break

    return fcn_storage, fcn_err_storage, fcn_time_storage


def integrate_over_time_clement( *args ) -> list:
    # Unpack input arguments
    fcn     = args[0][0]
    mu      = args[0][1]
    beta    = args[0][2]
    pos     = args[0][3]

    # Loop over time to integrate function
    fcn_storage         = np.zeros( ( beta.shape[0], ) )
    fcn_err_storage     = np.zeros( ( beta.shape[0], ) )
    fcn_time_storage    = np.zeros( ( beta.shape[0], ) )
    hist_np             = 100
    for i in range( beta.shape[0] ):
        t0          = time.perf_counter( )
        G_int       = integrate_adapt_fcn( lambda ul: fcn( ul, mu, beta[i], 9.81 ), 0.0, 30e3, ref_step=2.0 )
        t1          = time.perf_counter( )

        # Storage results
        fcn_storage[i]      = G_int[0]
        fcn_err_storage[i]  = G_int[1]
        fcn_time_storage[i] = t1 - t0

        # Check if latest results are zero so the function is damped out
        if i > hist_np:
            if np.abs( fcn_storage[i-hist_np:i] ).max( ) < 1e-6:
                break

    return fcn_storage, fcn_err_storage, fcn_time_storage, pos


def raw_fcn_0( k: np.ndarray, h: float, ep: float, z: float, R: float ) -> np.ndarray:
    return ( 
                np.exp( - k * h ) / np.cosh( k * h )
                *
                np.cosh( k * ( ep + h ) ) * np.cosh( k * ( z + h ) )
                *
                sp.special.jv( 0, k * R )
            )


def raw_fcn_1( k: np.ndarray, h: float, ep: float, z: float, R: float, g: float, t: float ) -> np.ndarray:
    return ( 
                np.sqrt( g * k * np.tanh( k * h ) ) / ( np.cosh( k * h ) * np.sinh( k * h ) )
                *
                np.sin( t * np.sqrt( g * k * np.tanh( k * h ) ) )
                *
                np.cosh( k * ( ep + h ) ) * np.cosh( k * ( z + h ) )
                *
                sp.special.jv( 0, k * R )
            )


def main( ) -> None:
    # Configure case
    g   = 9.81
    R   = 250
    ep  = -1
    z   = -2
    t   = 10.0
    h   = 1000.0

    # Define wave numbers
    k   = np.linspace( 0, 10.0, int( 3e5 ) ) # kmax=200 for high oscillation cases

    # Define non dimensional wave numbers
    u   = k * h

    # Define non-dimensional variables
    A   = R / h
    T   = t * np.sqrt( g / h )
    b1  = np.abs( z - ep ) / h
    b2  = ( z + ep + 2 * h ) / h
    b3  = np.abs( z + ep ) / h

    # # Calculate raw function values
    # rf0 = 2.0 * raw_fcn_0( k, h, ep, z, R )
    # rf1 = raw_fcn_1( k, h, ep, z, R, g, t )

    # # Calculate equations with a variable change
    # cf0 = cvar_fcn_0( u, b1, b2, A )
    # cf1 = cvar_fcn_1( u, T, 0.0, 2.0, 0.0, g ) / np.sqrt( h )
    # cf1_rep = np.zeros_like( cf1 )


    Al = 20.0
    Bl = 1.99
    Tl = 0.1

    T           = [ ]
    G           = [ ]
    G_err       = [ ]
    G_inf       = [ ]
    G_inf_err   = [ ]
    t_init      = time.perf_counter( )
    for i in range( 100 ):
        ti          = (i+1) * 1.0
        t0          = time.perf_counter( )
        G_int       = integrate_adapt_fcn( lambda ul: G_fcn( ul, Al, Bl, ti, g ), 0.0, u[-1], ref_step=10 )
        G_inf_int   = [ 0.0, 0.0 ]
        # G_inf_int   = integrate_adapt_fcn( lambda ul: G_inf_star_fcn( ul, Al, Bl, ti, g ), 0.0, u[-1], ref_step=1 )
        t1          = time.perf_counter( )
        print( "i: ", i, t1 - t0, G_int, G_inf_int )

        T.append( ti )
        G.append( G_int[0] )
        G_err.append( G_int[1] )
        G_inf.append( G_inf_int[0] )
        G_inf_err.append( G_inf_int[1] )
    t_end       = time.perf_counter( )

    print( "Total time: ", t_end - t_init )

    T           = np.array( T )
    G           = np.array( G )
    G_err       = np.array( G_err )
    G_inf       = np.array( G_inf )
    G_inf_err   = np.array( G_inf_err )

    res_finame = f"results_A_{Al:0.6f}_B_{Bl:0.6f}.h5"
    res_fopath = r"E:\sergio\0050_OASIS_SM\TimeDomain"
    res_fipath = os.path.join( res_fopath, res_finame )
    with h5py.File( res_fipath, "w" ) as fid:
        fid.create_dataset( "T", data=T )
        fid.create_dataset( "G", data=G )
        fid.create_dataset( "G_err", data=G_err )
        fid.create_dataset( "G_inf", data=G_inf )
        fid.create_dataset( "G_inf_err", data=G_inf_err )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )

    ax0.plot( T, G )
    ax1.plot( T, G_err )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )

    ax0.plot( T, G_inf )
    ax1.plot( T, G_inf_err )

    plt.show( )

    # t0      = time.perf_counter( )
    # int_0   = integrate_fcn( lambda ul: G_inf_star_fcn( ul, Al, Bl, Tl, g ), 0.0, u[-1], 0.01*h )
    # t1      = time.perf_counter( )
    # int_1   = integrate_fcn( lambda ul: G_inf_star_fcn( ul, Al, Bl, Tl, g ), 0.0, u[-1], 0.001*h )
    # t2      = time.perf_counter( )
    # int_2   = integrate_fcn( lambda ul: G_inf_star_fcn( ul, Al, Bl, Tl, g ), 0.0, u[-1], 0.0001*h )
    # t3      = time.perf_counter( )

    # print( "Du 0.01*h:      ", int_0[0], " - Error: ", int_0[1], " - Exec Time: ", t1 - t0 )
    # print( "Du 0.001*h:     ", int_1[0], " - Error: ", int_1[1], " - Exec Time: ", t2 - t1 )
    # print( "Du 0.0001*h:    ", int_2[0], " - Error: ", int_2[1], " - Exec Time: ", t3 - t2 )

    f0 = G_fcn( u, Al, Bl, Tl, g )
    f1 = G_inf_star_fcn( u, Al, Bl, Tl, g )
    f2 = f1 - f0
    f3 = G_inf_star2_fcn( u, Al, Bl, Tl, g )
    f4 = f1 - G_inf_star2_fcn( u, Al, Bl, Tl, g )

    fig = plt.figure( )
    ax0 = fig.add_subplot( 511 )
    ax1 = fig.add_subplot( 512 )
    ax2 = fig.add_subplot( 513 )
    ax3 = fig.add_subplot( 514 )
    ax4 = fig.add_subplot( 515 )
    ax0.plot( u, f0, label="f0" )
    ax1.plot( u, f1, label="f1" )
    ax2.plot( u, f2, label="f2" )
    ax3.plot( u, f3, label="f3" )
    ax4.plot( u, f4, label="f4" )
    plt.show( )
    raise ValueError( "Stop by user" )

    # for i in range( cf1.shape[0] ):
    #     cf1_rep[i] = cvar_fcn_1( u[i], T, 0.0, 2.0, 0.0, g ) / np.sqrt( h )

    # plt.plot( u, cf1 )
    # plt.plot( u, cf1_rep )
    # plt.show( )
    # raise ValueError( "Stop by user" )

    # print( integrate_fcn( lambda u: cvar_fcn_1( u, T, 0.0, 2.0, 0.0, g ), h ) )

    t_mat       = [ ]
    f_b1_mat    = [ ]
    f_b2_mat    = [ ]
    f_inf_mat   = [ ]
    f_infs_mat  = [ ]
    f_infss_mat = [ ]
    A           = 1.0
    # A_vec       = [ 0.0, 0.25, 0.5, 0.75, 1.0 ]
    # B_vec       = [ 0.0, 0.25, 0.5, 0.75, 0.99, 1.0 ]
    B_vec       = [ 0.999 ]
    B           = 0.1
    for bi in B_vec:
        t_vec       = [ ]
        f_b1_vec    = [ ]
        f_b2_vec    = [ ]
        f_inf_vec   = [ ]
        f_infs_vec  = [ ]
        f_infss_vec = [ ]
        for i in range( 300 ):
            # A   = ai
            T   = i * np.sqrt( g / h )
            print( f"A: {A} - T: {T}" )
            
            # int_f       = integrate_fcn( lambda u: cvar_fcn_1( u, T, b1, b2, ai, g ), h )
            int_f_b1    = integrate_fcn( lambda u: G_fcn( u, A, bi, T, g ), h )
            int_f_b2    = integrate_fcn( lambda u: G_fcn( u, A, bi+1, T, g ), h )
            int_f_inf   = integrate_fcn( lambda u: G_inf_star_fcn( u, A, bi+1, T, g ), h )
            int_f_infs  = integrate_fcn( lambda u: G_inf_star2_fcn( u, A, bi+1, T, g ), h )
            int_f_infss = integrate_fcn( lambda u: G_inf_star3_fcn( u, A, bi+1, T, g ), h )
            t_vec.append( T )
            f_b1_vec.append( int_f_b1 )
            f_b2_vec.append( int_f_b2 )
            f_inf_vec.append( int_f_inf )
            f_infs_vec.append( int_f_infs )
            f_infss_vec.append( int_f_infss )

        t_mat.append( t_vec )
        f_b1_mat.append( f_b1_vec )
        f_b2_mat.append( f_b2_vec )
        f_inf_mat.append( f_inf_vec )
        f_infs_mat.append( f_infs_vec )
        f_infss_mat.append( f_infss_vec )

    fig = plt.figure( )
    fig.suptitle( "G( B1 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        ax.plot( t_mat[i], f_b1_mat[i], label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        ax.plot( t_mat[i], f_b2_mat[i], label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G( B1 ) + G( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        v0  = np.array( f_b1_mat[i] )
        v1  = np.array( f_b2_mat[i] )
        ax.plot( t_mat[i], v0 + v1, label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G_inf_star( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        ax.plot( t_mat[i], f_inf_mat[i], label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G_inf_star2( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        ax.plot( t_mat[i], f_infs_mat[i], label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G_inf_star3( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        ax.plot( t_mat[i], f_infss_mat[i], label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G ( B2 ) - G_inf_star( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        v0  = np.array( f_b2_mat[i] )
        v1  = np.array( f_inf_mat[i] )
        ax.plot( t_mat[i], v0 - v1, label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G_inf_star( B2 ) - G_inf_star2( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        v0  = np.array( f_inf_mat[i] )
        v1  = np.array( f_infs_mat[i] )
        ax.plot( t_mat[i], v0 - v1, label=f"A: {B_vec[i]}" )
        ax.legend( )

    fig = plt.figure( )
    fig.suptitle( "G_inf_star( B2 ) - G_inf_star3( B2 )" )
    for i in range( len( B_vec ) ):
        ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
        v0  = np.array( f_inf_mat[i] )
        v1  = np.array( f_infss_mat[i] )
        ax.plot( t_mat[i], v0 - v1, label=f"A: {B_vec[i]}" )
        ax.legend( )
    
    # fig = plt.figure( )
    # fig.suptitle( "G( B1 ) + G( B2 ) - G( B3 )" )
    # for i in range( len( B_vec ) ):
    #     ax  = fig.add_subplot( len( B_vec ), 1, i+1 )
    #     v0  = np.array( f_b1_mat[i] )
    #     v1  = np.array( f_b2_mat[i] )
    #     v2  = np.array( f_inf_mat[i] )
    #     ax.plot( t_mat[i], v0 + v1 - v2, label=f"A: {B_vec[i]}" )
    #     ax.legend( )
    
    plt.show( )

    # plt.plot( u, np.cosh( u ) * np.sinh( u ) )
    # plt.plot( u, np.sinh( 2.0 * u ) / 2.0 )
    # plt.show( )

    # Plot raw functions
    plot_raw_fcn( k, rf0, rf1 )

    # Plot comparison in between raw and variable composed
    plot_comparison_raw_comp( u, rf0, rf1, cf0, cf1 )

    plt.show( )


def plot_comparison_raw_comp( u: np.ndarray, rf0: np.ndarray, rf1: np.ndarray, cf0: np.ndarray, cf1: np.ndarray, show=False ) -> None:
    fig = plt.figure( )
    ax0 = fig.add_subplot( 211 )
    ax1 = fig.add_subplot( 212 )

    ax0.plot( u, rf0, label="Raw0" )
    ax0.plot( u, cf0, label="CVar0" )
    ax0.set_xlabel( "u [-]" )
    ax0.legend( )

    ax1.plot( u, rf1, label="Raw1" )
    ax1.plot( u, cf1, label="CVar1" )
    ax1.set_xlabel( "u [-]" )
    ax1.legend( )

    if show:
        plt.show( )


def plot_raw_fcn( k: np.ndarray, f0: np.ndarray, f1: np.ndarray, show=False ) -> None:
    plt.plot( k, f0, label="f0" )
    plt.plot( k, f1, label="f1" )
    plt.plot( k, f0 + f1, label="Total" )
    plt.xlabel( "Wave Number [-]" )
    plt.legend( )

    if show:
        plt.show( )


def read_input_file( fipath: str ) -> list:
    # Read data from disk
    with open( fipath, "r" ) as fid:
        lines = fid.readlines( )

    # Read target function
    fcn_spec = lines[0].strip( ).split( )
    if fcn_spec[0] != "FCN":
        raise ValueError( "First channel sould be the FCN spec" )
    fcn_spec = eval( fcn_spec[1] )

    # Read A values
    a_spec = read_input_spec( lines[1].strip( ) )
    if a_spec[0] != "A":
        raise ValueError( "Second channel should be the A spec." )
    
    # Read B values
    b_spec = read_input_spec( lines[2].strip( ) )
    if b_spec[0] != "B":
        raise ValueError( "Third channel should be the B spec." )
    
    # Read T values
    t_spec = read_input_spec( lines[3].strip( ) )
    if t_spec[0] != "T":
        raise ValueError( "Fourth channel should be the T spec." )
    
    return fcn_spec, a_spec[1], b_spec[1], t_spec[1]


def read_input_spec( spec: str ) -> list:
    # Read channel specification
    spec    = spec.strip( ).split( )
    name    = spec[0]
    spec[1] = int( spec[1] )

    # Check for channel specification
    if spec[1] == 0:
        if len( spec ) != 3:
            print( "Error in the specification of the input channels" )
            print( "The channel should have three components. The channel was specified like: " )
            print( spec )
            raise ValueError( "" )
    
    elif spec[1] == 1:
        if len( spec ) < 4:
            print( "Error in the specification of the input channels" )
            print( "The channel should have at least three components. The channel was specified like: " )
            print(  spec )
            raise ValueError( "" )
        
    # Expand compact list if any
    if spec[1] == 0:
        s   = spec[2].strip( )
        ss  = s.strip( "(" ).strip( ")" ).strip( "[" ).strip( "]" ).split( ":" )
        
        start   = float( ss[0] )
        dx      = float( ss[1] )
        end     = float( ss[2] )
        if s[0] == "(":
            start += dx

        if s[-1] == "]":
            end += dx
        
        values = np.arange( start, end, dx )

    else:
        values = np.array( spec[2:], dtype=float )
        
    return name, values


def tanhs( x: np.ndarray ) -> np.ndarray:
    # Compute function argument
    args    = 2 * x

    # Check positions where the function argument is lower than 300
    pos     = args < 300.0

    if type( x ) is np.ndarray:
        # Allocate space for the resulting function
        value   = np.ones( ( x.shape[0],  ) )

        # Calculate function value for each interval
        value[ pos ] = ( np.exp( 2 * x[ pos ] ) - 1.0 ) / ( np.exp( 2 * x[ pos ] ) + 1.0 )
    
    else:
        value = 1.0
        if pos:
            value = ( np.exp( 2 * x ) - 1.0 ) / ( np.exp( 2 * x ) + 1.0 )

    return value


if __name__ == "__main__":
    parser = argparse.ArgumentParser( "Calculate transient finite depth Green function coefficients" )
    parser.add_argument( "file_path", type=str )
    parser.add_argument( "-np", type=int, nargs="?", default=1 )

    args = parser.parse_args( )

    integrate_full_matrix_clement( args.np, args.file_path )

    # main( )

    # u = np.linspace( 0, 1000, int( 1e5 ) )
    # f = G_fcn_clement( u, 0.01, 10.0, 9.81 )
    # plt.plot( u, f )
    # plt.show( )
    # f = sp.special.jv( 0, u )
    # f1 = besselj0( u )
    # f2 = besselj0_fit( u )

    # fig = plt.figure( )
    # ax0 = fig.add_subplot( 211 )
    # ax1 = fig.add_subplot( 212 )

    # ax0.plot( u, f )
    # ax0.plot( u, f1 )
    # ax0.plot( u, f2 )

    # ax1.plot( u, f - f1 )
    # ax1.plot( u, f - f2 )

    # plt.show( )