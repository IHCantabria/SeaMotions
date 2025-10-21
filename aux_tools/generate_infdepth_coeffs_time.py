
# Import general usage libraries
import argparse
import h5py
import os
import time

# Import general usage scientific libraries
import numpy as np
import scipy as sp
from scipy.special import jv as besselj

# Import general usage scientific plotting libraries
import matplotlib.pyplot as plt

# Import local modules
from time_domain import integrate_adapt_fcn, read_input_spec


#######################################################################
################## Define F_beta time derivatives #####################
#######################################################################
def fas_lt( mu ) -> float:
    return 1.0

def fas( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fast( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return np.sqrt( lam ) * np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fastt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam * np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fasttt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(3/2) * np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


#######################################################################
################## Define F_beta time derivatives #####################
#######################################################################
def fac_lt( mu ) -> float:
    return 1.0

def fac( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fact( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return  - np.sqrt( lam ) * np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def factt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam * np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def facttt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(3/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


#######################################################################
################## Define F_beta time derivatives #####################
#######################################################################
def fb_lt( mu ) -> float:
    return 1.0

def fb( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return np.sqrt( lam ) * np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fbt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam * np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fbtt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(3/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


def fbttt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_b( lam, mu )


#######################################################################
############### Define F_beta_mu_0 time derivatives ###################
#######################################################################
def fbm0_lt( mu ) -> float:
    return -1.0


def fbm0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(3/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bm0( lam, mu )


def fbmt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bm0( lam, mu )


def fbmtt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bm0( lam, mu )


def fbmttt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bm0( lam, mu )


#######################################################################
############### Define F_beta_mu_1 time derivatives ###################
#######################################################################
def fbm1_lt( mu: float ) -> float:
    return mu / feta( mu )


def fbm1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(3/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bm1( lam, mu )


def fbmt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bm1( lam, mu )


def fbmtt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bm1( lam, mu )


def fbmttt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bm1( lam, mu )


#######################################################################
############## Define F_beta_mu_mu_0 time derivatives #################
#######################################################################
def fbmm0_lt( mu ) -> float:
    return 1.5


def fbmm0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm0( lam, mu )


def fbmmt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm0( lam, mu )


def fbmmtt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm0( lam, mu )


def fbmmttt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm0( lam, mu )


#######################################################################
############## Define F_beta_mu_mu_1 time derivatives #################
#######################################################################
def fbmm1_lt( mu ) -> float:
    return - 2.0 * mu / feta( mu )


def fbmm1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm1( lam, mu )


def fbmmt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm1( lam, mu )


def fbmmtt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm1( lam, mu )


def fbmmttt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm1( lam, mu )


#######################################################################
############## Define F_beta_mu_mu_2 time derivatives #################
#######################################################################
def fbmm2_lt( mu ) -> float:
    return 0.5 * ( mu**2.0 + 1.0 ) / feta( mu )**2.0


def fbmm2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm2( lam, mu )


def fbmmt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm2( lam, mu )


def fbmmtt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bmm2( lam, mu )


def fbmmttt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bmm2( lam, mu )


#######################################################################
############### Define F_beta_beta time derivatives ###################
#######################################################################
def fbb_lt( mu ) -> float:
    return 1.0

def fbb( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam * np.cos( np.sqrt( lam ) * beta ) * fkernel_bb( lam, mu )


def fbbt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(3/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bb( lam, mu )


def fbbtt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bb( lam, mu )


def fbbttt( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bb( lam, mu )


#######################################################################
############### Define F_beta_beta time derivatives ###################
#######################################################################
def fbbm0_lt( mu ) -> float:
    return -1.0


def fbbm0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbm0( lam, mu )


def fbbmt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbm0( lam, mu )


def fbbmtt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbm0( lam, mu )


def fbbmttt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbm0( lam, mu )


#######################################################################
############# Define F_beta_beta_mu time derivatives ##################
#######################################################################
def fbbm1_lt( mu ) -> float:
    return mu / feta( mu )


def fbbm1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**2.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbm1( lam, mu )


def fbbmt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(5/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbm1( lam, mu )


def fbbmtt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbm1( lam, mu )


def fbbmttt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbm1( lam, mu )


#######################################################################
############ Define F_beta_beta_mu_mu time derivatives ################
#######################################################################
def fbbmm0_lt( mu ) -> float:
    return 1.5


def fbbmm0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm0( lam, mu )


def fbbmmt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm0( lam, mu )


def fbbmmtt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm0( lam, mu )


def fbbmmttt0( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(9/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm0( lam, mu )


#######################################################################
############ Define F_beta_beta_mu_mu time derivatives ################
#######################################################################
def fbbmm1_lt( mu ) -> float:
    return -2.0 * mu / feta( mu )


def fbbmm1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm1( lam, mu )


def fbbmmt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm1( lam, mu )


def fbbmmtt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm1( lam, mu )


def fbbmmttt1( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(9/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm1( lam, mu )


#######################################################################
############ Define F_beta_beta_mu_mu time derivatives ################
#######################################################################
def fbbmm2_lt( mu ) -> float:
    return 0.5 * ( ( mu**2.0 + 1.0 ) / feta( mu )**2.0 )


def fbbmm2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**3.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm2( lam, mu )


def fbbmmt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**(7/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm2( lam, mu )


def fbbmmtt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return - lam**4.0 * np.cos( np.sqrt( lam ) * beta ) * fkernel_bbmm2( lam, mu )


def fbbmmttt2( lam: np.ndarray, mu: float, beta: float ) -> np.ndarray:
    return lam**(9/2) * np.sin( np.sqrt( lam ) * beta ) * fkernel_bbmm2( lam, mu )


def feta( mu: float ) -> float:
    return np.sqrt( 1 - mu**2.0 )


def fkernel_b( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bb( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bbm0( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bbm1( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 1, lam * feta( mu ) )


def fkernel_bbmm0( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bbmm1( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 1, lam * feta( mu ) )


def fkernel_bbmm2( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 2, lam * feta( mu ) )


def fkernel_bm0( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bm1( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 1, lam * feta( mu ) )


def fkernel_bmm0( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 0, lam * feta( mu ) )


def fkernel_bmm1( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 1, lam * feta( mu ) )


def fkernel_bmm2( lam: np.ndarray, mu: float ) -> np.ndarray:
    return np.exp( -lam * mu ) * besselj( 2, lam * feta( mu ) )


def calculate_Gas( fopath: str ) -> None:
    print( "Calculating Gas" )
    # Calculate Fb
    f, ft, ftt, fttt    = integrate_over_time( 
                                                calculate_initial_conditions_integration_fas,
                                                fas_lt,
                                                0,
                                                0
                                            )
    
    # Get calculation matrix to pass storage system
    beta, mu = get_calculation_matrix( )
    
    # Save data
    finame = "Gas.h5"
    fipath = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gas: DONE" )


def calculate_Gac( fopath: str ) -> None:
    print( "Calculating Gac" )
    # Calculate Fb
    f, ft, ftt, fttt    = integrate_over_time( 
                                                calculate_initial_conditions_integration_fac,
                                                fac_lt,
                                                0,
                                                0
                                            )
    
    # Get calculation matrix to pass storage system
    beta, mu = get_calculation_matrix( )
    
    # Save data
    finame = "Gac.h5"
    fipath = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gac: DONE" )


def calculate_Gt( fopath: str ) -> None:
    print( "Calculating Gt" )
    # Calculate Fb
    nu0                     = 0
    l0                      = 1/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fb,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f, ft, ftt, fttt        = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, True ),
                                                        fb_lt,
                                                        l0,
                                                        nu0
                                                    )
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gt.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gt: DONE" )


def calculate_Gtx( fopath: str ) -> None:
    print( "Calculating Gtx" )
    # Calculate Fbm0
    nu0                     = 0
    l0                      = 3/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbm0,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f0, ft0, ftt0, fttt0    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, True ),
                                                        fbm0_lt,
                                                        l0,
                                                        nu0
                                                    )
    
    # Calculate Fbm1
    nu1                     = 1
    l1                      = 3/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbm1,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f1, ft1, ftt1, fttt1    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, lt, True ),
                                                        fbm1_lt,
                                                        l1,
                                                        nu1
                                                    )
    
    # Add contributions
    f           = f0 + f1
    ft          = ft0 + ft1
    ftt         = ftt0 + ftt1
    fttt        = fttt0 + fttt1
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gtx.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gtx: DONE" )


def calculate_Gtxx( fopath: str ) -> None:
    print( "Calculating Gtxx" )
    # Calculate Fbmm0
    print( " -> Calculate Fbmm0" )
    nu0                     = 0
    l0                      = 5/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbmm0,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f0, ft0, ftt0, fttt0    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, True ),
                                                        fbmm0_lt,
                                                        l0,
                                                        nu0
                                                    )
    print( " -> Calculate Fbmm0 -> Done!" )

    # Calculate Fbmm1
    print( " -> Calculate Fbmm1" )
    nu1                     = 1
    l1                      = 5/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbmm1,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f1, ft1, ftt1, fttt1    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, lt, True ),
                                                        fbmm1_lt,
                                                        l1,
                                                        nu1
                                                    )
    print( " -> Calculate Fbmm1 -> Done!" )
    
    # Calculate Fbmm2
    print( " -> Calculate Fbmm2" )
    nu2                     = 2
    l2                      = 5/2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbmm2,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu2, l2, np.ones( ( 4,  ) ), True ),
                                                        True
                                                    )
    f2, ft2, ftt2, fttt2    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu2, l2, lt, True ),
                                                        fbmm2_lt,
                                                        l2,
                                                        nu2
                                                    )
    print( " -> Calculate Fbmm2 -> Done!" )
    
    # Add contributions
    f           = f0 + f1 + f2
    ft          = ft0 + ft1 + ft2
    ftt         = ftt0 + ftt1 + ftt2
    fttt        = fttt0 + fttt1 + fttt2
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gtxx.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gtxx: DONE" )


def calculate_Gtt( fopath: str ) -> None:
    print( "Calculating Gtt" )
    # Calculate Fb
    nu0                 = 0
    l0                  = 1
    lt                  = calculate_leading_term( 
                                                    calculate_initial_conditions_integration_fbb,
                                                    lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), False ),
                                                    False
                                                )
    f, ft, ftt, fttt    = integrate_over_time( 
                                                lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, False ),
                                                fbb_lt,
                                                l0,
                                                nu0
                                            )
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gtt.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gtt: DONE" )


def calculate_Gttx( fopath: str ) -> None:
    print( "Calculating Gttx" )
    # Calculate Fb
    nu0                     = 0
    l0                      = 2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbbm0,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), False ),
                                                        False
                                                    )
    f0, ft0, ftt0, fttt0    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, False ),
                                                        fbbm0_lt,
                                                        l0,
                                                        nu0
                                                    )
    
    # Calculate Fb
    nu1                     = 1
    l1                      = 2
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbbm1,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, np.ones( ( 4,  ) ), False ),
                                                        False
                                                    )
    f1, ft1, ftt1, fttt1    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, lt, False ),
                                                        fbbm1_lt,
                                                        l1,
                                                        nu1
                                                    )
    
    # Sum up contributions
    f           = f0 + f1
    ft          = ft0 + ft1
    ftt         = ftt0 + ftt1
    fttt        = fttt0 + fttt1
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gttx.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gttx: DONE" )


def calculate_Gttxx( fopath: str ) -> None:
    print( "Calculating Gttxx" )
    # Calculate Fb
    nu0                     = 0
    l0                      = 3
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbbmm0,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, np.ones( ( 4,  ) ), False ),
                                                        False
                                                    )
    f0, ft0, ftt0, fttt0    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu0, l0, lt, False ),
                                                        fbbmm0_lt,
                                                        l0,
                                                        nu0
                                                    )
    
    # Calculate Fb
    nu1                     = 1
    l1                      = 3
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbbmm1,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, np.ones( ( 4,  ) ), False ),
                                                        False
                                                    )
    f1, ft1, ftt1, fttt1    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu1, l1, lt, False ),
                                                        fbbmm1_lt,
                                                        l1,
                                                        nu1
                                                    )
    
    # Calculate Fb
    nu2                     = 2
    l2                      = 3
    lt                      = calculate_leading_term( 
                                                        calculate_initial_conditions_integration_fbbmm2,
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu2, l2, np.ones( ( 4,  ) ), False ),
                                                        False
                                                    )
    f2, ft2, ftt2, fttt2    = integrate_over_time( 
                                                        lambda mu: calculate_initial_conditions_integration_gen( mu, nu2, l2, lt, False ),
                                                        fbbmm2_lt,
                                                        l2,
                                                        nu2
                                                    )
    
    # Sum up contributions
    f           = f0 + f1 + f2
    ft          = ft0 + ft1 + ft2
    ftt         = ftt0 + ftt1 + ftt2
    fttt        = fttt0 + fttt1 + fttt2
    
    # Get calculation matrix to pass storage system
    beta, mu    = get_calculation_matrix( )
    
    # Save data
    finame      = "Gttxx.h5"
    fipath      = os.path.join( fopath, finame )
    save_data( fipath, beta, mu, f, ft, ftt, fttt )
    print( " -> Calculating Gttxx: DONE" )


def calculate_initial_conditions_integration_gen( mu, nu: int, l: float, lt: float, first_odd: bool ) -> list:
    if first_odd:   # Time dependence sin(beta) [ Gt, Gttt ]
        A       = 0.0
        Adt1    = lt[1] * init_conds_fcn( mu, nu, l, 0 )
        Adt2    = 0.0
        Adt3    = lt[3] * init_conds_fcn( mu, nu, l, 1 )
    
    else:           # Time dependence cos( beta ) [ Gtt, Gtttt ]
        A       = lt[0] * init_conds_fcn( mu, nu, l, 0 )
        Adt1    = 0.0
        Adt2    = lt[2] * init_conds_fcn( mu, nu, l, 1 )
        Adt3    = 0.0

    return A, Adt1, Adt2, Adt3


def calculate_initial_conditions_integration_fas( mu, is_ode=True ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 90e3 ] )

    # Define kernel functions
    kernel      = lambda ul: fas( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fast( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fastt( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fasttt( ul, mu, 0.0 )

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )

    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fac( mu, is_ode=True ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 90e3 ] )

    # Define kernel functions
    kernel      = lambda ul: fac( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fact( ul, mu, 0.0 )
    kerneldt2   = lambda ul: factt( ul, mu, 0.0 )
    kerneldt3   = lambda ul: facttt( ul, mu, 0.0 )

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )

    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fb( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fb( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbt( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbtt( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbttt( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fb" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )

    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbm0( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbm0( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbmt0( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbmtt0( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbmttt0( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbm0" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    print( "A0:     ", A0 )
    print( "A0dt1:  ", A0dt1 )
    print( "A0dt2:  ", A0dt2 )
    print( "A0dt3:  ", A0dt3 )

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbm1( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbm1( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbmt1( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbmtt1( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbmttt1( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbm1" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbmm0( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbmm0( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbmmt0( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbmmtt0( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbmmttt0( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbmm0" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbmm1( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbmm1( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbmmt1( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbmmtt1( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbmmttt1( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbmm1" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbmm2( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbmm2( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbmmt2( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbmmtt2( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbmmttt2( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbmm2" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbb( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbb( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbt( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbtt( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbttt( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbb" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbbm0( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbbm0( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbmt0( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbmtt0( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbmttt0( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbbm0" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbbm1( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbbm1( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbmt1( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbmtt1( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbmttt1( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbbm1" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbbmm0( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbbmm0( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbmmt0( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbmmtt0( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbmmttt0( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbbmm0" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbbmm1( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbbmm1( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbmmt1( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbmmtt1( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbmmttt1( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbbmm1" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_initial_conditions_integration_fbbmm2( mu, is_ode=True, is_show_fcn=False ) -> list:
    # Define time bounds
    timeb       = np.array( [ 0.0, 1e4 ] )

    # Define kernel functions
    kernel      = lambda ul: fbbmm2( ul, mu, 0.0 )
    kerneldt1   = lambda ul: fbbmmt2( ul, mu, 0.0 )
    kerneldt2   = lambda ul: fbbmmtt2( ul, mu, 0.0 )
    kerneldt3   = lambda ul: fbbmmttt2( ul, mu, 0.0 )

    # Show function if any
    if is_show_fcn:
        fig, ax = plt.subplots( 4, 1, sharex=True )

        fig.suptitle( "Target function: Fbbmm2" )
        
        time = np.linspace( 0, 1e7, int( 1e7 ) )

        ax[0].plot( time, kernel( time ) )
        ax[0].set_xlabel( "lambda" )
        ax[0].set_ylabel( "kernel" )

        ax[1].plot( time, kerneldt1( time ) )
        ax[1].set_xlabel( "lambda" )
        ax[1].set_ylabel( "kerneldt1" )

        ax[2].plot( time, kerneldt2( time ) )
        ax[2].set_xlabel( "lambda" )
        ax[2].set_ylabel( "kerneldt2" )

        ax[3].plot( time, kerneldt3( time ) )
        ax[3].set_xlabel( "lambda" )
        ax[3].set_ylabel( "kerneldt3" )

        plt.show( )
        return

    # Calculate initial values
    if is_ode:
        A0      = integrate_like_ode( kernel, timeb[0], timeb[1] )
        A0dt1   = integrate_like_ode( kerneldt1, timeb[0], timeb[1] )
        A0dt2   = integrate_like_ode( kerneldt2, timeb[0], timeb[1] )
        A0dt3   = integrate_like_ode( kerneldt3, timeb[0], timeb[1] )
		
    else:
        A0      = integrate_adapt_fcn( kernel, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt1   = integrate_adapt_fcn( kerneldt1, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt2   = integrate_adapt_fcn( kerneldt2, timeb[0], timeb[1], ref_step=10 )[0]
        A0dt3   = integrate_adapt_fcn( kerneldt3, timeb[0], timeb[1], ref_step=10 )[0]

    return A0, A0dt1, A0dt2, A0dt3


def calculate_leading_term( init_conds_ode, init_conds_leg, first_odd: bool ) -> float:
    # Define the point to calculate the initial conditions
    mu = 0.5

    # Calculate initial conditions by using ODE
    a0_ode = init_conds_ode( mu )

    # Calculate initial conditions by using legendre functions
    a0_leg = init_conds_leg( mu )

    # Calculate leading term ration
    lt = np.zeros( ( 4, ) )

    if first_odd:   # Time dependence sin(beta) [ Gt, Gttt ]
        lt[1] = a0_ode[1] / a0_leg[1]
        lt[3] = a0_ode[3] / a0_leg[3]
    
    else:           # Time dependence cos( beta ) [ Gtt, Gtttt ]
        lt[0] = a0_ode[0] / a0_leg[0]
        lt[2] = a0_ode[2] / a0_leg[2]

    for i in range( len( a0_ode ) ):
        print( f"A{i}: ", a0_ode[i], a0_leg[i], lt[i] )

    return lt
    


def clement_fcn( t, y, l, nu, mu ) -> np.ndarray:
    # Define ODE coefficients
    a0  = ( l + 1 )**2.0 - nu**2.0
    a1  = ( l + 5.0 / 4.0 ) * t
    a2  = ( t**2.0 / 4.0 + mu * ( 3.0 + 2.0 * l ) )
    a3  = mu * t
    a4  = 1.0

    # Define solution derivative vector
    y_new       = np.zeros_like( y )
    y_new[0]    = y[1]
    y_new[1]    = y[2]
    y_new[2]    = y[3]
    y_new[3]    = - ( 
                        a0 * y[0]
                        +
                        a1 * y[1]
                        +
                        a2 * y[2]
                        +
                        a3 * y[3]
                        ) / a4
    
    return y_new


def G0_Magee( beta: np.ndarray, mu ) -> np.ndarray:
    lt = np.pi * beta**3.0 / 8 / np.sqrt( 2 ) * np.exp( -beta**2.0 * mu / 4.0 )
    return lt * (
                    sp.special.jv( 0.25, beta**2.0 / 8.0 )
                    *
                    sp.special.jv( -0.25, beta**2.0 / 8.0 )
                    +
                    sp.special.jv( 0.75, beta**2.0 /8.0 )
                    *
                    sp.special.jv( -0.75, beta**2.0 / 8.0 )
                )


def get_calculation_matrix( ) -> list:
    t_start     = 0.0
    t_end       = 60.0
    dt          = 0.01
    beta        = np.arange( t_start, t_end, dt )
    mu          = 10**np.linspace( -4, np.log10( 0.9999 ), 100 )
    # mu          = 10**np.array( [ -4.0 ] )

    return beta, mu


def init_conds_fcn( mu: float, nu: int, l: float, k: int ) -> float:
    deg = int( l + k + 1/2 )
    fv  = ( 
                ( -1 ) ** ( k + nu ) 
                * 
                sp.special.gamma( l + ( 2.0 * k + 3.0 ) / 2.0 - nu )
                *
                sp.special.lpmn( nu, deg, mu )[0][nu, deg]
            )
    
    return fv


def integrate_like_ode( kernel_fcn, a, b ) -> None:
    def kernel_ode( t, y ):
        return kernel_fcn( t )
    
    print( a, b, kernel_fcn( a ) )
    result = sp.integrate.solve_ivp( 
                                        kernel_ode, 
                                        [ a, b ],
                                        [ kernel_fcn( a ) ],
                                        atol=1e-8, 
                                        rtol=1e-8,
                                        max_step=10.0
                                    )
    
    return result.y[0, -1] - kernel_fcn( a )


def integrate_over_time( init_cond_fcn, lt_fcn, l: float, nu: float ) -> list:
    # Define calculation matrix
    beta, mu    = get_calculation_matrix( )

    # Loop over mu domain
    f_val       = np.zeros( ( beta.shape[0], mu.shape[0] ) )
    f_beta1_val = np.zeros( ( beta.shape[0], mu.shape[0] ) )
    f_beta2_val = np.zeros( ( beta.shape[0], mu.shape[0] ) )
    f_beta3_val = np.zeros( ( beta.shape[0], mu.shape[0] ) )
    for i in range( mu.shape[0] ):
        print( " --> Integrating mu: ", mu[i] )
        # Calculate initial concitions
        init_conds  = init_cond_fcn( mu[i] )

        # Integrate over time ODE
        result      = sp.integrate.solve_ivp( 
                                                lambda t, y: clement_fcn( t, y, l, nu, mu[i] ), 
                                                                [beta[0], beta[-1]], 
                                                                init_conds,
                                                                t_eval=beta,
                                                                atol=1e-8, 
                                                                rtol=1e-8,
                                                                max_step=0.001
                                            )
        
        # print( init_conds )
        # plt.plot( result.t, result.y[0, :], label="IVP" )
        # plt.plot( result.t, G0_Magee( result.t, mu[i] ) / 2.0, label="Magee" )
        # plt.xlabel( "Time [s]" )
        # plt.ylabel( "Function Value" )
        # plt.legend( )
        # plt.show( )
        
        # Storage data
        f_val[:, i]         = lt_fcn( mu[i] ) * result.y[0, :]
        f_beta1_val[:, i]   = lt_fcn( mu[i] ) * result.y[1, :]
        f_beta2_val[:, i]   = lt_fcn( mu[i] ) * result.y[2, :]
        f_beta3_val[:, i]   = lt_fcn( mu[i] ) * result.y[3, :]
    
    return f_val, f_beta1_val, f_beta2_val, f_beta3_val


def save_data( fipath: str, beta: np.ndarray, mu: np.ndarray, f: np.ndarray, ft: np.ndarray, ftt: np.ndarray, fttt: np.ndarray ) -> None:
    with h5py.File( fipath, "w" ) as fid:
        fid.create_dataset( "mu",       data=mu     )
        fid.create_dataset( "beta",     data=beta   )
        fid.create_dataset( "fcn",      data=f      )
        fid.create_dataset( "fcn_b1",   data=ft     )
        fid.create_dataset( "fcn_b2",   data=ftt    )
        fid.create_dataset( "fcn_b3",   data=fttt   )


if __name__ == "__main__":
    folder_path = r"D:\sergio\0050_OASIS_SM\TimeDomain"
    # calculate_Gas( folder_path )
    # calculate_Gac( folder_path )
    calculate_Gt( folder_path )
    calculate_Gtx( folder_path )
    calculate_Gtxx( folder_path )
    calculate_Gtt( folder_path )
    calculate_Gttx( folder_path )
    calculate_Gttxx( folder_path )



    # mu  = np.array( [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5 ] )
    # P   = np.zeros( ( mu.shape[0], 4 ) )
    # P_dt   = np.zeros( mu.shape )
    # for i in range( mu.shape[0] ):
    #     print( "Calculating mu: ", i )
    #     P[i, 0], P[i, 1], P[i, 2], P[i, 3] = calculate_initial_conditions_integration_fbmm0( mu[i] )

    # plt.plot( np.log10( mu ), P[:, 0], label="A" )
    # plt.plot( np.log10( mu ), P[:, 1], label="A_dt" )
    # plt.plot( np.log10( mu ), P[:, 2], label="A_dtt" )
    # plt.plot( np.log10( mu ), P[:, 3], label="A_dttt" )
    # plt.legend( )
    # plt.show( )



    # with h5py.File( r"D:\sergio\0050_OASIS_SM\TimeDomain\Gas.h5", "r" ) as fid:
    #     beta   = fid[ "beta" ][:]
    #     mu     = fid[ "mu" ][:]
    #     fcns   = fid[ "fcn" ][:]

    # with h5py.File( r"D:\sergio\0050_OASIS_SM\TimeDomain\Gt.h5", "r" ) as fid:
    #     beta   = fid[ "beta" ][:]
    #     mu     = fid[ "mu" ][:]
    #     fcnc   = fid[ "fcn" ][:]

    # B, M = np.meshgrid( beta, mu, indexing="ij" )
    
    # fig, ax = plt.subplots( 2, 2 )

    # cnf = ax[0, 0].contourf( B, np.log10( M ), fcns )
    # plt.colorbar( cnf, ax=ax[0, 0] )

    # cnf = ax[0, 1].contourf( B, np.log10( M ), fcnc )
    # plt.colorbar( cnf, ax=ax[0, 1] )

    # cnf = ax[1, 0].contourf( B, np.log10( M ), fcns - fcnc )
    # plt.colorbar( cnf, ax=ax[1, 0] )

    # ax[1, 1].plot( beta, fcns[:, 0] )
    # ax[1, 1].plot( beta, fcnc[:, 0] )

    # plt.show( )




    # idx = 85
    # with h5py.File( r"D:\sergio\0050_OASIS_SM\TimeDomain\Gas.h5", "r" ) as fid:
    #     betas   = fid[ "beta" ][:]
    #     ss      = fid[ "fcn" ][:, idx]

    # with h5py.File( r"D:\sergio\0050_OASIS_SM\TimeDomain\Gac.h5", "r" ) as fid:
    #     betac   = fid[ "beta" ][:]
    #     sc      = fid[ "fcn" ][:, idx]

    # beta1   = np.arange( 0, 20.0, 0.1 )
    # sols     = np.zeros( beta1.shape[0] )
    # solc     = np.zeros( beta1.shape[0] )
    # for i in range( beta1.shape[0] ):
    #     print( i )
    #     si = integrate_adapt_fcn( lambda ul: fas( ul, 0.521356061745598, beta1[i] ), 0.0, 90e3, ref_step=10 )
    #     sols[i] = si[0]

    #     si = integrate_adapt_fcn( lambda ul: fac( ul, 0.521356061745598, beta1[i] ), 0.0, 90e3, ref_step=10 )
    #     solc[i] = si[0]

    # plt.plot( beta1, sols )
    # plt.plot( beta1, solc )
    # plt.plot( betas, ss )
    # plt.plot( betac, sc )
    # plt.show( )