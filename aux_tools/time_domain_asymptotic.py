
import numpy as np
import scipy as sp


def get_alpha( n ):
    return sp.special.factorial( 2*n + 2 ) / sp.special.factorial( n )


def get_gamma_half( n ):
    value = 1.0
    for i in range( 2, 2*n ):
        value *= 2*i - 1

    return value * np.sqrt( np.pi ) / 2**n


def get_cn( n ):
    return get_gamma_half( n )**2.0 / np.pi / 2**n / sp.special.factorial( n )


def get_dmn( m, n ):
    if m == 0 and n == 0:
        return 1.0
    elif n == 0:
        return 0.0
    else:
        return ( 
                    sp.special.factorial( 2*m + 2*n - 2 ) 
                    / 
                    ( 
                        sp.special.factorial( 2*n - 2 ) 
                        * 
                        2**( 2*m )
                        *
                        sp.special.factorial( m )
                    ) 
                )


def dGdt_asymptotic( beta, mu ):
    """
    Compute the asymptotic approximation of the Green function time derivative dG/dt for large time differences.

    Parameters
    ----------
    beta : float
        The non-dimensional time difference (t - tau) * g / U^2.
    mu : float
        The non-dimensional frequency parameter omega * U / g.
    Returns
    -------
    dGdt : float
        The asymptotic approximation of dG/dt.
    """
    nmax = 15
    
    eta = np.sqrt( 1 - mu**2.0 )

    # Calculate f0
    f0 = 0.0
    for n in range( nmax ):
        f0 += get_alpha( n ) * beta**( -2*n - 3 ) * sp.special.eval_legendre( n, mu )
    f0 *= -4.0

    # Calculate f1
    f2a = ( 
                ( - 1j ) * 4.0
                *
                np.exp( -1.0 / 4.0 * beta**2.0 * mu )
                *
                np.exp( 
                            1j
                            *
                            ( 
                                0.25 * beta**2.0 * eta
                                -
                                0.5 * np.arccos( mu )
                                +
                                np.pi / 4.0
                            )
                        )
                /
                np.sqrt( 2.0 * eta )
            )
    
    f2b = 0.0
    for n in range( nmax ):
        r0 = ( 1j / eta )**n
        r1 = 0.0
        for m in range( n + 1 ):
            r1  += ( 
                        get_dmn( m, n )
                        *
                        ( 0.5 * beta * ( 1j * mu + eta ) )**( 1.0 - 2.0 * m - 2.0 * n )
                        *
                        ( mu - 1j * eta )**m
                    )
        f2b += r0 * r1

    return (f0 + f2a * f2b).real
    # return (f2a * f2b).real


def dGdtt_asymptotic( beta, mu ):
    """
    Compute the asymptotic approximation of d²G/dt², i.e. the derivative of
    dGdt_asymptotic with respect to beta, obtained by exact analytical
    differentiation.

    Parameters
    ----------
    beta : float or array_like
        Non-dimensional time parameter.
    mu : float or array_like
        Non-dimensional frequency parameter (cos θ).

    Returns
    -------
    dGdtt : float
        Asymptotic approximation of d²G/dt² at (beta, mu).

    Derivation
    ----------
    Writing G' = f0 + Re[f2A · f2B]:

      df0/dβ  = 4 Σ alpha(n) (2n+3) β^{-2n-4} P_n(μ)

      f2A = (-4i/√(2η)) exp(β²(-μ+iη)/4) · [β-independent phase]
        → df2A/dβ = f2A · (β/2)(−μ + iη)

      Each term in f2B ∝ [β(iμ+η)/2]^{1-2m-2n}
        → df2B/dβ = (1/β) Σ_n (i/η)^n Σ_{m=0}^n (1−2m−2n) d_nm
                                [β(iμ+η)/2]^{1-2m-2n} (μ−iη)^m

      d(f2A·f2B)/dβ = df2A · f2B + f2A · df2B
    """
    nmax = 15

    eta = np.sqrt( 1.0 - mu**2.0 )
    j   = 1j

    # ------------------------------------------------------------------
    # Derivative of f0:  df0/dβ = 4 Σ alpha(n)*(2n+3)*β^(-2n-4)*P_n(μ)
    # ------------------------------------------------------------------
    df0 = 0.0
    for n in range( nmax ):
        df0 += get_alpha( n ) * ( 2*n + 3 ) * beta**( -2*n - 4 ) * sp.special.eval_legendre( n, mu )
    df0 *= 4.0

    # ------------------------------------------------------------------
    # f2A  (identical to dGdt_asymptotic)
    # ------------------------------------------------------------------
    f2a = (
                ( - 1j ) * 4.0
                *
                np.exp( -0.25 * beta**2.0 * mu )
                *
                np.exp(
                            1j
                            *
                            (
                                0.25 * beta**2.0 * eta
                                -
                                0.5 * np.arccos( mu )
                                +
                                np.pi / 4.0
                            )
                        )
                /
                np.sqrt( 2.0 * eta )
            )

    # df2A/dβ = f2A · (β/2)(−μ + iη)
    df2a = f2a * ( 0.5 * beta * ( -mu + j * eta ) )

    # ------------------------------------------------------------------
    # f2B and df2B/dβ computed in a single shared loop.
    # Each term is multiplied by (1−2m−2n)/β for the derivative.
    # ------------------------------------------------------------------
    f2b  = 0.0
    df2b = 0.0
    for n in range( nmax ):
        r0  = ( j / eta )**n
        r1  = 0.0
        dr1 = 0.0
        for m in range( n + 1 ):
            base = (
                        get_dmn( m, n )
                        *
                        ( 0.5 * beta * ( j * mu + eta ) )**( 1.0 - 2.0*m - 2.0*n )
                        *
                        ( mu - j * eta )**m
                    )
            r1  += base
            dr1 += ( 1.0 - 2.0*m - 2.0*n ) * base
        f2b  += r0 * r1
        df2b += r0 * dr1
    df2b /= beta

    # ------------------------------------------------------------------
    # Product rule: d(f2A · f2B)/dβ = df2A · f2B + f2A · df2B
    # ------------------------------------------------------------------
    df2 = df2a * f2b + f2a * df2b

    # # return ( df0 + df2 ).real
    return df2.real


def dGdttt_asymptotic( beta, mu ):
    """
    Compute the asymptotic approximation of d³G/dt³, i.e. the second derivative
    of dGdt_asymptotic with respect to beta, obtained by exact analytical
    differentiation.

    Parameters
    ----------
    beta : float or array_like
        Non-dimensional time parameter.
    mu : float or array_like
        Non-dimensional frequency parameter (cos θ).

    Returns
    -------
    dGdttt : float
        Asymptotic approximation of d³G/dt³ at (beta, mu).

    Derivation
    ----------
    Differentiating dGdtt_asymptotic once more with respect to beta.

      d²f0/dβ²  = -4 Σ alpha(n) (2n+3)(2n+4) β^{-2n-5} P_n(μ)

      Writing h(β) = β²(-μ+iη)/4  so that  f2A = C·exp(h):
        h'  = (β/2)(−μ+iη),   h'' = (−μ+iη)/2
        df2A/dβ  = f2A · h'
        d²f2A/dβ² = f2A · (h'² + h'')

      Define three sums computed in one shared inner loop:
        S0 = Σ_n (i/η)^n Σ_{m=0}^n   1      · d_nm [β(iμ+η)/2]^p (μ-iη)^m = f2B
        S1 = Σ_n (i/η)^n Σ_{m=0}^n   p      · d_nm [β(iμ+η)/2]^p (μ-iη)^m
        S2 = Σ_n (i/η)^n Σ_{m=0}^n p(p-1)  · d_nm [β(iμ+η)/2]^p (μ-iη)^m
        where p = 1−2m−2n.

        df2B/dβ   = S1 / β
        d²f2B/dβ² = S2 / β²   (from d(S1/β)/dβ = S2/β² − S1/β²)

      Leibniz rule:
        d²(f2A·f2B)/dβ² = d²f2A · S0 + 2 df2A · S1/β + f2A · S2/β²
    """
    nmax = 15

    eta = np.sqrt( 1.0 - mu**2.0 )
    j   = 1j

    # ------------------------------------------------------------------
    # d²f0/dβ²  = -4 Σ alpha(n)*(2n+3)*(2n+4)*β^(-2n-5)*P_n(μ)
    # ------------------------------------------------------------------
    d2f0 = 0.0
    for n in range( nmax ):
        d2f0 += get_alpha( n ) * ( 2*n + 3 ) * ( 2*n + 4 ) * beta**( -2*n - 5 ) * sp.special.eval_legendre( n, mu )
    d2f0 *= -4.0

    # ------------------------------------------------------------------
    # f2A  (identical to dGdt_asymptotic)
    # ------------------------------------------------------------------
    f2a = (
                ( - 1j ) * 4.0
                *
                np.exp( -0.25 * beta**2.0 * mu )
                *
                np.exp(
                            1j
                            *
                            (
                                0.25 * beta**2.0 * eta
                                -
                                0.5 * np.arccos( mu )
                                +
                                np.pi / 4.0
                            )
                        )
                /
                np.sqrt( 2.0 * eta )
            )

    # h'(β) = (β/2)(−μ+iη),   h''(β) = (1/2)(−μ+iη)
    hp  = 0.5 * beta * ( -mu + j * eta )
    hpp = 0.5 *         ( -mu + j * eta )

    df2a  = f2a * hp                  # df2A/dβ  = f2A · h'
    d2f2a = f2a * ( hp**2 + hpp )     # d²f2A/dβ² = f2A · (h'² + h'')

    # ------------------------------------------------------------------
    # S0, S1, S2 built in a single shared double loop.
    #   S0  → f2B
    #   S1  → β · df2B/dβ
    #   S2  → β² · d²f2B/dβ²
    # ------------------------------------------------------------------
    S0 = 0.0
    S1 = 0.0
    S2 = 0.0
    for n in range( nmax ):
        r0 = ( j / eta )**n
        s0 = 0.0
        s1 = 0.0
        s2 = 0.0
        for m in range( n + 1 ):
            p    = 1.0 - 2.0*m - 2.0*n
            base = (
                        get_dmn( m, n )
                        *
                        ( 0.5 * beta * ( j * mu + eta ) )**p
                        *
                        ( mu - j * eta )**m
                    )
            s0 += base
            s1 += p           * base
            s2 += p * (p - 1.0) * base
        S0 += r0 * s0
        S1 += r0 * s1
        S2 += r0 * s2

    # ------------------------------------------------------------------
    # Leibniz rule:
    #   d²(f2A·f2B)/dβ² = d²f2A·f2B + 2·df2A·(S1/β) + f2A·(S2/β²)
    # ------------------------------------------------------------------
    d2f2 = d2f2a * S0 + 2.0 * df2a * ( S1 / beta ) + f2a * ( S2 / beta**2 )

    # # return ( d2f0 + d2f2 ).real
    return d2f2.real


def dGdtx_asymptotic( beta, mu ):
    """
    Compute the asymptotic approximation of ∂(dG/dt)/∂μ, i.e. the derivative
    of dGdt_asymptotic with respect to mu, obtained by exact analytical
    differentiation.

    Parameters
    ----------
    beta : float or array_like
        Non-dimensional time parameter.
    mu : float or array_like
        Non-dimensional frequency parameter (cos θ).

    Returns
    -------
    dGdtx : float
        Asymptotic approximation of ∂(dG/dt)/∂μ at (beta, mu).

    Derivation
    ----------
    With η = sqrt(1−μ²),  z1 = iμ+η,  z2 = μ−iη:

      ∂f0/∂μ  = -4 Σ alpha(n) β^{-2n-3} P_n'(μ)
        where P_n'(μ) = n(P_{n-1}(μ) − μ P_n(μ)) / η²  (standard recurrence)

      Log-differentiating f2A w.r.t. μ  (η'=−μ/η, (arccos μ)'=−1/η):
        ∂f2A/∂μ = f2A · [ μ/(2η²) − β²/4 + i(2−β²μ)/(4η) ]

      For each term T_mn = (i/η)^n d_nm [β z1/2]^p z2^m in f2B:
        dz1/dμ = i − μ/η = −z2/η,   dz2/dμ = 1 + iμ/η = z1/η
        ∂T_mn/∂μ = T_mn · [ n·μ/η² − p·z2/(η·z1) + m·z1/(η·z2) ]
      where p = 1−2m−2n.
      Splitting r0=(i/η)^n from base in the loop:
        ∂(r0·base)/∂μ = r0 · [ n·μ/η² · r1 + Σ_m base·(−p·z2/(η·z1) + m·z1/(η·z2)) ]

      Product rule:  ∂(f2A·f2B)/∂μ = ∂f2A/∂μ · f2B + f2A · ∂f2B/∂μ
    """
    nmax = 15

    eta = np.sqrt( 1.0 - mu**2.0 )
    j   = complex( 0.0, 1.0 )

    z1 = j * mu + eta          # iμ + η
    z2 = mu - j * eta          # μ − iη

    # ------------------------------------------------------------------
    # ∂f0/∂μ = -4 Σ alpha(n) β^{-2n-3} P_n'(μ)
    # P_0' = 0;  P_n'(μ) = n(P_{n-1}(μ) − μ P_n(μ)) / η²  for n ≥ 1
    # ------------------------------------------------------------------
    df0_dmu = 0.0
    for n in range( 1, nmax ):
        pn   = sp.special.eval_legendre( n,     mu )
        pnm1 = sp.special.eval_legendre( n - 1, mu )
        dPn  = n * ( pnm1 - mu * pn ) / eta**2
        df0_dmu += get_alpha( n ) * beta**( -2*n - 3 ) * dPn
    df0_dmu *= -4.0

    # ------------------------------------------------------------------
    # f2A  (identical to dGdt_asymptotic)
    # ------------------------------------------------------------------
    f2a = (
                ( - 1j ) * 4.0
                *
                np.exp( -0.25 * beta**2.0 * mu )
                *
                np.exp(
                            j
                            *
                            (
                                0.25 * beta**2.0 * eta
                                -
                                0.5 * np.arccos( mu )
                                +
                                np.pi / 4.0
                            )
                        )
                /
                np.sqrt( 2.0 * eta )
            )

    # ∂f2A/∂μ = f2A · [ μ/(2η²) − β²/4 + i(2 − β²μ)/(4η) ]
    df2a_dmu = f2a * ( mu / ( 2.0 * eta**2 ) - 0.25 * beta**2 + j * ( 2.0 - beta**2 * mu ) / ( 4.0 * eta ) )

    # ------------------------------------------------------------------
    # f2B and ∂f2B/∂μ in a single shared double loop.
    #
    # For each term the log-derivative (excluding the r0 = (i/η)^n factor)
    # contributes:  base · ( −p·z2/(η·z1) + m·z1/(η·z2) )
    # The r0 factor contributes:  r0 · n·μ/η² · r1
    # Together:  ∂(r0·r1)/∂μ = r0 · ( n·μ/η² · r1 + dr1_mu )
    # ------------------------------------------------------------------
    f2b      = 0.0
    df2b_dmu = 0.0
    for n in range( nmax ):
        r0 = ( j / eta )**n
        r1      = 0.0
        dr1_mu  = 0.0
        for m in range( n + 1 ):
            p    = 1.0 - 2.0*m - 2.0*n
            base = (
                        get_dmn( m, n )
                        *
                        ( 0.5 * beta * z1 )**p
                        *
                        z2**m
                    )
            r1 += base
            # log-deriv of base:  -p·z2/(η·z1)  +  m·z1/(η·z2)
            c = -p * z2 / ( eta * z1 )
            if m > 0:
                c += m * z1 / ( eta * z2 )
            dr1_mu += base * c
        f2b      += r0 * r1
        df2b_dmu += r0 * ( n * mu / eta**2 * r1 + dr1_mu )

    # ------------------------------------------------------------------
    # Product rule:  ∂(f2A·f2B)/∂μ = ∂f2A/∂μ·f2B + f2A·∂f2B/∂μ
    # ------------------------------------------------------------------
    df2_dmu = df2a_dmu * f2b + f2a * df2b_dmu

    # # return ( df0_dmu + df2_dmu ).real
    return df2_dmu.real


def dGdttx_asymptotic( beta, mu ):
    """
    Compute the asymptotic approximation of ∂²(dG/dt)/(∂β ∂μ), i.e. the mixed
    derivative of dGdt_asymptotic: first with respect to beta, then with respect
    to mu.  Obtained by exact analytical differentiation.

    Parameters
    ----------
    beta : float or array_like
        Non-dimensional time parameter.
    mu : float or array_like
        Non-dimensional frequency parameter (cos θ).

    Returns
    -------
    dGdttx : float
        Asymptotic approximation of ∂²(dG/dt)/(∂β ∂μ) at (beta, mu).

    Derivation
    ----------
    Differentiate d(f2A·f2B)/dβ  w.r.t. μ via the four-term product rule:

      ∂/∂μ [ f2A_β·f2B + f2A·f2B_β ]
        = f2A_βμ·f2B + f2A_β·f2B_μ + f2A_μ·f2B_β + f2A·f2B_βμ

    Building blocks (η = sqrt(1−μ²), z1 = iμ+η, z2 = μ−iη):

      f2A_β  = f2A · h'         with h' = (β/2)(−μ+iη)
      f2A_μ  = f2A · c_μ        with c_μ = μ/(2η²) − β²/4 + i(2−β²μ)/(4η)
      f2A_βμ = f2A · (c_μ·h' + h'_μ)   with h'_μ = ∂h'/∂μ = −β z1/(2η)

    Four sums built in one shared inner loop (p = 1−2m−2n,
    c_mn = −p z2/(η z1) + m z1/(η z2)):

      S0      = Σ r0·r1                         (= f2B)
      dS0_dmu = Σ r0·(n μ/η² · r1 + Σ_m base·c_mn)    (= ∂f2B/∂μ)
      S1      = Σ r0·(Σ_m p·base)               (= β · f2B_β)
      dS1_dmu = Σ r0·(n μ/η² · s1 + Σ_m p·base·c_mn)  (= β · ∂f2B_β/∂μ)

      f2B_β   = S1/β,   ∂f2B_β/∂μ = dS1_dmu/β

    Mixed partial of f2A·f2B:
      f2A_βμ·S0 + f2A_β·dS0_dmu + f2A_μ·S1/β + f2A·dS1_dmu/β
    """
    nmax = 15

    eta = np.sqrt( 1.0 - mu**2.0 )
    j   = complex( 0.0, 1.0 )

    z1 = j * mu + eta          # iμ + η
    z2 = mu - j * eta          # μ − iη

    # ------------------------------------------------------------------
    # Mixed derivative of f0:
    #   ∂²f0/(∂β∂μ) = 4 Σ_{n≥1} alpha(n)·(2n+3)·β^{-2n-4}·P_n'(μ)
    # ------------------------------------------------------------------
    d2f0_dxdt = 0.0
    for n in range( 1, nmax ):
        pn   = sp.special.eval_legendre( n,     mu )
        pnm1 = sp.special.eval_legendre( n - 1, mu )
        dPn  = n * ( pnm1 - mu * pn ) / eta**2
        d2f0_dxdt += get_alpha( n ) * ( 2*n + 3 ) * beta**( -2*n - 4 ) * dPn
    d2f0_dxdt *= 4.0

    # ------------------------------------------------------------------
    # f2A
    # ------------------------------------------------------------------
    f2a = (
                complex( 0.0, -4.0 )
                *
                np.exp( -0.25 * beta**2.0 * mu )
                *
                np.exp(
                            j
                            *
                            (
                                0.25 * beta**2.0 * eta
                                -
                                0.5 * np.arccos( mu )
                                +
                                np.pi / 4.0
                            )
                        )
                /
                np.sqrt( 2.0 * eta )
            )

    # β-derivative of f2A
    hp       = 0.5 * beta * ( -mu + j * eta )         # h'(β) = dh/dβ
    f2a_beta = f2a * hp

    # μ-derivative of f2A
    c_mu   = mu / ( 2.0 * eta**2 ) - 0.25 * beta**2 + j * ( 2.0 - beta**2 * mu ) / ( 4.0 * eta )
    f2a_mu = f2a * c_mu

    # mixed (β,μ)-derivative of f2A:
    #   h'_μ = ∂/∂μ [(β/2)(−μ+iη)] = (β/2)(−1 + i·η') = −β(η+iμ)/(2η) = −β z1/(2η)
    hp_mu       = -0.5 * beta * z1 / eta
    f2a_beta_mu = f2a * ( c_mu * hp + hp_mu )

    # ------------------------------------------------------------------
    # Four sums in a single shared double loop.
    # ------------------------------------------------------------------
    S0      = 0.0
    dS0_dmu = 0.0
    S1      = 0.0
    dS1_dmu = 0.0
    for n in range( nmax ):
        r0 = ( j / eta )**n
        r1      = 0.0
        dr1_mu  = 0.0
        s1      = 0.0
        ds1_mu  = 0.0
        for m in range( n + 1 ):
            p    = 1.0 - 2.0*m - 2.0*n
            base = (
                        get_dmn( m, n )
                        *
                        ( 0.5 * beta * z1 )**p
                        *
                        z2**m
                    )
            # log-deriv of base w.r.t. μ (z1^p and z2^m parts only)
            c = -p * z2 / ( eta * z1 )
            if m > 0:
                c += m * z1 / ( eta * z2 )
            r1     += base
            dr1_mu += base * c
            s1     += p * base
            ds1_mu += p * base * c
        # accumulate including the r0 log-deriv  n·μ/η²
        nm_eta2   = n * mu / eta**2
        S0      += r0 * r1
        dS0_dmu += r0 * ( nm_eta2 * r1 + dr1_mu )
        S1      += r0 * s1
        dS1_dmu += r0 * ( nm_eta2 * s1 + ds1_mu )

    # ------------------------------------------------------------------
    # Four-term product rule:
    #   ∂²(f2A·f2B)/(∂β∂μ) = f2A_βμ·S0 + f2A_β·dS0_dmu
    #                        + f2A_μ·S1/β + f2A·dS1_dmu/β
    # ------------------------------------------------------------------
    d2f2_mixed = ( f2a_beta_mu * S0
                   + f2a_beta   * dS0_dmu
                   + f2a_mu     * ( S1      / beta )
                   + f2a        * ( dS1_dmu / beta ) )

    # # return ( d2f0_dxdt + d2f2_mixed ).real
    return d2f2_mixed.real


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import h5py
    fipath = r"D:\sergio\developments\SeaMotionsTimeDev\aux_data\0_integrals_database\1_time_domain\Gt.h5"
    with h5py.File( fipath, 'r' ) as f:
        beta = f['beta'][:]
        mu   = f['mu'][:]
        fcn  = f['fcn'][:]
    
    pos = beta < 8.0
    beta_sub = beta[pos]
    pos2 = beta > 7.0
    beta_sub_2 = beta[pos2]
    fcn_sub  = fcn[pos, :]
    fcn_sub_2  = fcn[pos2, :]
    # beta = 7.5
    mu = 1e-4

    f = 0.0
    for i in range( 100 ):
        f += ( 
                ( 1j * beta_sub )**i 
                / 
                sp.special.factorial( i ) 
                * 
                sp.special.gamma( 0.5 * i + 1.5 )
                *
                sp.special.lpmv( 0, 0.5*i + 0.5, mu )
            )

    F = - 2.0 * complex( 0.0, 1 ) * f
    G = dGdt_asymptotic( beta_sub_2, mu )
    plt.plot( beta_sub, F.real, label="F" )
    plt.plot( beta_sub, 2*fcn_sub[:, 0], label="fcn" )
    plt.plot( beta_sub_2, G, label="G" )
    plt.plot( beta_sub_2, 2*fcn_sub_2[:, 0], label="fcn_2" )
    plt.legend( )
    plt.show( )
    # print( F )
