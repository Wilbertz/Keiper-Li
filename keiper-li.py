# keiper_li.py
from mpmath import mp, zetazero, power, mpc, mpf

def set_precision(dps):
    mp.dps = dps

def rho_from_index(k):
    """Return the k-th (positive imaginary) nontrivial zero as a complex mp.mpc:
       rho_k = 1/2 + i * gamma_k where gamma_k = zetazero(k).
       Note: zetazero(1) is the first positive imaginary part.
    """
    gamma_k = zetazero(k)
    return mpc(mpf(1)/2, gamma_k)

def lambda_by_zeros(n, maxzeros=2000, tol=None, verbose=False):
    """
    Compute Li coefficient lambda_n by summing contributions of zeta zeros in conjugate pairs.
    Parameters:
      n        : index of Li coefficient (positive integer)
      maxzeros : maximum number of positive zeros to use (will try to stop earlier if converged)
      tol      : stopping tolerance (mp.mpf). If None, auto-chosen from mp.dps.
      verbose  : print convergence progress
    Returns:
      (lambda_n, used_zeros)
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if tol is None:
        # choose tolerance a few orders below precision
        tol = mp.mpf(10) ** (-(mp.dps - 6))

    total = mp.mpf('0')
    last_total = None

    # Sum pairwise: for each positive imaginary zero rho, add contribution of rho and conj(rho)
    for k in range(1, maxzeros + 1):
        rho = rho_from_index(k)
        # contribution from the pair (rho and conj(rho)):
        #  (1 - (1 - 1/rho)^n) + (1 - (1 - 1/rho_conj)^n)
        # = 2 - (1 - 1/rho)^n - (1 - 1/rho_conj)^n
        term = mp.mpf(2) - power(1 - 1/rho, n) - power(1 - 1/mp.conj(rho), n)
        total += term

        # convergence check every few zeros
        if k % 10 == 0:
            if last_total is None:
                last_total = total
            else:
                diff = abs(total - last_total)
                if verbose:
                    print(f"zeros used={k}, partial lambda_{n} ≈ {total}, diff={diff}")
                if diff < tol:
                    return (mp.nstr(total, mp.dps), k)
                last_total = total

    # if we reach here, didn't converge within maxzeros
    return (mp.nstr(total, mp.dps), maxzeros)


if __name__ == "__main__":
    # Example usage:
    # compute first 10 Li coefficients with ~80 decimal digits precision
    set_precision(80)
    print("mp.dps =", mp.dps)
    results = {}
    for n in range(1, 2):
        lam, used = lambda_by_zeros(n, maxzeros=200, tol=mp.mpf(10) ** (-70), verbose=True)
        print(f"lambda_{n} (approx) =\n{lam}\n  zeros used: {used}\n")
        results[n] = lam
