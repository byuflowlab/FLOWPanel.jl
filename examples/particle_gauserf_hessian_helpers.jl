module ParticleGauserfHessian

export HAVE_SF, ZETA0, zeta_sigma, green_profile, phi_gauserf,
       phi_prime_gauserf, hess_green_times_gamma,
       hess_green_times_gamma_components

const HAVE_SF = let
    try
        @eval import SpecialFunctions
        true
    catch
        false
    end
end

if HAVE_SF
    _erf(x) = SpecialFunctions.erf(x)
else
    function _erf(x::T) where {T}
        z = abs(x)
        if z < 4
            s = z
            term = z
            n = 0
            while abs(term) > 1e-18 * max(abs(s), one(T)) && n < 200
                n += 1
                term *= -z^2 / n
                s += term / (2n + 1)
            end
            e = 2 / sqrt(T(pi)) * s
        else
            e = one(T) - exp(-z^2) / (z * sqrt(T(pi)))
        end
        return x < 0 ? -e : e
    end
end

const ZETA0 = 1 / (2pi)^(3 / 2)

g_gauserf(rho) = _erf(rho / sqrt(2)) - sqrt(2 / pi) * rho * exp(-rho^2 / 2)

function phi_gauserf(rho2)
    if rho2 < 1e-4
        return sqrt(2 / pi) * (1 / 3 - rho2 / 10 + rho2^2 / 56 - rho2^3 / 432)
    else
        rho = sqrt(rho2)
        return g_gauserf(rho) / (rho2 * rho)
    end
end

function phi_prime_gauserf(rho2)
    if rho2 < 1e-4
        return sqrt(2 / pi) * (-1 / 10 + rho2 / 28 - rho2^2 / 144 + rho2^3 / 864)
    else
        rho = sqrt(rho2)
        return (sqrt(2 / pi) * rho^3 * exp(-rho2 / 2) - 3 * g_gauserf(rho)) /
               (2 * rho2^2 * rho)
    end
end

zeta_sigma(rho2, sigma) = ZETA0 * exp(-rho2 / 2) / sigma^3

function green_profile(rho2)
    if rho2 < 1e-4
        return sqrt(2 / pi) *
               (1 - rho2 / 6 + rho2^2 / 40 - rho2^3 / 336 + rho2^4 / 3456)
    else
        rho = sqrt(rho2)
        return _erf(rho / sqrt(2)) / rho
    end
end

@inline function hess_green_times_gamma_components(rx, ry, rz, sigma, g1, g2, g3)
    s = (rx * rx + ry * ry + rz * rz) / (sigma * sigma)
    A = -phi_gauserf(s) / (4pi * sigma^3)
    B = -phi_prime_gauserf(s) / (2pi * sigma^5)
    yg = rx * g1 + ry * g2 + rz * g3
    return A * g1 + B * yg * rx,
           A * g2 + B * yg * ry,
           A * g3 + B * yg * rz
end

function hess_green_times_gamma(y, sigma, Gamma)
    q1, q2, q3 = hess_green_times_gamma_components(
        y[1], y[2], y[3], sigma, Gamma[1], Gamma[2], Gamma[3])
    return [q1, q2, q3]
end

end # module
