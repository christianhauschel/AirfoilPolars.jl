import Xfoil
using AirfoilFast
using .AirfoilPolars: Polar

name_af = ""

function init(af::Airfoil)
    Xfoil.set_coordinates(af.x, af.y)
    global name_af
    name_af = af.name
end

function solve(
    alpha::Vector,
    re::Number;
    n_iter::Int = 100,
    n_crit::Int = 9,
    xtrip::Tuple = (1.0, 1.0),
    mach::Number = 0.0,
)
    n_alpha = length(alpha)
    cl = zeros(n_alpha)
    cd = zeros(n_alpha)
    cdp = zeros(n_alpha)
    cm = zeros(n_alpha)
    converged = zeros(n_alpha)
    for i = 1:n_alpha
        cl[i], cd[i], cdp[i], cm[i], converged[i] = Xfoil.solve_alpha(
            alpha[i],
            re;
            iter = n_iter,
            reinit = i == 1,
            ncrit = n_crit,
            xtrip = xtrip,
            mach = mach,
        )
    end
    global name_af
    return Polar(re, alpha, cl, cd, cm, mach, n_crit, xtrip, name_af)
end