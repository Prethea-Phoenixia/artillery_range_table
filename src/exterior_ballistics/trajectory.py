from __future__ import annotations

from bisect import bisect
from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import cached_property, total_ordering
from math import atan2, cos, sin, inf

from .dekker import dekker

R_e = 6356766  # nominal earth's radius
g0 = 9.80665  # nominal earth's gravitational acceleration
R = 287.05287  # J/(K · kg), R / M
kappa = 1.40  # adiabatic index of atmospheric gas


@dataclass(frozen=True)
class AtmosphericCondition:
    """
    Attributes
    ----------
    c: float
        Speed of sound in m/s
    rho : float
        Density in kg m^-3
    g : float
        Acceleration due gravity in m/s^2
    """

    c: float
    rho: float
    g: float


def delta_t(H: float) -> float:
    dT = 0
    gradients = (
        (0.00e3, -6.50e-3),
        (11.00e3, 0.00),
        (20.00e3, 1.00e-3),
        (32.00e3, 2.80e-3),
        (47.00e3, 0.00e-3),
        (51.00e3, -2.80e-3),
        (71.00e3, -2.00e-3),
        (80.00e3, 0),
        (inf, 0),
    )
    for i in range(len(gradients) - 1):
        Hb, beta = gradients[i]
        Hc, _ = gradients[i + 1]
        if H > Hc:
            dT = dT + beta * (Hc - Hb)
        else:
            dT = dT + beta * (H - Hb)
            break

    return dT


def atmosphere(h: float, p: float, ambient_temperature: float) -> AtmosphericCondition:
    """
    ICAO standard atmosphere

    Parameters
    ----------
    h: float
        Geometrical altitude
    p: float
        pressure
    ambient_temperature: float
        reference temperature
    """
    g = g0 * (R_e / (R_e + h)) ** 2  # approximate local gravitational acceleration
    H = R_e * h / (R_e + h)  # H: approximate geopotential height in kilometer
    T = delta_t(H) + ambient_temperature
    c = (kappa * R * T) ** 0.5
    rho = p / (R * T)

    return AtmosphericCondition(c=c, rho=rho, g=g)


def cd43(Ma: float) -> float:
    """
    Curve fit for the 1943 drag curve, a.la[1]. Maximum error in the domain (< Ma.5) is less than 3.27%

    Parameters
    ----------
    Ma: float
        Projectile's Mach number

    Returns
    -------
    cd43: float
        The projectile's drag coefficient according to the 1943 drag curve.


    References
    ----------
    - **[1]** NI Qingle，WANG Yushi，WEN Quan，ZHANG Zhibiao, Empirical Formulas of Projectile Air Resistance
     Law in Whole Definition Domain in Analytic Function, Journal of Projectiles, Rockets, Missiles, and
     Guidance, Vol.36, No.6, Dec. 2016.

    """
    if Ma < 0.79:
        return 0.157
    elif Ma > 3.95:
        return 0.260
    Cd = 0.15494 - 0.38989 * Ma**2 + 0.38100 * Ma**4 - 0.23778 * Ma**6 + 0.11412 * Ma**8
    Cd /= 1 - 2.6458 * Ma**2 + 2.7429 * Ma**4 - 1.4990 * Ma**6 + 0.47141 * Ma**8
    return Cd


@dataclass(frozen=True)
class Increment:
    v_gr: float
    v_h: float
    gr: float
    h: float
    p: float

    def __mul__(self, value: float) -> Increment:
        return Increment(
            v_gr=self.v_gr * value,
            v_h=self.v_h * value,
            gr=self.gr * value,
            h=self.h * value,
            p=self.p * value,
            # m=self.m * value,
        )

    def __add__(self, other: Increment) -> Increment:
        return Increment(
            v_gr=self.v_gr + other.v_gr,
            v_h=self.v_h + other.v_h,
            gr=self.gr + other.gr,
            h=self.h + other.h,
            p=self.p + other.p,
            # m=self.m + other.m,
        )

    def __truediv__(self, value: float) -> Increment:
        return Increment(
            v_gr=self.v_gr / value,
            v_h=self.v_h / value,
            gr=self.gr / value,
            h=self.h / value,
            p=self.p / value,
            # m=self.m / value,
        )

    def __rmul__(self, value: float) -> Increment:
        return self * value


@dataclass(frozen=True)
class Gradient:
    dv_gr: float
    dv_h: float
    dgr: float
    dh: float
    dp: float

    def __mul__(self, dt: float) -> Increment:
        return Increment(
            v_gr=self.dv_gr * dt,
            v_h=self.dv_h * dt,
            gr=self.dgr * dt,
            h=self.dh * dt,
            p=self.dp * dt,  # m=self.dm * dt
        )

    def __rmul__(self, dt: float) -> Increment:  # make the operation commutative
        return self * dt


@total_ordering
@dataclass(frozen=True)
class Point:
    """
    class

    Attributes
    ----------
    t: float
        time of flight
    v_gr, v_h: float
        the horizontal and vertical components of velocity
    gr: float
        ground range (as projected onto an imaginary sphere of radius R_e + h_0)
    h: float
        height of shell (above reference height h_0)
    t: float
        time of flight
    """

    trajectory: BaseTrajectory = field(repr=False)
    t: float
    v_gr: float
    v_h: float
    gr: float
    h: float
    p: float

    @cached_property
    def v(self) -> float:
        return (self.v_gr**2 + self.v_h**2) ** 0.5

    @cached_property
    def theta(self) -> float:
        return atan2(self.v_h, self.v_gr)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Point) and self.trajectory == other.trajectory:
            return self.t == other.t
        elif isinstance(other, float) or isinstance(other, int):
            return self.t == other
        return False

    def __lt__(self, other: Point | float) -> bool:
        if isinstance(other, Point) and self.trajectory == other.trajectory:
            return self.t < other.t
        else:
            return self.t < other

    def increment(self, delta: Increment, dt: float) -> Point:
        return Point(
            trajectory=self.trajectory,
            t=self.t + dt,
            v_gr=self.v_gr + delta.v_gr,
            v_h=self.v_h + delta.v_h,
            gr=self.gr + delta.gr,
            h=self.h + delta.h,
            p=self.p + delta.p,
            # m=self.m + delta.m,
        )


@dataclass(frozen=True)
class BaseTrajectory:

    sectional_density: float  # kg/m^2
    velocity: float
    elevation: float  # initial shell elevation, in radian
    i_43: float
    gun_alt: float  # altitude of the gun site above mean sea level
    tgt_alt: float
    ambient_pressure: float  # initial pressure
    ambient_temperature: float  # change in virtual temperature
    time_increment: float

    @cached_property
    def h_t(self) -> float:
        return self.tgt_alt - self.gun_alt

    def __iter__(self) -> Iterable[Point]:
        return iter(self.points)

    @property
    def points(self) -> list[Point]:
        raise NotImplementedError

    def get_gradient(self, point: Point) -> Gradient:
        raise NotImplementedError("")

    def propagate_rk4(self, point: Point, dt: float) -> Point:
        """in particular, the mid-point flavour"""
        k1 = self.get_gradient(point)
        k2 = self.get_gradient(point.increment(delta=0.5 * k1 * dt, dt=0.5 * dt))
        k3 = self.get_gradient(point.increment(delta=0.5 * k2 * dt, dt=0.5 * dt))
        k4 = self.get_gradient(point.increment(delta=k3 * dt, dt=dt))

        return point.increment(delta=k1 * dt / 6 + k2 * dt / 3 + k3 * dt / 3 + k4 * dt / 6, dt=dt)

    def find_intercept(self, p_i: Point, p_j: Point) -> Point:
        h_t = self.tgt_alt - self.gun_alt

        def f(t: float) -> float:
            p = self.at_time(t)
            return p.h - h_t

        if h_t == p_i.h:
            return p_i
        elif h_t == p_j.h:
            return p_j
        else:
            t, _ = dekker(f=f, x_0=p_i.t, x_1=p_j.t, tol=1e-6)
            return self.at_time(t)

    def find_descending_intercept(self) -> Point:
        return self.find_intercept(p_i=self.crest, p_j=self.points[-1])

    def find_ascending_intercept(self) -> Point:
        return self.find_intercept(p_i=self.points[0], p_j=self.crest)

    def at_time(self, t: float) -> Point:
        """
        calculate the Point object corresponding to supplied time.
        """
        if t < 0 or t > self.points[-1].t:
            raise ValueError(f"trajectory not calculated to {t} s")
        i = max(min(bisect(self.points, t) - 1, len(self.points) - 1), 0)
        p = self.points[i]

        if p.t == t:
            return p
        else:
            return self.propagate_rk4(p, t - p.t)

    @cached_property
    def crest(self) -> Point:
        """use a bisection based algorithm to iterating
        find and return the point on trajectory where product of velocity and
        theta changes sign. The product has the same sign as the vertical component
        of velocity, namely sin(theta) * velocity.
        relies on walking the point records calculated up until this point,
        and is done in a left inclusive, right exclusive manner"""

        def f(t: float) -> float:
            p = self.at_time(t)
            return p.v_h

        for p_i, p_j in zip(self.points[:-1], self.points[1:]):
            if p_i.v_h * p_j.v_h <= 0:
                if p_i.v_h == 0.0:
                    return p_i
                elif p_j.v_h == 0.0:
                    return p_j
                else:
                    t, _ = dekker(f=f, x_0=p_i.t, x_1=p_j.t, tol=1e-6)
                    return self.at_time(t)

        return self.points[0]


@dataclass(frozen=True)
class ArtilleryTrajectory(BaseTrajectory):

    @cached_property
    def points(self) -> list[Point]:
        # populate the entire trajectory to the class-defined end conditions

        p = Point(
            trajectory=self,
            t=0.0,
            v_gr=self.velocity * cos(self.elevation),
            v_h=self.velocity * sin(self.elevation),
            gr=0.0,
            h=0.0,
            p=self.ambient_pressure,
        )
        p_spi = self.get_spi(p)
        while self.get_spi(self.propagate_rk4(point=p, dt=self.time_increment)) > p_spi:
            self.time_increment *= 0.5

        points = [p]

        while p.h >= min(self.h_t, 0):
            points.append(p := self.propagate_rk4(p, self.time_increment))

        return points

    def get_spi(self, p: Point) -> float:
        """get the specific energy of a point"""
        return 0.5 * p.v**2 - (g0 * R_e**2) / (R_e + self.gun_alt + p.h)

    def get_gradient(self, point: Point) -> Gradient:
        i_43, h_0, sd = self.i_43, self.gun_alt, self.sectional_density
        v_gr, v_h, v, p, h = point.v_gr, point.v_h, point.v, point.p, point.h

        ac = atmosphere(h=h + h_0, p=p, ambient_temperature=self.ambient_temperature)
        g, rho, c = ac.g, ac.rho, ac.c

        cd = cd43(v / c)

        dv_gr = -0.5 * rho * v_gr * v * cd * i_43 / sd - v_gr * v_h / (h + h_0 + R_e)
        dv_h = -0.5 * rho * v_h * v * cd * i_43 / sd + v_gr**2 / (h + h_0 + R_e) - g
        dgr, dh, dp = v_gr / (1 + h / (h_0 + R_e)), v_h, -rho * g * v_h

        return Gradient(dv_gr=dv_gr, dv_h=dv_h, dgr=dgr, dh=dh, dp=dp)


if __name__ == "__main__":
    pass
