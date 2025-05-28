from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from math import pi, ceil
from multiprocessing import Pool
from typing import TYPE_CHECKING

from tabulate import tabulate

from .dekker import dekker
from .gss import gss_max
from .trajectory import ArtilleryTrajectory, BaseTrajectory

if TYPE_CHECKING:
    from .trajectory import Point


ONE_ARCMIN = pi / 180 / 60


def rad_to_deg(rad: float) -> float:
    return rad * 180 / pi


def deg_to_rad(deg: float) -> float:
    return deg * pi / 180


def dd_to_dms(dd: float) -> tuple[float, float, float]:
    mult = -1 if dd < 0 else 1
    mnt, sec = divmod(abs(dd) * 3600, 60)
    deg, mnt = divmod(mnt, 60)
    return mult * deg, mult * mnt, mult * sec


@dataclass(frozen=True)
class RangeTableEntry:
    """Class implementing a Soviet style artillery table entry"""

    ground_range: float  # in meters
    gun_elevation: float  # stored as radians, convert to value when necessary
    impact_velocity: float  # m/s
    impact_angle: float  # same as above
    time_of_flight: float  # seconds
    apex_height: float  # m

    marker: str = ""  # M for maximum range shots

    def format(self) -> tuple[str, str, str, str, str, str, str]:

        d, m, s = dd_to_dms(rad_to_deg(self.gun_elevation))

        return (
            f"{self.marker: ^1}",
            f"{self.ground_range: >6,.0f}",
            f"{d: >3.0f}° {abs(m): >2.0f}'",
            f"{rad_to_deg(self.impact_angle): >3.0f}°",
            f"{self.impact_velocity: >4.0f}",
            f"{self.time_of_flight: >4.1f}",
            f"{self.apex_height: >.0f}",
        )

    def __str__(self) -> str:
        return " ".join(self.format())


@dataclass(frozen=True)
class RangeTable:
    sectional_density: float
    velocity: float
    i_43: float = 1.0
    gun_alt: float = 0.0
    ambient_pressure: float = 101.325e3
    ambient_temperature: float = 288.15
    tgt_alt: float = 0.0
    time_increment: float = 1.0
    range_increment: float = 1000
    max_elev: float = deg_to_rad(90)
    min_elev: float = deg_to_rad(-90)

    def at_theta(self, theta: float) -> BaseTrajectory:
        """Return a Trajectory object with the current load at a given gun elevation"""

        return ArtilleryTrajectory(
            sectional_density=self.sectional_density,
            velocity=self.velocity,
            elevation=theta,
            i_43=self.i_43,
            gun_alt=self.gun_alt,
            tgt_alt=self.tgt_alt,
            ambient_pressure=self.ambient_pressure,
            ambient_temperature=self.ambient_temperature,
            time_increment=self.time_increment,
        )

    def crest_at_theta(self, theta: float) -> Point:
        return self.at_theta(theta).crest

    def ascending_intercept_at_theta(self, theta: float) -> Point:
        return self.at_theta(theta).find_ascending_intercept()  # self.tgt_alt)

    def descending_intercept_at_theta(self, theta: float) -> Point:
        return self.at_theta(theta).find_descending_intercept()  # self.tgt_alt)

    @cached_property
    def min_air_range(self) -> float:
        return self.ascending_intercept_at_theta(self.max_elev).gr

    @cached_property
    def max_air_range(self) -> float:
        return self.ascending_intercept_at_theta(self.crit_elev).gr

    @cached_property
    def low_min_ground_range(self) -> float:
        return self.descending_intercept_at_theta(self.crit_elev).gr

    @cached_property
    def max_ground_range(self) -> float:
        return self.descending_intercept_at_theta(self.opt_elev).gr

    @cached_property
    def high_min_ground_range(self) -> float:
        return self.descending_intercept_at_theta(self.max_elev).gr

    @cached_property
    def crit_elev(self) -> float:
        h_t = self.tgt_alt - self.gun_alt

        def f(theta: float) -> float:
            return self.crest_at_theta(theta).h - h_t

        if f(0) >= 0:
            return self.min_elev
        else:
            if f(self.max_elev) < 0:
                raise ValueError("gun cannot be brought to bear.")

            a, b = dekker(f=f, x_0=self.min_elev, x_1=self.max_elev, tol=ONE_ARCMIN)
            return max(a, b)

    @cached_property
    def opt_elev(self) -> float:
        """Find the gun elevation angle to achieve maximum range with the projectile."""

        def f(theta: float) -> float:
            p = self.descending_intercept_at_theta(theta)
            return abs(p.gr)

        return 0.5 * sum(gss_max(f=f, x_0=self.crit_elev, x_1=self.max_elev, tol=ONE_ARCMIN))

    def ground_entry_at_range(self, gr: float, theta_min: float, theta_max: float, marker: str = "") -> RangeTableEntry:

        def f(theta: float) -> float:
            p = self.descending_intercept_at_theta(theta)
            return p.gr - gr

        if f(theta_min) == 0:
            theta = theta_min
        elif f(theta_max) == 0:
            theta = theta_max
        else:
            theta, _ = dekker(f=f, x_0=theta_min, x_1=theta_max, tol=ONE_ARCMIN)
        return self.ground_entry_at_theta(theta=theta, marker=marker)

    def ground_entry_at_theta(self, theta: float, marker: str = "") -> RangeTableEntry:
        p = self.descending_intercept_at_theta(theta)

        return RangeTableEntry(
            marker=marker,
            ground_range=p.gr,
            time_of_flight=p.t,
            gun_elevation=theta,
            impact_angle=p.theta,
            impact_velocity=p.v,
            apex_height=p.trajectory.crest.h + p.trajectory.gun_alt,
        )

    def low_ground_entry_at_range(self, gr: float, marker: str = "") -> RangeTableEntry:
        return self.ground_entry_at_range(gr, theta_min=self.crit_elev, theta_max=self.opt_elev, marker=marker)

    def high_ground_entry_at_range(self, gr: float, marker: str = "") -> RangeTableEntry:
        return self.ground_entry_at_range(gr, theta_min=self.opt_elev, theta_max=self.max_elev, marker=marker)

    def air_entry_at_theta(self, theta: float, marker: str = "") -> RangeTableEntry:
        p = self.ascending_intercept_at_theta(theta)
        return RangeTableEntry(
            marker=marker,
            ground_range=p.gr,
            time_of_flight=p.t,
            gun_elevation=theta,
            impact_angle=p.theta,
            impact_velocity=p.v,
            apex_height=p.trajectory.crest.h + p.trajectory.gun_alt,
        )

    def air_entry_at_range(self, gr: float, marker: str = "") -> RangeTableEntry:
        def f(theta: float) -> float:
            p = self.ascending_intercept_at_theta(theta)
            if p:
                return p.gr - gr
            else:
                raise ValueError(f"ascending intercept is invalid at {theta}")

        theta_min = self.crit_elev
        theta_max = self.max_elev

        if f(theta_min) == 0:
            theta = theta_min
        elif f(theta_max) == 0:
            theta = theta_max
        else:
            theta, _ = dekker(f=f, x_0=theta_min, x_1=theta_max, tol=ONE_ARCMIN)

        return self.air_entry_at_theta(theta=theta, marker=marker)

    def generate_notches(self, min_range: float, max_range: float) -> list[int]:
        notches = []
        min_range, max_range = min((min_range, max_range)), max((min_range, max_range))

        range_notch = ceil(min_range / self.range_increment) * self.range_increment

        while range_notch < max_range:
            notches.append(range_notch)
            range_notch += self.range_increment

        # if 0 in notches:
        #     notches.remove(0)

        return notches

    @cached_property
    def air_entries(self) -> list[RangeTableEntry]:
        if self.tgt_alt < self.gun_alt:
            raise ValueError("anti-air solutions does not exist for gun sited above the target.")

        grs = self.generate_notches(min_range=self.min_air_range, max_range=self.max_air_range)

        with Pool() as p:
            return [
                *p.map(self.air_entry_at_range, grs),
                self.air_entry_at_theta(theta=self.crit_elev, marker="M"),
            ]

    @cached_property
    def ground_entries(self) -> list[RangeTableEntry]:

        grs_low = self.generate_notches(self.low_min_ground_range, self.max_ground_range)
        grs_high = self.generate_notches(self.high_min_ground_range, self.max_ground_range)

        with Pool() as p:
            return [
                *p.map(self.low_ground_entry_at_range, grs_low),
                self.ground_entry_at_theta(theta=self.opt_elev, marker="M"),
                *p.map(self.high_ground_entry_at_range, grs_high[::-1]),
            ]

        # return (
        #     [self.low_ground_entry_at_range(gr) for gr in grs_low]
        #     + [self.ground_entry_at_theta(theta=self.opt_elev, marker="M")]
        #     + [self.high_ground_entry_at_range(gr) for gr in grs_high[::-1]]
        # )

    @staticmethod
    def prettyprint(entries: list[RangeTableEntry]) -> str:
        return tabulate(
            (v.format() for v in entries),
            headers=[
                "Range\n\n\n(m)",
                "Gun\nElev.",
                "Impact\nAngle",
                "Impact\nVelocity\n\n(m/s)",
                "Time\nto\nTarget\n(s)",
            ],
            disable_numparse=True,
            stralign=None,
        )


if __name__ == "__main__":
    """
    Calculates an example range table for the 122mm howitzer M-30.
    """
    rt = RangeTable(
        velocity=510,
        sectional_density=21.76 / (0.25 * pi * 121.92e-3**2),
        i_43=1.0,
        gun_alt=0,
        tgt_alt=0,
        ambient_pressure=101.325e3,
        ambient_temperature=288.15,
    )

    print(RangeTable.prettyprint(rt.ground_entries))
    print(RangeTable.prettyprint(rt.air_entries))
