"""Classical orbital elements and Cartesian <-> Kepler conversions."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import WGS84
from .errors import ConvergenceError, FieldRangeError, NonPhysicalOrbitError
from .time_utils import Epoch

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class ClassicalElements:
    """J2 mean Keplerian elements at ``epoch``.

    Angles are radians; frame is the service inertial frame
    (see docs/conventions.md).
    """

    semi_major_axis: float
    eccentricity: float
    inclination: float
    raan: float
    argument_of_periapsis: float
    mean_anomaly: float
    epoch: Epoch


def solve_kepler(mean_anomaly: float, eccentricity: float) -> float:
    """Eccentric anomaly E from M and e (radians)."""
    if not (math.isfinite(mean_anomaly) and math.isfinite(eccentricity)):
        raise NonPhysicalOrbitError("non-finite Kepler inputs")
    if eccentricity < 0.0 or eccentricity >= 1.0:
        raise NonPhysicalOrbitError(
            f"eccentricity {eccentricity} outside [0, 1)", location="eccentricity"
        )
    m = math.fmod(mean_anomaly, TWO_PI)
    if m > math.pi:
        m -= TWO_PI
    elif m < -math.pi:
        m += TWO_PI
    if eccentricity < 0.8:
        anomaly = m
    else:
        anomaly = math.copysign(math.pi, m) if abs(m) > math.pi / 2.0 else m
    for _ in range(100):
        f = anomaly - eccentricity * math.sin(anomaly) - m
        fp = 1.0 - eccentricity * math.cos(anomaly)
        step = f / fp
        anomaly -= step
        if abs(step) <= 1e-14 * max(1.0, abs(anomaly)):
            return anomaly
    raise ConvergenceError("Kepler equation did not converge", location="mean_anomaly")


def _validate_elements(elements: ClassicalElements) -> None:
    for name in (
        "semi_major_axis",
        "eccentricity",
        "inclination",
        "raan",
        "argument_of_periapsis",
        "mean_anomaly",
    ):
        if not math.isfinite(getattr(elements, name)):
            raise NonPhysicalOrbitError(f"{name} is not finite", location=name)
    if elements.semi_major_axis <= 0.0:
        raise NonPhysicalOrbitError(
            f"semi-major axis must be positive, got {elements.semi_major_axis}",
            location="semi_major_axis",
        )
    if not 0.0 <= elements.eccentricity < 1.0:
        raise NonPhysicalOrbitError(
            f"eccentricity {elements.eccentricity} outside [0, 1)",
            location="eccentricity",
        )
    if not 0.0 <= elements.inclination <= math.pi:
        raise FieldRangeError(
            f"inclination {elements.inclination} outside [0, pi]",
            location="inclination",
        )
    min_radius = elements.semi_major_axis * (1.0 - elements.eccentricity)
    if min_radius < WGS84.radius_equator:
        raise NonPhysicalOrbitError(
            f"perigee radius {min_radius:.1f} m is inside the Earth "
            f"(equatorial radius {WGS84.radius_equator:.0f} m)",
            location="semi_major_axis",
        )


def coe_to_eci(elements: ClassicalElements, mu: float = WGS84.mu):
    """Osculating inertial (position, velocity) at the element epoch."""
    _validate_elements(elements)
    a = elements.semi_major_axis
    e = elements.eccentricity
    eccentric = solve_kepler(elements.mean_anomaly, e)
    cos_e, sin_e = math.cos(eccentric), math.sin(eccentric)

    radius = a * (1.0 - e * cos_e)
    perigee_x = a * (cos_e - e)
    perigee_y = a * math.sqrt(max(1.0 - e * e, 0.0)) * sin_e
    speed = math.sqrt(mu * a) / radius
    perigee_vx = -speed * sin_e
    perigee_vy = speed * math.sqrt(max(1.0 - e * e, 0.0)) * cos_e

    cos_o, sin_o = math.cos(elements.raan), math.sin(elements.raan)
    cos_i, sin_i = math.cos(elements.inclination), math.sin(elements.inclination)
    cos_w, sin_w = math.cos(elements.argument_of_periapsis), math.sin(
        elements.argument_of_periapsis
    )
    r11 = cos_o * cos_w - sin_o * sin_w * cos_i
    r12 = -cos_o * sin_w - sin_o * cos_w * cos_i
    r21 = sin_o * cos_w + cos_o * sin_w * cos_i
    r22 = -sin_o * sin_w + cos_o * cos_w * cos_i
    r31 = sin_w * sin_i
    r32 = cos_w * sin_i

    position = (
        r11 * perigee_x + r12 * perigee_y,
        r21 * perigee_x + r22 * perigee_y,
        r31 * perigee_x + r32 * perigee_y,
    )
    velocity = (
        r11 * perigee_vx + r12 * perigee_vy,
        r21 * perigee_vx + r22 * perigee_vy,
        r31 * perigee_vx + r32 * perigee_vy,
    )
    return position, velocity


def eci_to_coe(position, velocity, epoch: Epoch, mu: float = WGS84.mu) -> ClassicalElements:
    rx, ry, rz = position
    vx, vy, vz = velocity
    radius = math.sqrt(rx * rx + ry * ry + rz * rz)
    if radius <= 0.0:
        raise NonPhysicalOrbitError("zero position radius")
    speed2 = vx * vx + vy * vy + vz * vz
    hx, hy, hz = (
        ry * vz - rz * vy,
        rz * vx - rx * vz,
        rx * vy - ry * vx,
    )
    h_mag = math.sqrt(hx * hx + hy * hy + hz * hz)
    if h_mag <= 0.0:
        raise NonPhysicalOrbitError("degenerate orbit (zero angular momentum)")

    energy = 0.5 * speed2 - mu / radius
    if energy >= 0.0:
        raise NonPhysicalOrbitError(
            "state is not on a bound elliptic orbit", location="semi_major_axis"
        )
    a = -mu / (2.0 * energy)
    rv = rx * vx + ry * vy + rz * vz
    ex = (speed2 - mu / radius) * rx / mu - rv * vx / mu
    ey = (speed2 - mu / radius) * ry / mu - rv * vy / mu
    ez = (speed2 - mu / radius) * rz / mu - rv * vz / mu
    e_mag = math.sqrt(ex * ex + ey * ey + ez * ez)

    node_x, node_y = -hy, hx
    node_mag = math.sqrt(node_x * node_x + node_y * node_y)
    circular = e_mag < 1e-11
    equatorial = node_mag < 1e-11

    inclination = math.acos(max(-1.0, min(1.0, hz / h_mag)))
    raan = 0.0 if equatorial else math.atan2(node_y, node_x) % TWO_PI
    if circular:
        argument_of_periapsis = 0.0
    elif equatorial:
        argument_of_periapsis = math.atan2(ey, ex) % TWO_PI
    else:
        cos_argp = (node_x * ex + node_y * ey) / (node_mag * e_mag)
        sin_argp = (
            hx * node_y * ez
            - hy * node_x * ez
            + hz * (node_x * ey - node_y * ex)
        ) / (h_mag * node_mag * e_mag)
        argument_of_periapsis = math.atan2(sin_argp, cos_argp) % TWO_PI

    if circular and equatorial:
        true_anomaly = math.atan2(ry, rx) % TWO_PI
    elif circular:
        cos_nu = (node_x * rx + node_y * ry) / (node_mag * radius)
        sin_nu = (node_x * ry - node_y * rx) / (node_mag * radius)
        true_anomaly = math.atan2(sin_nu, cos_nu) % TWO_PI
    else:
        cos_nu = (ex * rx + ey * ry + ez * rz) / (e_mag * radius)
        true_anomaly = math.acos(max(-1.0, min(1.0, cos_nu)))
        if rv < 0.0:
            true_anomaly = TWO_PI - true_anomaly

    eccentric = 2.0 * math.atan2(
        math.sqrt(max(1.0 - e_mag, 0.0)) * math.sin(true_anomaly / 2.0),
        math.sqrt(1.0 + e_mag) * math.cos(true_anomaly / 2.0),
    )
    mean_anomaly = (eccentric - e_mag * math.sin(eccentric)) % TWO_PI

    return ClassicalElements(
        semi_major_axis=a,
        eccentricity=e_mag,
        inclination=inclination,
        raan=raan,
        argument_of_periapsis=argument_of_periapsis,
        mean_anomaly=mean_anomaly,
        epoch=epoch,
    )


def specific_energy(position, velocity, mu: float = WGS84.mu) -> float:
    vx, vy, vz = velocity
    speed2 = vx * vx + vy * vy + vz * vz
    radius = math.sqrt(sum(component * component for component in position))
    return 0.5 * speed2 - mu / radius


def angular_momentum(position, velocity):
    rx, ry, rz = position
    vx, vy, vz = velocity
    return (
        ry * vz - rz * vy,
        rz * vx - rx * vz,
        rx * vy - ry * vx,
    )
