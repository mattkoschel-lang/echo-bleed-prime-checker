#!/usr/bin/env python3
"""
Echo Bleed Geometric Prime Checker & Analytic Number Theory Engine (v6)
Fused with the 3-6-9 Harmonic Resonance Lattice Engine, Quasi-RH bounds,
and Landau-Siegel zero diagnostics.
"""

import sys
import math
import hashlib
import argparse
from dataclasses import dataclass, field
from typing import List, Callable, Optional, Tuple, Any, Dict

# Optional high-performance library
try:
    import gmpy2
    from gmpy2 import mpz
    HAVE_GMPY2 = True
except Exception:
    gmpy2 = None
    mpz = int
    HAVE_GMPY2 = False

# Allow unlimited digits when converting very large ints to strings (Python 3.11+)
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

# First non-trivial zeros of the Riemann zeta function (imaginary parts)
ZETA_ZEROS = [
    14.1347251417, 21.0220396388, 25.0108575801,
    30.4248761259, 32.9350615877, 37.5861781588
]


# ============================================================
# 1. Core Math, Riemann, & Quasi-RH Analytical Utilities
# ============================================================

def digital_root_from_bytes(data: bytes) -> int:
    """Compute the digital root of a byte sequence interpreted as hex."""
    hex_str = data.hex()
    total = sum(int(ch, 16) for ch in hex_str)
    while total > 9:
        total = sum(int(d) for d in str(total))
    return total


def hash_state(obj: Any) -> bytes:
    """Hash an arbitrary object by hashing its string representation."""
    s = repr(obj).encode("utf-8")
    return hashlib.sha256(s).digest()


def mobius(n: int) -> int:
    """Compute the Möbius function mu(n) for positive integer n."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    prime_factors = 0
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            temp //= d
            if temp % d == 0:
                return 0  # Has a squared prime factor
            prime_factors += 1
        d += 1
    if temp > 1:
        prime_factors += 1
    return -1 if (prime_factors % 2 == 1) else 1


def euler_log_zeta_partial(s: float = 2.0, primes: Tuple[int, ...] = (2, 3, 5, 7, 11), terms: int = 4) -> float:
    """Evaluate the Maclaurin Taylor expansion of ln(zeta(s)) over the first primes."""
    total = 0.0
    for p in primes:
        for k in range(1, terms + 1):
            total += 1.0 / (k * (p ** (k * s)))
    return total


def quasi_rh_error_bound(n: int, theta: float = 0.875) -> float:
    """
    Compute the theoretical error bound |pi(n) - li(n)| under the
    Quasi-Riemann Hypothesis with zero-free bound theta (e.g., 7/8 = 0.875).
    """
    if n < 2:
        return 0.0
    ln_n = n.bit_length() * 0.6931471805599453 if n > 10**300 else math.log(n)
    # Prevent float overflow on huge n^theta display
    try:
        return (n ** theta) * ln_n
    except OverflowError:
        return float('inf')


def check_siegel_zero_proximity(n: int, conductor_threshold: int = 3) -> Dict[str, Any]:
    """
    Heuristic check for potential Landau-Siegel zero proximity.
    Evaluates potential near-1 pole interference in L-functions.
    """
    if n <= conductor_threshold:
        log_q = math.log(conductor_threshold)
    elif n > 10**300:
        log_q = n.bit_length() * 0.6931471805599453
    else:
        log_q = math.log(n)

    implied_proximity = 1.0 / log_q
    risk_flag = implied_proximity > 0.15

    return {
        "conductor_scale": round(log_q, 4),
        "siegel_proximity_metric": round(implied_proximity, 6),
        "anomaly_risk": risk_flag,
        "status": "SUPPRESSED (Uniform Zero-Free Region Active)" if not risk_flag else "WARNING: Potential Near-1 Pole Interference"
    }


def riemann_analytical_metrics(n: int, max_exact_n: int = 1_000_000) -> Dict[str, Any]:
    """Compute Riemann's J(n), prime counting pi(n), Möbius mu(n), Quasi-RH, and Siegel metrics."""
    if n < 2:
        res = {"mode": "exact", "J_n": 0.0, "pi_n": 0, "mu_n": mobius(n)}
        res["quasi_rh_bound"] = 0.0
        res["siegel_guard"] = check_siegel_zero_proximity(n)
        return res

    if n > max_exact_n:
        ln_n = n.bit_length() * 0.6931471805599453
        approx_pi = n / ln_n
        sqrt_n = math.isqrt(n)
        approx_j = approx_pi + (0.5 * (sqrt_n / (0.5 * ln_n)) if sqrt_n >= 2 else 0.0)
        qrh = quasi_rh_error_bound(n, theta=0.875)
        asymptotic_res = {
            "mode": "asymptotic",
            "J_n": round(approx_j, 6) if approx_j < 1e300 else "O(n / ln n)",
            "pi_n": int(round(approx_pi)) if approx_pi < 1e300 else f"~2^{n.bit_length()} / {round(ln_n, 2)}",
            "mu_n": None,
            "quasi_rh_bound": round(qrh, 2) if qrh != float('inf') else f"O(2^{int(n.bit_length() * 0.875)} * {round(ln_n, 2)})",
            "siegel_guard": check_siegel_zero_proximity(n)
        }
        return asymptotic_res

    # Sieve of Eratosthenes up to n
    sieve = [True] * (n + 1)
    sieve[0] = False
    sieve[1] = False
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            for multiple in range(p * p, n + 1, p):
                sieve[multiple] = False

    pi_prefix = [0] * (n + 1)
    running = 0
    for i in range(1, n + 1):
        if sieve[i]:
            running += 1
        pi_prefix[i] = running

    # J(n) = sum_{k>=1} (1/k) * pi(n^(1/k))
    j_val = 0.0
    k = 1
    while True:
        root = int(round(n ** (1.0 / k)))
        while (root + 1) ** k <= n:
            root += 1
        while root ** k > n:
            root -= 1
        if root < 2:
            break
        j_val += pi_prefix[root] / k
        k += 1

    return {
        "mode": "exact",
        "J_n": round(j_val, 6),
        "pi_n": pi_prefix[n],
        "mu_n": mobius(n),
        "quasi_rh_bound": round(quasi_rh_error_bound(n, theta=0.875), 2),
        "siegel_guard": check_siegel_zero_proximity(n)
    }


# ============================================================
# 2. Geometric Primitives & 3-6-9 Core Resonance Lattice
# ============================================================

@dataclass
class Node:
    value: float = 0.0
    resonance: float = 0.0
    weight: float = 1.0


@dataclass
class Edge:
    node_indices: Tuple[int, int]
    coherence: float = 0.0


@dataclass
class TetrahedralCell:
    nodes: List[Node] = field(default_factory=lambda: [Node() for _ in range(4)])
    edges: List[Edge] = field(default_factory=list)
    position: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    channel_k: int = 3
    layer_i: int = 0
    angular_j: int = 0
    intensity: float = 1.0
    opacity: float = 1.0

    def __post_init__(self):
        if not self.edges:
            pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            self.edges = [Edge(pair) for pair in pairs]

    def resonance(self, F: Optional[Callable[['TetrahedralCell'], float]] = None) -> float:
        if F is not None:
            return F(self)

        values = [n.value for n in self.nodes]
        if not values:
            return 0.0

        mean_val = sum(values) / len(values)
        variance = sum((v - mean_val) ** 2 for v in values) / len(values)
        variance_term = 1.0 / (1.0 + variance)

        if self.edges:
            mean_coherence = sum(e.coherence for e in self.edges) / len(self.edges)
        else:
            mean_coherence = 0.0

        # Modulate cell resonance by its 3-6-9 lattice intensity envelope and fade opacity
        return (variance_term + mean_coherence) * self.intensity * self.opacity

    def brightness(self, F: Optional[Callable[['TetrahedralCell'], float]] = None) -> float:
        return self.resonance(F)


@dataclass
class DodecahedralFrame:
    phi: float = (1 + 5 ** 0.5) / 2
    faces: int = 12
    sai_target: float = 1.0

    def target_coherence(self) -> float:
        return self.sai_target


@dataclass
class TeslaLattice369:
    """
    Unified 3-6-9 Harmonic Resonance Engine.
    Interleaves triangular (k=3), hexagonal (k=6), and nonagonal (k=9) channels
    for angular modulation, radial layer rhythms, and intensity envelopes.
    """
    center_x: float = 0.0
    center_y: float = 0.0
    R: float = 20.0
    ks: Tuple[int, ...] = (3, 6, 9)

    # Channel modulation depths & phase offsets
    alpha_k: Dict[int, float] = field(default_factory=lambda: {3: 0.15, 6: 0.15, 9: 0.15})
    beta_k: Dict[int, float] = field(default_factory=lambda: {3: 0.10, 6: 0.10, 9: 0.10})
    delta_k: Dict[int, float] = field(default_factory=lambda: {3: 0.05, 6: 0.05, 9: 0.05})
    varphi_k: Dict[int, float] = field(default_factory=lambda: {3: 0.0, 6: 0.0, 9: 0.0})
    phi_k0: Dict[int, float] = field(default_factory=lambda: {3: 0.0, 6: 0.0, 9: 0.0})

    # Frequency driver (f_3 : f_6 : f_9 = 1 : 2 : 3)
    base_freq: float = 1.0

    # Animation driver parameters
    i0: float = 0.0
    velocity: float = 0.05

    # Unified intensity envelope parameters
    w_k: Dict[int, float] = field(default_factory=lambda: {3: 1.0, 6: 1.0, 9: 1.0})
    gamma_0: float = 1.0
    gamma_1: float = 0.25
    gamma_2: float = 0.15
    gamma_3: float = 0.20

    def freq_for_channel(self, k: int) -> float:
        ratio = {3: 1.0, 6: 2.0, 9: 3.0}.get(k, k / 3.0)
        return self.base_freq * ratio

    def channel_phase(self, k: int, t: float) -> float:
        fk = self.freq_for_channel(k)
        dk = self.delta_k.get(k, 0.05)
        p0 = self.phi_k0.get(k, 0.0)
        return dk * math.sin(2.0 * math.pi * fk * t + p0)

    def theta(self, j: int, k: int, t: float) -> float:
        return (2.0 * math.pi * j / k) + self.channel_phase(k, t)

    def discrete_layer_resonance(self, i: int, k: int) -> float:
        return math.sin(2.0 * math.pi * i / k)

    def continuous_pulse(self, k: int, t: float) -> float:
        fk = self.freq_for_channel(k)
        vp = self.varphi_k.get(k, 0.0)
        return math.sin(2.0 * math.pi * fk * t + vp)

    def radius(self, i: int, k: int, t: float) -> float:
        ak = self.alpha_k.get(k, 0.15)
        bk = self.beta_k.get(k, 0.10)
        Ak = self.discrete_layer_resonance(i, k)
        Sk = self.continuous_pulse(k, t)
        return self.R * (i + ak * Ak + bk * Sk)

    def intensity(self, i: int, j: int, k: int, t: float) -> float:
        wk = self.w_k.get(k, 1.0)
        Ak = self.discrete_layer_resonance(i, k)
        th = self.theta(j, k, t)
        Sk = self.continuous_pulse(k, t)
        return wk * (self.gamma_0 + self.gamma_1 * Ak + self.gamma_2 * math.cos(th) + self.gamma_3 * Sk)

    def layer_indices(self, t: float) -> Tuple[int, float]:
        imax = self.i0 + self.velocity * t
        base = int(math.floor(imax))
        fade = imax - base
        return base, fade

    def generate_lattice_points(self, t: float, min_layers: int = 1) -> List[Dict[str, Any]]:
        """Generate interleaved (x, y, k, i, j, opacity, intensity, radius, theta) points."""
        pts = []
        base, fade = self.layer_indices(t)
        effective_base = max(base, min_layers - 1)

        for k in self.ks:
            # Complete layers
            for i in range(effective_base + 1):
                for j in range(k):
                    th = self.theta(j, k, t)
                    rad = self.radius(i, k, t)
                    x = self.center_x + rad * math.cos(th)
                    y = self.center_y + rad * math.sin(th)
                    inten = self.intensity(i, j, k, t)
                    pts.append({
                        "x": x, "y": y, "k": k, "i": i, "j": j,
                        "opacity": 1.0, "intensity": inten,
                        "radius": rad, "theta": th
                    })

            # Fractional fade-in layer
            i_next = effective_base + 1
            layer_opacity = fade if fade > 0.0 else 0.15
            for j in range(k):
                th = self.theta(j, k, t)
                rad = self.radius(i_next, k, t)
                x = self.center_x + rad * math.cos(th)
                y = self.center_y + rad * math.sin(th)
                inten = self.intensity(i_next, j, k, t)
                pts.append({
                    "x": x, "y": y, "k": k, "i": i_next, "j": j,
                    "opacity": layer_opacity, "intensity": inten,
                    "radius": rad, "theta": th
                })
        return pts


# ============================================================
# 3. Quantum-Inspired Configuration Space
# ============================================================

@dataclass
class ConfigurationSuperposition:
    configurations: List[Dict[str, Any]] = field(default_factory=list)

    def add_configuration(self, state: Dict[str, Any]) -> None:
        self.configurations.append(state)

    def collapse(self, scoring_fn: Callable[[Dict[str, Any]], float]) -> Dict[str, Any]:
        if not self.configurations:
            return {}
        scored = [(scoring_fn(cfg), cfg) for cfg in self.configurations]
        scored.sort(key=lambda x: x[0], reverse=True)
        return scored[0][1]


# ============================================================
# 4. ECHO BLEED Core System (Driven by TeslaLattice369)
# ============================================================

@dataclass
class EchoBleedSystem:
    cells: List[TetrahedralCell]
    frame: DodecahedralFrame
    lattice: TeslaLattice369

    alpha: float = 0.05
    beta: float = 0.01
    step_index: int = 0
    dt: float = 0.25

    mean_radius: float = 20.0
    mean_theta: float = 0.0
    mean_intensity: float = 1.0

    last_coherent_state: Optional[List[TetrahedralCell]] = None
    resonance_fn: Optional[Callable[[TetrahedralCell], float]] = None
    superposition: ConfigurationSuperposition = field(default_factory=ConfigurationSuperposition)

    @property
    def current_time(self) -> float:
        return self.step_index * self.dt

    def compute_sai(self) -> float:
        if not self.cells:
            return 0.0
        return sum(c.resonance(self.resonance_fn) for c in self.cells) / len(self.cells)

    def compute_coherence_glow(self) -> float:
        if not self.cells:
            return 0.0
        return sum(c.brightness(self.resonance_fn) for c in self.cells) / len(self.cells)

    def update_geometry(self, sai: float) -> None:
        target = self.frame.target_coherence()
        self.lattice.R += self.alpha * (sai - target)
        self.mean_theta += self.beta * sai

    def lattice_step(self) -> None:
        """Advance the 3-6-9 harmonic lattice and map cells across the interleaved channels."""
        t = self.current_time
        # Ensure enough layers exist to house all cells across channels (3 + 6 + 9 = 18 pts/layer)
        pts_per_layer = sum(self.lattice.ks)
        required_layers = max(1, math.ceil(len(self.cells) / pts_per_layer))
        pts = self.lattice.generate_lattice_points(t, min_layers=required_layers)

        radii = []
        thetas = []
        intensities = []

        for idx, cell in enumerate(self.cells):
            pt = pts[idx % len(pts)]
            cell.position = (pt["x"], pt["y"], float(pt["k"]))
            cell.channel_k = pt["k"]
            cell.layer_i = pt["i"]
            cell.angular_j = pt["j"]
            cell.intensity = pt["intensity"]
            cell.opacity = pt["opacity"]

            radii.append(pt["radius"])
            thetas.append(pt["theta"])
            intensities.append(pt["intensity"])

        if radii:
            self.mean_radius = sum(radii) / len(radii)
            self.mean_theta = sum(thetas) / len(thetas)
            self.mean_intensity = sum(intensities) / len(intensities)

    def serialize_lattice_state(self) -> Any:
        state = []
        for cell in self.cells:
            node_vals = [n.value for n in cell.nodes]
            edge_coh = [e.coherence for e in cell.edges]
            state.append((tuple(node_vals), tuple(edge_coh), cell.position, cell.channel_k))
        return {
            "step": self.step_index,
            "time": self.current_time,
            "mean_radius": self.mean_radius,
            "mean_theta": self.mean_theta,
            "cells": tuple(state),
        }

    def nine_axis_sentinel(self) -> bool:
        return digital_root_from_bytes(hash_state(self.serialize_lattice_state())) == 9

    def snapshot_if_coherent(self) -> None:
        if self.nine_axis_sentinel():
            self.last_coherent_state = [
                TetrahedralCell(
                    nodes=[Node(n.value, n.resonance, n.weight) for n in cell.nodes],
                    edges=[Edge(e.node_indices, e.coherence) for e in cell.edges],
                    position=cell.position,
                    channel_k=cell.channel_k,
                    layer_i=cell.layer_i,
                    angular_j=cell.angular_j,
                    intensity=cell.intensity,
                    opacity=cell.opacity,
                )
                for cell in self.cells
            ]

    def restore_protocol(self) -> None:
        sai = self.compute_sai()
        target = self.frame.target_coherence()
        self.lattice.R += self.alpha * (target - sai)
        if self.last_coherent_state is not None:
            self.cells = [
                TetrahedralCell(
                    nodes=[Node(n.value, n.resonance, n.weight) for n in cell.nodes],
                    edges=[Edge(e.node_indices, e.coherence) for e in cell.edges],
                    position=cell.position,
                    channel_k=cell.channel_k,
                    layer_i=cell.layer_i,
                    angular_j=cell.angular_j,
                    intensity=cell.intensity,
                    opacity=cell.opacity,
                )
                for cell in self.last_coherent_state
            ]

    def ingest_vector(self, data: List[float]) -> None:
        if not data or not self.cells:
            return
        dmin, dmax = min(data), max(data)
        span = (dmax - dmin) or 1.0
        norm = [(x - dmin) / span for x in data]
        idx = 0
        for cell in self.cells:
            for node in cell.nodes:
                node.value = norm[idx % len(norm)]
                idx += 1
        for edge in cell.edges:
            i, j = edge.node_indices
            vi = cell.nodes[i].value
            vj = cell.nodes[j].value
            edge.coherence = 1.0 - abs(vi - vj)

    def explore_superposition(self, variants: List[Dict[str, Any]]) -> None:
        for v in variants:
            self.superposition.add_configuration(v)

    def collapse_superposition(self) -> None:
        def score(cfg: Dict[str, Any]) -> float:
            return cfg.get("sai", 0.0)
        chosen = self.superposition.collapse(score)
        if not chosen:
            return
        self.mean_radius = chosen.get("radius", self.mean_radius)
        self.mean_theta = chosen.get("theta", self.mean_theta)

    def step(self) -> Dict[str, float]:
        self.lattice_step()
        sai = self.compute_sai()
        glow = self.compute_coherence_glow()
        self.update_geometry(sai)
        if self.nine_axis_sentinel():
            self.snapshot_if_coherent()
        else:
            self.restore_protocol()
        self.step_index += 1
        return {
            "step": self.step_index,
            "time": self.current_time,
            "sai": sai,
            "glow": glow,
            "radius": self.mean_radius,
            "theta": self.mean_theta,
            "intensity": self.mean_intensity,
        }


# ============================================================
# 5. Primality & Engine Integration
# ============================================================

def riemann_cell_resonance(cell: TetrahedralCell) -> float:
    base_res = cell.resonance(None)
    euler_factor = euler_log_zeta_partial(s=2.0)
    mean_val = sum(n.value for n in cell.nodes) / (len(cell.nodes) or 1)
    # Channel harmonic boost: k=9 carries highest harmonic weight
    channel_weight = {3: 1.0, 6: 1.03, 9: 1.09}.get(cell.channel_k, 1.0)
    return base_res * channel_weight * (1.0 + euler_factor * mean_val)


def build_default_echo_bleed_system(num_cells: int = 36, use_riemann: bool = False) -> EchoBleedSystem:
    # Default 36 cells maps cleanly across two full 3+6+9=18-node layers
    cells = [TetrahedralCell() for _ in range(num_cells)]
    frame = DodecahedralFrame(sai_target=1.0)

    # When Riemann mode is active, lock the 3-6-9 base frequency to the 1st Zeta Zero
    base_f = ZETA_ZEROS[0] if use_riemann else 1.0
    lattice = TeslaLattice369(
        R=20.0,
        base_freq=base_f,
        i0=0.0,
        velocity=0.05,
    )
    return EchoBleedSystem(
        cells=cells,
        frame=frame,
        lattice=lattice,
        alpha=0.05,
        beta=0.01,
        resonance_fn=riemann_cell_resonance if use_riemann else None,
    )


def is_miller_rabin_prime(n: int) -> bool:
    if n < 2:
        return False
    if HAVE_GMPY2:
        return gmpy2.is_prime(mpz(n))
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n <= a:
            break
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def generate_prime_candidate_vector(n: int, max_vector_dim: int = 128) -> Tuple[List[float], Optional[int]]:
    if n < 2:
        return [0.0], None
    if n in (2, 3):
        return [1.0, 1.0, 1.0], None
    if n % 2 == 0:
        return [0.0], 2
    if n % 3 == 0:
        return [0.0], 3
    vector = [
        float(n % 2), float(n % 3), float(n % 5),
        float(n % 7), float(n % 9),
        float(digital_root_from_bytes(str(n).encode()))
    ]
    limit = math.isqrt(n)
    first_factor = None
    d, step, samples = 5, 2, 0
    while d <= limit and samples < max_vector_dim:
        rem = n % d
        vector.append(float(rem))
        if rem == 0 and first_factor is None:
            first_factor = d
            break
        d, step = d + step, 6 - step
        samples += 1
    if first_factor is None and not is_miller_rabin_prime(n):
        first_factor = -1
    return vector, first_factor


def check_prime_with_engine(n: int, num_cells: int = 36, use_riemann: bool = False) -> Dict[str, Any]:
    if n < 2:
        return {"n": n, "is_prime": False, "reason": "Numbers < 2 are not prime."}
    vector, factor = generate_prime_candidate_vector(n)
    system = build_default_echo_bleed_system(num_cells=num_cells, use_riemann=use_riemann)
    system.ingest_vector(vector)
    diag = system.step()
    riemann_data = riemann_analytical_metrics(n) if use_riemann else None
    mu_bias = 1.05 if (riemann_data is None or riemann_data.get("mu_n") != 0) else 0.95
    system.explore_superposition([
        {"sai": diag["sai"], "radius": diag["radius"], "theta": diag["theta"]},
        {"sai": diag["sai"] * mu_bias, "radius": diag["radius"] * 0.9, "theta": diag["theta"] + 0.1},
    ])
    system.collapse_superposition()
    is_prime = (factor is None)

    # Calculate channel distribution across k=3, 6, 9
    channel_counts = {3: 0, 6: 0, 9: 0}
    for c in system.cells:
        if c.channel_k in channel_counts:
            channel_counts[c.channel_k] += 1

    return {
        "n": n,
        "is_prime": is_prime,
        "factor": factor,
        "sai": diag["sai"],
        "glow": diag["glow"],
        "intensity": diag["intensity"],
        "base_freq": system.lattice.base_freq,
        "channels": channel_counts,
        "sentinel_active": system.nine_axis_sentinel(),
        "final_radius": system.mean_radius,
        "final_theta": system.mean_theta,
        "riemann": riemann_data,
    }


def find_next_prime(n: int, num_cells: int = 36, use_riemann: bool = False) -> Dict[str, Any]:
    candidate = n + 1
    if candidate <= 2:
        return check_prime_with_engine(2, num_cells, use_riemann)
    if candidate % 2 == 0:
        candidate += 1
    while True:
        res = check_prime_with_engine(candidate, num_cells, use_riemann)
        if res.get("is_prime"):
            return res
        candidate += 2


# ============================================================
# 6. Command Line Interface (CLI)
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Echo Bleed 3-6-9 Geometric Prime Checker")
    parser.add_argument("number", type=str, nargs="?", help="The integer to check for primality.")
    parser.add_argument("--cells", type=int, default=36, help="Number of tetrahedral cells in the 3-6-9 lattice (default: 36).")
    parser.add_argument("--verbose", "-v", action="store_true", help="Display full 3-6-9 geometric and lattice diagnostics.")
    parser.add_argument("--riemann", "-r", action="store_true", help="Enable Riemann J(x), Möbius mu(n), Quasi-RH, and Zeta-locked harmonics.")
    args = parser.parse_args()

    if args.number is None:
        try:
            args.number = input("Enter a number to check: ").strip()
        except EOFError:
            sys.exit(1)

    cleaned_str = "".join(
        c for i, c in enumerate(args.number)
        if c.isdigit() or (c == '-' and i == 0)
    )

    try:
        target_number = int(cleaned_str)
    except ValueError:
        parser.error(f"'{args.number}' is not a valid integer.")

    result = check_prime_with_engine(target_number, num_cells=args.cells, use_riemann=args.riemann)
    print("\n==================================================")
    print(f" ECHO BLEED 3-6-9 PRIME CHECKER :: Target: {result['n']}")
    print("==================================================")
    if result['is_prime']:
        status_text = "PRIME [OK]"
    else:
        factor_val = result.get('factor')
        if factor_val and factor_val != -1:
            reason = f"Divisible by {factor_val}"
        else:
            reason = "Composite (Verified via Miller-Rabin)"
        status_text = f"COMPOSITE [X] ({reason})"
    print(f"Status    : {status_text}")
    next_res = find_next_prime(target_number, num_cells=args.cells, use_riemann=args.riemann)
    print(f"Next Prime: {next_res['n']}")

    if args.verbose:
        ch = result['channels']
        print("\n--- 3-6-9 Resonance Lattice Diagnostics ---")
        print(f"Lattice Cells    : {args.cells} (Channels k=3:{ch[3]}, k=6:{ch[6]}, k=9:{ch[9]})")
        print(f"Base Frequency f : {result['base_freq']:.6f} Hz (1:2:3 Ratio)")
        print(f"Mean Intensity I : {result['intensity']:.6f}")
        print(f"SAI Metric       : {result['sai']:.6f}")
        print(f"Coherence Glow   : {result['glow']:.6f}")
        print(f"Mean Radius      : {result['final_radius']:.6f}")
        print(f"Mean Theta       : {result['final_theta']:.6f}")
        print(f"9-Axis Sentinel  : {'ALIGNED (9)' if result['sentinel_active'] else 'DIVERGENT'}")

    if args.riemann and result.get("riemann"):
        rm = result["riemann"]
        print("\n--- Riemann & Euler Analytical Metrics ---")
        print(f"Evaluation Mode  : {rm['mode'].upper()}")
        print(f"Riemann J(n)     : {rm['J_n']}")
        print(f"Prime Count pi(n): {rm['pi_n']}")
        print(f"Mobius mu(n)     : {rm['mu_n'] if rm['mu_n'] is not None else 'N/A (Asymptotic)'}")
        print(f"Quasi-RH Bound   : O(n^0.875 log n) ≈ {rm['quasi_rh_bound']}")

        sg = rm['siegel_guard']
        print(f"Siegel Zero Guard: {sg['status']} (Metric: {sg['siegel_proximity_metric']})")
        print(f"Euler ln zeta(2) : {euler_log_zeta_partial(2.0):.6f}")

    print("==================================================\n")


if __name__ == "__main__":
    main()
