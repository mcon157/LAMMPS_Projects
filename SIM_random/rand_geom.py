"""rand_geom.py - shared geometry for the random-orientation SIM figures.

Same maths as in_rand_W.step1 / in_rand_WC.step1: lattice (+ random origin
shift) -> rotate about the box centre by R = Rz(alpha) Rx(beta) Rz(gamma) ->
keep r <= R, z <= 0.  Surface normal in the crystal frame:
n = R^T z = (sin b sin g, sin b cos g, cos b).

Euler angles of surface s, first found of:
  1. rand_<MAT>_surfaces.dat from step 1            exact
  2. rand_angles_lammps.csv (LAMMPS random stream)   exact
  3. the LAMMPS python module                        exact
  4. numpy                                           NOT the runs' surfaces (flagged)
"""
import math
import os
from functools import lru_cache
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

A_WC, C_WC = 2.92949147531862, 2.82420965326672   # in_rand_WC.step1
A_W = 3.165                                         # in_rand_W.step1 (CHECK value at T)

COL_W, COL_C = "#8c8c8c", "#3b6fd6"


# ---------------------------------------------------------------- angles --
def _angles_from_u(u):
    return (360.0 * u[0], math.degrees(math.acos(2.0 * u[1] - 1.0)), 360.0 * u[2],
            0.999 * u[3], 0.999 * u[4], 0.999 * u[5])


@lru_cache(maxsize=None)
def _csv_table():
    p = os.path.join(HERE, "rand_angles_lammps.csv")
    if not os.path.isfile(p):
        return {}
    out = {}
    with open(p) as fh:
        for line in fh:
            if line.startswith(("#", "s")) or not line.strip():
                continue
            v = line.replace(",", " ").split()
            out[int(v[0])] = tuple(float(x) for x in v[1:7])
    return out


@lru_cache(maxsize=None)
def _surfaces_dat(path):
    out = {}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                v = line.split()
                if len(v) >= 11 and not line.startswith("#"):
                    # tag s alpha beta gamma nx ny nz ox oy oz ...
                    out[(v[0], int(v[1]))] = tuple(float(x) for x in (v[2], v[3], v[4], v[8], v[9], v[10]))
    return out


def get_angles(s, surfaces_dat=None, tag="r"):
    """(alpha, beta, gamma, ox, oy, oz), source-string."""
    d = _surfaces_dat(surfaces_dat)
    if (tag, s) in d:
        return d[(tag, s)], "surfaces.dat"
    t = _csv_table()
    if s in t:
        return t[s], "LAMMPS stream"
    try:
        import lammps
        l = lammps.lammps(cmdargs=["-log", "none", "-screen", "none"])
        l.command(f"variable seed0 equal 1000003*{s}+17")
        u = []
        for k in range(6):
            l.command(f"variable u{k} equal $(random(0.0,1.0,v_seed0))")
            u.append(l.extract_variable(f"u{k}"))
        l.close()
        return _angles_from_u(u), "LAMMPS stream"
    except Exception:
        rng = np.random.default_rng(1000003 * s + 17)
        return _angles_from_u(rng.random(6)), "numpy (not the runs' surfaces)"


# -------------------------------------------------------------- rotations --
def rot_z(deg):
    t = math.radians(deg); c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rot_x(deg):
    t = math.radians(deg); c, s = math.cos(t), math.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def rotation(alpha, beta, gamma):
    """R = Rz(alpha) Rx(beta) Rz(gamma): the deck's three displace_atoms rotate calls, gamma first."""
    return rot_z(alpha) @ rot_x(beta) @ rot_z(gamma)


def normal(beta, gamma):
    b, g = math.radians(beta), math.radians(gamma)
    return np.array([math.sin(b) * math.sin(g), math.sin(b) * math.cos(g), math.cos(b)])


# ---------------------------------------------------------------- lattices --
def lattice_W(lo, hi, origin=(0, 0, 0), a=A_W):
    """bcc W atoms inside the box [lo, hi]^3 (Angstrom); type 1."""
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    n0, n1 = int(math.floor(lo / a)) - 1, int(math.ceil(hi / a)) + 1
    ijk = np.array(list(product(range(n0, n1 + 1), repeat=3)), float)
    pts = ((ijk[:, None, :] + basis[None, :, :] + np.asarray(origin)) * a).reshape(-1, 3)
    keep = np.all((pts >= lo) & (pts <= hi), axis=1)
    return pts[keep], np.ones(keep.sum(), int)


def lattice_WC(lo, hi, origin=(0, 0, 0), a=A_WC, c=C_WC):
    """WC (orthorhombic 4-atom cell of in_rand_WC.step1) inside [lo, hi]^3; W = 1, C = 2."""
    cell = np.array([a, math.sqrt(3) * a, c])
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0, 2 / 3, 0.5], [0.5, 1 / 6, 0.5]])
    types = np.array([1, 1, 2, 2])
    rng = [range(int(math.floor(lo / L)) - 1, int(math.ceil(hi / L)) + 2) for L in cell]
    ijk = np.array(list(product(*rng)), float)
    frac = ijk[:, None, :] + basis[None, :, :] + np.asarray(origin)
    pts = (frac * cell).reshape(-1, 3)
    typ = np.tile(types, len(ijk))
    keep = np.all((pts >= lo) & (pts <= hi), axis=1)
    return pts[keep], typ[keep]


def lattice(mat, lo, hi, origin=(0, 0, 0)):
    return lattice_W(lo, hi, origin) if mat == "W" else lattice_WC(lo, hi, origin)


def hemisphere(mat, R, angles):
    """Atoms of the (unrelaxed) hemisphere in the LAB frame, plus their crystal-frame positions."""
    al, be, ga, ox, oy, oz = angles
    p0, t = lattice(mat, -R - 2, R + 2, (ox, oy, oz))
    keep0 = np.linalg.norm(p0, axis=1) <= R
    p0, t = p0[keep0], t[keep0]
    p = p0 @ rotation(al, be, ga).T
    keep = p[:, 2] <= 0.0
    return p[keep], t[keep], p0[keep]


# ----------------------------------------------------------- Miller labels --
def _ints(nmax):
    for h, k, l in product(range(-nmax, nmax + 1), repeat=3):
        if (h, k, l) != (0, 0, 0) and math.gcd(math.gcd(abs(h), abs(k)), abs(l)) == 1:
            yield h, k, l


def miller_cubic(n, nmax=6):
    """Closest (hkl) with |h|,|k|,|l| <= nmax to the unit normal n; returns (hkl), deviation (deg)."""
    best, bd = None, 1e9
    for hkl in _ints(nmax):
        v = np.array(hkl, float)
        d = math.degrees(math.acos(min(1.0, float(np.dot(v, n)) / np.linalg.norm(v))))
        # prefer the simplest index within 0.5 deg of the best
        if d < bd - 0.5 or (abs(d - bd) <= 0.5 and sum(map(abs, hkl)) < sum(map(abs, best))):
            best, bd = hkl, d
    return best, bd


def miller_hex(n, nmax=6, a=A_WC, c=C_WC):
    """Closest (h k i l) plane to the unit normal n (crystal frame of the WC deck:
    x = a[2-1-10], y = [01-10], z = c[0001]); returns (h,k,i,l), deviation (deg)."""
    a1 = np.array([a, 0, 0]); a2 = np.array([-a / 2, math.sqrt(3) * a / 2, 0]); a3 = np.array([0, 0, c])
    V = np.dot(a1, np.cross(a2, a3))
    b1, b2, b3 = np.cross(a2, a3) / V, np.cross(a3, a1) / V, np.cross(a1, a2) / V
    best, bd = None, 1e9
    for h, k, l in _ints(nmax):
        g = h * b1 + k * b2 + l * b3
        d = math.degrees(math.acos(min(1.0, float(np.dot(g, n)) / np.linalg.norm(g))))
        if d < bd - 0.5 or (abs(d - bd) <= 0.5 and abs(h) + abs(k) + abs(l) < sum(map(abs, best[:2])) + abs(best[3])):
            best, bd = (h, k, -(h + k), l), d
    return best, bd


def miller_tex(idx):
    """Mathtext string with overbars for negative indices, e.g. (1 0 \\bar{1} 0)."""
    parts = [r"\bar{%d}" % -i if i < 0 else "%d" % i for i in idx]
    return r"$(" + r"\,".join(parts) + r")$"


def plane_label(mat, n, nmax=6):
    if mat == "W":
        idx, dev = miller_cubic(n, nmax)
    else:
        idx, dev = miller_hex(n, nmax)
    return miller_tex(idx), dev


# -------------------------------------------------------- plane in a cube --
def plane_cube_polygon(n, d, half):
    """Vertices of the polygon where the plane n.x = d cuts the cube [-half, half]^3, ordered."""
    corners = np.array(list(product((-half, half), repeat=3)), float)
    edges = [(i, j) for i in range(8) for j in range(i + 1, 8)
             if np.sum(corners[i] != corners[j]) == 1]
    pts = []
    for i, j in edges:
        fi, fj = np.dot(n, corners[i]) - d, np.dot(n, corners[j]) - d
        if fi == fj:
            continue
        t = fi / (fi - fj)
        if 0 <= t <= 1:
            pts.append(corners[i] + t * (corners[j] - corners[i]))
    if len(pts) < 3:
        return np.empty((0, 3))
    pts = np.unique(np.round(pts, 9), axis=0)
    c = pts.mean(axis=0)
    u = pts[0] - c; u /= np.linalg.norm(u)
    w = np.cross(n, u)
    ang = np.arctan2((pts - c) @ w, (pts - c) @ u)
    return pts[np.argsort(ang)]
