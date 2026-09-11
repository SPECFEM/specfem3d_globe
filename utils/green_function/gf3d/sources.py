"""Seismic sources, and how to read them from specfem's own files.

Parsing happens here, in Python, and not in the library. That is deliberate:
reading a CMTSOLUTION in Fortran means ``get_cmt()``, which ends the process
on malformed input, and a ``stop`` inside a shared object takes the
interpreter with it. So the C API takes numbers, and turning a file into
numbers is this module's job -- where a bad file is a ``ValueError`` with a
line number.

Both readers follow the format the solver reads: one ``name: value`` per
line, in a fixed order, with the value being whatever follows the last
colon. Fortran ``d`` exponents (``1.0d15``) are accepted.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from ._lib import CSource, GF_SRC_CMT, GF_SRC_FORCE

__all__ = ["CMTSource", "ForceSource"]


def _to_float(text: str, where: str) -> float:
    """A float, accepting Fortran's d exponent."""
    cleaned = text.strip().replace("d", "e").replace("D", "E")
    try:
        return float(cleaned)
    except ValueError:
        raise ValueError(f"{where}: expected a number, got {text.strip()!r}") from None


def _value(line: str, where: str) -> str:
    """Everything after the last colon."""
    if ":" not in line:
        raise ValueError(f"{where}: expected 'name: value', got {line.strip()!r}")
    return line.rsplit(":", 1)[1]


def _lines(path, nexpected: int, kind: str) -> list[str]:
    path = Path(path)
    try:
        raw = path.read_text().splitlines()
    except OSError as exc:
        raise ValueError(f"cannot read {kind} {path}: {exc}") from None
    if len(raw) < nexpected:
        raise ValueError(
            f"{path}: a {kind} has {nexpected} lines, this one has {len(raw)}"
        )
    return raw


@dataclass
class CMTSource:
    """A moment-tensor (centroid) source.

    Field names and units are the CMTSOLUTION's own: degrees, kilometres,
    seconds, and dyne-cm for the moment tensor. ``hdur`` is the *triangle*
    half duration as the file writes it; the Gaussian width specfem actually
    uses is ``hdur/1.628``, and the library applies that itself.

    ``time_shift`` is metadata. The traces returned by this package have
    ``t = 0`` at the centroid time, so the shift belongs in an absolute
    start time and never in the trace.
    """

    latitude: float = 0.0
    longitude: float = 0.0
    depth: float = 0.0
    """Centroid depth, km."""
    hdur: float = 0.0
    time_shift: float = 0.0
    Mrr: float = 0.0
    Mtt: float = 0.0
    Mpp: float = 0.0
    Mrt: float = 0.0
    Mrp: float = 0.0
    Mtp: float = 0.0
    origin_time: datetime | None = None
    """The PDE line's time, if the file carried one."""
    event_name: str = ""
    region: str = ""

    @property
    def tensor(self) -> tuple:
        """``(Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)`` in dyne-cm, the order the
        moment-tensor partials come back in."""
        return (self.Mrr, self.Mtt, self.Mpp, self.Mrt, self.Mrp, self.Mtp)

    @property
    def centroid_time(self) -> datetime | None:
        """``origin_time + time_shift``: the instant ``t = 0`` refers to."""
        if self.origin_time is None:
            return None
        from datetime import timedelta

        return self.origin_time + timedelta(seconds=self.time_shift)

    @classmethod
    def read(cls, path) -> "CMTSource":
        """Read a 13-line CMTSOLUTION."""
        raw = _lines(path, 13, "CMTSOLUTION")
        src = cls()

        # the PDE header: fixed columns, not name: value
        pde = raw[0].split()
        if len(pde) >= 7:
            try:
                year, month, day = int(pde[1]), int(pde[2]), int(pde[3])
                hour, minute = int(pde[4]), int(pde[5])
                second = float(pde[6])
                whole = int(second)
                src.origin_time = datetime(
                    year, month, day, hour, minute, whole,
                    int(round((second - whole) * 1e6)), tzinfo=timezone.utc,
                )
            except (ValueError, IndexError):
                src.origin_time = None  # a header we do not recognise is not fatal
            src.region = " ".join(pde[12:]) if len(pde) > 12 else ""

        src.event_name = _value(raw[1], f"{path}:2").strip()
        src.time_shift = _to_float(_value(raw[2], f"{path}:3"), f"{path}:3")
        src.hdur = _to_float(_value(raw[3], f"{path}:4"), f"{path}:4")
        src.latitude = _to_float(_value(raw[4], f"{path}:5"), f"{path}:5")
        src.longitude = _to_float(_value(raw[5], f"{path}:6"), f"{path}:6")
        src.depth = _to_float(_value(raw[6], f"{path}:7"), f"{path}:7")
        for i, name in enumerate(("Mrr", "Mtt", "Mpp", "Mrt", "Mrp", "Mtp")):
            setattr(src, name, _to_float(_value(raw[7 + i], f"{path}:{8+i}"), f"{path}:{8+i}"))

        return src

    def _to_struct(self) -> CSource:
        c = CSource()
        c.source_type = GF_SRC_CMT
        c.force_stf = 0
        c.latitude = self.latitude
        c.longitude = self.longitude
        c.depth_km = self.depth
        c.hdur = self.hdur
        c.time_shift = self.time_shift
        for i, v in enumerate(self.tensor):
            c.moment[i] = v
        c.force_factor = 0.0
        for i in range(3):
            c.force_dir[i] = 0.0
        return c


@dataclass
class ForceSource:
    """A point-force source.

    ``f0`` is the FORCESOLUTION's own field: a dominant frequency for a
    Ricker (``stf = 1``), a width otherwise. ``factor`` is in Newtons and
    ``direction`` is ``(E, N, Z_up)`` of arbitrary length -- the magnitude
    lives entirely in ``factor``.
    """

    latitude: float = 0.0
    longitude: float = 0.0
    depth: float = 0.0
    """Source depth, km."""
    f0: float = 0.0
    time_shift: float = 0.0
    stf: int = 0
    """0 Gaussian, 1 Ricker, 2 step, 3 monochromatic, 4 Gaussian (Meschede)."""
    factor: float = 1.0e15
    """Newtons."""
    direction: tuple = (0.0, 0.0, 1.0)
    """(E, N, Z up)."""
    label: str = ""

    @classmethod
    def read(cls, path) -> "ForceSource":
        """Read an 11-line FORCESOLUTION.

        Accepts both the globe's spelling (``f0:``, ``comp dir vect source
        Z:``) and the Cartesian package's (``half duration:``, ``...
        Z_UP:``): every line is read by position, and only the value after
        the last colon is used, so the two differ in nothing that matters
        here.
        """
        raw = _lines(path, 11, "FORCESOLUTION")
        src = cls()

        first = raw[0].split()
        src.label = first[-1] if len(first) > 1 else ""

        src.time_shift = _to_float(_value(raw[1], f"{path}:2"), f"{path}:2")
        src.f0 = _to_float(_value(raw[2], f"{path}:3"), f"{path}:3")
        src.latitude = _to_float(_value(raw[3], f"{path}:4"), f"{path}:4")
        src.longitude = _to_float(_value(raw[4], f"{path}:5"), f"{path}:5")
        src.depth = _to_float(_value(raw[5], f"{path}:6"), f"{path}:6")
        src.stf = int(_to_float(_value(raw[6], f"{path}:7"), f"{path}:7"))
        src.factor = _to_float(_value(raw[7], f"{path}:8"), f"{path}:8")
        src.direction = (
            _to_float(_value(raw[8], f"{path}:9"), f"{path}:9"),
            _to_float(_value(raw[9], f"{path}:10"), f"{path}:10"),
            _to_float(_value(raw[10], f"{path}:11"), f"{path}:11"),
        )
        return src

    def _to_struct(self) -> CSource:
        c = CSource()
        c.source_type = GF_SRC_FORCE
        c.force_stf = int(self.stf)
        c.latitude = self.latitude
        c.longitude = self.longitude
        c.depth_km = self.depth
        c.hdur = self.f0
        c.time_shift = self.time_shift
        for i in range(6):
            c.moment[i] = 0.0
        c.force_factor = self.factor
        for i in range(3):
            c.force_dir[i] = self.direction[i]
        return c
