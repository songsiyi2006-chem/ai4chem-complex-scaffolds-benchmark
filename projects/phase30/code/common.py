"""Shared contracts, SI constants, bounded configuration and portable evidence."""
from dataclasses import dataclass, asdict, replace
import hashlib
import importlib.util
import json
import math
from pathlib import Path

F = 96485.33212                 # C mol-1
R = 8.314462618                 # J mol-1 K-1
KB_EV = 8.617333262145e-5       # eV K-1
KB_H = 2.083661912e10           # K-1 s-1
PRODUCTS = ("CO", "HCOOH", "C2H4", "H2")
ELECTRONS = dict(CO=2, HCOOH=2, C2H4=12, H2=2)
CARBON = dict(CO=1, HCOOH=1, C2H4=2, H2=0)
MW_KG = dict(CO=.028010, HCOOH=.046025, C2H4=.028054, H2=.002016)
LABEL = "UNCALIBRATED_REDUCED_ORDER_SCENARIO"


@dataclass(frozen=True)
class Config:
    seed: int = 30
    material: str = "Cu"
    target: str = "C2H4"
    potential_RHE_V: float = -.75
    ligand_shift_eV: float = 0.
    temperature_K: float = 298.15
    bulk_pH: float = 7.
    ionic_strength_mol_m3: float = 100.
    site_density_mol_m2: float = 2e-5
    co2_equilibrium_mol_m3: float = 33.
    co2_inlet_mol_m3: float = 33.
    gas_feed_mol_s: float = 2e-6
    length_m: float = .05
    width_m: float = .002
    height_m: float = .0005
    velocity_m_s: float = .05
    kla_s: float = 5.
    buffer_capacity_mol_m3_pH: float = 100.
    liquid_film_m_s: float = 1e-4
    carbonate_loss_s: float = .03
    area_resistance_ohm_m2: float = 5e-5
    heat_fraction: float = .25
    heat_transfer_W_m2_K: float = 100.
    max_temperature_K: float = 333.15
    max_current_A_m2: float = 20000.
    max_pressure_drop_Pa: float = 50000.
    cells: int = 12
    stack_kW: float = 100.
    stacks: int = 10
    years: int = 10
    capacity_factor: float = .85
    discount_rate: float = .08
    degradation_per_operating_h: float = 2e-5
    replacement_threshold: float = .8
    replacement_downtime_h: float = 24.
    capex_USD_kW: float = 1500.
    replacement_fraction: float = .25
    fixed_opex_fraction_year: float = .03
    electricity_USD_kWh: float = .05
    grid_kgCO2e_kWh: float = .1
    co2_USD_kg: float = .08
    feed_kgCO2e_kg: float = .2
    embodied_kgCO2e_kW: float = 200.
    separation_kWh_kg: float = 1.
    recovery_fraction: float = .9
    initial_samples: int = 6
    rounds: int = 3
    batch_size: int = 3

    def validate(self):
        for key, value in asdict(self).items():
            if isinstance(value, (int, float)) and (isinstance(value, bool) or not math.isfinite(value)):
                raise ValueError(f"Nonfinite/boolean numerical configuration: {key}")
        if self.material not in ("Cu", "Ag") or self.target not in ("CO", "HCOOH", "C2H4"):
            raise ValueError("Unsupported material or sale product")
        for key in ("seed", "cells", "stacks", "years", "initial_samples", "rounds", "batch_size"):
            if type(getattr(self, key)) is not int:
                raise ValueError(f"{key} must be an integer")
        if not (4 <= self.cells <= 192 and 1 <= self.years <= 30 and 1 <= self.stacks <= 100):
            raise ValueError("Discretization or planning horizon outside bounded pilot")
        if not (2 <= self.initial_samples <= 12 and 1 <= self.rounds <= 6 and 1 <= self.batch_size <= 3):
            raise ValueError("Bounded active-learning budget exceeded")
        if not (0 <= self.seed < 2**32 and -.95 <= self.potential_RHE_V <= -.4
                and abs(self.ligand_shift_eV) <= .08 and 4 <= self.bulk_pH <= 10):
            raise ValueError("Outside microkinetic pilot envelope")
        positive = ("temperature_K", "ionic_strength_mol_m3", "site_density_mol_m2",
                    "co2_equilibrium_mol_m3", "gas_feed_mol_s", "length_m", "width_m",
                    "height_m", "velocity_m_s", "kla_s", "buffer_capacity_mol_m3_pH",
                    "liquid_film_m_s", "heat_transfer_W_m2_K", "max_current_A_m2",
                    "max_pressure_drop_Pa", "stack_kW", "capex_USD_kW")
        if any(getattr(self, k) <= 0 for k in positive):
            raise ValueError("Positive physical parameter required")
        if not (273.15 <= self.temperature_K < self.max_temperature_K <= 373.15):
            raise ValueError("Unsupported temperature envelope")
        if not (0 <= self.co2_inlet_mol_m3 <= self.co2_equilibrium_mol_m3):
            raise ValueError("Inlet CO2 must be in [0, solubility reference]")
        for k in ("capacity_factor", "heat_fraction", "recovery_fraction",
                  "replacement_threshold", "replacement_fraction"):
            if not 0 < getattr(self, k) <= 1:
                raise ValueError(f"Invalid fraction: {k}")
        for k in ("carbonate_loss_s", "area_resistance_ohm_m2", "discount_rate",
                  "degradation_per_operating_h", "replacement_downtime_h",
                  "fixed_opex_fraction_year", "electricity_USD_kWh", "grid_kgCO2e_kWh",
                  "co2_USD_kg", "feed_kgCO2e_kg", "embodied_kgCO2e_kW", "separation_kWh_kg"):
            if getattr(self, k) < 0:
                raise ValueError(f"Negative parameter: {k}")
        if self.replacement_threshold >= 1 or self.discount_rate > 1:
            raise ValueError("Replacement threshold must be <1 and discount rate <=1")
        if self.height_m > self.width_m:
            raise ValueError("Rectangular laminar closure requires height <= width")
        return self

    def with_changes(self, **kwargs):
        return replace(self, **kwargs).validate()


def load_config(path=None):
    values = json.loads(Path(path).read_text(encoding="utf-8")) if path else {}
    return Config(**values).validate()  # unknown keys fail, never silently ignored


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(raw, encoding="utf-8")
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(filename):
    path = Path(__file__).parent / filename
    spec = importlib.util.spec_from_file_location("phase30_" + path.stem, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def validate_micro_artifact(path):
    """Audit an external input packet. This does not validate its chemistry."""
    path = Path(path).resolve()
    data = json.loads(path.read_text(encoding="utf-8"))
    if data["schema_version"] != "phase30.micro.v1":
        raise ValueError("Unsupported micro schema")
    if data["evidence_type"] not in ("computed", "literature"):
        raise ValueError("External micro packet must identify computed/literature evidence")
    required = ("method", "functional", "solvation", "potential_reference", "temperature_K",
                "surface_id", "geometry_sha256", "convergence", "records", "raw_files")
    if any(k not in data for k in required) or data["potential_reference"] != "SHE":
        raise ValueError("Missing method metadata or ambiguous electrode reference")
    if not data["records"] or not data["raw_files"]:
        raise ValueError("No rates/barriers or supporting raw files")
    for row in data["records"]:
        if row["unit"] != "eV" or not math.isfinite(row["barrier_eV"]) or row["barrier_eV"] < 0:
            raise ValueError("Invalid barrier or units")
        if not row.get("reaction_id") or not math.isfinite(row["potential_SHE_V"]):
            raise ValueError("Missing reaction/potential")
    for entry in data["raw_files"]:
        raw = (path.parent / entry["path"]).resolve()
        if not raw.is_relative_to(path.parent) or not raw.is_file() or sha(raw) != entry["sha256"]:
            raise ValueError("Raw evidence traversal, absence or hash mismatch")
    return dict(schema_valid=True, raw_hashes_verified=True,
                chemistry_validated=False, packet_sha256=sha(path),
                production_adapter_implemented=False)
