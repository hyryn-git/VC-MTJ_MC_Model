"""Public interface for the VC-MTJ compact model."""
from .model import AP, P, DevicePopulation, MTJModel, ResistanceResult, SimulationResult
from .parameters import ModelParameters

__all__ = ["MTJModel", "DevicePopulation", "SimulationResult", "ResistanceResult", "ModelParameters", "P", "AP"]
