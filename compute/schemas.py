from pydantic import BaseModel, ConfigDict, Field
from typing import Optional, Any

class BaseComputeParams(BaseModel):
    """
    Base validation schema for compute tasks.
    We allow extra fields but strictly limit array allocations and resolutions
    to prevent memory exhaustion (DoS attacks).
    """
    model_config = ConfigDict(extra='allow', arbitrary_types_allowed=True)
    
    # Common parameters passed to almost all methods
    pump_wl: Optional[float] = Field(default=None, ge=100e-9, le=10000e-9)
    temperature: Optional[float] = Field(default=None, ge=-273.15, le=5000.0)
    poling_period: Optional[float] = Field(default=None, ge=1e-9, le=1e-2)
    crystal_length: Optional[float] = Field(default=None, ge=1e-9, le=1.0)
    qpm_order: Optional[int] = Field(default=None, ge=-100, le=100)
    pulsewidth: Optional[float] = Field(default=None, ge=1e-15, le=1e-3)
    pump_shape: Optional[str] = None
    crystal_material: Optional[str] = None
    
    # ----------------------------------------------------------------
    # SECURITY BOUNDS: Prevent DoS via Memory Exhaustion
    # Limit max points for any resolution array to reasonable sizes
    # ----------------------------------------------------------------
    resolution: Optional[int] = Field(default=None, ge=2, le=2000)
    jsi_resolution: Optional[int] = Field(default=None, ge=2, le=2000)
    hom_resolution: Optional[int] = Field(default=None, ge=2, le=5000)
    fwhm_resolution: Optional[int] = Field(default=None, ge=2, le=5000)
    wl_resolution: Optional[int] = Field(default=None, ge=2, le=5000)
    tau_resolution: Optional[int] = Field(default=None, ge=2, le=5000)
    L_resolution: Optional[int] = Field(default=None, ge=2, le=5000)
