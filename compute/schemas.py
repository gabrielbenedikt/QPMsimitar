from pydantic import BaseModel, ConfigDict
from typing import Optional

class BaseComputeParams(BaseModel):
    """
    Base validation schema for compute tasks.
    We allow extra fields because the GUI sends many dynamic parameters.
    """
    model_config = ConfigDict(extra='allow', arbitrary_types_allowed=True)
    
    # Common parameters passed to almost all methods
    pump_wl: Optional[float] = None
    temperature: Optional[float] = None
    poling_period: Optional[float] = None
    crystal_length: Optional[float] = None
    qpm_order: Optional[int] = None
    pulsewidth: Optional[float] = None
    pump_shape: Optional[str] = None
    crystal_material: Optional[str] = None
