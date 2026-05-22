import msgpack
import numpy as np

def encode_numpy(obj):
    if isinstance(obj, np.ndarray):
        return {
            '__nd__': True,
            'dtype': obj.dtype.str,
            'shape': obj.shape,
            'data': obj.tobytes()
        }
    elif isinstance(obj, (np.float32, np.float64, np.int32, np.int64)):
        return obj.item()
    return obj

def decode_numpy(obj):
    if '__nd__' in obj:
        return np.frombuffer(obj['data'], dtype=np.dtype(obj['dtype'])).reshape(obj['shape'])
    return obj

def serialize(data: dict) -> bytes:
    """Serialize a dictionary containing numpy arrays to msgpack bytes."""
    return msgpack.packb(data, default=encode_numpy, use_bin_type=True)

def deserialize(data: bytes) -> dict:
    """Deserialize msgpack bytes to a dictionary containing numpy arrays."""
    return msgpack.unpackb(data, object_hook=decode_numpy, raw=False)
