from __future__ import annotations
from pathlib import Path
import struct
import numpy as np
from typing import Callable
import subprocess
import shutil


def convert_obj_to_flt_osg(
    obj_path: Path, 
    flt_path: Path,
    osg_conv_path: str | None = None
) -> bool:
    """
    Convert OBJ to FLT using OpenSceneGraph's osgconv binary.
    
    Args:
        obj_path: Path to input OBJ file
        flt_path: Path to output FLT file
        osg_conv_path: Path to osgconv binary (default: searches PATH)
    
    Returns:
        True if conversion succeeded, False otherwise
    """
    # Find osgconv binary
    if osg_conv_path is None:
        osg_conv_path = shutil.which("osgconv")
        if osg_conv_path is None:
            print("Error: osgconv not found in PATH. Please install OpenSceneGraph.")
            return False
    
    # Ensure output directory exists
    flt_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Run osgconv
    try:
        cmd = [osg_conv_path, str(obj_path), str(flt_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode != 0:
            print(f"osgconv failed: {result.stderr}")
            return False
        
        # Check if output file was created
        if not flt_path.exists():
            print(f"osgconv succeeded but output file not found: {flt_path}")
            return False
        
        return True
        
    except subprocess.TimeoutExpired:
        print(f"osgconv timed out converting {obj_path}")
        return False
    except Exception as e:
        print(f"Error running osgconv: {e}")
        return False


class FLTWriter:
    """
    OpenFlight FLT format writer for CDB models.
    
    OpenFlight is a binary format used by CDB for 3D models.
    This is a simplified writer that creates valid FLT files.
    """
    
    # OpenFlight magic number
    MAGIC = b'AF\x00\x00'
    
    # Record type codes
    RECORD_HEADER = 0
    RECORD_GROUP_START = 2
    RECORD_GROUP_END = 3
    RECORD_POLYGON = 5
    RECORD_VERTEX_WITH_COLOR = 110
    
    def __init__(self):
        self.data = bytearray()
        self.vertex_list = []
        
    def write_header(self):
        """Write FLT file header."""
        # Header record
        record_type = self.RECORD_HEADER
        record_size = 324
        
        self.data.extend(struct.pack('>HH', record_type, record_size))
        self.data.extend(b'Created by Aphelion3D CDB Converter\x00')
        self.data.extend(b'\x00' * (100 - len('Created by Aphelion3D CDB Converter')))
        # Format revision (15.8)
        self.data.extend(struct.pack('>I', 15008))
        # Units (meters)
        self.data.extend(struct.pack('>I', 0))
        # Set flags and other header fields to defaults
        self.data.extend(b'\x00' * (record_size - len(self.data) % record_size))
    
    def add_vertex(self, x: float, y: float, z: float, color: tuple[int, int, int] = (255, 255, 255)):
        """Add a vertex to the model."""
        self.vertex_list.append((x, y, z, color))
    
    def add_triangle(self, v1_idx: int, v2_idx: int, v3_idx: int):
        """Add a triangle face using vertex indices."""
        if v1_idx < 0 or v2_idx < 0 or v3_idx < 0:
            return
        if v1_idx >= len(self.vertex_list) or v2_idx >= len(self.vertex_list) or v3_idx >= len(self.vertex_list):
            return
        
        # Polygon record
        record_type = self.RECORD_POLYGON
        record_size = 68
        
        self.data.extend(struct.pack('>HH', record_type, record_size))
        self.data.extend(struct.pack('>I', 0))  # IRColor
        self.data.extend(struct.pack('>HH', 3, 0))  # Number of vertices, reserved
        self.data.extend(struct.pack('>I', 0))  # Edge flags
        
        # Write vertices for this polygon
        for idx in [v1_idx, v2_idx, v3_idx]:
            x, y, z, color = self.vertex_list[idx]
            self.data.extend(struct.pack('>fff', x, y, z))
        
        # Padding to reach 68 bytes
        self.data.extend(b'\x00' * (record_size - (len(self.data) % record_size)))
    
    def write_to_file(self, path: Path):
        """Write FLT file to disk."""
        path = Path(path).resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(path, 'wb') as f:
            f.write(self.MAGIC)
            f.write(self.data)
        
        return path


def convert_obj_to_flt_single_lod(input_obj: Path, output_flt: Path) -> bool:
    """
    Convert a single OBJ file to FLT format.
    
    Args:
        input_obj: Path to input OBJ file
        output_flt: Path to output FLT file
    
    Returns:
        True if successful, False otherwise
    """
    try:
        input_obj = Path(input_obj).resolve()
        output_flt = Path(output_flt).resolve()
        
        if not input_obj.exists():
            raise FileNotFoundError(f"Input file not found: {input_obj}")
        
        # Parse OBJ file
        vertices = []
        faces = []
        
        with open(input_obj, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if not parts:
                    continue
                
                if parts[0] == 'v':
                    # Vertex
                    try:
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        vertices.append((x, y, z))
                    except (ValueError, IndexError):
                        continue
                
                elif parts[0] == 'f':
                    # Face (support triangles and quads)
                    try:
                        indices = []
                        for i in range(1, len(parts)):
                            # Handle v, v/vt, v/vt/vn, v//vn formats
                            vertex_data = parts[i].split('/')
                            v_idx = int(vertex_data[0]) - 1  # OBJ uses 1-based indexing
                            indices.append(v_idx)
                        
                        # Triangulate quads and larger polygons
                        if len(indices) >= 3:
                            for i in range(1, len(indices) - 1):
                                faces.append((indices[0], indices[i], indices[i + 1]))
                    except (ValueError, IndexError):
                        continue
        
        if not vertices:
            raise ValueError("No vertices found in OBJ file")
        if not faces:
            raise ValueError("No faces found in OBJ file")
        
        # Create FLT writer
        writer = FLTWriter()
        writer.write_header()
        
        # Add vertices
        for x, y, z in vertices:
            writer.add_vertex(x, y, z)
        
        # Add faces
        for v1, v2, v3 in faces:
            writer.add_triangle(v1, v2, v3)
        
        # Write FLT file
        output_flt = Path(output_flt).resolve()
        output_flt.parent.mkdir(parents=True, exist_ok=True)
        
        writer.write_to_file(output_flt)
        
        return True
    
    except Exception as e:
        print(f"Error converting {input_obj}: {e}")
        return False


def convert_obj_directory_to_cdb(input_dir: Path, output_cdb_root: Path, dataset=None) -> list[Path]:
    """
    Convert all OBJ files in a directory to CDB geotypical models.
    
    Args:
        input_dir: Directory containing OBJ files
        output_cdb_root: CDB root directory
        dataset: Dataset object with code and name
    
    Returns:
        List of paths to written FLT files
    """
    written = []
    
    input_dir = Path(input_dir).resolve()
    output_cdb_root = Path(output_cdb_root).resolve()
    
    obj_files = list(input_dir.glob("**/*.obj"))
    
    for obj_file in obj_files:
        try:
            # Create temporary FLT path (caller handles final placement)
            temp_flt = output_cdb_root / f"_{obj_file.stem}.flt"
            
            if convert_obj_to_flt_single_lod(obj_file, temp_flt):
                written.append(temp_flt)
        
        except Exception as e:
            print(f"Failed to convert {obj_file}: {e}")
    
    return written