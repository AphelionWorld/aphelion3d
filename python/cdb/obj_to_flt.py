"""
Convert OBJ files in ENU coordinates to OpenFlight (FLT) format for CDB models.

OpenFlight is a binary format used by CDB for 3D models. This module provides
a simplified writer that outputs valid FLT files for static geometry with LOD support.
"""
from __future__ import annotations
import struct
from pathlib import Path
from dataclasses import dataclass
from typing import BinaryIO
import numpy as np


# OpenFlight record opcodes
OPCODE_HEADER = 1
OPCODE_GROUP = 2
OPCODE_OBJECT = 4
OPCODE_FACE = 5
OPCODE_PUSH_LEVEL = 10
OPCODE_POP_LEVEL = 11
OPCODE_LOD = 73
OPCODE_VERTEX_PALETTE = 67
OPCODE_VERTEX_WITH_COLOR = 68
OPCODE_VERTEX_WITH_NORMAL = 69
OPCODE_VERTEX_WITH_NORMAL_UV = 70
OPCODE_VERTEX_WITH_COLOR_NORMAL = 71
OPCODE_MATERIAL_PALETTE = 113
OPCODE_VERTEX_LIST = 72
OPCODE_LONG_ID = 33


@dataclass
class OBJVertex:
    """Vertex data from OBJ file"""
    position: tuple[float, float, float]
    normal: tuple[float, float, float] | None = None
    texcoord: tuple[float, float] | None = None


@dataclass
class OBJFace:
    """Face data from OBJ file"""
    vertex_indices: list[int]  # 0-based indices into vertex list


@dataclass
class LODLevel:
    """LOD level with geometry and distance thresholds"""
    vertices: list[OBJVertex]
    faces: list[OBJFace]
    switch_in_distance: float  # Distance to switch TO this LOD
    switch_out_distance: float  # Distance to switch FROM this LOD
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)


class OpenFlightWriter:
    """Writer for OpenFlight FLT format files with LOD support"""
    
    def __init__(self, filepath: Path):
        self.filepath = Path(filepath)
        self.lod_levels: list[LODLevel] = []
        
    def write_header(self, f: BinaryIO):
        """Write OpenFlight header record"""
        # Header is 320 bytes
        opcode = OPCODE_HEADER
        length = 320
        
        # Write opcode and length
        f.write(struct.pack('>HH', opcode, length))
        
        # Format revision level (e.g., 1620 for v16.2)
        f.write(struct.pack('>i', 1620))
        
        # Database origin (ENU = 0)
        f.write(struct.pack('>i', 0))
        
        # Coordinate units (1 = meters)
        f.write(struct.pack('>i', 1))
        
        # Flags
        f.write(struct.pack('>i', 0))
        
        # Reserved fields
        f.write(b'\x00' * 24)
        
        # Vertex storage type (0 = normal)
        f.write(struct.pack('>i', 0))
        
        # Database origin (again)
        f.write(struct.pack('>i', 0))
        
        # SW corner lat/lon (for georeferenced models, 0 for ENU)
        f.write(struct.pack('>dd', 0.0, 0.0))
        
        # NE corner lat/lon
        f.write(struct.pack('>dd', 0.0, 0.0))
        
        # Origin lat/lon/height
        f.write(struct.pack('>ddd', 0.0, 0.0, 0.0))
        
        # Lambert projection parameters (unused for ENU)
        f.write(b'\x00' * 80)
        
        # Next group ID, face ID, etc.
        f.write(struct.pack('>HHHH', 0, 0, 0, 0))
        
        # Earth model (0 = WGS84)
        f.write(struct.pack('>i', 0))
        
        # Padding to 320 bytes
        current = 4 + 4 + 4*6 + 24 + 4*2 + 8*6 + 80 + 8 + 4
        padding = 320 - current
        f.write(b'\x00' * padding)
        
    def write_vertex_palette(self, f: BinaryIO, vertices: list[OBJVertex]):
        """Write vertex palette record"""
        if not vertices:
            return
            
        opcode = OPCODE_VERTEX_PALETTE
        vertex_data_size = 0
        
        for v in vertices:
            if v.normal is not None:
                vertex_data_size += 64
            else:
                vertex_data_size += 36
        
        length = 8 + vertex_data_size
        f.write(struct.pack('>HH', opcode, length))
        f.write(struct.pack('>i', 0))
        
        for idx, v in enumerate(vertices):
            if v.normal is not None:
                self._write_vertex_with_normal(f, idx, v)
            else:
                self._write_basic_vertex(f, idx, v)
    
    def _write_vertex_with_normal(self, f: BinaryIO, index: int, vertex: OBJVertex):
        """Write a vertex with normal (opcode 70)"""
        opcode = OPCODE_VERTEX_WITH_NORMAL
        length = 64
        
        f.write(struct.pack('>HH', opcode, length))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>H', 0))
        
        f.write(struct.pack('>ddd', 
                           vertex.position[0],
                           vertex.position[1], 
                           vertex.position[2]))
        
        if vertex.normal:
            f.write(struct.pack('>fff',
                               vertex.normal[0],
                               vertex.normal[1],
                               vertex.normal[2]))
        else:
            f.write(struct.pack('>fff', 0.0, 0.0, 1.0))
        
        if vertex.texcoord:
            f.write(struct.pack('>ff', vertex.texcoord[0], vertex.texcoord[1]))
        else:
            f.write(struct.pack('>ff', 0.0, 0.0))
        
        f.write(struct.pack('>I', 0xFFFFFFFF))
        f.write(struct.pack('>i', index))
        
    def _write_basic_vertex(self, f: BinaryIO, index: int, vertex: OBJVertex):
        """Write a basic vertex without normal"""
        opcode = OPCODE_VERTEX_WITH_COLOR
        length = 36
        
        f.write(struct.pack('>HH', opcode, length))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>H', 0))
        
        f.write(struct.pack('>ddd',
                           vertex.position[0],
                           vertex.position[1],
                           vertex.position[2]))
        
        f.write(struct.pack('>I', 0xFFFFFFFF))
        f.write(struct.pack('>i', index))
        
    def write_lod_node(self, f: BinaryIO, lod_level: LODLevel):
        """Write LOD node with switch distances"""
        opcode = OPCODE_LOD
        length = 52
        
        f.write(struct.pack('>HH', opcode, length))
        
        # LOD ID (8 bytes)
        f.write(b'lod\x00\x00\x00\x00\x00')
        
        # Reserved
        f.write(struct.pack('>i', 0))
        
        # Switch-in distance
        f.write(struct.pack('>d', lod_level.switch_in_distance))
        
        # Switch-out distance
        f.write(struct.pack('>d', lod_level.switch_out_distance))
        
        # Special effect ID1, ID2
        f.write(struct.pack('>hh', 0, 0))
        
        # Flags (0 = use previous slant range)
        f.write(struct.pack('>I', 0))
        
        # Center X, Y, Z
        f.write(struct.pack('>ddd', 
                           lod_level.center[0],
                           lod_level.center[1],
                           lod_level.center[2]))
        
    def write_group(self, f: BinaryIO, name: str = "model"):
        """Write a group record"""
        opcode = OPCODE_GROUP
        length = 44
        
        f.write(struct.pack('>HH', opcode, length))
        
        group_id = name[:7].ljust(7, '\x00')
        f.write(group_id.encode('ascii'))
        
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>i', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>b', 0))
        f.write(b'\x00' * 5)
        f.write(struct.pack('>ii', 0, 0))
        f.write(struct.pack('>ii', 0, 0))
        
    def write_object(self, f: BinaryIO, name: str = "geometry"):
        """Write an object record"""
        opcode = OPCODE_OBJECT
        length = 28
        
        f.write(struct.pack('>HH', opcode, length))
        
        obj_id = name[:7].ljust(7, '\x00')
        f.write(obj_id.encode('ascii'))
        
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>I', 0x80000000))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>H', 0))
        
    def write_face(self, f: BinaryIO, face: OBJFace):
        """Write a face record"""
        opcode = OPCODE_FACE
        length = 80
        
        f.write(struct.pack('>HH', opcode, length))
        f.write(b'face\x00\x00\x00\x00')
        f.write(struct.pack('>i', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>b', 2))
        f.write(struct.pack('>b', 1))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>b', 0))
        f.write(struct.pack('>b', 0))
        f.write(struct.pack('>h', -1))
        f.write(struct.pack('>h', -1))
        f.write(struct.pack('>h', -1))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>i', 0))
        f.write(struct.pack('>H', 0))
        f.write(struct.pack('>B', 0))
        f.write(struct.pack('>B', 0))
        f.write(struct.pack('>I', 0))
        f.write(struct.pack('>B', 2))
        f.write(b'\x00' * 7)
        f.write(struct.pack('>I', 0xFFFFFFFF))
        f.write(struct.pack('>I', 0xFFFFFFFF))
        f.write(struct.pack('>h', -1))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>i', 0))
        f.write(struct.pack('>i', 0))
        f.write(struct.pack('>h', 0))
        f.write(struct.pack('>h', -1))
        
    def write_vertex_list(self, f: BinaryIO, face: OBJFace):
        """Write vertex list for a face"""
        opcode = OPCODE_VERTEX_LIST
        num_verts = len(face.vertex_indices)
        length = 4 + 4 * num_verts
        
        f.write(struct.pack('>HH', opcode, length))
        
        for idx in face.vertex_indices:
            f.write(struct.pack('>i', idx))
    
    def write_push_level(self, f: BinaryIO):
        """Write push level (start hierarchy)"""
        f.write(struct.pack('>HH', OPCODE_PUSH_LEVEL, 4))
        
    def write_pop_level(self, f: BinaryIO):
        """Write pop level (end hierarchy)"""
        f.write(struct.pack('>HH', OPCODE_POP_LEVEL, 4))
    
    def write(self):
        """Write complete FLT file with LOD support"""
        with open(self.filepath, 'wb') as f:
            # Write header
            self.write_header(f)
            
            if not self.lod_levels:
                return
            
            # Write vertex palettes for all LOD levels
            for lod in self.lod_levels:
                self.write_vertex_palette(f, lod.vertices)
            
            # Write hierarchy
            self.write_group(f, "model")
            self.write_push_level(f)
            
            # Write each LOD level
            for lod_idx, lod in enumerate(self.lod_levels):
                self.write_lod_node(f, lod)
                self.write_push_level(f)
                
                self.write_object(f, f"lod{lod_idx}")
                self.write_push_level(f)
                
                # Write faces for this LOD
                for face in lod.faces:
                    self.write_face(f, face)
                    self.write_push_level(f)
                    self.write_vertex_list(f, face)
                    self.write_pop_level(f)
                
                self.write_pop_level(f)  # Close object
                self.write_pop_level(f)  # Close LOD
            
            self.write_pop_level(f)  # Close group


def load_obj(filepath: Path) -> tuple[list[OBJVertex], list[OBJFace]]:
    """Load OBJ file and return vertices and faces (assumes ENU coordinates)."""
    vertices: list[OBJVertex] = []
    normals: list[tuple[float, float, float]] = []
    texcoords: list[tuple[float, float]] = []
    faces: list[OBJFace] = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if not parts:
                continue
            
            cmd = parts[0]
            
            if cmd == 'v':
                x, y, z = map(float, parts[1:4])
                vertices.append(OBJVertex(position=(x, y, z)))
                
            elif cmd == 'vn':
                nx, ny, nz = map(float, parts[1:4])
                normals.append((nx, ny, nz))
                
            elif cmd == 'vt':
                u, v = map(float, parts[1:3])
                texcoords.append((u, v))
                
            elif cmd == 'f':
                indices = []
                for vert_str in parts[1:]:
                    components = vert_str.split('/')
                    v_idx = int(components[0]) - 1
                    
                    if len(components) >= 3 and components[2]:
                        n_idx = int(components[2]) - 1
                        if n_idx < len(normals):
                            vertices[v_idx].normal = normals[n_idx]
                    
                    if len(components) >= 2 and components[1]:
                        t_idx = int(components[1]) - 1
                        if t_idx < len(texcoords):
                            vertices[v_idx].texcoord = texcoords[t_idx]
                    
                    indices.append(v_idx)
                
                # Triangulate
                if len(indices) == 3:
                    faces.append(OBJFace(vertex_indices=indices))
                elif len(indices) == 4:
                    faces.append(OBJFace(vertex_indices=[indices[0], indices[1], indices[2]]))
                    faces.append(OBJFace(vertex_indices=[indices[0], indices[2], indices[3]]))
                elif len(indices) > 4:
                    for i in range(1, len(indices) - 1):
                        faces.append(OBJFace(vertex_indices=[indices[0], indices[i], indices[i+1]]))
    
    return vertices, faces


def calculate_bounding_center(vertices: list[OBJVertex]) -> tuple[float, float, float]:
    """Calculate the center of bounding box for vertices"""
    if not vertices:
        return (0.0, 0.0, 0.0)
    
    positions = np.array([v.position for v in vertices])
    min_bounds = positions.min(axis=0)
    max_bounds = positions.max(axis=0)
    center = (min_bounds + max_bounds) / 2.0
    
    return tuple(center)


def convert_obj_to_flt_single_lod(obj_path: Path, flt_path: Path) -> bool:
    """Convert a single OBJ file to FLT (single LOD, no switching)."""
    try:
        vertices, faces = load_obj(obj_path)
        
        if not vertices or not faces:
            print(f"Warning: No geometry found in {obj_path}")
            return False
        
        writer = OpenFlightWriter(flt_path)
        
        # Create single LOD with no distance switching
        center = calculate_bounding_center(vertices)
        lod = LODLevel(
            vertices=vertices,
            faces=faces,
            switch_in_distance=0.0,
            switch_out_distance=1e10,
            center=center
        )
        writer.lod_levels = [lod]
        writer.write()
        
        print(f"Converted {obj_path.name} -> {flt_path.name}")
        print(f"  Vertices: {len(vertices)}, Faces: {len(faces)}")
        
        return True
        
    except Exception as e:
        print(f"Error converting {obj_path}: {e}")
        return False


def convert_obj_directory_to_cdb(
    obj_dir: Path,
    cdb_root: Path,
    lod_distances: dict[int, tuple[float, float]] | None = None
) -> list[Path]:
    """
    Convert directory of OBJ files organized by LOD to CDB FLT models.
    
    Expected directory structure:
        obj_dir/
            lod0/
                model1.obj
                model2.obj
            lod1/
                model1.obj
                model2.obj
            ...
    
    Or flat structure (all models at single LOD):
        obj_dir/
            model1.obj
            model2.obj
    
    Args:
        obj_dir: Directory containing OBJ files (organized by lod0/, lod1/, etc. or flat)
        cdb_root: CDB root directory
        lod_distances: Optional dict mapping LOD index to (switch_in, switch_out) distances
                      Default: LOD0=(0, 10000), LOD1=(10000, 5000), LOD2=(5000, 1000), etc.
    
    Returns:
        List of created FLT file paths
    """
    from .paths import CDBPaths
    from .datasets import Dataset
    
    if lod_distances is None:
        # Default LOD switching distances (in meters)
        lod_distances = {
            0: (0.0, 10000.0),
            1: (10000.0, 5000.0),
            2: (5000.0, 2000.0),
            3: (2000.0, 500.0),
            4: (500.0, 100.0),
            5: (100.0, 0.0),
        }
    
    cdb = CDBPaths(cdb_root)
    models_dataset = Dataset(300, "Models", "flt")
    
    obj_dir = Path(obj_dir)
    written = []
    
    # Check if directory has LOD subdirectories
    lod_dirs = [d for d in obj_dir.iterdir() if d.is_dir() and d.name.lower().startswith('lod')]
    
    if lod_dirs:
        # Multi-LOD structure: group by model name
        print(f"Found LOD subdirectories: {[d.name for d in lod_dirs]}")
        
        # Collect all model names across LODs
        model_names = set()
        for lod_dir in lod_dirs:
            for obj_file in lod_dir.glob("*.obj"):
                model_names.add(obj_file.stem)
        
        # Process each model
        for model_name in sorted(model_names):
            print(f"\nProcessing model: {model_name}")
            
            writer = OpenFlightWriter(obj_dir / f"{model_name}.flt")  # Temp path
            
            # Load all LOD levels for this model
            lod_levels = []
            for lod_dir in sorted(lod_dirs):
                lod_num_str = lod_dir.name.lower().replace('lod', '')
                try:
                    lod_num = int(lod_num_str)
                except ValueError:
                    continue
                
                obj_file = lod_dir / f"{model_name}.obj"
                if not obj_file.exists():
                    continue
                
                vertices, faces = load_obj(obj_file)
                if not vertices or not faces:
                    continue
                
                center = calculate_bounding_center(vertices)
                switch_in, switch_out = lod_distances.get(lod_num, (0.0, 1e10))
                
                lod_level = LODLevel(
                    vertices=vertices,
                    faces=faces,
                    switch_in_distance=switch_in,
                    switch_out_distance=switch_out,
                    center=center
                )
                lod_levels.append((lod_num, lod_level))
                print(f"  Loaded {lod_dir.name}: {len(vertices)} vertices, {len(faces)} faces")
            
            if not lod_levels:
                print(f"  No geometry found for {model_name}")
                continue
            
            # Use the highest LOD number to determine CDB LOD placement
            max_lod_num = max(lod_num for lod_num, _ in lod_levels)
            
            # Place in CDB structure (use geocell 0,0 and tile 0,0 for geospecific models)
            flt_path = cdb.tile_path(models_dataset, max_lod_num, 0, 0, 0, 0)
            flt_dir = flt_path.parent / model_name
            flt_dir.mkdir(parents=True, exist_ok=True)
            flt_file = flt_dir / f"{model_name}.flt"
            
            # Write FLT with all LOD levels
            writer.filepath = flt_file
            writer.lod_levels = [lod_level for _, lod_level in sorted(lod_levels)]
            writer.write()
            
            print(f"  Created {flt_file.relative_to(cdb_root)}")
            written.append(flt_file)
    
    else:
        # Flat structure: single LOD per model
        print(f"Processing flat directory (single LOD)")
        
        for obj_file in obj_dir.glob("*.obj"):
            model_name = obj_file.stem
            
            vertices, faces = load_obj(obj_file)
            if not vertices or not faces:
                continue
            
            # Place at LOD 5 by default for flat structure
            lod_num = 5
            flt_path = cdb.tile_path(models_dataset, lod_num, 0, 0, 0, 0)
            flt_dir = flt_path.parent / model_name
            flt_dir.mkdir(parents=True, exist_ok=True)
            flt_file = flt_dir / f"{model_name}.flt"
            
            writer = OpenFlightWriter(flt_file)
            center = calculate_bounding_center(vertices)
            lod_level = LODLevel(
                vertices=vertices,
                faces=faces,
                switch_in_distance=0.0,
                switch_out_distance=1e10,
                center=center
            )
            writer.lod_levels = [lod_level]
            writer.write()
            
            print(f"Created {flt_file.relative_to(cdb_root)}")
            written.append(flt_file)
    
    return written


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage:")
        print("  Single file: python obj_to_flt.py <input.obj> <output.flt>")
        print("  Directory:   python obj_to_flt.py <obj_dir> <cdb_root>")
        sys.exit(1)
    
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    
    if input_path.is_file():
        convert_obj_to_flt_single_lod(input_path, output_path)
    elif input_path.is_dir():
        written = convert_obj_directory_to_cdb(input_path, output_path)
        print(f"\nTotal: {len(written)} models converted to CDB")
    else:
        print(f"Error: {input_path} not found")
        sys.exit(1)