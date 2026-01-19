from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from .datasets import Dataset

@dataclass
class CDBLayout:
    use_folders_per_lod: bool = True
    use_geocell_dirs: bool = True

class CDBPaths:
    def __init__(self, root: Path, layout: CDBLayout | None = None):
        self.root = Path(root)
        self.layout = layout or CDBLayout()

    @staticmethod
    def geocell_dir(lat_cell: int, lon_cell: int) -> str:
        ns = "N" if lat_cell >= 0 else "S"
        ew = "E" if lon_cell >= 0 else "W"
        return f"{ns}{abs(lat_cell):02d}{ew}{abs(lon_cell):03d}"

    def dataset_root(self, dataset: Dataset) -> Path:
        return self.root / f"D{dataset.id}_{dataset.name}"

    def tile_path(self, dataset: Dataset, lod: int, lat_cell: int, lon_cell: int, u: int, v: int) -> Path:
        base = self.dataset_root(dataset)
        if self.layout.use_folders_per_lod:
            base = base / f"LOD{lod}"
        if self.layout.use_geocell_dirs:
            base = base / self.geocell_dir(lat_cell, lon_cell)
        base.mkdir(parents=True, exist_ok=True)
        return base / f"U{u:02d}_V{v:02d}.{dataset.default_ext}"