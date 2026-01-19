from __future__ import annotations
from pathlib import Path

class ModelWriter:
    def __init__(self, root: Path):
        self.root = Path(root)

    def write_model_placeholder(self) -> Path:
        p = self.root / "D300_Models" / "README.txt"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("Models dataset placeholder. Wire model export here.")
        return p