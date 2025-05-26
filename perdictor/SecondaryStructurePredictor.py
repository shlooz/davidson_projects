import os
from typing import Dict, List, Tuple

class StructurePredictor:
    def __init__(self, propensities_path: str):
        self.propensities = self._load_propensities(propensities_path)

    def _load_propensities(self, path: str) -> Dict[str, Dict[str, float]]:
        propensities = {}
        with open(path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.strip() and not line.startswith('-') and not line.startswith('AA'):
                    parts = line.split()
                    if len(parts) == 4:
                        aa, helix, sheet, coil = parts
                        propensities[aa.upper()] = {
                            'H': float(helix),
                            'E': float(sheet),
                            'C': float(coil)
                        }
        return propensities

    def predict(self,
                sequence: str,
                window_size: int = 6,
                min_region_length: int = 4,
                min_propensity: float = 1.0) -> List[Tuple[int, int, str]]:
        predictions = []
        seq_len = len(sequence)
        half_window = window_size // 2

        # Initialize scores for each position
        scores = [{'H': 0.0, 'E': 0.0, 'C': 0.0} for _ in range(seq_len)]

        # Calculate average propensities using sliding window
        for i in range(seq_len):
            window_start = max(0, i - half_window)
            window_end = min(seq_len, i + half_window + 1)
            window = sequence[window_start:window_end]
            counts = {'H': 0.0, 'E': 0.0, 'C': 0.0}
            for aa in window:
                aa_props = self.propensities.get(aa.upper())
                if aa_props:
                    for key in counts:
                        counts[key] += aa_props[key]
            window_length = len(window)
            for key in counts:
                scores[i][key] = counts[key] / window_length

        # Assign secondary structure based on highest average propensity
        assignments = []
        for score in scores:
            max_struct = max(score, key=score.get)
            assignments.append(max_struct)

        # Identify continuous regions
        regions = []
        start = 0
        current_struct = assignments[0]
        for i in range(1, seq_len):
            if assignments[i] != current_struct:
                if (i - start) >= min_region_length:
                    region_avg = sum([scores[j][current_struct] for j in range(start, i)]) / (i - start)
                    if region_avg >= min_propensity:
                        regions.append((start, i - 1, current_struct))
                start = i
                current_struct = assignments[i]
        # Check last region
        if (seq_len - start) >= min_region_length:
            region_avg = sum([scores[j][current_struct] for j in range(start, seq_len)]) / (seq_len - start)
            if region_avg >= min_propensity:
                regions.append((start, seq_len - 1, current_struct))

        return regions
