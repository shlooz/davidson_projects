from FileMenager import AminoAcidSecondaryStructureCounter  # import your original class
import os

class AminoAcidPropensityCalculator:
    # Mapping from 3-letter to 1-letter amino acid codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    def __init__(self, pdbFolderPath):
        self.pdbFolderPath = pdbFolderPath
        self.totalCounts = {
            "helix": {},
            "sheet": {},
            "coil": {}
        }
        # Use same amino acids set to keep consistent
        self.standardAminoAcids = {
            'ALA' , 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }

    def addCounts(self, counts):
        for structure in ["helix", "sheet", "coil"]:
            for aa, count in counts[structure].items():
                if aa not in self.totalCounts[structure]:
                    self.totalCounts[structure][aa] = 0
                self.totalCounts[structure][aa] += count

    def calculatePropensities(self):
        propensities = {}
        for aa in sorted(self.standardAminoAcids):
            helix_count = self.totalCounts["helix"].get(aa, 0)
            sheet_count = self.totalCounts["sheet"].get(aa, 0)
            coil_count = self.totalCounts["coil"].get(aa, 0)
            total = helix_count + sheet_count + coil_count

            if total > 0:
                propensities[aa] = [
                    helix_count / total,
                    sheet_count / total,
                    coil_count / total
                ]
            else:
                propensities[aa] = [0.0, 0.0, 0.0]
        return propensities

    def processAllPDBFiles(self):
        # Go through all .pdb files in the folder and accumulate counts
        for filename in os.listdir(self.pdbFolderPath):
            if filename.lower().endswith(".pdb"):
                filepath = os.path.join(self.pdbFolderPath, filename)
                print(f"Processing {filepath}...")
                counter = AminoAcidSecondaryStructureCounter(filepath)
                counter.parseSecondaryStructureRegions()
                counter.processStructure()
                self.addCounts(counter.aminoAcidCounts)

    def savePropensitiesToFile(self, outputFilePath):
        propensities = self.calculatePropensities()
        with open(outputFilePath, "w") as f:
            f.write("AA  Helix   Sheet    Coil\n")
            f.write("-------------------------\n")
            for aa3 in sorted(self.standardAminoAcids):
                aa1 = self.three_to_one.get(aa3, "?")
                helix_prop, sheet_prop, coil_prop = propensities[aa3]
                f.write(f"{aa1:<3}{helix_prop:7.3f}{sheet_prop:8.3f}{coil_prop:8.3f}\n")

if __name__ == "__main__":
    pdb_folder = r"C:\Users\eitha\programing\Waitzman\final_project\dataFiles"
    output_file = os.path.join(r"C:\Users\eitha\programing\Waitzman\final_project", "propensities.txt")

    calculator = AminoAcidPropensityCalculator(pdb_folder)
    calculator.processAllPDBFiles()
    calculator.savePropensitiesToFile(output_file)
    print(f"Propensities saved to {output_file}")
