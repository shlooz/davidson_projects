from Bio.PDB import PDBParser

class AminoAcidSecondaryStructureCounter:
    def __init__(self, pdbFilePath):
        self.pdbFilePath = pdbFilePath
        self.helixRegions = []
        self.sheetRegions = []
        self.structure = None

        self.aminoAcidCounts = {
            "helix": {},
            "sheet": {},
            "coil": {}
        }

        self.standardAminoAcids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }

    def parseSecondaryStructureRegions(self):
        pdbFile = open(self.pdbFilePath, "r")
        for line in pdbFile:
            if line.startswith("HELIX "):
                chainId = line[19].strip()
                startIndex = int(line[21:25].strip())
                endIndex = int(line[33:37].strip())
                self.helixRegions.append((chainId, startIndex, endIndex))
            elif line.startswith("SHEET "):
                chainId = line[21].strip()
                startIndex = int(line[22:26].strip())
                endIndex = int(line[33:37].strip())
                self.sheetRegions.append((chainId, startIndex, endIndex))
        pdbFile.close()

    def isResidueInRegion(self, chainId, residueIndex, regionList):
        for regionChainId, startIndex, endIndex in regionList:
            if chainId == regionChainId and startIndex <= residueIndex <= endIndex:
                return True
        return False

    def getStructureTypeForResidue(self, chainId, residueIndex):
        if self.isResidueInRegion(chainId, residueIndex, self.helixRegions):
            return "helix"
        elif self.isResidueInRegion(chainId, residueIndex, self.sheetRegions):
            return "sheet"
        else:
            return "coil"

    def incrementAminoAcidCount(self, structureType, aminoAcidName):
        if aminoAcidName not in self.aminoAcidCounts[structureType]:
            self.aminoAcidCounts[structureType][aminoAcidName] = 0
        self.aminoAcidCounts[structureType][aminoAcidName] += 1

    def processResidue(self, chain, residue):
        aminoAcidName = residue.get_resname()
        if aminoAcidName not in self.standardAminoAcids:
            return

        chainId = chain.id
        residueIndex = residue.get_id()[1]
        structureType = self.getStructureTypeForResidue(chainId, residueIndex)
        self.incrementAminoAcidCount(structureType, aminoAcidName)

    def processStructure(self):
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("pdbStructure", self.pdbFilePath)

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    self.processResidue(chain, residue)

    def displayAminoAcidCounts(self):
        for structureType in ["helix", "sheet", "coil"]:
            print("\nAmino acid counts in", structureType + ":")
            for aminoAcidName in sorted(self.aminoAcidCounts[structureType].keys()):
                print(aminoAcidName + ":", self.aminoAcidCounts[structureType][aminoAcidName])

    def run(self):
        self.parseSecondaryStructureRegions()
        self.processStructure()
        self.displayAminoAcidCounts()

