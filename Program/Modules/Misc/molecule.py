from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import SetDihedralDeg
from rdkit.Chem.rdmolfiles import MolToXYZBlock
from rdkit.Chem.AllChem import ConstrainedEmbed


from convertFile import convertFile
    

class Molecule():
    def __init__(self, molecule):
        if type(molecule) == str:
            if molecule.endswith(".xyz"):
                self.mol = Chem.MolFromXYZFile(molecule)
            else:
                self.mol = Chem.MolFromXYZBlock(molecule)
        else:
            self.mol = molecule

        try:
            if len(self.mol.GetAtoms()) > 100:
                raise ValueError
            from rdkit.Chem.rdDetermineBonds import DetermineBonds
            DetermineBonds(self.mol,charge = 0)
        except:
            from rdkit.Chem.rdDetermineBonds import DetermineConnectivity
            DetermineConnectivity(self.mol,charge = 0)

            self.mol.UpdatePropertyCache() #Black magic to make the structure work
            Chem.GetSymmSSSR(self.mol)
    
    def getMol(self):
        return self.mol
    
    def removeAtom(self, atomNr):
        edFrag = Chem.EditableMol(self.mol)
        edFrag.RemoveAtom(atomNr)
        self.mol = edFrag.GetMol()

    def setDihedralAgle(self, angle, dihedralIndex):
        SetDihedralDeg(self.mol.GetConformer(), int(dihedralIndex[0]), int(dihedralIndex[1]), int(dihedralIndex[2]), int(dihedralIndex[3]), angle)

    def getXYZBlock(self):
        return MolToXYZBlock(self.mol)
    
    def getCOMBlock(self):
        molBlock =Chem.MolToMolBlock(self.mol)
        return convertFile(inFile=molBlock, inFormat="mol", outFormat="com")
    
    def getGROBlock(self):
        molBlock = Chem.MolToMolBlock(self.mol)
        return convertFile(inFile=molBlock, inFormat="mol", outFormat="gro")

    @classmethod
    def combineMolecules(cls, mol1, mol2, index:tuple[int,int]):  #Takes the rdkit mols and combines them
        """Function to combine two rdkit-mol objects at a specific point

        Args:
            mol1 (_mol_): rdkit structure of molecule Nr.1
            mol2 (_mol_): rdkit Structure of molecule Nr.2
            index (_tuple_): A tuple with the index of the atoms to connect (eg. (12, 24)). IMPORTANT: The order has to match the order of the given Molecules

        Returns:
            _mol_: A New Molecule
        """
        combo = Chem.CombineMols(mol1, mol2)  #Combine
        edcombo = Chem.EditableMol(combo)     #Make editabel
        a1, a2 = index
        
        edcombo.AddBond(a1, mol1.GetNumAtoms()+ a2, order=Chem.rdchem.BondType.SINGLE)  #Combine the molecules with a bond
        combinedMol = edcombo.GetMol()  #Get the finished structure

        combinedMol.UpdatePropertyCache() #Black magic to make the structure work
        Chem.GetSymmSSSR(combinedMol)

        ConstrainedEmbed(combinedMol, mol2)#Somewhat relax the structure to make a belivable Molecule

        return Molecule(combinedMol)
