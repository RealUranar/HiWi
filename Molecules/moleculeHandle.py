from rdkit import Chem
from rdkit.Chem.rdDetermineBonds import DetermineBonds
from rdkit.Chem import AllChem

class Molecule():
    def __init__(self, molecule):
        if type(molecule) == str:
            self.mol = Chem.MolFromXYZFile(molecule)
            DetermineBonds(self.mol,charge = 0)
        else:
            self.mol = molecule
    
    def getMol(self):
        return self.mol
    
    def removeAtom(self, atomNr):
        edFrag = Chem.EditableMol(self.mol)
        edFrag.RemoveAtom(atomNr)
        self.mol = edFrag.GetMol()
    
    def combineMolecules(mol1, mol2, index):  #Takes the rdkit mols and combines them
        """Function to combine two rdkit-mol objects at a specific point

        Args:
            mol1 (_type_): rdkit structure of molecule Nr.1
            mol2 (_type_): rdkit Structure of molecule Nr.2
            index (_type_): A tuple with the index of the atoms to connect (eg. (12, 24)). IMPORTANT: The order has to match the order of the given Molecules

        Returns:
            _type_: A New Molecule
        """
        combo = Chem.CombineMols(mol1, mol2)  #Combine
        edcombo = Chem.EditableMol(combo)     #Make editabel
        a1, a2 = index
        edcombo.AddBond(a1, mol1.GetNumAtoms()+ a2, order=Chem.rdchem.BondType.SINGLE)  #Combine the molecules with a bond
        combinedMol = edcombo.GetMol()  #Get the finished structure
        AllChem.EmbedMolecule(combinedMol)  #Somewhat relax the structure to make a belivable Molecule
        return Molecule(combinedMol)
    
    def Mol2COM(self):
        from openbabel import openbabel
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "com")

        molFile =Chem.MolToMolBlock(self.mol)
        comFile = openbabel.OBMol()
        obConversion.ReadString(comFile, molFile)
        return obConversion.WriteString(comFile)
    
    def writeCOMFile(self, fileName, head,coords,tail = "",comment="COMMENT"):
        with open(fileName, "w") as file:
            file.write(head)
            file.write(f"\n\n{comment}\n\n")
            for line in coords:
                file.write(f"{line}\n")
            file.write(tail)
            