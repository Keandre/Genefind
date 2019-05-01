import numpy
import copy
from collections import namedtuple
from typing import List
Vector = List[float]

def build_acid_table(pdbfile):
    if type(pdbfile) == str:
        pdbfile = open(pdbfile)
    lines = pdbfile.readlines()
    acid_table = {}
    PDBAtom = namedtuple('PDBAtom', ['record_name', 'serial', 'name',
                                     'residue_name', 'chain_id', 'residue_sequence',
                                     'x', 'y', 'z',
                                     'occupancy', 'temp_factor', 'charge'])
    model = False # Whether we're currently reading lines from a model
    previous_tail = Atom('origin', AminoAcid('placeholder'),
                         [0,0,0], 'H')
    acid_ranges = []
    current_range = []
    for line in enumerate(lines):
        if line[1].split()[0] == 'MODEL':
            model = True
            continue
        elif line[1].split()[0] == 'ENDMDL':
            model = False
            continue
        
        if model:
            try:
                atom = PDBAtom._make(line[1].split())
            except TypeError:
                if line[1].split()[0] == 'TER':
                    continue
                
        if model and (atom.charge == 'N') and atom.residue_name not in acid_table:
            current_range.append(line[0])
            if len(current_range) == 2:
                acid_ranges.append(current_range)
                current_range = []

    for line_range in acid_ranges:
        # TODO: Skip ones we've already created to speed this up
        acid_lines = lines[line_range[0]:line_range[1]]
        pdb_atoms = [PDBAtom._make(line.split()) for line in acid_lines
                     if line.split()[0] == 'ATOM']
        acid = build_reference_acid(pdb_atoms, previous_tail)
        previous_tail = acid.tail()
        acid_table[acid.name] = acid

    return acid_table
        
def build_reference_acid(pdb_atoms, previous_tail):
    acid = AminoAcid(pdb_atoms[0].residue_name)
    atoms = [Atom(atom.name, acid,
                  [float(atom.x), float(atom.y), float(atom.z)],
                  atom.charge)
             for atom in pdb_atoms]
    for atom in atoms:
        atom.update_pos(previous_tail)
        acid.push_atom(atom)
    return acid

class Gene:
    """Represents a sequence of RNA."""
    def __init__(self, codons):
        self._acids = []
        # Bla bla stuff that parses codons and turns them into acids

    def push_acid(self, new_acid):
        """Append a new amino acid to the gene sequence."""
        for acid in self._acids:
            acid.accomodate(new_acid)

    def to_gromacs(self):
        """Convert the gene to a gromacs text format string."""
        gromacs = ""
        for acid in self._acids:
            gromacs += acid.to_gromacs()
        return gromacs
        
class AminoAcid:
    """Represents a series of atoms making up an amino acid."""
    def __init__(self, name):
        self.name = name
        self._atoms = []
        
    def head(self):
        """Return the 'head' atom of the acid."""
        return self._atoms[0]

    def tail(self):
        """Return the 'tail' atom of the acid."""
        return self._atoms[-1]

    def push_atom(self, atom):
        self._atoms.append(atom)
    
    def accomodate(self, acid):
        """Moves atoms in the acid to accomodate the addition of a new acid to 
        the chain."""
        for atom in self._atoms:
            atom.update_pos(acid.head())

    def to_gromacs(self):
        """Convert the acid to a series of gromacs text format lines."""
        for atom in self._atoms:
            pass # TODO: Actually make this output to gromacs
            
class Atom:
    """Represents an atom."""
    def __init__(self, name: str, amino_acid: AminoAcid,
                 xyz: Vector, charge: str) -> None:
        """

        bond1:  Not sure what this is yet, denotes some kind of connection
        between atoms.
        amino_acid:  The parent amino acid to which this atom belongs.
        xyz:  The xyz coordinates of the atom.
        charge: The atoms charge.
        """
        self._name = name
        self.amino_acid = amino_acid
        self._coordinates = numpy.array(xyz)
        self._charge = charge
        
    def update_pos(self, atom):
        self._coordinates = self._coordinates - atom._coordinates

class GeneToModelConverter:
    """Conversion class to turn genefind codon strings into .gro atomic model files.

    The class currently requires a .pdb containing all the amino acids to compile 
    a table from."""
    def __init__(self, filepath):
        pdb_file = open(filepath)

        self.acid_table = self._extract_acid_table(pdb_file)
        
    def _extract_acid_table(self, infile):
        """Extract the amino acid to atomic position table from a provided pdb file.

        It is important that the pdb file contain all the amino acids for this 
        to work."""
        pass
        
    def gene_convert(self, codons):
        """Convert a given codon string(?) to a .gro position model."""
        
# Testing code

# table = build_acid_table("1akp.pdb")
