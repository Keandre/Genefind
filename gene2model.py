import numpy
import copy
from collections import namedtuple
from typing import List
from pathlib import Path
Vector = List[float]

class Gene:
    """Represents a sequence of RNA."""
    def __init__(self, acids):
        self._acids = []
        for acid in acids:
            self.push_acid(acid)

    def push_acid(self, new_acid):
        """Append a new amino acid to the gene sequence."""
        if self._acids:
            new_acid.accomodate(self._acids[-1])
        self._acids.append(new_acid)

    def to_gromacs(self, title):
        """Convert the gene to a gromacs text format string."""

        # Get number of atoms and box vector information
        num_atoms = 0
        xmax = ymax = zmax = 0
        xmin = ymin = zmin = 1000000
        for acid in self._acids:
            num_atoms += acid.num_atoms()
            xmax = max(xmax, acid.max_x())
            ymax = max(ymax, acid.max_y())
            zmax = max(zmax, acid.max_z())
            xmin = min(xmin, acid.min_x())
            ymin = min(ymin, acid.min_y())
            zmin = min(zmin, acid.min_z())
        box_size_x = str((2 * (xmax-xmin)))
        box_size_y = str((2 * (ymax-ymin)))
        box_size_z = str((2 * (zmax-zmin)))

        # Get column lengths so that we can format it in fixed width format
        column_lengths = [5 + len(str(len(self._acids)))] #first column: acid number and name
        column_lengths.append(7) #second column: atom name
        column_lengths.append(len(str(num_atoms))) #third column: atom sequence number
        column_lengths.append((3 + len(str(max(abs(xmax), abs(xmin)))))) #fourth column: atom x coordinate
        column_lengths.append((3 + len(str(max(abs(ymax), abs(ymin)))))) #fifth column: atom y coordinate
        column_lengths.append((3 + len(str(max(abs(zmax), abs(zmin)))))) #sixth column: atom z coordinate

        #add title to .gro string
        gromacs = title + "\n "
        #add number of atoms to .gro string
        gromacs += str(num_atoms) + "\n"
        #add all acid atom lines to .gro string
        acid_num = 0
        atom_num = 1
        for acid in self._acids:
            acid_num += 1
            gromacs += acid.to_gromacs(column_lengths, acid_num, atom_num)
            atom_num += acid.num_atoms()
        #add box vector to .gro string
        gromacs += "   " + box_size_x + "   " + box_size_y + "   " + box_size_z

        return gromacs
        
class AminoAcid:
    """Represents a series of atoms making up an amino acid."""
    def __init__(self, name):
        self.name = name
        self._atoms = []

    def __repr__(self):
        return "<{} Acid>".format(self.name)
                
    def head(self):
        """Return the 'head' atom of the acid."""
        return self._atoms[0]

    def tail(self):
        """Return the 'tail' atom of the acid."""
        return self._atoms[-1]

    def get_name(self):
        return self.name

    def change_name(self, name):
        self.name = name
        for atom in self._atoms:
            atom.change_acid(name)

    def atom_at_index(self, index):
        """return atom at indexed position."""
        return self._atoms[index]

    def num_atoms(self):
        """Return the number of atoms in the acid."""
        return len(self._atoms)

    def max_x(self):
        """Return the largest x coordinate"""
        return max((atom.get_coordinates()[0]) for atom in self._atoms)

    def max_y(self):
        """Return the largest y coordinate"""
        return max((atom.get_coordinates()[1]) for atom in self._atoms)

    def max_z(self):
        """Return the largest z coordinate"""
        return max((atom.get_coordinates()[2]) for atom in self._atoms)

    def min_x(self):
        """Return the smallest x coordinate"""
        return min((atom.get_coordinates()[0]) for atom in self._atoms)

    def min_y(self):
        """Return the smallest y coordinate"""
        return min((atom.get_coordinates()[1]) for atom in self._atoms)

    def min_z(self):
        """Return the smallest z coordinate"""
        return min((atom.get_coordinates()[2]) for atom in self._atoms)

    def push_atom(self, atom):
        self._atoms.append(atom)

    def convert_to_first_acid(self, beginning_atoms):
        """First acid in protein chain has a different amino group structure."""
        beginning_atoms.change_name(self.name)
        self._atoms[1] = beginning_atoms.atom_at_index(0)
        self._atoms.insert(2, beginning_atoms.atom_at_index(1))
        self._atoms.insert(3, beginning_atoms.atom_at_index(2))
        self._atoms[1].update_pos(self.head())
        self._atoms[2].update_pos(self.head())
        self._atoms[3].update_pos(self.head())

    def convert_to_last_acid(self, ending_atoms):
        """Last acid in protein chain has a different carboxyl group structure."""
        ending_atoms.change_name(self.name)
        self._atoms[-1] = ending_atoms.atom_at_index(-2)
        self._atoms.append(ending_atoms.atom_at_index(-1))
        self._atoms[-2].update_pos(self.head())
        self._atoms[-1].update_pos(self.head())
    
    def accomodate(self, acid):
        """Moves atoms in the acid so that it appears at the end of the chain."""
        for atom in self._atoms:
            atom.update_pos(acid.tail())

    def to_gromacs(self, column_lengths, acid_num, atom_num):
        """Convert the acid to a series of gromacs text format lines."""
        gromacs = ""
        for atom in self._atoms:
            gromacs += ' ' * (column_lengths[0] - len(self.name) - len(str(acid_num)))
            gromacs += str(acid_num) + self.name #residue name and number
            gromacs += ' ' * (column_lengths[1] - len(atom.get_name()))
            gromacs += atom.get_name() #atom name
            gromacs += ' ' * (column_lengths[2] - len(str(atom_num)))
            gromacs += str(atom_num) #atom number
            gromacs += ' ' * (column_lengths[3] - len(str(atom.get_coordinates()[0])))
            gromacs += str(atom.get_coordinates()[0]) #x coordinate
            gromacs += ' ' * (column_lengths[4] - len(str(atom.get_coordinates()[1])))
            gromacs += str(atom.get_coordinates()[1]) #y coordinate
            gromacs += ' ' * (column_lengths[5] - len(str(atom.get_coordinates()[2])))
            gromacs += str(atom.get_coordinates()[2]) + '\n' #z coordinate
            atom_num += 1
        return gromacs
            
class Atom:
    """Represents an atom."""
    def __init__(self, name: str, amino_acid: AminoAcid,
                 xyz: Vector) -> None:
        """

        name:  Not sure what this is yet, denotes some kind of connection
        between atoms.
        amino_acid:  The parent amino acid to which this atom belongs.
        xyz:  The xyz coordinates of the atom.
        """
        self._name = name
        self.amino_acid = amino_acid
        self._coordinates = numpy.array(xyz)

    def __repr__(self):
        if self.amino_acid:
            return "<{} Atom inside {} Acid>".format(self._name,
                                                      self.amino_acid.name)
        else:
            return "<{} Atom>".format(self._name)
    
    def update_pos(self, atom):
        self._coordinates = self._coordinates + atom._coordinates

    def update_pos_pdb(self, atom):
        self._coordinates = self._coordinates - atom._coordinates

    def get_name(self):
        return self._name

    def get_coordinates(self):
        return self._coordinates

    def change_acid(self, acid):
        self.amino_acid = acid

class GeneToModelConverter:
    """Conversion class to turn genefind codon strings into .gro atomic model files.

    The class requires a reference file containing all the amino acids to compile
    a table from."""
    def __init__(self, filepath):
        reference_file = open(filepath)

        self.acid_table = self._extract_acid_table(reference_file)
        

    def _extract_acid_table(self, reference_file):
        """Extract the amino acid to atomic position table from a provided reference file."""
        if type(reference_file) == str:
            reference_file = open(reference_file)
        lines = reference_file.readlines()
        acid_table = {}
        reference_atom = namedtuple('reference_atom', ['acid', 'name', 'x', 'y', 'z',])

        acid_ranges = []
        current_range = []
        last_acid = ""
        # Loop through each line of the file and construct a range of atoms for each acid
        for line in lines: 
            atom = reference_atom._make(line.split())
            if atom.acid != last_acid:
                acid_ranges.append(current_range)
                current_range = []
            current_range.append(atom)
            last_acid = atom.acid

        # Append the final acid range
        acid_ranges.append(current_range)
        # hacky: Delete errant empty range appearing in beginning
        del acid_ranges[0]
        # build the acids and add them to the acid table
        for acid_range in acid_ranges:
            acid = self._build_reference_acid(acid_range)
            acid_table[acid.name] = acid

        return acid_table

    def _build_reference_acid(self, atoms):
        acid = AminoAcid(atoms[0].acid)
        atoms = [Atom(atom.name, acid,
                      [float(atom.x), float(atom.y), float(atom.z)])
                 for atom in atoms]
        for atom in atoms:
            acid.push_atom(atom)
        return acid

    def gene_convert(self, codons, filepath):
        """Convert a given codon string to a .gro position model."""
        codon_to_acid = {
            "AAA":"LYS", "AAT":"ASN", "AAG":"LYS",
            "AAC":"ASN", "ATA":"LLE", "ATT":"ILE",
            "ATG":"MET", "ATC":"ILE", "AGA":"ARG",
            "AGT":"SER", "AGG":"ARG", "AGC":"SER",
            "ACA":"THR", "ACT":"THR", "ACG":"THR",
            "ACC":"THR", "TAA":"END", "TAT":"TYR",
            "TAG":"END", "TAC":"TYR", "TTA":"LEU",
            "TTT":"PHE", "TTG":"LEU", "TTC":"PHE",
            "TGA":"END", "TGT":"CYS", "TGG":"TRP",
            "TGC":"CYS", "TCA":"SER", "TCT":"SER",
            "TCG":"SER", "TCC":"SER", "GAA":"GLU",
            "GAT":"ASP", "GAG":"GLU", "GAC":"ASP",
            "GTA":"VAL", "GTT":"VAL", "GTG":"VAL",
            "GTC":"VAL", "GGA":"GLY", "GGT":"GLY",
            "GGG":"GLY", "GGC":"GLY", "GCA":"ALA",
            "GCT":"ALA", "GCG":"ALA", "GCC":"ALA",
            "CAA":"GLN", "CAT":"HIS", "CAG":"GLN",
            "CAC":"HIS", "CTA":"LEU", "CTT":"LEU",
            "CTG":"LEU", "CTC":"LEU", "CGA":"ARG",
            "CGT":"ARG", "CGG":"ARG", "CGC":"ARG",
            "CCA":"PRO", "CCT":"PRO", "CCG":"PRO",
            "CCC":"PRO"}
        
        codons_list = [codons[i:i+3] for i in range(0, len(codons), 3)]
        acids = [self.acid_table[codon_to_acid[codons_list[i]]] for i in range(0, len(codons_list))]

        end_index = -1
        for i in range(0, len(acids)):
            if acids[i].get_name() == "END":
                end_index = i
                break
        if end_index >= 0:
            del acids[end_index:len(acids)]

        acids[0].convert_to_first_acid(self.acid_table['BEG'])
        acids[-1].convert_to_last_acid(self.acid_table['END'])

        gene = Gene(acids)
        path = Path(filepath)
        gromacs = gene.to_gromacs(path.stem)
        file = open(filepath, 'w')
        file.write(gromacs)
        file.close()
    
# old pdb methods for making acid table

    def _extract_acid_table_from_pdb(self, pdbfile):
        """Extract the amino acid to atomic position table from a provided pdb file.

        It is important that the pdb file contain all the amino acids for this
        to work."""
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
                             [0,0,0])
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
            acid = self._build_reference_acid_from_pdb(pdb_atoms, previous_tail)
            previous_tail = acid.tail()
            acid_table[acid.name] = acid

    def _build_reference_acid_from_pdb(self, pdb_atoms, previous_tail):
        acid = AminoAcid(pdb_atoms[0].acid)
        atoms = [Atom(atom.atom, acid,
                      [float(atom.x), float(atom.y), float(atom.z)])
                 for atom in pdb_atoms]
        for atom in atoms:
            atom.update_pos_pdb(previous_tail)
            acid.push_atom(atom)
        return acid

# Testing code

table = GeneToModelConverter("amino_acid_reference.txt")
table.gene_convert("AAAGTCCGG", "test1.gro")