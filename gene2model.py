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
            new_acid.add_offset(0.13, 0, 0)
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
        box_size_x = "{:.3f}".format((2 * (xmax-xmin)))
        box_size_y = "{:.3f}".format((2 * (ymax-ymin)))
        box_size_z = "{:.3f}".format((2 * (zmax-zmin)))

        # Get column lengths so that we can format it in fixed width format
        column_lengths = [7 + len(str(len(self._acids)))] #first column: acid number and name
        column_lengths.append(7) #second column: atom name
        column_lengths.append((3 + len(str(num_atoms)))) #third column: atom sequence number
        column_lengths.append((3 + len("{:.3f}".format(max(abs(xmax), abs(xmin)))))) #fourth column: atom x coordinate
        column_lengths.append((3 + len("{:.3f}".format(max(abs(ymax), abs(ymin)))))) #fifth column: atom y coordinate
        column_lengths.append((3 + len("{:.3f}".format(max(abs(zmax), abs(zmin)))))) #sixth column: atom z coordinate

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
        beginning_atoms.accomodate_atom(self.head())
        self._atoms[1] = beginning_atoms.atom_at_index(0)
        self._atoms.insert(2, beginning_atoms.atom_at_index(1))
        self._atoms.insert(3, beginning_atoms.atom_at_index(2))
        
    def convert_to_last_acid(self, ending_atoms):
        """Last acid in protein chain has a different carboxyl group structure."""
        ending_atoms.change_name(self.name)
        ending_atoms.accomodate_atom(self.atom_at_index(-3))
        self._atoms[-1] = ending_atoms.atom_at_index(-2)
        self._atoms.append(ending_atoms.atom_at_index(-1))

    def accomodate_atom(self, atom):
        """Moves atoms in the acid so that it appears next to an atom."""
        for a in self._atoms:
            a.update_pos(atom)

    def accomodate(self, acid):
        """Moves atoms in the acid so that it appears at the end of the chain."""
        self.accomodate_atom(acid.atom_at_index(-2))

    def add_offset(self, x, y, z):
        """Adds a positional offset to each atom in the acid"""
        for atom in self._atoms:
            atom.add_offset([x, y, z])

    def to_gromacs(self, column_lengths, acid_num, atom_num):
        """Convert the acid to a series of gromacs text format lines."""
        gromacs = ""
        for atom in self._atoms:
            coords = atom.get_coordinates()
            xcoord = "{:.3f}".format(coords[0])
            ycoord = "{:.3f}".format(coords[1])
            zcoord = "{:.3f}".format(coords[2])
            gromacs += ' ' * (column_lengths[0] - len(self.name) - len(str(acid_num)))
            gromacs += str(acid_num) + self.name #residue name and number
            gromacs += ' ' * (column_lengths[1] - len(atom.get_name()))
            gromacs += atom.get_name() #atom name
            gromacs += ' ' * (column_lengths[2] - len(str(atom_num)))
            gromacs += str(atom_num) #atom number
            gromacs += ' ' * (column_lengths[3] - len(xcoord))
            gromacs += xcoord #x coordinate
            gromacs += ' ' * (column_lengths[4] - len(ycoord))
            gromacs += ycoord #y coordinate
            gromacs += ' ' * (column_lengths[5] - len(zcoord))
            gromacs += zcoord + '\n' #z coordinate
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
        """adds an offset to the atom's position based on another atom's position."""
        self.add_offset(atom._coordinates)

    def add_offset(self, xyz):
        """Adds a vector offset to the atom's position."""
        self._coordinates = self._coordinates + xyz

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

    def convert_acids_to_gene(self, acid_names, filepath):
        """Convert a given acid string to a .gro position model."""
        acids_list = [acid_names[i:i+3] for i in range(0, len(acid_names), 3)]
        acids = []
        for acid_name in acids_list:
            acids.append(copy.deepcopy(self.acid_table[acid_name]))
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

    def convert_codons_to_gene(self, codons, filepath):
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
        acid_names_list = [codon_to_acid[codon] for codon in codons_list]
        acid_names = ''.join(acid_names_list)
        self.convert_acids_to_gene(acid_names, filepath)

# Testing code

table = GeneToModelConverter("amino_acid_reference.txt")
table.convert_codons_to_gene("AAAGCCATG", "test1.gro")
table.convert_acids_to_gene("LYSVALPHEGLYARGCYSGLULEUALAALAALAMETLYSARGHISGLYLEUASPASNTYRARGGLYTYRSERLEUGLYASNTRPVALCYSALAALALYSPHEGLUSERASNPHEASNTHRGLNALATHRASNARGASNTHRASPGLYSERTHRASPTYRGLYILELEUGLNILEASNSERARGTRPTRPCYSASNASPGLYARGTHRPROGLYSERARGASNLEUCYSASNILEPROCYSSERALALEULEUSERSERASPILETHRALASERVALASNCYSALALYSLYSILEVALSERASPGLYASNGLYMETASNALATRPVALALATRPARGASNARGCYSLYSGLYTHRASPVALGLNALATRPILEARGGLYCYSARGLEU", "LYSOZYME.gro")