
INTRODUCTION
------------

P-MOST-PDB-Map is software that uses Protein Data Bank structure files to 
calculate optimal backmapping parameters for the Protein Model with Oriented 
SiTes (P-MOST), compute the root-mean-square displacement for each atom type,
and create molecular image files describing the mapping and backmapping.

LICENSE
-------

P-MOST-PDB-Map is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

P-MOST-PDB-Map is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
P-MOST-PDB-Map.  If not, see <http://www.gnu.org/licenses/>.

PUBLICATION
-----------

Details of the oriented coarse-grained protein model appear in my paper:

[1] T. K. Haxton, "High-resolution coarse-grained modeling using oriented 
coarse-grained sites," submitted (available online at arxiv.org/abs/1409.8658).

AUTHORS, CITATION, AND CONTACT INFORMATION
-------------------------------------------

P-MOST-PDB-Map was developed by Tom Haxton at the Molecular Foundry, Lawrence 
Berkeley National Laboratory.

I ask that users please cite my publication [1] in any publication presenting 
results produced with P-MOST-PDB-Map.

The latest version of P-MOST-PDB-Map can be found at ...

Contact: Tom Haxton (tomhaxton@gmail.com)

PLATFORM
--------

P-MOST-PDB-Map has been tested on Linux and OS X.

REQUIREMENTS
------------

1. C compiler (e.g. gcc) to compile source code
2. PDB structure files (available from www.rcsb.org)
3. VMD to view molecular files (optional)

INSTALLATION
------------

Run "make" in the directory "code".

CONTENTS/USAGE
--------------

roundtrip

   Performs a least-squares fit to minimize the rmsd between original and
   backmapped atomic configurations, given a set of Protein Data Bank (PDB)
   files, and outputs the optimal model parameters and root-mean-square-
   displacement by atom type.

   Options:
      Several parameter filenames. Defaults are set so that the code needs no 
         arguments if run from any second-level directory (e.g. /code or 
         /scripts).
      -top_directory: Top of the directory tree containing PDB files to be
         analyzed.
      -base: Base of filename (including relative filename) for output files.
      -discard: Boolean for whether to discard (1 or default) or not discard
         (0) atoms with irregular geometries

   Outputs:
      <base>.rmsd: Average position in the frame of the site (before applying
         any linear correction), component-wise root-mean-square displacement
         (after applying linear corrections), root-mean-square displacement,
         count, count if all residues were complete, and non-Gaussian parameter,
         listed for each atom type
      <base>.corrections.<RESIDUE>.<ATOM1>.<ATOM2>: Model parameters for ATOM1 
         in RESIDUE whose backmapping function is a linear function of the
         predicted position of ATOM2

      <base>.Nterminus_error: List of N-terminal hydrogen atoms that were
         removed because only one hydrogen was found 
      <base>.entry_average_file: Heavy-atom rmsd, non-Gaussian parameter, and
         count, and hydrogen rmsd, non-Gaussian parameter, and count, listed by
         PBD entry
      <base>.error_files.ensemble: List of PDB entries not used because they
         contained an ensemble of models
      <base>.error_files.no_DBREF: List of PDB entries not used because they
         contained no DBREF entry
      <base>.error_files.non_amino_acid: List of PDB entries not used because 
         they contained residues other than the 20 natural amino acids
      <base>.error_files.residue_indicies_out_of_order: List of PDB entries not 
         used because they contained residue indicies out of order within a 
         single chain
      <base>.error_files.terminal_atom_in_non_terminal_residue: List of PDB 
         entries not used because they contained terminal atom types in non-
         terminal positions
      <base>.error_files.unknown_atom: List of PDB entries not used because they
         contained an unrecognized atom type
      <base>.file_statistics: List of number of entries excluded by reason
      <base>.indistinguishable atoms: List of groups of indistinguishable 
         hydrogen atoms that could not be renamed according to their dihedral 
         angles
      <base>.large_disagreement_file: List of atoms with roundtrip displacements
         exceeding a threshold (default 1 angstrom)
      <base>.moments: Average position in the frame of the site, component-wise 
         root-mean-square displacement, root-mean-square displacement, and count
         for each atom type, before applying any linear corrections and listed
         for each atom type
      <base>.overlapping_atoms_error: List of pairs of atoms excluded because
         both atoms were found at the same position

   Outputs when discard option is on:
      <base>.bondangle_error: Atoms discarded due to irregular bond angles
      <base>.bondlength_error: Atoms discarded due to irregular bond lengths
      <base>.tetrahedra_error: List of groups of atoms excluded because of an
         irregular molecular geometry around a tetrahedral center

   Outputs (auxillary calculations on the PDB structures):
      <base>.Cterminus: Average and std dev bond lengths and angles for C-
         terminus C, O, OXT atoms
      <base>.radial: Average longitudinal and radial positions in the frame of
         three adjacent heavy atoms, for groups of terminal hydrogens

create_molecular_files

   Creates molecular image files describing the mapping and backmapping 
   functions

   Options:
      Several parameter filenames. Defaults are set so that the code needs no 
         arguments if run from any second-level directory (e.g. /code or 
         /scripts).
      -pdbfile: Input PDB structure file
      -input_base: Directory containing parameter files (outputs of roundtrip)
      -directory: Output directory
      -filebase: Base of output files

   Outputs:
      <directory>/<filebase>.pdb: PDB file matching the input file except that
         atom positions are replaced with the backmapped atom positions
      <directory>/<filebase>.<RESIDUE>.<TYPE>.tcl: TCL script to be used with 
         VMD (in TkConsole run "source <TCLFILE>") to render an image of the
         atoms and bonds, colored by coarse-grained site, in the first residue 
         of type RESIDUE in the PDB file.  TYPE denotes whether the atoms are
         placed in their original positions or their backmapped positions.
      <directory>/<filebase>.<RESIDUE>.<TYPE>.tcl: As above, and including
         rendering of vectors representing the coarse-grained sites 
      <directory>/<filebase>.<RESIDUE>.<TYPE>.xyz: Coordinates used in the
         above two TCL scripts

SCRIPTS
-------

Example scripts are in the directory "scripts."  They should be modified to
reflect the location of PDB files on your computer.

PRE-GENERATED OUTPUTS
---------------------

output/: Outputs generated by running "roundtrip" on the entire PDB, both with
and without the "discard" flag

molecular_files/: Outputs generated by running "create_molecular_files" on two 
large PDB structures with hydrogens

