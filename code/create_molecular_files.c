#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <err.h>
#include <fts.h>
#include "mystdlib.h"
#include "mymath.h"
#include "protein_functions.h"

int main(int argc, char *argv[]){
    int Nresidues, i, j, covariancecount=0, type, nrecords=0, terminiierror, errorfilefirst=0, N_sidechain_corrections, read, N_backbone_atom_types, backbone_Oindex, backbone_Nindex, backbone_Hindex, backbone_Cindex, discard, do_amino_acid_examples, N_overlapping_atoms;
    char **backbone_atom_types, **Cterminusatomnames;
	FILE *inp;
    amino_acid_struct *amino_acid_list;
    residuedata *residuearray;
	cgresidue *cgresiduearray;
    declare_array(char, amino_acid_directory, maxstringlength);
    declare_array(char, site_director_atoms_file, maxstringlength);
    declare_array(char, site_atoms_file, maxstringlength);
    declare_array(char, elements_file, maxstringlength);
    declare_array(char, Cterminus_atoms_file, maxstringlength);
    declare_array(char, change_atom_names_file, maxstringlength);
    declare_array(char, indistinguishable_atoms_file, maxstringlength);
    declare_array(char, radial_atoms_file, maxstringlength);
    declare_array(char, aromatic_hydrogens_file, maxstringlength);
    declare_array(char, corrections_parameter_file, maxstringlength);
    declare_array(char, tetrahedral_file, maxstringlength);
    declare_array(char, rings_file, maxstringlength);
    declare_array(char, covalent_bond_lengths_file, maxstringlength);
    declare_array(char, bond_angles_file, maxstringlength);
    
    declare_array(char, pdbfile, maxstringlength);
    declare_array(char, outputfile, maxstringlength);
    declare_array(char, input_base, maxstringlength);
    declare_array(char, directory, maxstringlength);
    declare_array(char, filebase, maxstringlength);
    declare_array(char, tetrahedral_error_file, maxstringlength);
    declare_array(char, rings_error_file, maxstringlength);
    declare_array(char, overlapping_atoms_error_file, maxstringlength);
    declare_array(char, bondlength_error_file, maxstringlength);
    declare_array(char, bondangle_error_file, maxstringlength);
    declare_array(char, Nterminus_error_file, maxstringlength);
    declare_array(char, indistinguishable_error_file, maxstringlength);
    declare_array(char, aromatic_hydrogens_error_file, maxstringlength);
    declare_array(char, large_disagreement_file_specific, maxstringlength);
    declare_array(char, entry_average_file, maxstringlength);
    char **sidechain_correction_identifiers;
    declare_array(int, Ntype, 9);
    declare_array_nozero(char, linechar, maxlinelength);
    
    //  Options
    
    discard=loadparam(argc, "-discard", argv, "1").i;
    do_amino_acid_examples=loadparam(argc, "-do_amino_acid_examples", argv, "1").i;

    //  Inputs
    
    amino_acid_directory=loadparam(argc, "-amino_acid_directory", argv, "../amino_acid_pdb_files").s;
    site_director_atoms_file=loadparam(argc, "-site_director_atoms_file", argv, "../parameter_files/site_director_atoms.txt").s;
    site_atoms_file=loadparam(argc, "-site_atoms_file", argv, "../parameter_files/site_atoms.txt").s;
    elements_file=loadparam(argc, "-elements_file", argv, "../parameter_files/elements.txt").s;
    Cterminus_atoms_file=loadparam(argc, "-Cterminus_atoms_file", argv, "../parameter_files/Cterminus_atoms.txt").s;
    change_atom_names_file=loadparam(argc, "-change_atom_names_file", argv, "../parameter_files/change_atom_names.txt").s;
    indistinguishable_atoms_file=loadparam(argc, "-indistinguishable_atoms_file", argv, "../parameter_files/indistinguishable_atoms.txt").s;
    radial_atoms_file=loadparam(argc, "-radial_atoms_file", argv, "../parameter_files/radial_atoms.txt").s;
    aromatic_hydrogens_file=loadparam(argc, "-aromatic_hydrogens_file", argv, "../parameter_files/aromatic_hydrogens.txt").s;
    corrections_parameter_file=loadparam(argc, "-corrections_parameter_file", argv, "../parameter_files/corrections_atoms.txt").s;
    if(discard==1){
        tetrahedral_file=loadparam(argc, "-tetrahedral_file", argv, "../parameter_files/tetrahedral_atoms.txt").s;
        rings_file=loadparam(argc, "-rings_file", argv, "../parameter_files/rings.txt").s;
        covalent_bond_lengths_file=loadparam(argc, "-covalent_bond_lengths_file", argv, "../parameter_files/bondlengths.txt").s;
        bond_angles_file=loadparam(argc, "-bond_angles_file", argv, "../parameter_files/bondangles.txt").s;
    }

    input_base=loadparam(argc, "-input_base", argv, "0").s;
    pdbfile=loadparam(argc, "-pdbfile", argv, "0").s;
    
    //  Outputs
    
    directory=loadparam(argc, "-directory", argv, "0").s;
    filebase=loadparam(argc, "-filebase", argv, "0").s;
    
    sprintf(tetrahedral_error_file, "%s/%s.tetrahedra_error", directory, filebase);
    delete_file_if_exists(tetrahedral_error_file);
    sprintf(rings_error_file, "%s/%s.ring_error_file", directory, filebase);
    delete_file_if_exists(rings_error_file);
    sprintf(overlapping_atoms_error_file, "%s.overlapping_atoms_error", filebase);
    delete_file_if_exists(overlapping_atoms_error_file);
    sprintf(bondlength_error_file, "%s/%s.bondlength_error", directory, filebase);
    delete_file_if_exists(bondlength_error_file);
    sprintf(bondangle_error_file, "%s/%s.bondangle_error", directory, filebase);
    delete_file_if_exists(bondangle_error_file);
    sprintf(Nterminus_error_file, "%s/%s.Nterminus_error", directory, filebase);
    delete_file_if_exists(Nterminus_error_file);
    sprintf(indistinguishable_error_file, "%s/%s.indistinguishable_atoms", directory, filebase);
    delete_file_if_exists(indistinguishable_error_file);
    sprintf(aromatic_hydrogens_error_file, "%s/%s.aromatic_hydrogens", directory, filebase);
    delete_file_if_exists(aromatic_hydrogens_error_file);
    sprintf(large_disagreement_file_specific, "%s/%s.large_disagreement_file_specific", directory, filebase);
    delete_file_if_exists(large_disagreement_file_specific);
    sprintf(entry_average_file, "%s/%s.entry_average_file", directory, filebase);
    delete_file_if_exists(entry_average_file);
    
    //  Create an array of amino acid structures characterizing the atoms in each amino acid
    
    create_amino_acid_struct(amino_acid_directory, site_director_atoms_file, site_atoms_file, corrections_parameter_file, elements_file, Cterminus_atoms_file, indistinguishable_atoms_file, radial_atoms_file, &amino_acid_list, &N_sidechain_corrections, &sidechain_correction_identifiers, &backbone_atom_types, &N_backbone_atom_types, &backbone_Oindex, &backbone_Nindex, &backbone_Hindex, &backbone_Cindex, &Cterminusatomnames, aromatic_hydrogens_file);
    change_atom_names_in_amino_acid_struct(amino_acid_list, change_atom_names_file, &N_backbone_atom_types, backbone_atom_types);
    
    //  Create and initialize structures for calculating corrections from covariances
    
    declare_array_nozero(correction_structure, carbonyl_correction_array, Naminoacids);
    declare_array_nozero(correction_structure, backboneH_correction_array, Naminoacids);
    for(i=0;i<Naminoacids;i++){
        initialize_correction_structure(&(carbonyl_correction_array[i]));
        initialize_correction_structure(&(backboneH_correction_array[i]));
    }
    declare_array_nozero(correction_structure, sidechain_corrections, N_sidechain_corrections);
    for(i=0;i<N_sidechain_corrections;i++){
        initialize_correction_structure(&(sidechain_corrections[i]));
    }
    
    //  Create structures for checking molecular geometries

    declare_array_nozero(atomlookup, overlapping_atoms_list, max_overlapping_atoms);
    for(i=0;i<max_overlapping_atoms;i++) allocate_array(char, (overlapping_atoms_list[i].filename), maxstringlength);
    int **foundbond, *Nangles, ***angle_index, *Ntetrahedra, ***tetrahedraindices, *Nrings, **ringsize, **number_dependent_atoms, ***ringatomindex;
    double **maxbondlength, **minbondlength, **angles;
    setup_bond_length_structure(covalent_bond_lengths_file, &foundbond, &maxbondlength, &minbondlength);
    setup_bond_angles_structure(bond_angles_file, amino_acid_list, &Nangles, &angles, &angle_index);
    setup_tetrahedra_structure(tetrahedral_file, amino_acid_list, &Ntetrahedra, &tetrahedraindices);
    setup_rings_structure(rings_file, amino_acid_list, &Nrings, &ringsize, &number_dependent_atoms, &ringatomindex);

    //  Input best-fit model
    
    sprintf(outputfile, "%s.moments", input_base);
    input_cg_coords(outputfile, amino_acid_list, Cterminusatomnames);
    for(i=0;i<Naminoacids;i++){
        sprintf(outputfile, "%s.corrections.%s. O  . N  .residue_specific", input_base, amino_acid_list[i].resname);
        input_correction(&(carbonyl_correction_array[i]), outputfile);
        sprintf(outputfile, "%s.corrections.%s. H  . C  .residue_specific", input_base, amino_acid_list[i].resname);
        input_correction(&(backboneH_correction_array[i]), outputfile);
    }
    for(i=0;i<N_sidechain_corrections;i++){
        sprintf(outputfile, "%s.corrections.%s", input_base, sidechain_correction_identifiers[i]);
        input_correction(&(sidechain_corrections[i]), outputfile);
    }
    
    type=read_pdb_to_residuearray(pdbfile, amino_acid_list, &residuearray, &Nresidues, 1);
    if(type>0) print_error(type);
    if(type==0){
        reassign_Nterminus_atoms(residuearray, Nresidues, amino_acid_list, Nterminus_error_file, pdbfile, 1, 1);
        terminiierror=assign_atomfound_terminii(Nresidues, residuearray, amino_acid_list);
        if(terminiierror==1) print_error(7);
        if(terminiierror==0){
            discard_overlapping_atoms(residuearray, Nresidues, amino_acid_list, overlapping_atoms_error_file, pdbfile, 1, overlapping_atoms_list, &N_overlapping_atoms);
            if(discard==1){
                check_covalent_bond_lengths(residuearray, Nresidues, amino_acid_list, foundbond, maxbondlength, minbondlength, bondlength_error_file, pdbfile, 1);
                check_bond_angles(residuearray, Nresidues, amino_acid_list, Nangles, angles, angle_index, bondangle_error_file, pdbfile, 1);
                discard_bad_tetrahedra(residuearray, Nresidues, Ntetrahedra, tetrahedraindices, amino_acid_list, tetrahedral_error_file, pdbfile, 1);
                check_rings(residuearray, Nresidues, Nrings, ringsize, number_dependent_atoms, ringatomindex, amino_acid_list, rings_error_file, pdbfile, 1);
            }
            resolve_indistinguishable_atoms(residuearray, Nresidues, amino_acid_list, 1, indistinguishable_error_file, pdbfile, 1);
            resolve_aromatic_hydrogens(residuearray, Nresidues, amino_acid_list, 1, aromatic_hydrogens_error_file, pdbfile, 1);
            map_to_cgmodel(Nresidues, residuearray, amino_acid_list, &cgresiduearray, 1);
            
            create_molecular_files_separate_corrections(pdbfile, directory, filebase, amino_acid_list, Nresidues, residuearray, cgresiduearray, carbonyl_correction_array, backboneH_correction_array, sidechain_corrections, do_amino_acid_examples);
            
            //  Calculate root mean square displacement between original all-atom configuration and all-atom configuration generated from coarse-grained representation
            
            calculate_rmsd_separate_corrections(amino_acid_list, residuearray, cgresiduearray, Nresidues, carbonyl_correction_array, backboneH_correction_array, sidechain_corrections, pdbfile, large_disagreement_file_specific, entry_average_file);
            
            free_cgresiduearray(&cgresiduearray, Nresidues);
        }
        //free_residuearray(&residuearray, Nresidues, 1);
    }
    sprintf(outputfile, "%s/%s.rmsd", directory, filebase);
    output_rmsd_onlyspecific(outputfile, amino_acid_list, nrecords, backbone_atom_types, N_backbone_atom_types, Cterminusatomnames);
	
    return 0;
}

