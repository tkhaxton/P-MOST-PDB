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
    int Nresidues, i, j, covariancecount=0, type, nrecords=0, terminiierror, errorfilefirst=0, N_sidechain_corrections, N_backbone_atom_types, backbone_Oindex, backbone_Nindex, backbone_Hindex, backbone_Cindex, first_indistinguishable_xyz=1, discard, print_indistinguishable, printboth, N_overlapping_atoms=0, Cterminus_count=0;
    char **backbone_atom_types, **Cterminusatomnames;
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

    declare_array(char, top_directory, maxstringlength);
    declare_array(char, myfile, maxstringlength);
    declare_array(char, outputfile, maxstringlength);
    declare_array(char, base, maxstringlength);
    declare_array(char, error_filename_base, maxstringlength);
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
    declare_array(char, ext, maxstringlength);
    char **sidechain_correction_identifiers;
    declare_array(int, Ntype, 9);
    
    //  Options
    
    discard=loadparam(argc, "-discard", argv, "1").i;
    print_indistinguishable=loadparam(argc, "-print_indistinguishable", argv, "0").i;
    printboth=loadparam(argc, "-printboth", argv, "0").i;

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

    top_directory=loadparam(argc, "-top_directory", argv, "0").s;
    
    //  Outputs
    
    base=loadparam(argc, "-base", argv, "0").s;
    sprintf(error_filename_base, "%s.error_files", base);
    
    sprintf(tetrahedral_error_file, "%s.tetrahedra_error", base);
    delete_file_if_exists(tetrahedral_error_file);
    sprintf(rings_error_file, "%s.ring_error_file", base);
    delete_file_if_exists(rings_error_file);
    sprintf(overlapping_atoms_error_file, "%s.overlapping_atoms_error", base);
    delete_file_if_exists(overlapping_atoms_error_file);
    sprintf(bondlength_error_file, "%s.bondlength_error", base);
    delete_file_if_exists(bondlength_error_file);
    sprintf(bondangle_error_file, "%s.bondangle_error", base);
    delete_file_if_exists(bondangle_error_file);
    sprintf(Nterminus_error_file, "%s.Nterminus_error", base);
    delete_file_if_exists(Nterminus_error_file);
    sprintf(indistinguishable_error_file, "%s.indistinguishable_atoms", base);
    delete_file_if_exists(indistinguishable_error_file);
    sprintf(aromatic_hydrogens_error_file, "%s.aromatic_hydrogens", base);
    delete_file_if_exists(aromatic_hydrogens_error_file);
    sprintf(large_disagreement_file_specific, "%s.large_disagreement_file_specific", base);
    delete_file_if_exists(large_disagreement_file_specific);
    sprintf(entry_average_file, "%s.entry_average_file", base);
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
    if(discard==1){
        setup_bond_length_structure(covalent_bond_lengths_file, &foundbond, &maxbondlength, &minbondlength);
        setup_bond_angles_structure(bond_angles_file, amino_acid_list, &Nangles, &angles, &angle_index);
        setup_tetrahedra_structure(tetrahedral_file, amino_acid_list, &Ntetrahedra, &tetrahedraindices);
        setup_rings_structure(rings_file, amino_acid_list, &Nrings, &ringsize, &number_dependent_atoms, &ringatomindex);
    }

    //  Create structures for calculating direct properties of PDB

    declare_matrix(double, Cterminus_avs, 3, 2);
    for(i=0;i<3;i++) for(j=0;j<2;j++) Cterminus_avs[i][j]=0;

    //  Look through file hierarchy downstream from top_directory
    //  Using code from Giorgios Keramidas, http://keramida.wordpress.com/2009/07/05/fts3-or-avoiding-to-reinvent-the-wheel/
    
    FTS *ftsp;
    FTSENT *p, *chp;
    int fts_options = FTS_COMFOLLOW | FTS_LOGICAL | FTS_NOCHDIR;
	char *paths[]={top_directory, NULL};
	if ((ftsp = fts_open(paths, fts_options, NULL)) == NULL) {
        warn("fts_open");
        return -1;
    }
    chp = fts_children(ftsp, 0);
    if (chp == NULL) {
        printf("No files found\n");
        return 0;               /* no files to traverse */
    }

    //  Loop through files once to calculate first and second moments

    while ((p = fts_read(ftsp)) != NULL) {
        switch (p->fts_info) {
            case FTS_D:
                break;
            case FTS_F:
                myfile=(*p).fts_path;
                ext = strrchr(myfile, '.');
                if (!ext) {
                    /* no extension */
                }
                else {
                    if((strcmp(ext+1, "pdb")==0)||(strcmp(ext+1, "ent")==0)||(strcmp(ext+1, "brk")==0)){
                        printf("1: %s\n", myfile);
                        
                        //  Read data from a .pdb file and convert to an array of residue structures containing atomic coordinates and connectivity information
                        //  Returns code indicating the type of molecule (1: protein consisting entirely of natural amino acids)
                        
                        type=read_pdb_to_residuearray(myfile, amino_acid_list, &residuearray, &Nresidues, 0);
                        Ntype[type]++;
                        if(type>0) output_error_filename(myfile, error_filename_base, type, &errorfilefirst);
                        if(type==0){
                            nrecords++;
                            
                            //  Check that right- and left-terminal atom types are found only where expected
                            
                            reassign_Nterminus_atoms(residuearray, Nresidues, amino_acid_list, Nterminus_error_file, myfile, 1, 0);
                            terminiierror=assign_atomfound_terminii(Nresidues, residuearray, amino_acid_list);
                            if(terminiierror==1){
                                Ntype[0]--;
                                Ntype[7]++;
                                output_error_filename(myfile, error_filename_base, 7, &errorfilefirst);
                            }
                            else{
                                discard_overlapping_atoms(residuearray, Nresidues, amino_acid_list, overlapping_atoms_error_file, myfile, 1, overlapping_atoms_list, &N_overlapping_atoms);
                                if(discard==1){
                                    check_covalent_bond_lengths(residuearray, Nresidues, amino_acid_list, foundbond, maxbondlength, minbondlength, bondlength_error_file, myfile, 1);
                                    check_bond_angles(residuearray, Nresidues, amino_acid_list, Nangles, angles, angle_index, bondangle_error_file, myfile, 1);
                                    discard_bad_tetrahedra(residuearray, Nresidues, Ntetrahedra, tetrahedraindices, amino_acid_list, tetrahedral_error_file, myfile, 1);
                                    check_rings(residuearray, Nresidues, Nrings, ringsize, number_dependent_atoms, ringatomindex, amino_acid_list, rings_error_file, myfile, 1);
                                }
                                if(print_indistinguishable==1) print_indistinguishable_xyz_files(amino_acid_list, residuearray, Nresidues, base, &first_indistinguishable_xyz);
                                resolve_indistinguishable_atoms(residuearray, Nresidues, amino_acid_list, 1, indistinguishable_error_file, myfile, 0);
                                resolve_aromatic_hydrogens(residuearray, Nresidues, amino_acid_list, 1, aromatic_hydrogens_error_file, myfile, 0);
                                
                                //  Directly calculate properties of PDB configurations
                                
                                calculate_average_radial_distances(amino_acid_list, Nresidues, residuearray);
                                calculate_Cterminus_bond_lengths_and_angles(amino_acid_list, Nresidues, residuearray, Cterminus_avs, &Cterminus_count);

                                //  Map all-atom configuration (residuearray) to coarse-grained configuration (cgresiduearray)
                                
                                map_to_cgmodel(Nresidues, residuearray, amino_acid_list, &cgresiduearray, 0);
                                
                                //  Average first and second moments of atomic coordinates relative to oriented coarse-grained configuration
                                
                                record_moments_relative_to_oriented_cg_sites(Nresidues, residuearray, amino_acid_list, cgresiduearray);
                                free_cgresiduearray(&cgresiduearray, Nresidues);
                            }
                            free_residuearray(&residuearray, Nresidues, 0);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    fts_close(ftsp);

    if(print_indistinguishable==1){
        add_header_to_xyz_files(base, amino_acid_list);
        return;
    }

    sprintf(outputfile, "%s.radial", base);
    average_and_output_radial_positions(outputfile, amino_acid_list);
    sprintf(outputfile, "%s.Cterminus", base);
    output_Cterminus_bond_lengths_and_angles(outputfile, Cterminus_avs, Cterminus_count);

    //  Calculate and output averages and standard deviations from first and second moments (moments stored in amino_acid_list)
    
    sprintf(outputfile, "%s.moments", base);
    average_and_output_moments_relative_to_oriented_cg_sites(outputfile, amino_acid_list, backbone_atom_types, N_backbone_atom_types, Cterminusatomnames);
    
    //  Loop through files a second time to calculate covariance, relative to average positions and orientations of sites
        
    ftsp = fts_open(paths, fts_options, NULL);
    while ((p = fts_read(ftsp)) != NULL) {
        switch (p->fts_info) {
            case FTS_D:
                break;
            case FTS_F:
                myfile=(*p).fts_path;
                ext = strrchr(myfile, '.');
                if (!ext) {
                    /* no extension */
                }
                else {
                    if((strcmp(ext+1, "pdb")==0)||(strcmp(ext+1, "ent")==0)||(strcmp(ext+1, "brk")==0)){
                        printf("2: %s\n", myfile);
                        type=read_pdb_to_residuearray(myfile, amino_acid_list, &residuearray, &Nresidues, 0);
                        if(type==0){
                            reassign_Nterminus_atoms(residuearray, Nresidues, amino_acid_list, Nterminus_error_file, myfile, 0, 0);
                            terminiierror=assign_atomfound_terminii(Nresidues, residuearray, amino_acid_list);
                            if(terminiierror==0){
                                discard_overlapping_atoms(residuearray, Nresidues, amino_acid_list, overlapping_atoms_error_file, myfile, 0, overlapping_atoms_list, &N_overlapping_atoms);
                                if(discard==1){
                                    check_covalent_bond_lengths(residuearray, Nresidues, amino_acid_list, foundbond, maxbondlength, minbondlength, bondlength_error_file, myfile, 0);
                                    check_bond_angles(residuearray, Nresidues, amino_acid_list, Nangles, angles, angle_index, bondangle_error_file, myfile, 0);
                                    discard_bad_tetrahedra(residuearray, Nresidues, Ntetrahedra, tetrahedraindices, amino_acid_list, tetrahedral_error_file, myfile, 0);
                                    check_rings(residuearray, Nresidues, Nrings, ringsize, number_dependent_atoms, ringatomindex, amino_acid_list, rings_error_file, myfile, 0);
                                }
                                resolve_indistinguishable_atoms(residuearray, Nresidues, amino_acid_list, 0, indistinguishable_error_file, myfile, 0);
                                resolve_aromatic_hydrogens(residuearray, Nresidues, amino_acid_list, 0, aromatic_hydrogens_error_file, myfile, 0);
                                map_to_cgmodel(Nresidues, residuearray, amino_acid_list, &cgresiduearray, 0);
                                
                                //  Average covariance between carbonyl Oxygen positions and predicted Nitrogen position of next residue
                                
                                //calculate_covariance(Nresidues, residuearray, amino_acid_list, cgresiduearray, &carbonyl_correction, &backboneH_correction);
                                calculate_covariance_separate(Nresidues, residuearray, amino_acid_list, cgresiduearray, carbonyl_correction_array, backboneH_correction_array);
                                calculate_covariance_sidechain(Nresidues, residuearray, amino_acid_list, cgresiduearray, sidechain_corrections);

                                free_cgresiduearray(&cgresiduearray, Nresidues);
                            }
                            free_residuearray(&residuearray, Nresidues, 0);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    fts_close(ftsp);

    //  Calculate linear correction based on (carbonyl Oxygen)-(predicted next Nitrogen) covariance; this accounts for cis-trans dihedral rotation

    for(i=0;i<Naminoacids;i++){
        sprintf(outputfile, "%s.corrections.%s. O  . N  .residue_specific", base, amino_acid_list[i].resname);
        calculate_correction_from_covariance(&(carbonyl_correction_array[i]), outputfile);
        sprintf(outputfile, "%s.corrections.%s. H  . C  .residue_specific", base, amino_acid_list[i].resname);
        calculate_correction_from_covariance(&(backboneH_correction_array[i]), outputfile);
    }
    for(i=0;i<N_sidechain_corrections;i++){
        //printf("%i: %s\n", i, sidechain_correction_identifiers[i]);
        sprintf(outputfile, "%s.corrections.%s", base, sidechain_correction_identifiers[i]);
        calculate_correction_from_covariance(&(sidechain_corrections[i]), outputfile);
    }
    
    //  Loop through files a third time to calculate root mean square displacement between original all-atom configuration and all-atom configuration generated from coarse-grained representation
    
    ftsp = fts_open(paths, fts_options, NULL);
    while ((p = fts_read(ftsp)) != NULL) {
        switch (p->fts_info) {
            case FTS_D:
                break;
            case FTS_F:
                myfile=(*p).fts_path;
                ext = strrchr(myfile, '.');
                if (!ext) {
                    /* no extension */
                }
                else {
                    if((strcmp(ext+1, "pdb")==0)||(strcmp(ext+1, "ent")==0)||(strcmp(ext+1, "brk")==0)){
                        printf("3: %s\n", myfile);
                        type=read_pdb_to_residuearray(myfile, amino_acid_list, &residuearray, &Nresidues, 0);
                        if(type==0){
                            reassign_Nterminus_atoms(residuearray, Nresidues, amino_acid_list, Nterminus_error_file, myfile, 0, 0);
                            terminiierror=assign_atomfound_terminii(Nresidues, residuearray, amino_acid_list);
                            if(terminiierror==0){
                                discard_overlapping_atoms(residuearray, Nresidues, amino_acid_list, overlapping_atoms_error_file, myfile, 0, overlapping_atoms_list, &N_overlapping_atoms);
                                if(discard==1){
                                    check_covalent_bond_lengths(residuearray, Nresidues, amino_acid_list, foundbond, maxbondlength, minbondlength, bondlength_error_file, myfile, 0);
                                    check_bond_angles(residuearray, Nresidues, amino_acid_list, Nangles, angles, angle_index, bondangle_error_file, myfile, 0);
                                    discard_bad_tetrahedra(residuearray, Nresidues, Ntetrahedra, tetrahedraindices, amino_acid_list, tetrahedral_error_file, myfile, 0);
                                    check_rings(residuearray, Nresidues, Nrings, ringsize, number_dependent_atoms, ringatomindex, amino_acid_list, rings_error_file, myfile, 0);
                                }
                                resolve_indistinguishable_atoms(residuearray, Nresidues, amino_acid_list, 0, indistinguishable_error_file, myfile, 0);
                                resolve_aromatic_hydrogens(residuearray, Nresidues, amino_acid_list, 0, aromatic_hydrogens_error_file, myfile, 0);
                                map_to_cgmodel(Nresidues, residuearray, amino_acid_list, &cgresiduearray, 1);
                                
                                //  Calculate root mean square displacement between original all-atom configuration and all-atom configuration generated from coarse-grained representation
                                
                                calculate_rmsd_separate_corrections(amino_acid_list, residuearray, cgresiduearray, Nresidues, carbonyl_correction_array, backboneH_correction_array, sidechain_corrections, myfile, large_disagreement_file_specific, entry_average_file);
                                
                                free_cgresiduearray(&cgresiduearray, Nresidues);
                            }
                            free_residuearray(&residuearray, Nresidues, 0);
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    fts_close(ftsp);
    
    //  Output root mean square displacement
    
    sprintf(outputfile, "%s.rmsd", base);
    output_rmsd_onlyspecific(outputfile, amino_acid_list, nrecords, backbone_atom_types, N_backbone_atom_types, Cterminusatomnames);
    sprintf(outputfile, "%s.file_statistics", base);
    output_file_type_statistics(Ntype, outputfile);
	
    return 0;
}

