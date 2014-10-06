#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mystdlib.h"
#include "mymath.h"
#include "protein_functions.h"

void create_amino_acid_struct(char *directory, char *site_director_atoms_file, char *site_atoms_file, char *corrections_parameter_file, char *elements_file, char *Cterminus_atoms_file, char *indistinguishable_atoms_file, char *radial_atoms_file, amino_acid_struct **pamino_acid_list, int *pN_sidechain_corrections, char ***psidechain_correction_identifiers, char ***pbackbone_atom_types, int *pN_backbone_atom_types, int *pbackbone_Oindex, int *pbackbone_Nindex, int *pbackbone_Hindex, int *pbackbone_Cindex, char ***pCterminusatomnames, char *aromatic_hydrogens_file){
	
    //  Create an array of amino acid structures characterizing the atoms in each amino acid
	
    int counter, atomcounter, i, j, k, read, Nmappedatoms, scannedbytes, stringsize, sitecounter, natoms, startcolumn, Nindex, CAindex, Nelements=0, found, NCterminusatoms=0, setcounter, Nconectcounter;					//	must not use i, j for counters because allocate_... macros use them
    (*pN_sidechain_corrections)=0;
    allocate_matrix(char, (*psidechain_correction_identifiers), max_sidechain_corrections, 14);
    declare_array_nozero(char, label, 6);
    declare_array_nozero(char, field, 6);
    declare_array_nozero(char, atomname, 3);
    declare_array_nozero(char, readres, 3);
    allocate_array(amino_acid_struct, (*pamino_acid_list), Naminoacids);
	for(i=0;i<Naminoacids;i++){
		allocate_array(char, ((*pamino_acid_list)[i]).resname, 3);
	}
    declare_matrix_nozero(char, amino_H_names, number_amino_hydrogen_atom_types, 4);
    declare_matrix_nozero(char, amino_H_inputnames, number_amino_hydrogen_atom_types, 4);
	
	strcpy((*pamino_acid_list)[0].resname, "ALA");
	strcpy((*pamino_acid_list)[1].resname, "ARG");
	strcpy((*pamino_acid_list)[2].resname, "ASN");
	strcpy((*pamino_acid_list)[3].resname, "ASP");
	strcpy((*pamino_acid_list)[4].resname, "CYS");
	strcpy((*pamino_acid_list)[5].resname, "GLN");
	strcpy((*pamino_acid_list)[6].resname, "GLU");
	strcpy((*pamino_acid_list)[7].resname, "GLY");
	strcpy((*pamino_acid_list)[8].resname, "HIS");
	strcpy((*pamino_acid_list)[9].resname, "ILE");
	strcpy((*pamino_acid_list)[10].resname, "LEU");
	strcpy((*pamino_acid_list)[11].resname, "LYS");
	strcpy((*pamino_acid_list)[12].resname, "MET");
	strcpy((*pamino_acid_list)[13].resname, "PHE");
	strcpy((*pamino_acid_list)[14].resname, "PRO");
	strcpy((*pamino_acid_list)[15].resname, "SER");
	strcpy((*pamino_acid_list)[16].resname, "THR");
	strcpy((*pamino_acid_list)[17].resname, "TRP");
	strcpy((*pamino_acid_list)[18].resname, "TYR");
	strcpy((*pamino_acid_list)[19].resname, "VAL");
	
	FILE *inp, *directorinp, *atominp;
	declare_array_nozero(char, filename, maxlinelength);
    declare_array_nozero(char, linechar, maxlinelength);
    char *originallinechar;
    originallinechar=linechar;
    declare_3d_tensor(char, siteatomname, 3, maxsiteatoms, 3);        //  big enough to store site atoms for any residue
    declare_array(int, Nfoundsiteatoms, 3);
    declare_array(int, Nsiteatoms, 3);
    declare_matrix_nozero(char, sitemapatomname, maxatomsperresidue, 3);
    declare_array(int, sitemapsitenumber, maxatomsperresidue);
    declare_array(int, inverseserial, maxatomsperresidue);
    declare_matrix_nozero(char, atomtype, maxatomtypes, 4);
    declare_array(int, element, maxatomtypes);
    declare_matrix_nozero(char, Cterminusatoms, maxatomsperresidue, 3);
    
    //  Read parameter file indicating atomic number of atom types
    
    inp=fopen(elements_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        extract_piece_of_string(linechar, (atomtype[Nelements]), 1, 4);
        extract_piece_of_string(linechar, field, 5, 6);
        element[Nelements]=atoi(field);
        read=mygetline(linechar, maxlinelength, inp);
        Nelements++;
        if(Nelements==maxatomtypes) my_exit("too many atom types");
    }
    fclose(inp);
    
    //  Read parameter file indicating atoms in C-terminus site
    
    inp=fopen(Cterminus_atoms_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        extract_piece_of_string(linechar, (Cterminusatoms[NCterminusatoms]), 1, 4);
        NCterminusatoms++;
        read=mygetline(linechar, maxlinelength, inp);
    }
    fclose(inp);
    allocate_matrix(char, (*pCterminusatomnames), NCterminusatoms, 4);
    for(i=0;i<NCterminusatoms;i++){
        strcpy((*pCterminusatomnames)[i], Cterminusatoms[i]);
    }
    
	for(counter=0;counter<Naminoacids;counter++){
        linechar=originallinechar;
        Nmappedatoms=0;
		sprintf(filename, "%s/%s.pdb", directory, (*pamino_acid_list)[counter].resname);
		inp=fopen(filename, "r");
		if(inp==NULL){
			printf("%s not found!\n", filename);
            exit(1);
		}
        
        //  Read parameter file indicating which atoms describe the orientation of each site
        
        directorinp=fopen(site_director_atoms_file, "r");
        read=mygetline(linechar, maxlinelength, directorinp);
        while(read>0){
            extract_piece_of_string(linechar, readres, 1, 3);
			if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                extract_piece_of_string(linechar, field, 4, 4);
                (*pamino_acid_list)[counter].Nsites=atoi(field);
				allocate_array(int, ((*pamino_acid_list)[counter]).sitecode, (*pamino_acid_list)[counter].Nsites);
				allocate_array(int, ((*pamino_acid_list)[counter]).sitecount, (*pamino_acid_list)[counter].Nsites);
				allocate_array(int, ((*pamino_acid_list)[counter]).element, maxatomsperresidue);
				allocate_array(int, ((*pamino_acid_list)[counter]).nconect, maxatomsperresidue);
                for(i=0;i<maxatomsperresidue;i++) ((*pamino_acid_list)[counter]).nconect[i]=0;
                ((*pamino_acid_list)[counter]).conect=xcalloc(maxatomsperresidue, sizeof(int *));
				((*pamino_acid_list)[counter]).siteatomindex=xcalloc((*pamino_acid_list)[counter].Nsites, sizeof(int *));
                startcolumn=5;
                for(sitecounter=0;sitecounter<(*pamino_acid_list)[counter].Nsites;sitecounter++){
                    extract_piece_of_string(linechar, field, startcolumn, startcolumn);
                    (*pamino_acid_list)[counter].sitecode[sitecounter]=atoi(field);
                    startcolumn++;
                    if((*pamino_acid_list)[counter].sitecode[sitecounter]==0) Nsiteatoms[sitecounter]=3;
                    else{
                        printf("Need to write code for how many site atoms for site code %i!\n", (*pamino_acid_list)[counter].sitecode[sitecounter]);
                        exit(1);
                    }
					allocate_array(int, (((*pamino_acid_list)[counter]).siteatomindex[sitecounter]), Nsiteatoms[sitecounter]);
					for(j=0;j<Nsiteatoms[sitecounter];j++){
                        extract_piece_of_string(linechar, siteatomname[sitecounter][j], startcolumn, startcolumn+3);
                        startcolumn+=4;
                    }
                }
			}
			read=mygetline(linechar, maxlinelength, directorinp);
        }
		fclose(directorinp);
        allocate_array(int, ((*pamino_acid_list)[counter]).Cterminussiteatomindex, 3);
		
		//	Read parameter file indicating which site each atom belongs to
		
		atominp=fopen(site_atoms_file, "r");
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            stringsize=strlen(linechar);
            natoms=(stringsize-3)/5;
            extract_piece_of_string(linechar, readres, 1, 3);
			if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                Nmappedatoms+=natoms;
                for(i=0;i<natoms;i++){
                    extract_piece_of_string(linechar, sitemapatomname[i], 4+5*i, 4+5*i+3);
                    extract_piece_of_string(linechar, field, 4+5*i+4, 4+5*i+4);
                    sitemapsitenumber[i]=atoi(field);
                }
            }
			read=mygetline(linechar, maxlinelength, atominp);
		}
		fclose(atominp);

        //  Read residue PDB file once to count number of atoms
        
        (*pamino_acid_list)[counter].Natoms=0;
        read=mygetline(linechar, maxlinelength, inp);
        while(read>0){
            extract_piece_of_string(linechar, label, 1, 6);
            if(strcmp(label, "ATOM  ")==0){
                ((*pamino_acid_list)[counter].Natoms)++;
            }
            read=mygetline(linechar, maxlinelength, inp);
        }
        fclose(inp);
                
        //  Fifteen extra slots for amino terminus hydrogens:
        //  three categories (PRO, GLY, and other) and five atom types (2 neutral and 3 protonated)
        
		(*pamino_acid_list)[counter].Natoms+=number_amino_hydrogen_atom_types;
        
        allocate_matrix(char, (((*pamino_acid_list)[counter]).atomnames), (*pamino_acid_list)[counter].Natoms, 3);
        allocate_matrix(char, (((*pamino_acid_list)[counter]).input_atomnames), (*pamino_acid_list)[counter].Natoms, 3);
        allocate_array(int, (((*pamino_acid_list)[counter]).sitemap), (*pamino_acid_list)[counter].Natoms);
        allocate_array(double_triple, (((*pamino_acid_list)[counter]).avpos), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        allocate_array(double_triple, (((*pamino_acid_list)[counter]).varpos), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        allocate_array(double_triple, (((*pamino_acid_list)[counter]).msd_residuespecific), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        allocate_array(double, (((*pamino_acid_list)[counter]).fourthmoment_residuespecific), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        allocate_array(int, (((*pamino_acid_list)[counter]).count), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        allocate_array(int, (((*pamino_acid_list)[counter]).msdcount), (*pamino_acid_list)[counter].Natoms+NCterminusatoms);
        (*pamino_acid_list)[counter].Oindex=-1;
        (*pamino_acid_list)[counter].Nindex=-1;
        (*pamino_acid_list)[counter].Hindex=-1;
        (*pamino_acid_list)[counter].Cindex=-1;
        (*pamino_acid_list)[counter].CAindex=-1;
        (*pamino_acid_list)[counter].OXTindex=-1;
        //(*pamino_acid_list)[counter].HXTindex=-1;
        allocate_array(int, (((*pamino_acid_list)[counter]).Cterminusindex), (*pamino_acid_list)[counter].Natoms);
        for(i=0;i<(*pamino_acid_list)[counter].Natoms;i++) (*pamino_acid_list)[counter].Cterminusindex[i]=-1;
        (*pamino_acid_list)[counter].NCterminusatoms=NCterminusatoms;
		
        //  Read PBD file again to assign atom names
		
        inp=fopen(filename, "r");
        atomcounter=0;
        read=mygetline(linechar, maxlinelength, inp);
        for(i=0;i<(*pamino_acid_list)[counter].Natoms;i++) (*pamino_acid_list)[counter].sitemap[i]=-1;      //   to indicate no map to site found
        while(read>0){
            extract_piece_of_string(linechar, label, 1, 6);
            if(strcmp(label, "ATOM  ")==0){
                extract_piece_of_string(linechar, field, 7, 11);
                extract_piece_of_string(linechar, atomname, 13, 16);
                extract_piece_of_string(linechar, readres, 18, 20);
                if(strcmp(readres, (*pamino_acid_list)[counter].resname)!=0){
                    printf("Read residue name %s does not match filename %s!\n", readres, (*pamino_acid_list)[counter].resname);
                    exit(1);
                }
 				strcpy((*pamino_acid_list)[counter].atomnames[atomcounter], atomname);
 				strcpy((*pamino_acid_list)[counter].input_atomnames[atomcounter], atomname);
                if(atoi(field)>=maxatomsperresidue){
                    printf("serial=%i!\n", atoi(field));
                    exit(1);
                }
                inverseserial[atoi(field)]=atomcounter;
                found=0;
                for(i=0;i<Nelements;i++){
                    if(strcmp(atomname, atomtype[i])==0){
                        (*pamino_acid_list)[counter].element[atomcounter]=element[i];
                        found=1;
                        break;
                    }
                }
                if(found==0){
                    printf("Atomname %s in residue %s not mapped to an element\n", atomname, readres);
                    exit(1);
                }
                if(strcmp(atomname, " O  ")==0) (*pamino_acid_list)[counter].Oindex=atomcounter;
                if(strcmp(atomname, " N  ")==0) (*pamino_acid_list)[counter].Nindex=atomcounter;
                if(strcmp(atomname, " H  ")==0) (*pamino_acid_list)[counter].Hindex=atomcounter;
                if(strcmp(atomname, " C  ")==0) (*pamino_acid_list)[counter].Cindex=atomcounter;
                if(strcmp(atomname, " CA ")==0) (*pamino_acid_list)[counter].CAindex=atomcounter;
                if(strcmp(atomname, " OXT")==0) (*pamino_acid_list)[counter].OXTindex=atomcounter;
                //if(strcmp(atomname, " HXT")==0) (*pamino_acid_list)[counter].HXTindex=atomcounter;
                for(i=0;i<NCterminusatoms;i++){
                    if(strcmp(atomname, Cterminusatoms[i])==0){
                        (*pamino_acid_list)[counter].Cterminusindex[atomcounter]=i;
                        (*pamino_acid_list)[counter].element[(*pamino_acid_list)[counter].Natoms+i]=(*pamino_acid_list)[counter].element[atomcounter];
                        if(i<3){
                            (*pamino_acid_list)[counter].Cterminussiteatomindex[i]=atomcounter;
                        }
                        break;
                    }
                }
                for(i=0;i<(*pamino_acid_list)[counter].Nsites;i++){
                    for(j=0;j<Nsiteatoms[i];j++){
                        if(strcmp(atomname, siteatomname[i][j])==0){
                            (*pamino_acid_list)[counter].siteatomindex[i][j]=atomcounter;
                            Nfoundsiteatoms[i]++;
                        }
                    }
                }
                found=0;
                for(i=0;i<Nmappedatoms;i++){
                    if(strcmp(atomname, sitemapatomname[i])==0){
                        (*pamino_acid_list)[counter].sitemap[atomcounter]=sitemapsitenumber[i];
                        found=1;
                        break;
                    }
                }
                if(found==0){
                    printf("Atom %s in residue %s not mapped to a site!\n", atomname, readres);
                    printf("Nmappedatoms=%i, last %s\n", Nmappedatoms, sitemapatomname[Nmappedatoms-1]);
                    exit(1);
                }
                atomcounter++;
            }
            else if(strcmp(label, "CONECT")==0){
                stringsize=strlen(linechar);
                if((stringsize-12)%5!=0){
                    printf("CONECT line length %i!\n", stringsize-1);
                    exit(1);
                }
                extract_piece_of_string(linechar, field, 7, 11);
                atomcounter=inverseserial[atoi(field)];
                (*pamino_acid_list)[counter].nconect[atomcounter]=(stringsize-12)/5;
                
                //  Add fifteen covalent neighbors to N for the fifteen amino acid terminus hydrogen atom types
                
                if(strcmp((*pamino_acid_list)[counter].atomnames[atomcounter], " N  ")==0) (*pamino_acid_list)[counter].nconect[atomcounter]+=number_amino_hydrogen_atom_types;
                
                allocate_array(int, ((((*pamino_acid_list)[counter]).conect)[atomcounter]), (*pamino_acid_list)[counter].nconect[atomcounter]);
                for(i=0;i<(stringsize-12)/5;i++){
                    extract_piece_of_string(linechar, field, i*5+12, i*5+16);
                    (*pamino_acid_list)[counter].conect[atomcounter][i]=inverseserial[atoi(field)];
                }                
            }
            read=mygetline(linechar, maxlinelength, inp);
        }
		
		//	Extra atoms for neutral (HN1, HN2) or protonated HM1, HM2, HM3 atoms; assuming they all in backbone (0) site and not in site definitions
		
        (*pamino_acid_list)[counter].aminoHindex=(*pamino_acid_list)[counter].Natoms-number_amino_hydrogen_atom_types;
        
        atomcounter=(*pamino_acid_list)[counter].aminoHindex;
        Nindex=(*pamino_acid_list)[counter].Nindex;
        Nconectcounter=(*pamino_acid_list)[counter].nconect[Nindex]-number_amino_hydrogen_atom_types;
        
        strcpy(amino_H_names[0], " HN1");
        strcpy(amino_H_names[1], " HN2");
        strcpy(amino_H_names[2], " HM1");
        strcpy(amino_H_names[3], " HM2");
        strcpy(amino_H_names[4], " HM3");
        strcpy(amino_H_names[5], "HGN1");
        strcpy(amino_H_names[6], "HGN2");
        strcpy(amino_H_names[7], "HGM1");
        strcpy(amino_H_names[8], "HGM2");
        strcpy(amino_H_names[9], "HGM3");
        strcpy(amino_H_names[10], "HPN1");
        strcpy(amino_H_names[11], "HPN2");
        strcpy(amino_H_names[12], "HPM1");
        strcpy(amino_H_names[13], "HPM2");
        strcpy(amino_H_names[14], "HPM3");
        strcpy(amino_H_inputnames[0], " HN1");
        strcpy(amino_H_inputnames[1], " HN2");
        strcpy(amino_H_inputnames[2], " H1 ");
        strcpy(amino_H_inputnames[3], " H2 ");
        strcpy(amino_H_inputnames[4], " H3 ");
        strcpy(amino_H_inputnames[5], "HGN1");
        strcpy(amino_H_inputnames[6], "HGN2");
        strcpy(amino_H_inputnames[7], "HGM1");
        strcpy(amino_H_inputnames[8], "HGM2");
        strcpy(amino_H_inputnames[9], "HGM3");
        strcpy(amino_H_inputnames[10], "HPN1");
        strcpy(amino_H_inputnames[11], "HPN2");
        strcpy(amino_H_inputnames[12], "HPM1");
        strcpy(amino_H_inputnames[13], "HPM2");
        strcpy(amino_H_inputnames[14], "HPM3");
        
        for(i=0;i<number_amino_hydrogen_atom_types;i++){
            strcpy((*pamino_acid_list)[counter].atomnames[atomcounter], amino_H_names[i]);
            strcpy((*pamino_acid_list)[counter].input_atomnames[atomcounter], amino_H_inputnames[i]);
            (*pamino_acid_list)[counter].sitemap[atomcounter]=0;
            (*pamino_acid_list)[counter].element[atomcounter]=1;
            (*pamino_acid_list)[counter].nconect[atomcounter]=1;
            allocate_array(int, ((((*pamino_acid_list)[counter]).conect)[atomcounter]), (*pamino_acid_list)[counter].nconect[atomcounter]);
            (*pamino_acid_list)[counter].conect[atomcounter][0]=Nindex;
            (*pamino_acid_list)[counter].conect[Nindex][Nconectcounter]=atomcounter;
            atomcounter++;
            Nconectcounter++;
        }
        
		//	Check that there are no duplicates
		
        for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
            for(k=0;k<j;k++){
                if(strcmp((*pamino_acid_list)[counter].atomnames[j], (*pamino_acid_list)[counter].atomnames[k])==0){
                    printf("residue %s: atoms %i and %i both %s (%s)!\n", (*pamino_acid_list)[counter].resname, j, k, (*pamino_acid_list)[counter].atomnames[j], (*pamino_acid_list)[counter].atomnames[k]);
                    exit(1);
                }
            }
        }
        fclose(inp);        

		//	Read parameter file indicating atoms that need corrections
                
        (*pamino_acid_list)[counter].N_corrections=0;
        allocate_array(int, ((*pamino_acid_list)[counter]).atom_correction_index, (*pamino_acid_list)[counter].Natoms);
        allocate_array(int, ((*pamino_acid_list)[counter]).atom_correction_index_withinaa, (*pamino_acid_list)[counter].Natoms);
        for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
            ((*pamino_acid_list)[counter]).atom_correction_index[j]=-1;
            ((*pamino_acid_list)[counter]).atom_correction_index_withinaa[j]=-1;
        }
        atominp=fopen(corrections_parameter_file, "r");
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            extract_piece_of_string(linechar, readres, 1, 3);
			if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                ((*pamino_acid_list)[counter].N_corrections)++;
            }
            read=mygetline(linechar, maxlinelength, atominp);
        }
        fclose(atominp);
        allocate_array(int, ((*pamino_acid_list)[counter]).correction_index, (*pamino_acid_list)[counter].N_corrections);
        allocate_array(int, ((*pamino_acid_list)[counter]).corrected_atom_index, (*pamino_acid_list)[counter].N_corrections);
        allocate_array(int, ((*pamino_acid_list)[counter]).correcting_atom_index, (*pamino_acid_list)[counter].N_corrections);
        atominp=fopen(corrections_parameter_file, "r");
        setcounter=0;
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            stringsize=strlen(linechar);
            if(stringsize!=12){
                printf("line in %s of length %i\n", corrections_parameter_file, stringsize);
                exit(1);
            }
            extract_piece_of_string(linechar, readres, 1, 3);
			if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                (*pamino_acid_list)[counter].correction_index[setcounter]=(*pN_sidechain_corrections);
                strcpy((*psidechain_correction_identifiers)[(*pN_sidechain_corrections)], readres);
                extract_piece_of_string(linechar, field, 4, 7);
                for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
                    if(strcmp(field, (*pamino_acid_list)[counter].atomnames[j])==0){
                        (*pamino_acid_list)[counter].atom_correction_index[j]=(*pN_sidechain_corrections);
                        (*pamino_acid_list)[counter].atom_correction_index_withinaa[j]=setcounter;
                        (*pamino_acid_list)[counter].corrected_atom_index[setcounter]=j;
                        strcat((*psidechain_correction_identifiers)[(*pN_sidechain_corrections)], ".");
                        strcat((*psidechain_correction_identifiers)[(*pN_sidechain_corrections)], field);
                   }
                }
                extract_piece_of_string(linechar, field, 8, 11);
                for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
                    if(strcmp(field, (*pamino_acid_list)[counter].atomnames[j])==0){
                        (*pamino_acid_list)[counter].correcting_atom_index[setcounter]=j;
                        strcat((*psidechain_correction_identifiers)[(*pN_sidechain_corrections)], ".");
                        strcat((*psidechain_correction_identifiers)[(*pN_sidechain_corrections)], field);
                    }
                }
                (*pN_sidechain_corrections)++;
                setcounter++;
            }
            read=mygetline(linechar, maxlinelength, atominp);
        }
		fclose(atominp);
                
        //  Read parameter file indicating indistinguishable atoms
        //  Read through once to count number of sets of indistinguishable atoms
        
        (*pamino_acid_list)[counter].N_indistinguishable_sets=0;
        atominp=fopen(indistinguishable_atoms_file, "r");
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            extract_piece_of_string(linechar, readres, 1, 3);
            if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                ((*pamino_acid_list)[counter].N_indistinguishable_sets)++;
            }
            read=mygetline(linechar, maxlinelength, atominp);
        }
        fclose(atominp);
        allocate_array(int, (((*pamino_acid_list)[counter]).indistinguishable_set_size), (*pamino_acid_list)[counter].N_indistinguishable_sets);
        allocate_array(double, (((*pamino_acid_list)[counter]).first_dividing_angle), (*pamino_acid_list)[counter].N_indistinguishable_sets);
        allocate_matrix(int, (((*pamino_acid_list)[counter]).indistinguishable_atoms), (*pamino_acid_list)[counter].N_indistinguishable_sets, 6);

        //  Read through again to assign atoms
        
        setcounter=0;
        atominp=fopen(indistinguishable_atoms_file, "r");
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            stringsize=strlen(linechar);
            extract_piece_of_string(linechar, readres, 1, 3);
            if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                if(stringsize==27){
                    (*pamino_acid_list)[counter].indistinguishable_set_size[setcounter]=2;
                    extract_piece_of_string(linechar, field, 24, 26);
                    (*pamino_acid_list)[counter].first_dividing_angle[setcounter]=atof(field)/180.*M_PI;
                }
                else if(stringsize==31){
                    (*pamino_acid_list)[counter].indistinguishable_set_size[setcounter]=3;
                    extract_piece_of_string(linechar, field, 28, 30);
                    (*pamino_acid_list)[counter].first_dividing_angle[setcounter]=atof(field)/180.*M_PI;
                }
                else{
                    printf("stringsize %i in %s!\n", stringsize, indistinguishable_atoms_file);
                    exit(1);
                }
                for(i=0;i<(*pamino_acid_list)[counter].indistinguishable_set_size[setcounter]+3;i++){
                    extract_piece_of_string(linechar, field, 4+i*4, 7+i*4);
                    for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
                        if(strcmp(field, (*pamino_acid_list)[counter].atomnames[j])==0){
                            (*pamino_acid_list)[counter].indistinguishable_atoms[setcounter][i]=j;
                            break;
                        }
                    }
                    if(j==(*pamino_acid_list)[counter].Natoms){
                        printf("indistinguishable atom <%s> not found in residue <%s>!\n", field, (*pamino_acid_list)[counter].resname);
                        exit(1);
                    }
                }
                setcounter++;
            }
            read=mygetline(linechar, maxlinelength, atominp);
        }
        fclose(atominp);
        
        //  Read parameter file indicating aromatic hydrogens
        
        atominp=fopen(aromatic_hydrogens_file, "r");
        read=mygetline(linechar, maxlinelength, atominp);
        while(read>0){
            stringsize=strlen(linechar);
            if((stringsize-4)%8!=0){
                printf("line length %i in aromatic hydrogens file!\n", stringsize);
                exit(1);
            }
            extract_piece_of_string(linechar, readres, 1, 3);
            if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                (*pamino_acid_list)[counter].N_aromatic_hydrogens=(stringsize-4)/8;
                if((*pamino_acid_list)[counter].N_aromatic_hydrogens>max_aromatic_hydrogens){
                    my_exit("Too many aromatic hydrogens");
                }
                allocate_array(int, ((*pamino_acid_list)[counter].aromatic_hydrogen_list), 2*(*pamino_acid_list)[counter].N_aromatic_hydrogens);
                for(i=0;i<2*(*pamino_acid_list)[counter].N_aromatic_hydrogens;i++){
                    extract_piece_of_string(linechar, field, 4+i*4, 7+i*4);
                    for(j=0;j<(*pamino_acid_list)[counter].Natoms;j++){
                        if(strcmp(field, (*pamino_acid_list)[counter].atomnames[j])==0){
                            (*pamino_acid_list)[counter].aromatic_hydrogen_list[i]=j;
                            break;
                        }
                    }
                    if(j==(*pamino_acid_list)[counter].Natoms){
                        printf("indistinguishable atom <%s> not found in residue <%s>!\n", field, (*pamino_acid_list)[counter].resname);
                        exit(1);
                    }
                }
            }
            read=mygetline(linechar, maxlinelength, atominp);
        }
        fclose(atominp);
        
        //  Read parameter file indicating atoms to calculate average radial position
        
        atominp=fopen(radial_atoms_file, "r");
        if(atominp!=NULL){

            //  Read through once to count the number of radial sets

            read=mygetline(linechar, maxlinelength, atominp);
            while(read>0){
                stringsize=strlen(linechar);
                if(stringsize%4!=0){
                    printf("line length %i in radial atoms file!\n", stringsize);
                    exit(1);
                }
                extract_piece_of_string(linechar, readres, 1, 3);
                if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                    (*pamino_acid_list)[counter].N_radial_sets++;
                }
                read=mygetline(linechar, maxlinelength, atominp);
            }
            fclose(atominp);
            allocate_array(int, ((*pamino_acid_list)[counter].radial_set_size), (*pamino_acid_list)[counter].N_radial_sets);
            allocate_matrix(int, ((*pamino_acid_list)[counter].radial_atoms), (*pamino_acid_list)[counter].N_radial_sets, max_atoms_in_radial_set);
            allocate_array(int, ((*pamino_acid_list)[counter].radial_count), (*pamino_acid_list)[counter].N_radial_sets);
            allocate_array(double_double, ((*pamino_acid_list)[counter].radial_average), (*pamino_acid_list)[counter].N_radial_sets);
            for(i=0;i<(*pamino_acid_list)[counter].N_radial_sets;i++){
                ((*pamino_acid_list)[counter].radial_average)[i].x=0;
                ((*pamino_acid_list)[counter].radial_average)[i].y=0;
                (*pamino_acid_list)[counter].radial_count[i]=0;
            }
            atominp=fopen(radial_atoms_file, "r");
            i=0;
            read=mygetline(linechar, maxlinelength, atominp);
            while(read>0){
                stringsize=strlen(linechar);
                extract_piece_of_string(linechar, readres, 1, 3);
                if(strcmp(readres, (*pamino_acid_list)[counter].resname)==0){
                    (*pamino_acid_list)[counter].radial_set_size[i]=stringsize/4-1;
                    if((*pamino_acid_list)[counter].radial_set_size[i]>max_atoms_in_radial_set) my_exit("Big radial set size!");
                    for(j=0;j<(*pamino_acid_list)[counter].radial_set_size[i];j++){
                        extract_piece_of_string(linechar, field, 4+j*4, 7+j*4);
                        for(k=0;k<(*pamino_acid_list)[counter].Natoms;k++){
                            if(strcmp(field, (*pamino_acid_list)[counter].atomnames[k])==0){
                                (*pamino_acid_list)[counter].radial_atoms[i][j]=k;
                                break;
                            }
                        }
                        if(k==(*pamino_acid_list)[counter].Natoms){
                            printf("<%s> in radial atoms not found in %s!\n", field, (*pamino_acid_list)[counter].resname);
                            exit(1);
                        }
                    }
                    i++;
                }
                read=mygetline(linechar, maxlinelength, atominp);
            }
            fclose(atominp);
            
        }
        
    }
        
    //  Create array of backbone atom types
    
    (*pN_backbone_atom_types)=0;
    declare_matrix_nozero(char, temporary_backbone_atom_types, maxbackboneatomtypes, 4);
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<(*pamino_acid_list)[i].Natoms;j++){
            if((*pamino_acid_list)[i].sitemap[j]==0){
                found=0;
                for(counter=0;counter<(*pN_backbone_atom_types);counter++){
                    if(strcmp((*pamino_acid_list)[i].atomnames[j], temporary_backbone_atom_types[counter])==0){
                        found=1;
                        break;
                    }
                }
                if(found==0){
                    if((*pN_backbone_atom_types)==maxbackboneatomtypes){
                        my_exit("Found too many backbone atoms!\n");
                    }
                    strcpy(temporary_backbone_atom_types[(*pN_backbone_atom_types)], (*pamino_acid_list)[i].atomnames[j]);
                    if(strcmp((*pamino_acid_list)[i].atomnames[j], " O  ")==0) (*pbackbone_Oindex)=(*pN_backbone_atom_types);
                    if(strcmp((*pamino_acid_list)[i].atomnames[j], " N  ")==0) (*pbackbone_Nindex)=(*pN_backbone_atom_types);
                    if(strcmp((*pamino_acid_list)[i].atomnames[j], " H  ")==0) (*pbackbone_Hindex)=(*pN_backbone_atom_types);
                    if(strcmp((*pamino_acid_list)[i].atomnames[j], " C  ")==0) (*pbackbone_Cindex)=(*pN_backbone_atom_types);
                    (*pN_backbone_atom_types)++;
                }
            }
        }
    }
    allocate_matrix_nozero(char, (*pbackbone_atom_types), (*pN_backbone_atom_types+extra_backbone_atom_types), 4);
    for(i=0;i<(*pN_backbone_atom_types);i++){
        strcpy((*pbackbone_atom_types)[i], temporary_backbone_atom_types[i]);
    }
    free_matrix(temporary_backbone_atom_types, maxbackboneatomtypes);
        
    free(label);
    free(readres);
    free(atomname);
	free(filename);
	free(originallinechar);
    free_matrix(siteatomname, 3);
    free(Nsiteatoms);
    free(Nfoundsiteatoms);
    free_matrix(sitemapatomname, 20);
    free(sitemapsitenumber);
    free(field);
    free(inverseserial);
    free_matrix(atomtype, maxatomtypes);
    free(element);
    free_matrix(Cterminusatoms, 3);
    free_matrix(amino_H_names, number_amino_hydrogen_atom_types);
    free_matrix(amino_H_inputnames, number_amino_hydrogen_atom_types);
}

void change_atom_names_in_amino_acid_struct(amino_acid_struct *amino_acid_list, char *change_atom_names_file, int *pN_backbone_atom_types, char **backbone_atom_types){
    int i, j, k, read, new_backbone_atom_types=0;
    FILE *inp;
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, field, 4);
    for(i=0;i<Naminoacids;i++){
        
        //  Change atom names for certain residues according to change_atom_names_file
        
        for(j=0;j<amino_acid_list[i].Natoms;j++){
            inp=fopen(change_atom_names_file, "r");
            read=mygetline(linechar, maxlinelength, inp);
            while(read>0){
                if(read!=12){
                    printf("line length %i in change_atom_names!\n", read);
                    exit(1);
                }
                extract_piece_of_string(linechar, field, 1, 3);
                if(strcmp(field, amino_acid_list[i].resname)==0){
                    extract_piece_of_string(linechar, field, 4, 7);
                    if(strcmp(field, amino_acid_list[i].atomnames[j])==0){
                        extract_piece_of_string(linechar, field, 8, 11);
                        for(k=0;k<amino_acid_list[i].Natoms;k++){
                            if(k!=j){
                                if(strcmp(field, amino_acid_list[i].atomnames[k])==0){
                                    printf("new atom type %s already in amino acid %s!\n", field, amino_acid_list[i].resname);
                                    exit(1);
                                }
                            }
                        }
                        strcpy(amino_acid_list[i].atomnames[j], field);
                        for(k=0;k<(*pN_backbone_atom_types);k++){
                            if(strcmp(field, backbone_atom_types[k])==0){
                                printf("new atom type %s already in backbone atom types!\n", field);
                                exit(1);
                            }
                        }
                        if(new_backbone_atom_types==extra_backbone_atom_types){
                            my_exit("too many new backbone atom types!");
                        }
                        strcpy(backbone_atom_types[*pN_backbone_atom_types+new_backbone_atom_types], field);
                        new_backbone_atom_types++;
                    }
                }
                read=mygetline(linechar, maxlinelength, inp);
            }
        }
    }
    (*pN_backbone_atom_types)+=extra_backbone_atom_types;
}

void initialize_correction_structure(correction_structure *pmystructure){
    int i, j;
    allocate_array(double, ((*pmystructure).predicted_pos_av), 3);
    allocate_array(double, ((*pmystructure).predicted_pos_var), 3);
    allocate_array(double, ((*pmystructure).pos_difference_av), 3);
    allocate_matrix(double, ((*pmystructure).covariance_matrix_cross), 3, 3);
    allocate_matrix(double, ((*pmystructure).covariance_matrix_self), 3, 3);
    allocate_matrix(double, ((*pmystructure).slope), 3, 3);
    allocate_array(double, ((*pmystructure).intercept), 3);
    (*pmystructure).count=0;
}

int read_pdb_to_residuearray(char *filename, amino_acid_struct *amino_acid_list, residuedata **presiduearray, int *pNresidues, int record_all_data){
	
    //  Read data from a .pdb file and convert to an array of residue structures containing atomic coordinates and connectivity information
	
    int rescounter, chaincounter, i, j, k, read, first=1, Nchains=0, foundseqres=0, formats, residuesonline, thischainid, foundatom=0, resseq, lastresseq=0, founddbref=0;
    int *chainidarray;
    double_triple position;
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, label, 6);
    declare_array_nozero(char, field, 20);
    declare_array_nozero(char, chainid, 1);
    declare_array_nozero(char, lastchainid, 1);
    declare_array_nozero(char, icode, 1);
    declare_array_nozero(char, lasticode, 1);
    declare_array_nozero(char, atomname, 4);
    declare_array_nozero(char, readres, 3);
    declare_array_nozero(chaindata, mychaindata, maxNchains);
	(*pNresidues)=0;
    FILE *inp;
    inp=fopen(filename, "r");
	if(inp==NULL){
		printf("%s not found!\n", filename);
		exit(1);
	}
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        extract_piece_of_string(linechar, label, 1, 6);
        if(strcmp(label, "DBREF ")==0){
			founddbref=1;
            extract_piece_of_string(linechar, field, 15, 18);
            mychaindata[Nchains].seqbegin=atoi(field);
            extract_piece_of_string(linechar, field, 21, 24);
            mychaindata[Nchains].seqend=atoi(field);
            extract_piece_of_string(linechar, field, 63, 67);
            mychaindata[Nchains].dbseqend=atoi(field);
            Nchains++;
        }
        else if(strcmp(label, "DBREF1")==0){
			founddbref=1;
            extract_piece_of_string(linechar, field, 15, 18);
            mychaindata[Nchains].seqbegin=atoi(field);
            extract_piece_of_string(linechar, field, 21, 24);
            mychaindata[Nchains].seqend=atoi(field);
        }
        else if(strcmp(label, "DBREF2")==0){
			extract_piece_of_string(linechar, field, 58, 67);
			mychaindata[Nchains].dbseqend=atoi(field);
            Nchains++;
		}
        else if(strcmp(label, "MODEL ")==0){
            
            //  Ensemble of models found
            
            free(linechar);
            free(label);
            free(field);
            free(chainid);
            free(lastchainid);
            free(icode);
            free(lasticode);
            free(atomname);
            free(readres);
            free(mychaindata);
            fclose(inp);
            return 1;
        }
        else if(strcmp(label, "ATOM  ")==0){
            extract_piece_of_string(linechar, field, 17, 17);
            if((strcmp(field, " ")==0)){
                
                //  Only assign atom if its position is the only one (blank) or the most common among alternates (A)
                
                extract_piece_of_string(linechar, atomname, 13, 16);
                extract_piece_of_string(linechar, readres, 18, 20);
                extract_piece_of_string(linechar, chainid, 22, 22);
                extract_piece_of_string(linechar, field, 23, 26);
                resseq=atoi(field);
                extract_piece_of_string(linechar, icode, 27, 27);
                if(foundatom==0){
                    if(founddbref==0){
                        
                        //  No DBREF entries found
                        
						free(linechar);
						free(label);
						free(field);
						free(chainid);
						free(lastchainid);
						free(icode);
						free(lasticode);
						free(atomname);
						free(readres);
						free(mychaindata);
						fclose(inp);
						return 4;
					}
                    if(Nchains==0) my_exit("ATOM data found before any DBREF data");
                    for(i=0;i<Nchains;i++){
						if(mychaindata[i].seqend>mychaindata[i].dbseqend){
							(*pNresidues)+=(1+mychaindata[i].seqend);
						}
						else (*pNresidues)+=(1+mychaindata[i].dbseqend);
                        if(mychaindata[i].seqbegin<0){
                            (*pNresidues)-=mychaindata[i].seqbegin;
                        }
                    }
                    rescounter=0;
                    chaincounter=0;
                    j=0;
                    while(strcmp(readres, amino_acid_list[j].resname)!=0){
                        j++;
                        if(j==Naminoacids){
                            
                            //  Residue not among the known amino acids
                            
                            free(linechar);
                            free(label);
                            free(field);
                            free(chainid);
                            free(lastchainid);
                            free(icode);
                            free(lasticode);
                            free(atomname);
                            free(readres);
                            free(mychaindata);
                            fclose(inp);
                            return 2;
                        }
                    }
                    (*pNresidues)+=residuebuffer;
                    allocate_array(residuedata, (*presiduearray), (*pNresidues));
                    for(i=0;i<(*pNresidues);i++){
                        (*presiduearray)[i].residnumber=-1;
                        (*presiduearray)[i].chainid=-1;
                        (*presiduearray)[i].Natoms=0;
                        (*presiduearray)[i].Nterminusflag=0;
                    }
                    (*presiduearray)[rescounter].chainid=chaincounter;
                    (*presiduearray)[rescounter].residnumber=j;
                    (*presiduearray)[rescounter].Natoms=amino_acid_list[j].Natoms;
                    allocate_array(int, ((*presiduearray)[rescounter]).atomfound, (*presiduearray)[rescounter].Natoms);
                    for(i=0;i<(*presiduearray)[rescounter].Natoms;i++) (*presiduearray)[rescounter].atomfound[i]=0;
                    allocate_array(double_triple, ((*presiduearray)[rescounter]).atomposition, (*presiduearray)[rescounter].Natoms);
                    allocate_array(char, ((*presiduearray)[rescounter]).chainidchar, 1);
                    (*presiduearray)[rescounter].resseq=resseq;
                    strcpy(((*presiduearray)[rescounter].chainidchar), chainid);
                    if(record_all_data==1){
                        allocate_matrix(char, ((*presiduearray)[rescounter]).serialchar, (*presiduearray)[rescounter].Natoms, 5);
                        allocate_matrix(char, ((*presiduearray)[rescounter]).icodechar, (*presiduearray)[rescounter].Natoms, 1);
                        allocate_matrix(char, ((*presiduearray)[rescounter]).occupancychar, (*presiduearray)[rescounter].Natoms, 6);
                        allocate_matrix(char, ((*presiduearray)[rescounter]).tempfactorchar, (*presiduearray)[rescounter].Natoms, 6);
                        allocate_matrix(char, ((*presiduearray)[rescounter]).elementchar, (*presiduearray)[rescounter].Natoms, 2);
                        allocate_matrix(char, ((*presiduearray)[rescounter]).chargechar, (*presiduearray)[rescounter].Natoms, 2);
                    }
                }
                else{
                    if(strcmp(chainid, lastchainid)==0){
                        if(resseq<lastresseq){
                            
                            //  Residue indices out of order
                            
                            free_residuearray(&(*presiduearray), (*pNresidues), record_all_data);
							free(linechar);
							free(label);
							free(field);
							free(chainid);
							free(lastchainid);
							free(icode);
							free(lasticode);
							free(atomname);
							free(readres);
							free(mychaindata);
							fclose(inp);
							return 6;
						}
                        else if(resseq>lastresseq){
                            
                            //  New residue on same chain
                            
                            rescounter+=(resseq-lastresseq);
                            if(rescounter==(*pNresidues)){
                                my_exit("Found too many residues in ATOM!\n");
							}
                            j=0;
                            while(strcmp(readres, amino_acid_list[j].resname)!=0){
                                j++;
                                if(j==Naminoacids){
                                    
                                    //  Residue not among the known amino acids
                                    
                                    free_residuearray(&(*presiduearray), (*pNresidues), record_all_data);
                                    free(linechar);
                                    free(label);
                                    free(field);
                                    free(chainid);
                                    free(lastchainid);
                                    free(icode);
                                    free(lasticode);
                                    free(atomname);
                                    free(readres);
                                    free(mychaindata);
                                    fclose(inp);
                                    return 2;
                                }
                            }
                            (*presiduearray)[rescounter].chainid=chaincounter;
                            (*presiduearray)[rescounter].residnumber=j;
                            (*presiduearray)[rescounter].Natoms=amino_acid_list[j].Natoms;
                            allocate_array(int, ((*presiduearray)[rescounter]).atomfound, (*presiduearray)[rescounter].Natoms);
                            for(i=0;i<(*presiduearray)[rescounter].Natoms;i++) (*presiduearray)[rescounter].atomfound[i]=0;
                            allocate_array(double_triple, ((*presiduearray)[rescounter]).atomposition, (*presiduearray)[rescounter].Natoms);
                            allocate_array(char, ((*presiduearray)[rescounter]).chainidchar, 1);
                            (*presiduearray)[rescounter].resseq=resseq;
                            strcpy(((*presiduearray)[rescounter].chainidchar), chainid);
                            if(record_all_data==1){
                                allocate_matrix(char, ((*presiduearray)[rescounter]).serialchar, (*presiduearray)[rescounter].Natoms, 5);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).icodechar, (*presiduearray)[rescounter].Natoms, 1);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).occupancychar, (*presiduearray)[rescounter].Natoms, 6);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).tempfactorchar, (*presiduearray)[rescounter].Natoms, 6);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).elementchar, (*presiduearray)[rescounter].Natoms, 2);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).chargechar, (*presiduearray)[rescounter].Natoms, 2);
                            }
                        }
                        else if(strcmp(icode, lasticode)!=0){
                            rescounter++;
                            if(rescounter==(*pNresidues)){
                                my_exit("Found too many residues in ATOM!\n");
							}
                            j=0;
                            while(strcmp(readres, amino_acid_list[j].resname)!=0){
                                j++;
                                if(j==Naminoacids){
                                    
                                    //  Residue not among the known amino acids

                                    free_residuearray(&(*presiduearray), (*pNresidues), record_all_data);
                                    free(linechar);
                                    free(label);
                                    free(field);
                                    free(chainid);
                                    free(lastchainid);
                                    free(icode);
                                    free(lasticode);
                                    free(atomname);
                                    free(readres);
                                    free(mychaindata);
                                    fclose(inp);
                                    return 2;
                                }
                            }
                            (*presiduearray)[rescounter].chainid=chaincounter;
                            (*presiduearray)[rescounter].residnumber=j;
                            (*presiduearray)[rescounter].Natoms=amino_acid_list[j].Natoms;
                            allocate_array(int, ((*presiduearray)[rescounter]).atomfound, (*presiduearray)[rescounter].Natoms);
                            for(i=0;i<(*presiduearray)[rescounter].Natoms;i++) (*presiduearray)[rescounter].atomfound[i]=0;
                            allocate_array(double_triple, ((*presiduearray)[rescounter]).atomposition, (*presiduearray)[rescounter].Natoms);
                            allocate_array(char, ((*presiduearray)[rescounter]).chainidchar, 1);
                            (*presiduearray)[rescounter].resseq=resseq;
                            strcpy(((*presiduearray)[rescounter].chainidchar), chainid);
                            if(record_all_data==1){
                                allocate_matrix(char, ((*presiduearray)[rescounter]).serialchar, (*presiduearray)[rescounter].Natoms, 5);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).icodechar, (*presiduearray)[rescounter].Natoms, 1);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).occupancychar, (*presiduearray)[rescounter].Natoms, 6);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).tempfactorchar, (*presiduearray)[rescounter].Natoms, 6);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).elementchar, (*presiduearray)[rescounter].Natoms, 2);
                                allocate_matrix(char, ((*presiduearray)[rescounter]).chargechar, (*presiduearray)[rescounter].Natoms, 2);
                            }
                        }
                    }
                    else{
                        
                        //  New residue on new chain
                        
                        chaincounter++;
                        rescounter++;
						if(rescounter==(*pNresidues)){
                            my_exit("Found too many residues in ATOM!\n");
						}
                        j=0;
                        while(strcmp(readres, amino_acid_list[j].resname)!=0){
                            j++;
                            if(j==Naminoacids){
                                
                                //  Residue not among the known amino acids
                                
                                free_residuearray(&(*presiduearray), (*pNresidues), record_all_data);
                                free(linechar);
                                free(label);
                                free(field);
                                free(chainid);
                                free(lastchainid);
                                free(icode);
                                free(lasticode);
                                free(atomname);
                                free(readres);
                                free(mychaindata);
                                fclose(inp);
                                return 2;
                            }
                        }
                        (*presiduearray)[rescounter].Nterminusflag=1;
                        (*presiduearray)[rescounter].chainid=chaincounter;
                        (*presiduearray)[rescounter].residnumber=j;
                        (*presiduearray)[rescounter].Natoms=amino_acid_list[j].Natoms;
                        allocate_array(int, ((*presiduearray)[rescounter]).atomfound, (*presiduearray)[rescounter].Natoms);
                        for(i=0;i<(*presiduearray)[rescounter].Natoms;i++) (*presiduearray)[rescounter].atomfound[i]=0;
                        allocate_array(double_triple, ((*presiduearray)[rescounter]).atomposition, (*presiduearray)[rescounter].Natoms);
                        allocate_array(char, ((*presiduearray)[rescounter]).chainidchar, 1);
                        (*presiduearray)[rescounter].resseq=resseq;
                        strcpy(((*presiduearray)[rescounter].chainidchar), chainid);
                        if(record_all_data==1){
                            allocate_matrix(char, ((*presiduearray)[rescounter]).serialchar, (*presiduearray)[rescounter].Natoms, 5);
                            allocate_matrix(char, ((*presiduearray)[rescounter]).icodechar, (*presiduearray)[rescounter].Natoms, 1);
                            allocate_matrix(char, ((*presiduearray)[rescounter]).occupancychar, (*presiduearray)[rescounter].Natoms, 6);
                            allocate_matrix(char, ((*presiduearray)[rescounter]).tempfactorchar, (*presiduearray)[rescounter].Natoms, 6);
                            allocate_matrix(char, ((*presiduearray)[rescounter]).elementchar, (*presiduearray)[rescounter].Natoms, 2);
                            allocate_matrix(char, ((*presiduearray)[rescounter]).chargechar, (*presiduearray)[rescounter].Natoms, 2);
                        }
                    }
                }
                extract_piece_of_string(linechar, field, 31, 38);
                position.x=atof(field);
                extract_piece_of_string(linechar, field, 39, 46);
                position.y=atof(field);
                extract_piece_of_string(linechar, field, 47, 54);
                position.z=atof(field);
                j=0;
                while(strcmp(atomname, amino_acid_list[(*presiduearray)[rescounter].residnumber].input_atomnames[j])!=0){
                    j++;
                    if(j==(*presiduearray)[rescounter].Natoms){
                        free_residuearray(&(*presiduearray), (*pNresidues), record_all_data);
                        free(linechar);
                        free(label);
                        free(field);
                        free(chainid);
                        free(lastchainid);
                        free(icode);
                        free(lasticode);
                        free(atomname);
                        free(readres);
                        free(mychaindata);
                        fclose(inp);
                        return 3;
                    }
                }
                if((*presiduearray)[rescounter].atomfound[j]==1){
                    printf("Found two atoms of same type %s (%i) in ATOM, residue #%i (%s)!\n", amino_acid_list[(*presiduearray)[rescounter].residnumber].input_atomnames[j], j, rescounter, readres);
                    exit(1);
                }
                else{
                    (*presiduearray)[rescounter].atomfound[j]=1;
                    (*presiduearray)[rescounter].atomposition[j]=position;
                    if(record_all_data==1){
                        extract_piece_of_string(linechar, field, 7, 11);
                        strcpy(((*presiduearray)[rescounter]).serialchar[j], field);
                        extract_piece_of_string(linechar, field, 27, 27);
                        strcpy(((*presiduearray)[rescounter]).icodechar[j], field);
                        extract_piece_of_string(linechar, field, 55, 60);
                        strcpy(((*presiduearray)[rescounter]).occupancychar[j], field);
                        extract_piece_of_string(linechar, field, 61, 66);
                        strcpy(((*presiduearray)[rescounter]).tempfactorchar[j], field);
                        extract_piece_of_string(linechar, field, 77, 78);
                        strcpy(((*presiduearray)[rescounter]).elementchar[j], field);
                        extract_piece_of_string(linechar, field, 79, 80);
                        strcpy(((*presiduearray)[rescounter]).chargechar[j], field);
                    }
                }
                strcpy(lastchainid, chainid);
                strcpy(lasticode, icode);
                foundatom=1;
				lastresseq=resseq;
            }
        }
        read=mygetline(linechar, maxlinelength, inp);
    }
	free(linechar);
	free(label);
    free(field);
    free(chainid);
    free(lastchainid);
    free(icode);
    free(lasticode);
	free(atomname);
	free(readres);
	free(mychaindata);
    fclose(inp);
	if(foundatom==0){
		return 5;
	}
    return 0;
}

void calculate_Cterminus_bond_lengths_and_angles(amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray, double **Cterminus_avs, int *pCterminus_count){
    int i, aminoacidindex, Cindex, Oindex, OXTindex;
    double normCO, normCOXT, angle;
    double_triple CO, COXT;
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>-1){
            if((i==Nresidues-1)||(residuearray[i].chainid!=residuearray[i+1].chainid)){
                Cindex=amino_acid_list[aminoacidindex].Cindex;
                Oindex=amino_acid_list[aminoacidindex].Oindex;
                OXTindex=amino_acid_list[aminoacidindex].OXTindex;
                if((residuearray[i].atomfound[Cindex]==1)&&(residuearray[i].atomfound[Oindex]==1)&&(residuearray[i].atomfound[OXTindex]==1)){
                    CO=subtract_double_triple(residuearray[i].atomposition[Oindex], residuearray[i].atomposition[Cindex]);
                    COXT=subtract_double_triple(residuearray[i].atomposition[OXTindex], residuearray[i].atomposition[Cindex]);
                    normCO=norm(CO);
                    normCOXT=norm(COXT);
                    Cterminus_avs[0][0]+=normCO;
                    Cterminus_avs[0][1]+=normCO*normCO;
                    Cterminus_avs[1][0]+=normCOXT;
                    Cterminus_avs[1][1]+=normCOXT*normCOXT;
                    CO=scalar_multiply_double_triple(CO, 1./normCO);
                    COXT=scalar_multiply_double_triple(COXT, 1./normCOXT);
                    angle=acos(dot_product(CO, COXT));
                    Cterminus_avs[2][0]+=angle;
                    Cterminus_avs[2][1]+=angle*angle;
                    (*pCterminus_count)++;
                }
            }
        }
    }
}

void output_Cterminus_bond_lengths_and_angles(char *filename, double **Cterminus_avs, int Cterminus_count){
    int i;
    for(i=0;i<3;i++){
        Cterminus_avs[i][0]/=(1.*Cterminus_count);
        Cterminus_avs[i][1]/=(1.*Cterminus_count);
        Cterminus_avs[i][1]=sqrt(Cterminus_avs[i][1]-Cterminus_avs[i][0]*Cterminus_avs[i][0]);
    }
    FILE *outp;
    outp=fopen(filename, "w");
    fprintf(outp, "C-O bond\t%.8f\t%.8f\n", Cterminus_avs[0][0], Cterminus_avs[0][1]);
    fprintf(outp, "C-OXT bond\t%.8f\t%.8f\n", Cterminus_avs[1][0], Cterminus_avs[1][1]);
    fprintf(outp, "O-C-OXT angle\t%.8f\t%.8f\n", Cterminus_avs[2][0]*180./M_PI, Cterminus_avs[2][1]*180./M_PI);
    fclose(outp);
}

void calculate_average_radial_distances(amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray){
    int i, j, k, aminoacidindex, foundall, index;
    double radial, radialsum, longitudinal, longitudinalsum;
    double_triple first, second, director, vector;
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>-1){
            if(amino_acid_list[aminoacidindex].N_radial_sets>0){
                for(j=0;j<amino_acid_list[aminoacidindex].N_radial_sets;j++){
                    index=amino_acid_list[aminoacidindex].radial_atoms[j][0];
                    if(residuearray[i].atomfound[index]==1){
                        first=residuearray[i].atomposition[index];
                    }
                    else break;
                    index=amino_acid_list[aminoacidindex].radial_atoms[j][1];
                    if(residuearray[i].atomfound[index]==1){
                        second=residuearray[i].atomposition[index];
                        director=subtract_double_triple(second, first);
                    }
                    else break;
                    normalize(&director);
                    foundall=1;
                    for(k=2;k<amino_acid_list[aminoacidindex].radial_set_size[j];k++){
                        index=amino_acid_list[aminoacidindex].radial_atoms[j][k];
                        if(residuearray[i].atomfound[index]==1){
                            vector=subtract_double_triple(residuearray[i].atomposition[index], second);
                            longitudinal=dot_product(vector, director);
                            radial=norm(subtract_double_triple(vector, scalar_multiply_double_triple(director, longitudinal)));
                            if(k==2){
                                radialsum=radial;
                                longitudinalsum=longitudinal;
                            }
                            else{
                                radialsum+=radial;
                                longitudinalsum+=longitudinal;
                            }
                        }
                        else{
                            foundall=0;
                            break;
                        }
                    }
                    if(foundall==1){
                        amino_acid_list[aminoacidindex].radial_average[j].x+=longitudinalsum;
                        amino_acid_list[aminoacidindex].radial_average[j].y+=radialsum;
                        amino_acid_list[aminoacidindex].radial_count[j]+=(amino_acid_list[aminoacidindex].radial_set_size[j]-2);
                    }
                }
            }
        }
    }
}

void output_error_filename(char *file, char *base, int type, int *pfirst){
    int i;
    declare_array_nozero(char, outputfile, maxstringlength);
    if(*pfirst==0){
        sprintf(outputfile, "%s.ensemble", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.non_amino_acid", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.unknown_atom", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.no_DBREF", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.residue_indicies_out_of_order", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.terminal_atom_in_non_terminal_residue", base);
        delete_file_if_exists(outputfile);
        sprintf(outputfile, "%s.file_not_found", base);
        delete_file_if_exists(outputfile);
    }
    FILE *outp;
    (*pfirst)=1;
    if(type==1) sprintf(outputfile, "%s.ensemble", base);
    else if(type==2) sprintf(outputfile, "%s.non_amino_acid", base);
    else if(type==3) sprintf(outputfile, "%s.unknown_atom", base);
    else if(type==4) sprintf(outputfile, "%s.no_DBREF", base);
    else if(type==5) sprintf(outputfile, "%s.no_DBREF", base);
    else if(type==6) sprintf(outputfile, "%s.residue_indicies_out_of_order", base);
    else if(type==7) sprintf(outputfile, "%s.terminal_atom_in_non_terminal_residue", base);
    else if(type==8) sprintf(outputfile, "%s.file_not_found", base);
    else my_exit("Bad type in output_error_filename");
    outp=fopen(outputfile, "a");
    fprintf(outp, "%s\n", file);
    fclose(outp);
    free(outputfile);
}

void reassign_Nterminus_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite, int record_all_data){
    int i, j, aminoacidindex, N_Nterminal_Hs, newindex, Hcounter, offset;
    FILE *outp;
    declare_array_nozero(double_triple, new_positions, 3);
    char **new_serialchar, **new_icodechar, **new_occupancychar, **new_tempfactorchar, **new_elementchar, **new_chargechar;
    if(record_all_data==1){
        allocate_matrix(char, new_serialchar, 3, 5);
        allocate_matrix(char, new_icodechar, 3, 1);
        allocate_matrix(char, new_occupancychar, 3, 6);
        allocate_matrix(char, new_tempfactorchar, 3, 6);
        allocate_matrix(char, new_elementchar, 3, 2);
        allocate_matrix(char, new_chargechar, 3, 2);
    }
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if((i==0)||(residuearray[i].Nterminusflag==1)){
            N_Nterminal_Hs=0;
            for(j=0;j<amino_acid_list[aminoacidindex].aminoHindex;j++){
                if(residuearray[i].atomfound[j]==1){
                    if(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H  ")==0) N_Nterminal_Hs++;
                    if(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H2 ")==0) N_Nterminal_Hs++;
                }
            }
            for(j=amino_acid_list[aminoacidindex].aminoHindex;j<amino_acid_list[aminoacidindex].Natoms;j++){
                if(residuearray[i].atomfound[j]==1) N_Nterminal_Hs++;
            }
            if(N_Nterminal_Hs>0){
                if(N_Nterminal_Hs==1){
                    for(j=amino_acid_list[aminoacidindex].Natoms-number_amino_hydrogen_atom_types;j<amino_acid_list[aminoacidindex].Natoms;j++){
                        if(residuearray[i].atomfound[j]==1){
                            residuearray[i].atomfound[j]=0;
                            if(dowrite==1){
                                outp=fopen(output_file, "a");
                                fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%.8f\t%.8f\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[j], residuearray[i].atomposition[j].x, residuearray[i].atomposition[j].y, residuearray[i].atomposition[j].z);
                                fclose(outp);
                            }
                        }
                    }
                }
                else if(N_Nterminal_Hs>3){
                    printf("%i terminal Hs found!\n", N_Nterminal_Hs);
                    exit(1);
                }
                else{
                    Hcounter=0;
                    for(j=0;j<amino_acid_list[aminoacidindex].Natoms-number_amino_hydrogen_atom_types;j++){
                        if(residuearray[i].atomfound[j]==1){
                            if(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H  ")==0){
                                new_positions[Hcounter]=residuearray[i].atomposition[j];
                                if(record_all_data==1){
                                    strcpy(new_serialchar[Hcounter], residuearray[i].serialchar[j]);
                                    strcpy(new_icodechar[Hcounter], residuearray[i].icodechar[j]);
                                    strcpy(new_occupancychar[Hcounter], residuearray[i].occupancychar[j]);
                                    strcpy(new_tempfactorchar[Hcounter], residuearray[i].tempfactorchar[j]);
                                    strcpy(new_elementchar[Hcounter], residuearray[i].elementchar[j]);
                                    strcpy(new_chargechar[Hcounter], residuearray[i].chargechar[j]);
                                }
                                residuearray[i].atomfound[j]=0;
                                Hcounter++;
                            }
                            if(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H2 ")==0){
                                new_positions[Hcounter]=residuearray[i].atomposition[j];
                                if(record_all_data==1){
                                    strcpy(new_serialchar[Hcounter], residuearray[i].serialchar[j]);
                                    strcpy(new_icodechar[Hcounter], residuearray[i].icodechar[j]);
                                    strcpy(new_occupancychar[Hcounter], residuearray[i].occupancychar[j]);
                                    strcpy(new_tempfactorchar[Hcounter], residuearray[i].tempfactorchar[j]);
                                    strcpy(new_elementchar[Hcounter], residuearray[i].elementchar[j]);
                                    strcpy(new_chargechar[Hcounter], residuearray[i].chargechar[j]);
                                }
                                residuearray[i].atomfound[j]=0;
                                Hcounter++;
                            }
                        }
                    }
                    for(j=amino_acid_list[aminoacidindex].Natoms-number_amino_hydrogen_atom_types;j<amino_acid_list[aminoacidindex].Natoms;j++){
                        if(residuearray[i].atomfound[j]==1){
                            new_positions[Hcounter]=residuearray[i].atomposition[j];
                            if(record_all_data==1){
                                strcpy(new_serialchar[Hcounter], residuearray[i].serialchar[j]);
                                strcpy(new_icodechar[Hcounter], residuearray[i].icodechar[j]);
                                strcpy(new_occupancychar[Hcounter], residuearray[i].occupancychar[j]);
                                strcpy(new_tempfactorchar[Hcounter], residuearray[i].tempfactorchar[j]);
                                strcpy(new_elementchar[Hcounter], residuearray[i].elementchar[j]);
                                strcpy(new_chargechar[Hcounter], residuearray[i].chargechar[j]);
                            }
                            residuearray[i].atomfound[j]=0;
                            Hcounter++;
                        }
                    }
                    
                    offset=amino_acid_list[aminoacidindex].aminoHindex;
                    if(N_Nterminal_Hs==3) offset+=2;
                    if(strcmp(amino_acid_list[aminoacidindex].resname, "GLY")==0) offset+=GLYoffset;
                    if(strcmp(amino_acid_list[aminoacidindex].resname, "PRO")==0) offset+=PROoffset;
                    for(Hcounter=0;Hcounter<N_Nterminal_Hs;Hcounter++){
                        newindex=offset+Hcounter;
                        if(newindex==amino_acid_list[aminoacidindex].Natoms) my_exit("index too big in reassign!");
                        if(residuearray[i].atomfound[newindex]==1){
                            printf("In reassign %s (%i=%i+%i; %s) already filled!", amino_acid_list[aminoacidindex].atomnames[newindex], newindex, amino_acid_list[aminoacidindex].aminoHindex, Hcounter, amino_acid_list[aminoacidindex].resname);
                            exit(1);
                        }
                        residuearray[i].atomfound[newindex]=1;
                        residuearray[i].atomposition[newindex]=new_positions[Hcounter];
                        if(record_all_data==1){
                            strcpy((residuearray[i].serialchar[newindex]), new_serialchar[Hcounter]);
                            strcpy((residuearray[i].icodechar[newindex]), new_icodechar[Hcounter]);
                            strcpy((residuearray[i].occupancychar[newindex]), new_occupancychar[Hcounter]);
                            strcpy((residuearray[i].tempfactorchar[newindex]), new_tempfactorchar[Hcounter]);
                            strcpy((residuearray[i].elementchar[newindex]), new_elementchar[Hcounter]);
                            strcpy((residuearray[i].chargechar[newindex]), new_chargechar[Hcounter]);
                        }
                    }
                    
                }
            }
        }
    }
    if(record_all_data==1){
        free(new_positions);
        free_matrix(new_serialchar, 3);
        free_matrix(new_icodechar, 3);
        free_matrix(new_occupancychar, 3);
        free_matrix(new_tempfactorchar, 3);
        free_matrix(new_elementchar, 3);
        free_matrix(new_chargechar, 3);
    }
}

int assign_atomfound_terminii(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list){
	
    //  Check that right- and left-terminal atom types are found only where expected
	
    int i, j, k;
    declare_matrix_nozero(char, amino_H_names, number_amino_hydrogen_atom_types, 4);
    strcpy(amino_H_names[0], " HN1");
    strcpy(amino_H_names[1], " HN2");
    strcpy(amino_H_names[2], " HM1");
    strcpy(amino_H_names[3], " HM2");
    strcpy(amino_H_names[4], " HM3");
    strcpy(amino_H_names[5], "HGN1");
    strcpy(amino_H_names[6], "HGN2");
    strcpy(amino_H_names[7], "HGM1");
    strcpy(amino_H_names[8], "HGM2");
    strcpy(amino_H_names[9], "HGM3");
    strcpy(amino_H_names[10], "HPN1");
    strcpy(amino_H_names[11], "HPN2");
    strcpy(amino_H_names[12], "HPM1");
    strcpy(amino_H_names[13], "HPM2");
    strcpy(amino_H_names[14], "HPM3");
	for(i=0;i<Nresidues;i++){
		for(j=0;j<residuearray[i].Natoms;j++){
            if(strcmp(amino_acid_list[residuearray[i].residnumber].resname, "ASP")==0){
                if(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " HD2")==0){
                    if(residuearray[i].atomfound[j]==0){
                        residuearray[i].atomfound[j]=3;
                    }
                }
            }
            if(strcmp(amino_acid_list[residuearray[i].residnumber].resname, "GLU")==0){
                if(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " HE2")==0){
                    if(residuearray[i].atomfound[j]==0){
                        residuearray[i].atomfound[j]=3;
                    }
                }
            }
            if(i>0){
                if(residuearray[i-1].chainid==residuearray[i].chainid){
                    if(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " H2 ")==0){
                        if(residuearray[i].atomfound[j]==1){
                            return 1;
                        }
                        residuearray[i].atomfound[j]=2;
                    }
                    if(strcmp(amino_acid_list[residuearray[i].residnumber].resname, "PRO")==0){
                        if(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " H  ")==0){
                            if(residuearray[i].atomfound[j]==1){
                                return 1;
                            }
                            residuearray[i].atomfound[j]=2;
                        }
                    }
                    for(k=0;k<number_amino_hydrogen_atom_types;k++){
                        if(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], amino_H_names[k])==0){
                            if(residuearray[i].atomfound[j]==1){
                                return 1;
                            }
                            residuearray[i].atomfound[j]=2;
                        }
                    }
                }
            }
            if(i<Nresidues-1){
                if(residuearray[i].chainid==residuearray[i+1].chainid){
                    if((strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " OXT")==0)||(strcmp(amino_acid_list[residuearray[i].residnumber].atomnames[j], " HXT")==0)){
                        if(residuearray[i].atomfound[j]==1){
                            return 1;
                        }
                        residuearray[i].atomfound[j]=2;
                    }
                }
            }
		}
	}
    free_matrix(amino_H_names, number_amino_hydrogen_atom_types);
    return 0;
}

void discard_overlapping_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite, atomlookup *overlapping_atoms_list, int *pN_overlapping_atoms){
    int i, j, k, aminoacidindex, nextaminoacidindex, Nindex;
    double_triple pos;
    FILE *outp;
    if(dowrite==1){
        for(i=0;i<Nresidues;i++){
            if(residuearray[i].chainid>=0){
                aminoacidindex=residuearray[i].residnumber;
                for(j=0;j<residuearray[i].Natoms;j++){
                    if(residuearray[i].atomfound[j]==1){
                        pos=residuearray[i].atomposition[j];
                        
                        //  Check bonds within residue
                        
                        for(k=0;k<j;k++){
                            if(residuearray[i].atomfound[k]==1){
                                if(residuearray[i].atomposition[k].x==pos.x){
                                    if(residuearray[i].atomposition[k].y==pos.y){
                                        if(residuearray[i].atomposition[k].z==pos.z){
                                            residuearray[i].atomfound[j]=0;
                                            residuearray[i].atomfound[k]=0;
                                            if((*pN_overlapping_atoms)+2>max_overlapping_atoms) my_exit("Too many overlapping atoms");
                                            strcpy((overlapping_atoms_list[(*pN_overlapping_atoms)].filename), pdb_file);
                                            overlapping_atoms_list[(*pN_overlapping_atoms)].residue=i;
                                            overlapping_atoms_list[(*pN_overlapping_atoms)].atom=j;
                                            (*pN_overlapping_atoms)++;
                                            strcpy((overlapping_atoms_list[(*pN_overlapping_atoms)].filename), pdb_file);
                                            overlapping_atoms_list[(*pN_overlapping_atoms)].residue=i;
                                            overlapping_atoms_list[(*pN_overlapping_atoms)].atom=k;
                                            (*pN_overlapping_atoms)++;
                                            outp=fopen(output_file, "a");
                                            fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[j], amino_acid_list[aminoacidindex].atomnames[k], residuearray[i].atomposition[j].x, residuearray[i].atomposition[j].y, residuearray[i].atomposition[j].z, residuearray[i].atomposition[k].x, residuearray[i].atomposition[k].y, residuearray[i].atomposition[k].z);
                                            fclose(outp);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else{
        for(i=0;i<(*pN_overlapping_atoms);i++){
            if(strcmp(pdb_file, overlapping_atoms_list[i].filename)==0){
                residuearray[overlapping_atoms_list[i].residue].atomfound[overlapping_atoms_list[i].atom]=0;
            }
        }
    }
}

void setup_bond_length_structure(char *covalent_bond_lengths_file, int ***pfoundbond, double ***pmaxbondlength, double ***pminbondlength){
    int read, stringsize, i, j;
    double minlength, maxlength;
    FILE *inp;
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, field, 4);
    allocate_matrix(int, (*pfoundbond), maxelement, maxelement);
    allocate_matrix(double, (*pmaxbondlength), maxelement, maxelement);
    allocate_matrix(double, (*pminbondlength), maxelement, maxelement);
    inp=fopen(covalent_bond_lengths_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        stringsize=strlen(linechar);
        if(stringsize!=13){
            printf("linesize %i not 9 in %s!\n", stringsize, covalent_bond_lengths_file);
            exit(1);
        }
        extract_piece_of_string(linechar, field, 1, 2);
        i=atoi(field);
        if((i<0)||(i>=maxelement)){
            printf("element=%i in %s!\n", i, covalent_bond_lengths_file);
            exit(1);
        }
        extract_piece_of_string(linechar, field, 3, 4);
        j=atoi(field);
        if((j<0)||(j>=maxelement)){
            printf("element=%i in %s!\n", j, covalent_bond_lengths_file);
            exit(1);
        }
        extract_piece_of_string(linechar, field, 5, 8);
        minlength=atof(field);
        extract_piece_of_string(linechar, field, 9, 12);
        maxlength=atof(field);
        if(((*pfoundbond)[i][j]==1)||((*pfoundbond)[j][i]==1)){
            printf("duplicate bonds in %s!\n", covalent_bond_lengths_file);
        }
        (*pfoundbond)[i][j]=1;
        (*pfoundbond)[j][i]=1;
        (*pmaxbondlength)[i][j]=maxbondfactor*maxlength;
        (*pmaxbondlength)[j][i]=maxbondfactor*maxlength;
        (*pminbondlength)[i][j]=minbondfactor*minlength;
        (*pminbondlength)[j][i]=minbondfactor*minlength;
        read=mygetline(linechar, maxlinelength, inp);
    }
    fclose(inp);
    free(linechar);
    free(field);
}

void check_covalent_bond_lengths(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int **foundbond, double **maxbondlength, double **minbondlength, char *output_file, char *pdb_file, int dowrite){
    int i, j, k, aminoacidindex, nextaminoacidindex, neighbor, element1, element2, Nindex;
    double distance;
    FILE *outp;
    for(i=0;i<Nresidues;i++){
        if(residuearray[i].chainid>=0){
            aminoacidindex=residuearray[i].residnumber;
            for(j=0;j<residuearray[i].Natoms;j++){
                if(residuearray[i].atomfound[j]==1){
                    element1=amino_acid_list[aminoacidindex].element[j];
                    
                    //  Check bonds within residue
                    
                    for(k=0;k<amino_acid_list[aminoacidindex].nconect[j];k++){
                        neighbor=amino_acid_list[aminoacidindex].conect[j][k];
                        if(residuearray[i].atomfound[neighbor]==1){
                            element2=amino_acid_list[aminoacidindex].element[neighbor];
                            if(foundbond[element1][element2]==0){
                                printf("no bond found between elements %i and %i!\n", element1, element2);
                                printf("\tresidue %s, atoms %s and %s\n", amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[j], amino_acid_list[aminoacidindex].atomnames[neighbor]);
                                printf("\tamino_acid_list[%i].conect[%i][%i of %i]=%i\n", aminoacidindex, j, k, amino_acid_list[aminoacidindex].nconect[j], neighbor);
                                exit(1);
                            }
                            distance=norm(subtract_double_triple(residuearray[i].atomposition[j], residuearray[i].atomposition[neighbor]));
                            if((distance>maxbondlength[element1][element2])||(distance<minbondlength[element1][element2])){
                                residuearray[i].atomfound[j]=0;
                                residuearray[i].atomfound[neighbor]=0;
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[j], amino_acid_list[aminoacidindex].atomnames[neighbor], distance, residuearray[i].atomposition[j].x, residuearray[i].atomposition[j].y, residuearray[i].atomposition[j].z, residuearray[i].atomposition[neighbor].x, residuearray[i].atomposition[neighbor].y, residuearray[i].atomposition[neighbor].z);
                                    fclose(outp);
                                }
                            }
                        }
                    }
                    
                    //  Check bonds across residues
                    
                    if(j==amino_acid_list[aminoacidindex].Cindex){
                        if((i<Nresidues-1)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                            nextaminoacidindex=residuearray[i+1].residnumber;
                            Nindex=amino_acid_list[nextaminoacidindex].Nindex;
                            if(residuearray[i+1].atomfound[Nindex]==1){
                                element2=amino_acid_list[nextaminoacidindex].element[Nindex];
                                if(foundbond[element1][element2]==0){
                                    printf("no bond found between elements %i and %i!\n", element1, element2);
                                    exit(1);
                                }
                                distance=norm(subtract_double_triple(residuearray[i].atomposition[j], residuearray[i+1].atomposition[Nindex]));
                                if((distance>maxbondlength[element1][element2])||(distance<minbondlength[element1][element2])){
                                    residuearray[i].atomfound[j]=0;
                                    residuearray[i+1].atomfound[Nindex]=0;
                                    if(dowrite==1){
                                        outp=fopen(output_file, "a");
                                        fprintf(outp, "%s\t%s-%s\t%i-%i\t%s-%s\t%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i+1].chainidchar, residuearray[i].resseq, residuearray[i+1].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[nextaminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[j], amino_acid_list[nextaminoacidindex].atomnames[Nindex], distance, residuearray[i].atomposition[j].x, residuearray[i].atomposition[j].y, residuearray[i].atomposition[j].z, residuearray[i+1].atomposition[Nindex].x, residuearray[i+1].atomposition[Nindex].y, residuearray[i+1].atomposition[Nindex].z);
                                        fclose(outp);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void setup_bond_angles_structure(char *bond_angles_file, amino_acid_struct *amino_acid_list, int **pNangles, double ***pangles, int ****pangle_index){
    int read, i, j, k, l, stringsize, foundconect;
    FILE *inp;
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, field, 5);
    allocate_array(int, (*pNangles), Naminoacids);
    allocate_matrix(double, (*pangles), Naminoacids, max_angles_per_residue);
    allocate_3d_tensor(int, (*pangle_index), Naminoacids, max_angles_per_residue, 3);
    inp=fopen(bond_angles_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        stringsize=strlen(linechar);
        if(stringsize!=21){
            printf("linesize %i not 21 in %s!\n", stringsize, bond_angles_file);
            exit(1);
        }
        extract_piece_of_string(linechar, field, 1, 3);
        for(i=0;i<Naminoacids;i++){
            if(strcmp(field, amino_acid_list[i].resname)==0){
                break;
            }
        }
        if(i==Naminoacids){
            printf("Residue %s in %s!\n", field, bond_angles_file);
            exit(1);
        }
        for(j=0;j<3;j++){
            extract_piece_of_string(linechar, field, 4*(j+1), 4*(j+1)+3);
            for(k=0;k<amino_acid_list[i].Natoms;k++){
                if(strcmp(field, amino_acid_list[i].atomnames[k])==0){
                    break;
                }
            }
            if(k==amino_acid_list[i].Natoms){
                printf("Atom <%s> in residue %s in %s!\n", field, amino_acid_list[i].resname, bond_angles_file);
                exit(1);
            }
            (*pangle_index)[i][(*pNangles)[i]][j]=k;
        }
        for(k=0;k<3;k+=2){
            foundconect=0;
            for(l=0;l<amino_acid_list[i].nconect[(*pangle_index)[i][(*pNangles)[i]][1]];l++){
                if(amino_acid_list[i].conect[(*pangle_index)[i][(*pNangles)[i]][1]][l]==(*pangle_index)[i][(*pNangles)[i]][k]){
                    foundconect=1;
                }
            }
            if(foundconect==0){
                printf("Bond angle involving non-bonded atoms %s and %s in %s!\n", amino_acid_list[i].atomnames[(*pangle_index)[i][(*pNangles)[i]][1]], amino_acid_list[i].atomnames[(*pangle_index)[i][(*pNangles)[i]][k]], amino_acid_list[i].resname);
                printf("Line: %s", linechar);
                exit(1);
            }
        }
        extract_piece_of_string(linechar, field, 16, 20);
        (*pangles)[i][(*pNangles)[i]]=atof(field)*M_PI/180.;
        (*pNangles)[i]++;
        if((*pNangles)[i]==max_angles_per_residue){
            my_exit("Too many angles\n");
        }
        read=mygetline(linechar, maxlinelength, inp);
    }
    fclose(inp);
    free(linechar);
    free(field);
}

void check_bond_angles(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int *Nangles, double **angles, int ***angle_index, char *output_file, char *pdb_file, int dowrite){
    int i, j, k, aminoacidindex, nextaminoacidindex, foundall;
    double angle;
    FILE *outp;
    double_triple director1, director2;
    for(i=0;i<Nresidues;i++){
        if(residuearray[i].chainid>=0){
            aminoacidindex=residuearray[i].residnumber;
            
            //  Check angles between residues
            
            if((i<Nresidues-1)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                nextaminoacidindex=residuearray[i+1].residnumber;
                if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].Cindex]==1){
                    if(residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Nindex]==1){
                        director1=subtract_double_triple(residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex]);
                        normalize(&director1);
                        if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].Oindex]==1){
                            director2=subtract_double_triple(residuearray[i].atomposition[amino_acid_list[aminoacidindex].Oindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex]);
                            normalize(&director2);
                            angle=acos(dot_product(director1, director2));
                            if((isnan(angle)==1)||(fabs(angle-O_C_N_angle*M_PI/180.)>maxangledifference)){
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].Cindex]=0;
                                residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Nindex]=0;
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].Oindex]=0;
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s\t%.4f\t%.1f", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, angle/M_PI*180., O_C_N_angle);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].Cindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].x, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].y, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[nextaminoacidindex].atomnames[amino_acid_list[nextaminoacidindex].Nindex], residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].x, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].y, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].Oindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Oindex].x, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Oindex].y, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Oindex].z);
                                    fprintf(outp, "\n");
                                    fclose(outp);
                                }
                            }
                        }
                        if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].CAindex]==1){
                            director2=subtract_double_triple(residuearray[i].atomposition[amino_acid_list[aminoacidindex].CAindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex]);
                            normalize(&director2);
                            angle=acos(dot_product(director1, director2));
                            if((isnan(angle)==1)||(fabs(angle-CA_C_N_angle*M_PI/180.)>maxangledifference)){
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].Cindex]=0;
                                residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Nindex]=0;
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].CAindex]=0;
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s\t%.4f\t%.1f", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, angle/M_PI*180., CA_C_N_angle);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].Cindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].x, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].y, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[nextaminoacidindex].atomnames[amino_acid_list[nextaminoacidindex].Nindex], residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].x, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].y, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].CAindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].CAindex].x, residuearray[i].atomposition[amino_acid_list[aminoacidindex].CAindex].y, residuearray[i].atomposition[amino_acid_list[aminoacidindex].CAindex].z);
                                    fprintf(outp, "\n");
                                    fclose(outp);
                                }
                            }
                        }
                        if(residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Hindex]==1){
                            director2=subtract_double_triple(residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex], residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Hindex]);
                            normalize(&director2);
                            angle=acos(dot_product(director1, director2));
                            if((isnan(angle)==1)||(fabs(angle-C_N_H_angle*M_PI/180.)>maxangledifference)){
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].Cindex]=0;
                                residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Nindex]=0;
                                residuearray[i+1].atomfound[amino_acid_list[nextaminoacidindex].Hindex]=0;
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s\t%.4f\t%.1f", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, angle/M_PI*180., C_N_H_angle);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].Cindex], residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].x, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].y, residuearray[i].atomposition[amino_acid_list[aminoacidindex].Cindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[nextaminoacidindex].atomnames[amino_acid_list[nextaminoacidindex].Nindex], residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].x, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].y, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Nindex].z);
                                    fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[nextaminoacidindex].atomnames[amino_acid_list[nextaminoacidindex].Hindex], residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Hindex].x, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Hindex].y, residuearray[i+1].atomposition[amino_acid_list[nextaminoacidindex].Hindex].z);
                                    fprintf(outp, "\n");
                                    fclose(outp);
                                }
                            }
                        }
                    }
                }
            }
            
            //  Check angles within residue
            
            for(j=0;j<Nangles[aminoacidindex];j++){
                foundall=1;
                for(k=0;k<3;k++){
                    if(residuearray[i].atomfound[angle_index[aminoacidindex][j][k]]!=1){
                        foundall=0;
                        break;
                    }
                }
                if(foundall==1){
                    director1=subtract_double_triple(residuearray[i].atomposition[angle_index[aminoacidindex][j][0]], residuearray[i].atomposition[angle_index[aminoacidindex][j][1]]);
                    director2=subtract_double_triple(residuearray[i].atomposition[angle_index[aminoacidindex][j][2]], residuearray[i].atomposition[angle_index[aminoacidindex][j][1]]);
                    normalize(&director1);
                    normalize(&director2);
                    angle=acos(dot_product(director1, director2));
                    if((isnan(angle)==1)||(fabs(angle-angles[aminoacidindex][j])>maxangledifference)){
                        for(k=0;k<3;k++){
                            residuearray[i].atomfound[angle_index[aminoacidindex][j][k]]=0;
                        }
                        if(dowrite==1){
                            outp=fopen(output_file, "a");
                            fprintf(outp, "%s\t%s\t%i\t%s\t%.4f\t%.1f", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, angle/M_PI*180., angles[aminoacidindex][j]/M_PI*180.);
                            for(k=0;k<3;k++){
                                fprintf(outp, "\t%s\t%.8f\t%.8f\t%.8f", amino_acid_list[aminoacidindex].atomnames[angle_index[aminoacidindex][j][k]], residuearray[i].atomposition[angle_index[aminoacidindex][j][k]].x, residuearray[i].atomposition[angle_index[aminoacidindex][j][k]].y, residuearray[i].atomposition[angle_index[aminoacidindex][j][k]].z);
                            }
                            fprintf(outp, "\n");
                            fclose(outp);
                        }
                    }
                }
            }
        }
    }
}

void setup_tetrahedra_structure(char *tetrahedral_file, amino_acid_struct *amino_acid_list, int **pNtetrahedra, int ****ptetrahedraindices){
    FILE *inp;
    inp=fopen(tetrahedral_file, "r");
    declare_array_nozero(char, linechar, maxlinelength);
    int i, j, k, read, stringsize;
    char readres[3], atomname[4];
    allocate_array(int, (*pNtetrahedra), Naminoacids);
    allocate_3d_tensor(int, (*ptetrahedraindices), Naminoacids, max_tetrahedra_per_residue, 5);
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        stringsize=strlen(linechar);
        if(stringsize!=24){
            printf("line in %s of length %i\n", tetrahedral_file, stringsize);
            exit(1);
        }
        extract_piece_of_string(linechar, readres, 1, 3);
        for(i=0;i<Naminoacids;i++){
            if(strcmp(readres, amino_acid_list[i].resname)==0) break;
        }
        if(i==Naminoacids){
            printf("residue <%s> in %s not found among amino acids!\n", readres, tetrahedral_file);
            exit(1);
        }
        for(j=0;j<5;j++){
            extract_piece_of_string(linechar, atomname, 4*(j+1), 4*(j+1)+3);
            for(k=0;k<amino_acid_list[i].Natoms;k++){
                if(strcmp(atomname, amino_acid_list[i].atomnames[k])==0) break;
            }
            if(k==amino_acid_list[i].Natoms){
                printf("atom <%s> in %s not found among atom names for residue <%s>!\n", atomname, tetrahedral_file, readres);
                exit(1);
            }
            (*ptetrahedraindices)[i][(*pNtetrahedra)[i]][j]=k;
        }
        (*pNtetrahedra)[i]++;
        read=mygetline(linechar, maxlinelength, inp);
    }
	fclose(inp);
    free(linechar);
}

void discard_bad_tetrahedra(residuedata *residuearray, int Nresidues, int *Ntetrahedra, int ***tetrahedraindices, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite){
    FILE *outp;
    int i, j, k, aminoacidindex, foundall;
    double out_of_plane_angle;
    double_triple first, second, normal, third;
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>=0){
            if(Ntetrahedra[aminoacidindex]>0){
                for(j=0;j<Ntetrahedra[aminoacidindex];j++){
                    foundall=1;
                    for(k=0;k<4;k++){
                        if(residuearray[i].atomfound[tetrahedraindices[aminoacidindex][j][k]]==0) foundall=0;
                    }
                    if(foundall==1){
                        
                        //  If all heavy atoms are found, check to make sure third vector is at least 30 degrees out of the plane defined by the other two
                        //  And with correct chirality
                        
                        first=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][1]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][0]]);
                        normalize(&first);
                        second=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][2]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][0]]);
                        normalize(&second);
                        third=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][3]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][0]]);
                        normalize(&third);
                        normal=normed(cross_product(first, second));
                        out_of_plane_angle=asin(dot_product(normal, third));
                        if(out_of_plane_angle<minimum_tetrahedral_angle){
							if(dowrite==1){
 								outp=fopen(output_file, "a");
								fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][0]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][1]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][2]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][3]], out_of_plane_angle/M_PI*180.);
                                fclose(outp);
                            }
							for(k=0;k<5;k++){
                                residuearray[i].atomfound[tetrahedraindices[aminoacidindex][j][k]]=0;
                            }
                        }
                        else if(residuearray[i].atomfound[tetrahedraindices[aminoacidindex][j][4]]==1){
                            
                            //  If the hydrogen is also found, check to make sure it is on the opposite side of the chiral center as the heavy atoms
                            
                            first=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][2]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][1]]);
                            normalize(&first);
                            second=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][3]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][1]]);
                            normalize(&second);
                            normal=normed(cross_product(first, second));
                            second=subtract_double_triple(residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][4]], residuearray[i].atomposition[tetrahedraindices[aminoacidindex][j][0]]);
                            if(dot_product(second, normal)>0){
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s\t%.8f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][4]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][0]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][1]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][2]], amino_acid_list[aminoacidindex].atomnames[tetrahedraindices[aminoacidindex][j][3]], dot_product(second, normal));
                                    fclose(outp);
                                }
                                for(k=0;k<5;k++){
                                    residuearray[i].atomfound[tetrahedraindices[aminoacidindex][j][k]]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void setup_rings_structure(char *rings_file, amino_acid_struct *amino_acid_list, int **pNrings, int ***pringsize, int ***pnumber_dependent_atoms, int ****patomindex){
    FILE *inp, *outp;
    int i, j, k, read, stringsize;
    declare_array_nozero(char, linechar, maxlinelength);
    allocate_array(int, (*pNrings), Naminoacids);
    allocate_matrix(int, (*pringsize), Naminoacids, max_rings_per_residue);
    allocate_matrix(int, (*pnumber_dependent_atoms), Naminoacids, max_rings_per_residue);
    allocate_3d_tensor(int, (*patomindex), Naminoacids, max_rings_per_residue, maxatomsperresidue);
    declare_array(char, readres, 3);
    declare_array(char, field, 1);
    declare_array(char, atomname, 4);
    inp=fopen(rings_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        stringsize=strlen(linechar);
        if(stringsize%4!=1){
            printf("stringsize %i in check_rings!\n", stringsize);
            exit(1);
        }
        extract_piece_of_string(linechar, readres, 1, 3);
        for(i=0;i<Naminoacids;i++){
            if(strcmp(readres, amino_acid_list[i].resname)==0) break;
        }
        if(i==Naminoacids){
            printf("residue <%s> in %s not found among amino acids!\n", readres, rings_file);
            exit(1);
        }
        extract_piece_of_string(linechar, field, 4, 4);
        (*pringsize)[i][(*pNrings)[i]]=atoi(field);
        (*pnumber_dependent_atoms)[i][(*pNrings)[i]]=(stringsize-5)/4;
        for(j=0;j<(*pnumber_dependent_atoms)[i][(*pNrings)[i]];j++){
            extract_piece_of_string(linechar, atomname, 4*(j+1)+1, 4*(j+1)+4);
            for(k=0;k<amino_acid_list[i].Natoms;k++){
                if(strcmp(atomname, amino_acid_list[i].atomnames[k])==0) break;
            }
            if(k==amino_acid_list[i].Natoms){
                printf("atom <%s> in %s not found among atom names for residue <%s>!\n", atomname, rings_file, readres);
                exit(1);
            }
            (*patomindex)[i][(*pNrings)[i]][j]=k;
        }
        (*pNrings)[i]++;
        read=mygetline(linechar, maxlinelength, inp);
    }
	fclose(inp);
    free(readres);
    free(field);
    free(atomname);
}

void check_rings(residuedata *residuearray, int Nresidues, int *Nrings, int **ringsize, int **number_dependent_atoms, int ***atomindex, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite){
    FILE *inp, *outp;
    int i, j, k, aminoacidindex, foundall, index, lastindex, count, lastfound;
    double totalangle, difference;
    double_triple center, normal, trianglenormal, ex, ey;
    declare_array_nozero(double_triple, director, 6);
    declare_array(int, atomfound, 6);
    declare_array(double, angle, 6);
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>=0){
            if(Nrings[aminoacidindex]>0){
                for(j=0;j<Nrings[aminoacidindex];j++){
                    center.x=center.y=center.z=0;
                    count=0;
                    for(k=0;k<ringsize[aminoacidindex][j];k++){
                        index=atomindex[aminoacidindex][j][k];
                        if(residuearray[i].atomfound[index]==1){
                            center=add_double_triple(center, residuearray[i].atomposition[index]);
                            count++;
                        }
                    }
                    if(count>0){
                        center=scalar_multiply_double_triple(center, 1./count);
                        normal.x=normal.y=normal.z=0;
                        count=0;
                        for(k=0;k<ringsize[aminoacidindex][j];k++){
                            index=atomindex[aminoacidindex][j][k];
                            atomfound[k]=0;
                            if(residuearray[i].atomfound[index]==1){
                                atomfound[k]=1;
                                director[k]=subtract_double_triple(residuearray[i].atomposition[index], center);
                                normalize(&(director[k]));
                                lastfound=k;
                            }
                        }
                        for(k=0;k<ringsize[aminoacidindex][j];k++){
                            if((atomfound[k]==1)&&(atomfound[mod(k-1, ringsize[aminoacidindex][j])]==1)){
                                trianglenormal=cross_product(director[mod(k-1, ringsize[aminoacidindex][j])], director[k]);
                                normalize(&(trianglenormal));
                                normal=add_double_triple(normal, trianglenormal);
                                count++;
                            }
                        }
                        if(count>0){
                            normal=scalar_multiply_double_triple(normal, 1./count);
                            ex=director[lastfound];
                            normalize(&ex);
                            ex=subtract_double_triple(ex, scalar_multiply_double_triple(normal, dot_product(ex, normal)));
                            normalize(&ex);
                            ey=cross_product(normal, ex);
                            normalize(&ey);
                            totalangle=0;
                            for(k=0;k<ringsize[aminoacidindex][j];k++){
                                if(atomfound[k]==1){
                                    angle[k]=atan2(dot_product(director[k], ey), dot_product(director[k], ex));
                                    if((k>0)&&(atomfound[k-1]==1)){
                                        difference=angle[k]-angle[k-1];
                                        if(difference>M_PI) difference-=(2*M_PI);
                                        else if(difference<-M_PI) difference+=(2*M_PI);
                                        totalangle+=difference;
                                    }
                                }
                            }
                            if((atomfound[0]==1)&&(atomfound[ringsize[aminoacidindex][j]-1]==1)){
                                difference=angle[0]-angle[ringsize[aminoacidindex][j]-1];
                                if(difference>M_PI) difference-=(2*M_PI);
                                else if(difference<-M_PI) difference+=(2*M_PI);
                                totalangle+=difference;
                            }
                            if((totalangle>3*M_PI)||(totalangle<-3*M_PI)){
                                if(dowrite==1){
                                    outp=fopen(output_file, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%s", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname);
                                    for(k=0;k<number_dependent_atoms[aminoacidindex][j];k++){
                                        index=atomindex[aminoacidindex][j][k];
                                        if(residuearray[i].atomfound[index]==1){
                                            fprintf(outp, "\t%s", amino_acid_list[aminoacidindex].atomnames[index]);
                                        }
                                    }
                                    fprintf(outp, "\n");
                                    fclose(outp);
                                }
                                for(k=0;k<number_dependent_atoms[aminoacidindex][j];k++){
                                    index=atomindex[aminoacidindex][j][k];
                                    residuearray[i].atomfound[index]=0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    free(director);
    free(atomfound);
    free(angle);
}

void print_indistinguishable_xyz_files(amino_acid_struct *amino_acid_list, residuedata *residuearray, int Nresidues, char *base, int *pfirst){
    int i, j, k, l, foundall;
    double_triple position, position_frame, center, first, second, firstaxis, secondaxis, thirdaxis, seconddif;
    declare_array_nozero(char, outputfile, maxstringlength);
    declare_array_nozero(char, label, 1);
    FILE *outp;
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].N_indistinguishable_sets;j++){
            if(amino_acid_list[i].indistinguishable_set_size[j]==2){
                sprintf(outputfile, "%s.indistinguishable_%s_%s_%s.xyz", base, amino_acid_list[i].resname, amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][3]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][4]]);
            }
            if(amino_acid_list[i].indistinguishable_set_size[j]==3){
                sprintf(outputfile, "%s.indistinguishable_%s_%s_%s_%s.xyz", base, amino_acid_list[i].resname, amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][3]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][4]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][5]]);
            }
            if(*pfirst==1){
                delete_file_if_exists(outputfile);
                outp=fopen(outputfile, "w");
            }
            else{
                if(outp=fopen(outputfile, "r")){
                    fclose(outp);
                    outp=fopen(outputfile, "a");
                }
                else outp=fopen(outputfile, "w");
            }
            for(k=0;k<Nresidues;k++){
                if(residuearray[k].residnumber==i){
                    foundall=1;
                    for(l=0;l<amino_acid_list[i].indistinguishable_set_size[j]+3;l++){
                        if(residuearray[k].atomfound[amino_acid_list[i].indistinguishable_atoms[j][l]]==0) foundall=0;
                    }
                    if(foundall==1){
                        center=residuearray[k].atomposition[amino_acid_list[i].indistinguishable_atoms[j][1]];
                        first=residuearray[k].atomposition[amino_acid_list[i].indistinguishable_atoms[j][2]];
                        second=residuearray[k].atomposition[amino_acid_list[i].indistinguishable_atoms[j][0]];
                        firstaxis=subtract_double_triple(first, center);
                        normalize(&firstaxis);
                        seconddif=subtract_double_triple(second, center);
                        secondaxis=subtract_double_triple(seconddif, scalar_multiply_double_triple(firstaxis, dot_product(firstaxis, seconddif)));
                        normalize(&secondaxis);
                        thirdaxis=cross_product(firstaxis, secondaxis);
                        normalize(&thirdaxis);
                        for(l=0;l<amino_acid_list[i].indistinguishable_set_size[j]+3;l++){
                            position=subtract_double_triple(residuearray[k].atomposition[amino_acid_list[i].indistinguishable_atoms[j][l]], center);
                            position_frame.x=dot_product(position, firstaxis);
                            position_frame.y=dot_product(position, secondaxis);
                            position_frame.z=dot_product(position, thirdaxis);
                            if(l==0) strcpy(label, "C");
                            if(l==1) strcpy(label, "N");
                            if(l==2) strcpy(label, "O");
                            if(l==3) strcpy(label, "P");
                            if(l==4) strcpy(label, "S");
                            if(l==5) strcpy(label, "I");
                            fprintf(outp, "%s\t%.16f\t%.16f\t%.16f\n", label, position_frame.x, position_frame.y, position_frame.z);
                        }
                    }
                }
            }
            fclose(outp);
        }
    }
    (*pfirst)=0;
    free(outputfile);
    free(label);
}

void resolve_indistinguishable_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int dowrite, char *output_file, char *pdb_file, int record_all_data){
    int i, j, k, l, aminoacidindex, foundall, permutation_counter, foundfirst=-1, target;
    double secondprojection, thirdprojection;
    double_triple origin, firstatom, secondatom, seconddifference, firstaxis, secondaxis, thirdaxis, indistinguishable_position, indistinguishable_position_frame;
    declare_array(double, angle, 3);
    declare_array_nozero(double_triple, new_position, 3);
    declare_array(int, new_atomfound, 3);
    declare_array(int, permutation, 3);
    declare_array(int, taken, 3);
    char **new_serialchar, **new_icodechar, **new_occupancychar, **new_tempfactorchar, **new_elementchar, **new_chargechar;
    if(record_all_data==1){
        allocate_matrix(char, new_serialchar, 3, 5);
        allocate_matrix(char, new_icodechar, 3, 1);
        allocate_matrix(char, new_occupancychar, 3, 6);
        allocate_matrix(char, new_tempfactorchar, 3, 6);
        allocate_matrix(char, new_elementchar, 3, 2);
        allocate_matrix(char, new_chargechar, 3, 2);
    }
    FILE *outp;
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>=0){
            for(j=0;j<amino_acid_list[aminoacidindex].N_indistinguishable_sets;j++){
                foundfirst=-1;
                for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                    if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][k+3]]==1){
                        foundfirst=k;
                        break;
                    }
                }
                if(foundfirst>=0){
                    foundall=1;
                    for(k=0;k<3;k++){
                        if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][k]]==0) foundall=0;
                    }
                    if(foundall==0){
                        if(dowrite==1){
                            outp=fopen(output_file, "a");
                            fprintf(outp, "%s\t%s\t%i\t%s", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname);
                            for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                                if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]==1){
                                    fprintf(outp, "\t%s", amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                }
                            }
                            fprintf(outp, "\n");
                            fclose(outp);
                        }
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                            residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]=0;
                        }
                    }
                    else{
                        origin=residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][1]];
                        firstatom=residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][2]];
                        secondatom=residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][0]];
                        firstaxis=subtract_double_triple(firstatom, origin);
                        normalize(&firstaxis);
                        seconddifference=subtract_double_triple(secondatom, origin);
                        secondaxis=subtract_double_triple(seconddifference, scalar_multiply_double_triple(firstaxis, dot_product(firstaxis, seconddifference)));
                        normalize(&secondaxis);
                        thirdaxis=cross_product(firstaxis, secondaxis);
                        normalize(&thirdaxis);
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++) taken[k]=0;
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                            if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][k+3]]==1){
                                indistinguishable_position=residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][k+3]];
                                indistinguishable_position_frame=subtract_double_triple(indistinguishable_position, origin);
                                secondprojection=dot_product(indistinguishable_position_frame, secondaxis);
                                thirdprojection=dot_product(indistinguishable_position_frame, thirdaxis);
                                angle[k]=atan2(thirdprojection, secondprojection);
                                target=(int) floor((angle[k]-amino_acid_list[aminoacidindex].first_dividing_angle[j])*amino_acid_list[aminoacidindex].indistinguishable_set_size[j]/(2*M_PI));
                                target=mod(target, amino_acid_list[aminoacidindex].indistinguishable_set_size[j]);
                                if(taken[target]==0){
                                    permutation[k]=target;
                                    taken[target]=1;
                                }
                                else{
                                    for(l=0;l<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];l++){
                                        if(taken[l]==0){
                                            permutation[k]=l;
                                            taken[l]=1;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                            if(residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][k+3]]==1){
                                new_position[permutation[k]]=residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]];
                                if(record_all_data==1){
                                    strcpy(new_serialchar[permutation[k]], residuearray[i].serialchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                    strcpy(new_icodechar[permutation[k]], residuearray[i].icodechar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                    strcpy(new_occupancychar[permutation[k]], residuearray[i].occupancychar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                    strcpy(new_tempfactorchar[permutation[k]], residuearray[i].tempfactorchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                    strcpy(new_elementchar[permutation[k]], residuearray[i].elementchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                    strcpy(new_chargechar[permutation[k]], residuearray[i].chargechar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]);
                                }
                                new_atomfound[permutation[k]]=residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]];
                            }
                        }
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                            residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]=0;
                        }
                        for(k=0;k<amino_acid_list[aminoacidindex].indistinguishable_set_size[j];k++){
                            if(taken[k]==1){
                                residuearray[i].atomposition[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]=new_position[k];
                                if(record_all_data==1){
                                    strcpy(residuearray[i].serialchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_serialchar[k]);
                                    strcpy(residuearray[i].icodechar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_icodechar[k]);
                                    strcpy(residuearray[i].occupancychar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_occupancychar[k]);
                                    strcpy(residuearray[i].tempfactorchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_tempfactorchar[k]);
                                    strcpy(residuearray[i].elementchar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_elementchar[k]);
                                    strcpy(residuearray[i].chargechar[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]], new_chargechar[k]);
                                }
                                residuearray[i].atomfound[amino_acid_list[aminoacidindex].indistinguishable_atoms[j][3+k]]=new_atomfound[k];
                            }
                        }
                    }
                }
            }
        }
    }
    free(angle);
    free(new_position);
    if(record_all_data==1){
        free_matrix(new_serialchar, 3);
        free_matrix(new_icodechar, 3);
        free_matrix(new_occupancychar, 3);
        free_matrix(new_tempfactorchar, 3);
        free_matrix(new_elementchar, 3);
        free_matrix(new_chargechar, 3);
    }
    free(new_atomfound);
    free(permutation);
    free(taken);
}

void resolve_aromatic_hydrogens(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int dowrite, char *output_file, char *pdb_file, int record_all_data){
    int i, j, k, heavyindex, Hindex, conflict, different, aminoacidindex, missinghydrogens, totalhydrogens;
    double close, sep;
    FILE *outp;
    declare_array(int, closest, max_aromatic_hydrogens);
    declare_array(int, assigned, max_aromatic_hydrogens);
    declare_array(double, distance, max_aromatic_hydrogens);
    declare_array_nozero(double_triple, new_position, max_aromatic_hydrogens);
    char **new_serialchar, **new_icodechar, **new_occupancychar, **new_tempfactorchar, **new_elementchar, **new_chargechar;
    if(record_all_data==1){
        allocate_matrix(char, new_serialchar, max_aromatic_hydrogens, 5);
        allocate_matrix(char, new_icodechar, max_aromatic_hydrogens, 1);
        allocate_matrix(char, new_occupancychar, max_aromatic_hydrogens, 6);
        allocate_matrix(char, new_tempfactorchar, max_aromatic_hydrogens, 6);
        allocate_matrix(char, new_elementchar, max_aromatic_hydrogens, 2);
        allocate_matrix(char, new_chargechar, max_aromatic_hydrogens, 2);
    }
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>=0){
            if(amino_acid_list[aminoacidindex].N_aromatic_hydrogens>0){
                
                conflict=0;
                different=0;
                missinghydrogens=0;
                totalhydrogens=0;
                for(j=0;j<amino_acid_list[aminoacidindex].N_aromatic_hydrogens;j++){
                    closest[j]=-1;
                    assigned[j]=-1;
                    heavyindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2+1];
                    Hindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2];
                    if(residuearray[i].atomfound[heavyindex]==1) missinghydrogens++;
                    if(residuearray[i].atomfound[Hindex]==1) totalhydrogens++;
                }
                if(totalhydrogens>0){
                    missinghydrogens-=totalhydrogens;
                    if(missinghydrogens<0) missinghydrogens=0;
                    for(j=0;j<amino_acid_list[aminoacidindex].N_aromatic_hydrogens;j++){
                        heavyindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2+1];
                        if(residuearray[i].atomfound[heavyindex]==1){
                            for(k=0;k<amino_acid_list[aminoacidindex].N_aromatic_hydrogens;k++){
                                Hindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[k*2];
                                if(residuearray[i].atomfound[Hindex]==1){
                                    if(closest[j]==-1){
                                        closest[j]=k;
                                        close=norm(subtract_double_triple(residuearray[i].atomposition[Hindex], residuearray[i].atomposition[heavyindex]));
                                    }
                                    else{
                                        sep=norm(subtract_double_triple(residuearray[i].atomposition[Hindex], residuearray[i].atomposition[heavyindex]));
                                        if(sep<close){
                                            closest[j]=k;
                                            close=sep;
                                        }
                                    }
                                }
                            }
                            if(closest[j]>-1){
                                if(assigned[closest[j]]>-1){
                                    
                                    //  Already assigned hydrogen closest[j] to another heavy atom (assigned[closest[j]])
                                    
                                    if(missinghydrogens>0){
                                        
                                        //  Check which heavy atom the hydrogen is closest to
                                        
                                        if(close<distance[closest[j]]){
                                            //  Closer to current heavy atom
                                            
                                            closest[assigned[closest[j]]]=-1;
                                            assigned[closest[j]]=j;
                                            distance[closest[j]]=close;
                                            if(closest[j]!=j) different=1;
                                            
                                            //  Previous heavy atom will no longer be associated with a hydrogen; that hydrogen slot can no longer be filled
                                            
                                            missinghydrogens--;
                                        }
                                        else{
                                            
                                            //  Current heavy atom will not be associated with a hydrogen; that hydrogen slot can no longer be filled
                                            
                                            closest[j]=-1;
                                            missinghydrogens--;
                                        }
                                    }
                                    else{
                                        outp=fopen(output_file, "a");
                                        fprintf(outp, "%s\t%s\t%i\t%s\t%s\t%s\t%s\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, amino_acid_list[aminoacidindex].resname, amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].aromatic_hydrogen_list[closest[j]*2]], amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2+1]], amino_acid_list[aminoacidindex].atomnames[amino_acid_list[aminoacidindex].aromatic_hydrogen_list[assigned[closest[j]]*2+1]]);
                                        fclose(outp);
                                        conflict=1;
                                    }
                                }
                                else{
                                    assigned[closest[j]]=j;
                                    distance[closest[j]]=close;
                                    if(closest[j]!=j) different=1;
                                }
                            }
                        }
                    }
                    if((conflict==0)&&(different==1)){
                        for(j=0;j<amino_acid_list[aminoacidindex].N_aromatic_hydrogens;j++){
                            
                            if(closest[j]>-1){
                                
                                heavyindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2+1];
                                if(residuearray[i].atomfound[heavyindex]==1){
                                    Hindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[closest[j]*2];
                                    new_position[j]=residuearray[i].atomposition[Hindex];
                                    if(record_all_data==1){
                                        strcpy(new_serialchar[j], residuearray[i].serialchar[Hindex]);
                                        strcpy(new_icodechar[j], residuearray[i].icodechar[Hindex]);
                                        strcpy(new_occupancychar[j], residuearray[i].occupancychar[Hindex]);
                                        strcpy(new_tempfactorchar[j], residuearray[i].tempfactorchar[Hindex]);
                                        strcpy(new_elementchar[j], residuearray[i].elementchar[Hindex]);
                                        strcpy(new_chargechar[j], residuearray[i].chargechar[Hindex]);
                                    }
                                    //  Wipe out hydrogens currently associated with the heavy atoms
                                    //  If there is a hydrogen associated with a missing heavy atom, it will remain there

                                    residuearray[i].atomfound[amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2]]=0;
                                }
                                
                            }
                        }
                        for(j=0;j<amino_acid_list[aminoacidindex].N_aromatic_hydrogens;j++){
                            if(closest[j]>-1){
                                heavyindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2+1];
                                if(residuearray[i].atomfound[heavyindex]==1){
                                    Hindex=amino_acid_list[aminoacidindex].aromatic_hydrogen_list[j*2];
                                    residuearray[i].atomposition[Hindex]=new_position[j];
                                    if(record_all_data==1){
                                        strcpy((residuearray[i].serialchar[Hindex]), new_serialchar[j]);
                                        strcpy((residuearray[i].icodechar[Hindex]), new_icodechar[j]);
                                        strcpy((residuearray[i].occupancychar[Hindex]), new_occupancychar[j]);
                                        strcpy((residuearray[i].tempfactorchar[Hindex]), new_tempfactorchar[j]);
                                        strcpy((residuearray[i].elementchar[Hindex]), new_elementchar[j]);
                                        strcpy((residuearray[i].chargechar[Hindex]), new_chargechar[j]);
                                    }
                                    residuearray[i].atomfound[Hindex]=1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    free(closest);
    free(assigned);
    free(new_position);
    free(distance);
    if(record_all_data==1){
        free_matrix(new_serialchar, max_aromatic_hydrogens);
        free_matrix(new_icodechar, max_aromatic_hydrogens);
        free_matrix(new_occupancychar, max_aromatic_hydrogens);
        free_matrix(new_tempfactorchar, max_aromatic_hydrogens);
        free_matrix(new_elementchar, max_aromatic_hydrogens);
        free_matrix(new_chargechar, max_aromatic_hydrogens);
    }
}

void map_to_cgmodel(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue **pcgresiduearray, int docountsites){
	
    //  Map all-atom configuration (residuearray) to coarse-grained configuration (cgresiduearray)
	
    int rescounter, i, j, aminoacidindex, atomindex, secondatomindex;
    double mynorm;
	double_triple vector;
    allocate_array(cgresidue, (*pcgresiduearray), Nresidues);
	for(i=0;i<Nresidues;i++) (*pcgresiduearray)[i].Nsites=0;
    for(rescounter=0;rescounter<Nresidues;rescounter++){
        aminoacidindex=residuearray[rescounter].residnumber;
		if(aminoacidindex>=0){
			(*pcgresiduearray)[rescounter].Nsites=amino_acid_list[aminoacidindex].Nsites;
            if((*pcgresiduearray)[rescounter].Nsites>0){
                if((rescounter==Nresidues-1)||(residuearray[rescounter].chainid!=residuearray[rescounter+1].chainid)){
                    ((*pcgresiduearray)[rescounter]).Cterminussite=(*pcgresiduearray)[rescounter].Nsites;
                    allocate_array(cgcoord, (((*pcgresiduearray)[rescounter]).coord), ((*pcgresiduearray)[rescounter].Nsites+1));
                    allocate_array(int, (((*pcgresiduearray)[rescounter]).sitedefined), ((*pcgresiduearray)[rescounter].Nsites+1));
                }
                else{
                    ((*pcgresiduearray)[rescounter]).Cterminussite=-1;
                    allocate_array(cgcoord, (((*pcgresiduearray)[rescounter]).coord), ((*pcgresiduearray)[rescounter].Nsites));
                    allocate_array(int, (((*pcgresiduearray)[rescounter]).sitedefined), ((*pcgresiduearray)[rescounter].Nsites));
                }
				for(j=0;j<(*pcgresiduearray)[rescounter].Nsites;j++){
					
					(*pcgresiduearray)[rescounter].sitedefined[j]=0;
					atomindex=amino_acid_list[aminoacidindex].siteatomindex[j][0];         //  Central atom
					if(atomindex==-1){
						printf("Trying to map to cg model but site atom index not defined!\n");
						exit(1);
					}
					if(residuearray[rescounter].atomfound[atomindex]==1){
						(*pcgresiduearray)[rescounter].coord[j].pos=residuearray[rescounter].atomposition[atomindex];
						if(amino_acid_list[aminoacidindex].sitecode[j]==0){
							secondatomindex=amino_acid_list[aminoacidindex].siteatomindex[j][1];
							if(residuearray[rescounter].atomfound[secondatomindex]==1){
								(*pcgresiduearray)[rescounter].coord[j].ex=subtract_double_triple(residuearray[rescounter].atomposition[secondatomindex], residuearray[rescounter].atomposition[atomindex]);
                                mynorm=norm((*pcgresiduearray)[rescounter].coord[j].ex);
                                if(mynorm>0){
                                    (((*pcgresiduearray)[rescounter].coord[j].ex).x)/=mynorm;
                                    (((*pcgresiduearray)[rescounter].coord[j].ex).y)/=mynorm;
                                    (((*pcgresiduearray)[rescounter].coord[j].ex).z)/=mynorm;
                                    secondatomindex=amino_acid_list[aminoacidindex].siteatomindex[j][2];
                                    if(residuearray[rescounter].atomfound[secondatomindex]==1){
                                        
                                        //  Sufficient atoms found to define the site
                                        
                                        (*pcgresiduearray)[rescounter].sitedefined[j]=1;
                                        if(docountsites==1) amino_acid_list[aminoacidindex].sitecount[j]++;
                                        
                                        vector=subtract_double_triple(residuearray[rescounter].atomposition[secondatomindex], residuearray[rescounter].atomposition[atomindex]);		//	not yet orthogonal to ex
                                        vector=subtract_double_triple(vector, scalar_multiply_double_triple((*pcgresiduearray)[rescounter].coord[j].ex, dot_product(vector, (*pcgresiduearray)[rescounter].coord[j].ex)));	//	project off part parallel to ex
                                        mynorm=norm(vector);
                                        if(mynorm>0){
                                            ((*pcgresiduearray)[rescounter].coord[j].ey).x=vector.x/mynorm;
                                            ((*pcgresiduearray)[rescounter].coord[j].ey).y=vector.y/mynorm;
                                            ((*pcgresiduearray)[rescounter].coord[j].ey).z=vector.z/mynorm;
                                            (*pcgresiduearray)[rescounter].coord[j].ez=cross_product((*pcgresiduearray)[rescounter].coord[j].ex, (*pcgresiduearray)[rescounter].coord[j].ey);
                                        }
                                    }
                                }
                            }
                        }
						else{
							printf("Mapping to cg model for site code %i not defined!\n", amino_acid_list[aminoacidindex].sitecode[j]);
							exit(1);
						}
					}
				}
                if((*pcgresiduearray)[rescounter].Cterminussite>=0){
                    atomindex=amino_acid_list[aminoacidindex].Cterminussiteatomindex[0];
                    if(residuearray[rescounter].atomfound[atomindex]==1){
                        (*pcgresiduearray)[rescounter].coord[j].pos=residuearray[rescounter].atomposition[atomindex];
                        secondatomindex=amino_acid_list[aminoacidindex].Cterminussiteatomindex[1];
                        if(residuearray[rescounter].atomfound[secondatomindex]==1){
                            (*pcgresiduearray)[rescounter].coord[j].ex=subtract_double_triple(residuearray[rescounter].atomposition[secondatomindex], residuearray[rescounter].atomposition[atomindex]);
                            normalize(&((*pcgresiduearray)[rescounter].coord[j].ex));
                            secondatomindex=amino_acid_list[aminoacidindex].Cterminussiteatomindex[2];
                            if(residuearray[rescounter].atomfound[secondatomindex]==1){
                                (*pcgresiduearray)[rescounter].sitedefined[j]=1;
                                if(docountsites==1) amino_acid_list[aminoacidindex].sitecount[j]++;
                                vector=subtract_double_triple(residuearray[rescounter].atomposition[secondatomindex], residuearray[rescounter].atomposition[atomindex]);		//	not yet orthogonal to ex
                                vector=subtract_double_triple(vector, scalar_multiply_double_triple((*pcgresiduearray)[rescounter].coord[j].ex, dot_product(vector, (*pcgresiduearray)[rescounter].coord[j].ex)));	//	project off part parallel to ex
                                (*pcgresiduearray)[rescounter].coord[j].ey=normed(vector);
                                (*pcgresiduearray)[rescounter].coord[j].ez=cross_product((*pcgresiduearray)[rescounter].coord[j].ex, (*pcgresiduearray)[rescounter].coord[j].ey);
                            }
                        }
                    }
                    ((*pcgresiduearray)[rescounter].Nsites)++;
                }
			}
		}
    }
}

void record_moments_relative_to_oriented_cg_sites(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray){
	
    //  Average first and second moments of atomic coordinates relative to oriented coarse-grained configuration
	
    int i, j, k, aminoacidindex, site, centerindex, Cterminus, index;
    double_triple pos, centerpos, relativepos, orientedpos;
	for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        for(j=0;j<residuearray[i].Natoms;j++){
            if(residuearray[i].atomfound[j]==1){
                if((cgresiduearray[i].Cterminussite>-1)&&(amino_acid_list[aminoacidindex].Cterminusindex[j]>-1)){
                    site=cgresiduearray[i].Cterminussite;
                    Cterminus=1;
                }
                else{
                    site=amino_acid_list[aminoacidindex].sitemap[j];
                    Cterminus=0;
                }
                if(site>=0){
                    if(site>=cgresiduearray[i].Nsites){
                        printf("Mapping atom %s (res %s) to site %i greater than total number of sites %i!\n", amino_acid_list[aminoacidindex].atomnames[j], amino_acid_list[aminoacidindex].resname, site, cgresiduearray[i].Nsites);
                        exit(1);
                    }
                    if(cgresiduearray[i].sitedefined[site]==1){
                        pos=residuearray[i].atomposition[j];
                        if(Cterminus==1) centerindex=amino_acid_list[aminoacidindex].Cterminussiteatomindex[0];
                        else centerindex=amino_acid_list[aminoacidindex].siteatomindex[site][0];
                        centerpos=residuearray[i].atomposition[centerindex];
                        relativepos=subtract_double_triple(pos, centerpos);
                        orientedpos.x=dot_product(relativepos, cgresiduearray[i].coord[site].ex);
                        orientedpos.y=dot_product(relativepos, cgresiduearray[i].coord[site].ey);
                        orientedpos.z=dot_product(relativepos, cgresiduearray[i].coord[site].ez);
                        if(Cterminus==1) index=amino_acid_list[aminoacidindex].Natoms+amino_acid_list[aminoacidindex].Cterminusindex[j];
                        else index=j;
                        amino_acid_list[aminoacidindex].avpos[index].x+=orientedpos.x;
                        amino_acid_list[aminoacidindex].avpos[index].y+=orientedpos.y;
                        amino_acid_list[aminoacidindex].avpos[index].z+=orientedpos.z;
                        amino_acid_list[aminoacidindex].varpos[index].x+=(orientedpos.x*orientedpos.x);
                        amino_acid_list[aminoacidindex].varpos[index].y+=(orientedpos.y*orientedpos.y);
                        amino_acid_list[aminoacidindex].varpos[index].z+=(orientedpos.z*orientedpos.z);
                        amino_acid_list[aminoacidindex].count[index]++;
                    }
                }
            }
        }
    }
}

void free_cgresiduearray(cgresidue **pcgresiduearray, int Nresidues){
    int i;
    for(i=0;i<Nresidues;i++){
        if((*pcgresiduearray)[i].Nsites>0){
            free((*pcgresiduearray)[i].coord);
            free((*pcgresiduearray)[i].sitedefined);
        }
    }
    free((*pcgresiduearray));
}

void free_residuearray(residuedata **presiduearray, int Nresidues, int record_all_data){
    int i;
    for(i=0;i<Nresidues;i++){
        if((*presiduearray)[i].chainid>=0){
            free((*presiduearray)[i].atomfound);
            free((*presiduearray)[i].atomposition);
            free((*presiduearray)[i].chainidchar);
            if(record_all_data==1){
                free_matrix((*presiduearray)[i].serialchar, (*presiduearray)[i].Natoms);
                free_matrix((*presiduearray)[i].icodechar, (*presiduearray)[i].Natoms);
                free_matrix((*presiduearray)[i].occupancychar, (*presiduearray)[i].Natoms);
                free_matrix((*presiduearray)[i].tempfactorchar, (*presiduearray)[i].Natoms);
                free_matrix((*presiduearray)[i].elementchar, (*presiduearray)[i].Natoms);
                free_matrix((*presiduearray)[i].chargechar, (*presiduearray)[i].Natoms);
            }
            
        }
    }
    free((*presiduearray));
}

void add_header_to_xyz_files(char *base, amino_acid_struct *amino_acid_list){
    int i, j, k, read, counter;
    FILE *inp;
    declare_array_nozero(char, outputfile, maxstringlength);
    declare_array_nozero(char, linechar, maxlinelength);
    declare_matrix_nozero(char, filechar, 10000, maxlinelength);
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].N_indistinguishable_sets;j++){
            counter=0;
            if(amino_acid_list[i].indistinguishable_set_size[j]==2){
                sprintf(outputfile, "%s.indistinguishable_%s_%s_%s.xyz", base, amino_acid_list[i].resname, amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][3]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][4]]);
            }
            if(amino_acid_list[i].indistinguishable_set_size[j]==3){
                sprintf(outputfile, "%s.indistinguishable_%s_%s_%s_%s.xyz", base, amino_acid_list[i].resname, amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][3]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][4]], amino_acid_list[i].atomnames[amino_acid_list[i].indistinguishable_atoms[j][5]]);
            }
            inp=fopen(outputfile, "r");
            read=mygetline(linechar, maxlinelength, inp);
            while(read>0){
                strcpy(filechar[counter], linechar);
                counter++;
                read=mygetline(linechar, maxlinelength, inp);
                if(counter==10000){
                    printf("counter=%i!\n", counter);
                    exit(1);
                }
            }
            fclose(inp);
            inp=fopen(outputfile, "w");
            fprintf(inp, "%i\n\n", counter);
            for(k=0;k<counter;k++){
                fprintf(inp, "%s", filechar[k]);
            }
        }
    }
    free(outputfile);
    free(linechar);
    free_matrix(filechar, 10000);
}

void average_and_output_radial_positions(char *filename, amino_acid_struct *amino_acid_list){
    int i, j, k;
    FILE *outp;
    outp=fopen(filename, "w");
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].N_radial_sets;j++){
            fprintf(outp, "%s", amino_acid_list[i].resname);
            for(k=0;k<amino_acid_list[i].radial_set_size[j];k++){
                fprintf(outp, "\t%s", amino_acid_list[i].atomnames[amino_acid_list[i].radial_atoms[j][k]]);
            }
            if(amino_acid_list[i].radial_count[j]>0){
                fprintf(outp, "\t%.8f\t%.8f\t%i\n", amino_acid_list[i].radial_average[j].x/amino_acid_list[i].radial_count[j], amino_acid_list[i].radial_average[j].y/amino_acid_list[i].radial_count[j], amino_acid_list[i].radial_count[j]);
            }
            else{
                fprintf(outp, "\tN/A\tN/A\t%i\n", amino_acid_list[i].radial_count[j]);
            }
        }
    }
    fclose(outp);
}

void average_and_output_moments_relative_to_oriented_cg_sites(char *filename, amino_acid_struct *amino_acid_list, char **backbone_atom_types, int N_backbone_atom_types, char **Cterminusatomnames){
	
    //  Calculate and output averages and standard deviations from first and second moments (moments stored in amino_acid_list)
	
    int i, j, totalcountsum=0, k, found;
    double mag, totalvarsum=0;
    FILE *outp;
    outp=fopen(filename, "w");
    declare_array_nozero(double_triple, totalvarpos, N_backbone_atom_types);
    declare_array(int, totalcount, N_backbone_atom_types);
    declare_array_nozero(double_triple, Cterminusvarpos, N_backbone_atom_types);
    declare_array(int, Cterminuscount, N_backbone_atom_types);
    for(i=0;i<N_backbone_atom_types;i++){
        totalvarpos[i].x=totalvarpos[i].y=totalvarpos[i].z=0;
    }
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].Natoms;j++){
            if(amino_acid_list[i].count[j]>0){
                amino_acid_list[i].avpos[j]=scalar_multiply_double_triple(amino_acid_list[i].avpos[j], 1./amino_acid_list[i].count[j]);
                amino_acid_list[i].varpos[j].x=amino_acid_list[i].varpos[j].x/(1.*amino_acid_list[i].count[j])-pow(amino_acid_list[i].avpos[j].x, 2);
                amino_acid_list[i].varpos[j].y=amino_acid_list[i].varpos[j].y/(1.*amino_acid_list[i].count[j])-pow(amino_acid_list[i].avpos[j].y, 2);
                amino_acid_list[i].varpos[j].z=amino_acid_list[i].varpos[j].z/(1.*amino_acid_list[i].count[j])-pow(amino_acid_list[i].avpos[j].z, 2);
                mag=amino_acid_list[i].varpos[j].x+amino_acid_list[i].varpos[j].y+amino_acid_list[i].varpos[j].z;
                if(amino_acid_list[i].sitemap[j]==0){
                    
                    //  Calculate residue-average moments (average over amino acids) for backbone only
                    
                    found=0;
                    for(k=0;k<N_backbone_atom_types;k++){
                        if(strcmp(amino_acid_list[i].atomnames[j], backbone_atom_types[k])==0){
                            found=1;
                            break;
                        }
                    }
                    totalvarpos[k]=add_double_triple(totalvarpos[j], scalar_multiply_double_triple(amino_acid_list[i].varpos[j], amino_acid_list[i].count[j]));
                    totalcount[k]+=amino_acid_list[i].count[j];
                }
                mag=sqrt(mag);
                fprintf(outp, "%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%i\n", amino_acid_list[i].resname, amino_acid_list[i].atomnames[j], amino_acid_list[i].avpos[j].x, amino_acid_list[i].avpos[j].y, amino_acid_list[i].avpos[j].z, sqrt(amino_acid_list[i].varpos[j].x), sqrt(amino_acid_list[i].varpos[j].y), sqrt(amino_acid_list[i].varpos[j].z), mag, amino_acid_list[i].count[j]);
            }
        }
        for(j=0;j<amino_acid_list[i].NCterminusatoms;j++){
            k=amino_acid_list[i].Natoms+j;
            if(amino_acid_list[i].count[k]>0){
                amino_acid_list[i].avpos[k]=scalar_multiply_double_triple(amino_acid_list[i].avpos[k], 1./amino_acid_list[i].count[k]);
                amino_acid_list[i].varpos[k].x=amino_acid_list[i].varpos[k].x/(1.*amino_acid_list[i].count[k])-pow(amino_acid_list[i].avpos[k].x, 2);
                amino_acid_list[i].varpos[k].y=amino_acid_list[i].varpos[k].y/(1.*amino_acid_list[i].count[k])-pow(amino_acid_list[i].avpos[k].y, 2);
                amino_acid_list[i].varpos[k].z=amino_acid_list[i].varpos[k].z/(1.*amino_acid_list[i].count[k])-pow(amino_acid_list[i].avpos[k].z, 2);
                mag=amino_acid_list[i].varpos[k].x+amino_acid_list[i].varpos[k].y+amino_acid_list[i].varpos[k].z;
                Cterminusvarpos[j]=add_double_triple(Cterminusvarpos[j], scalar_multiply_double_triple(amino_acid_list[i].varpos[k], amino_acid_list[i].count[k]));
                Cterminuscount[j]+=amino_acid_list[i].count[k];
                mag=sqrt(mag);
                fprintf(outp, "%s\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%i\n", amino_acid_list[i].resname, Cterminusatomnames[j], amino_acid_list[i].avpos[k].x, amino_acid_list[i].avpos[k].y, amino_acid_list[i].avpos[k].z, sqrt(amino_acid_list[i].varpos[k].x), sqrt(amino_acid_list[i].varpos[k].y), sqrt(amino_acid_list[i].varpos[k].z), mag, amino_acid_list[i].count[k]);
            }
        }
    }
    fprintf(outp, "\n");
    for(i=0;i<N_backbone_atom_types;i++){
        if(totalcount[i]>0){
            totalvarpos[i]=scalar_multiply_double_triple(totalvarpos[i], 1./totalcount[i]);
            mag=totalvarpos[i].x+totalvarpos[i].y+totalvarpos[i].z;
            totalvarsum+=mag*totalcount[i];
            totalcountsum+=totalcount[i];
            totalvarpos[i].x=sqrt(totalvarpos[i].x);
            totalvarpos[i].y=sqrt(totalvarpos[i].y);
            totalvarpos[i].z=sqrt(totalvarpos[i].z);
            mag=sqrt(mag);
            fprintf(outp, "Varying\t%s\t\tN/A\t\tN/A\t\tN/A\t%.8f\t%.8f\t%.8f\t%.8f\t%i\n", backbone_atom_types[i], totalvarpos[i].x, totalvarpos[i].y, totalvarpos[i].z, mag, totalcount[i]);
        }
    }
    for(i=0;i<amino_acid_list[i].NCterminusatoms;i++){
        if(Cterminuscount[i]>0){
            Cterminusvarpos[i]=scalar_multiply_double_triple(Cterminusvarpos[i], 1./Cterminuscount[i]);
            mag=Cterminusvarpos[i].x+Cterminusvarpos[i].y+Cterminusvarpos[i].z;
            totalvarsum+=mag*Cterminuscount[i];
            totalcountsum+=Cterminuscount[i];
            Cterminusvarpos[i].x=sqrt(Cterminusvarpos[i].x);
            Cterminusvarpos[i].y=sqrt(Cterminusvarpos[i].y);
            Cterminusvarpos[i].z=sqrt(Cterminusvarpos[i].z);
            mag=sqrt(mag);
            fprintf(outp, "Varying\t%s\t\tN/A\t\tN/A\t\tN/A\t%.8f\t%.8f\t%.8f\t%.8f\t%i\n", Cterminusatomnames[i], Cterminusvarpos[i].x, Cterminusvarpos[i].y, Cterminusvarpos[i].z, mag, Cterminuscount[i]);
        }
    }
    fprintf(outp, "Varying\tOverall\t\tN/A\t\tN/A\t\tN/A\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\n", sqrt(totalvarsum/totalcountsum), totalcountsum);
    fprintf(outp, "\n");
    totalvarsum=0;
    totalcountsum=0;
    fclose(outp);
    free(totalvarpos);
    free(totalcount);
    free(Cterminusvarpos);
    free(Cterminuscount);
}

void calculate_covariance_separate(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array){
	
    //  Average covariance between carbonyl Oxygen positions and predicted Nitrogen position of next residue
	
    int i, j, k, selfaminoacidindex, nextaminoacidindex, Oindex, Nindex, Hindex, Cindex;
    double_triple Opos, Oposframe, Nposnextframe[3], Npos;
    double dOposframe[3], Nposframe[3];
    double_triple Hpos, Hposframe, Cpospreviousframe[3], Cpos;
    double dHposframe[3], Cposframe[3];
    for(i=0;i<Nresidues-1;i++){
        
        //  Stop at Nresidues-2, because right neighbor needs to be defined
        
        if(residuearray[i].chainid>=0){
            if(residuearray[i].chainid==residuearray[i+1].chainid){
                if((cgresiduearray[i].sitedefined[0]==1)&&(cgresiduearray[i+1].sitedefined[0]==1)){
                    selfaminoacidindex=residuearray[i].residnumber;
                    nextaminoacidindex=residuearray[i+1].residnumber;
                    Oindex=amino_acid_list[selfaminoacidindex].Oindex;
                    Nindex=amino_acid_list[nextaminoacidindex].Nindex;
                    if(Oindex==-1) my_exit("Oindex not assigned");
                    if(Nindex==-1) my_exit("Nindex not assigned");
                    Cindex=amino_acid_list[selfaminoacidindex].Cindex;
                    Hindex=amino_acid_list[nextaminoacidindex].Hindex;
                    if(Cindex==-1) my_exit("Cindex not assigned");
                    if(Hindex==-1) my_exit("Hindex not assigned");
                    if(residuearray[i].atomfound[Oindex]==1){
                        Opos=subtract_double_triple(residuearray[i].atomposition[Oindex], cgresiduearray[i].coord[0].pos);
                        Oposframe.x=dot_product(Opos, cgresiduearray[i].coord[0].ex);
                        Oposframe.y=dot_product(Opos, cgresiduearray[i].coord[0].ey);
                        Oposframe.z=dot_product(Opos, cgresiduearray[i].coord[0].ez);
                        dOposframe[0]=Oposframe.x-amino_acid_list[selfaminoacidindex].avpos[Oindex].x;
                        dOposframe[1]=Oposframe.y-amino_acid_list[selfaminoacidindex].avpos[Oindex].y;
                        dOposframe[2]=Oposframe.z-amino_acid_list[selfaminoacidindex].avpos[Oindex].z;
                        Nposnextframe[0]=scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ex, amino_acid_list[nextaminoacidindex].avpos[Nindex].x);
                        Nposnextframe[1]=scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ey, amino_acid_list[nextaminoacidindex].avpos[Nindex].y);
                        Nposnextframe[2]=scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ez, amino_acid_list[nextaminoacidindex].avpos[Nindex].z);
                        Npos=add_double_triple(cgresiduearray[i+1].coord[0].pos, add_double_triple(Nposnextframe[0], add_double_triple(Nposnextframe[1], Nposnextframe[2])));
                        Npos=subtract_double_triple(Npos, cgresiduearray[i].coord[0].pos);
                        Nposframe[0]=dot_product(Npos, cgresiduearray[i].coord[0].ex);
                        Nposframe[1]=dot_product(Npos, cgresiduearray[i].coord[0].ey);
                        Nposframe[2]=dot_product(Npos, cgresiduearray[i].coord[0].ez);
                        
                        for(j=0;j<3;j++){
                            (carbonyl_correction_array[selfaminoacidindex]).predicted_pos_av[j]+=Nposframe[j];
                            (carbonyl_correction_array[selfaminoacidindex]).predicted_pos_var[j]+=pow(Nposframe[j], 2);
                            (carbonyl_correction_array[selfaminoacidindex]).pos_difference_av[j]+=dOposframe[j];
                            for(k=0;k<3;k++){
                                (carbonyl_correction_array[selfaminoacidindex]).covariance_matrix_cross[j][k]+=Nposframe[k]*dOposframe[j];
                                (carbonyl_correction_array[selfaminoacidindex]).covariance_matrix_self[j][k]+=Nposframe[j]*Nposframe[k];
                            }
                        }
                        (carbonyl_correction_array[selfaminoacidindex]).count++;
                    }
                    if(residuearray[i+1].atomfound[Hindex]==1){
                        Hpos=subtract_double_triple(residuearray[i+1].atomposition[Hindex], cgresiduearray[i+1].coord[0].pos);
                        Hposframe.x=dot_product(Hpos, cgresiduearray[i+1].coord[0].ex);
                        Hposframe.y=dot_product(Hpos, cgresiduearray[i+1].coord[0].ey);
                        Hposframe.z=dot_product(Hpos, cgresiduearray[i+1].coord[0].ez);
                        dHposframe[0]=Hposframe.x-amino_acid_list[nextaminoacidindex].avpos[Hindex].x;
                        dHposframe[1]=Hposframe.y-amino_acid_list[nextaminoacidindex].avpos[Hindex].y;
                        dHposframe[2]=Hposframe.z-amino_acid_list[nextaminoacidindex].avpos[Hindex].z;
                        Cpospreviousframe[0]=scalar_multiply_double_triple(cgresiduearray[i].coord[0].ex, amino_acid_list[selfaminoacidindex].avpos[Cindex].x);
                        Cpospreviousframe[1]=scalar_multiply_double_triple(cgresiduearray[i].coord[0].ey, amino_acid_list[selfaminoacidindex].avpos[Cindex].y);
                        Cpospreviousframe[2]=scalar_multiply_double_triple(cgresiduearray[i].coord[0].ez, amino_acid_list[selfaminoacidindex].avpos[Cindex].z);
                        Cpos=add_double_triple(cgresiduearray[i].coord[0].pos, add_double_triple(Cpospreviousframe[0], add_double_triple(Cpospreviousframe[1], Cpospreviousframe[2])));
                        Cpos=subtract_double_triple(Cpos, cgresiduearray[i+1].coord[0].pos);
                        Cposframe[0]=dot_product(Cpos, cgresiduearray[i+1].coord[0].ex);
                        Cposframe[1]=dot_product(Cpos, cgresiduearray[i+1].coord[0].ey);
                        Cposframe[2]=dot_product(Cpos, cgresiduearray[i+1].coord[0].ez);
                        
                        for(j=0;j<3;j++){
                            (backboneH_correction_array)[selfaminoacidindex].predicted_pos_av[j]+=Cposframe[j];
                            (backboneH_correction_array)[selfaminoacidindex].predicted_pos_var[j]+=pow(Cposframe[j], 2);
                            (backboneH_correction_array)[selfaminoacidindex].pos_difference_av[j]+=dHposframe[j];
                            for(k=0;k<3;k++){
                                (backboneH_correction_array)[selfaminoacidindex].covariance_matrix_cross[j][k]+=Cposframe[k]*dHposframe[j];
                                (backboneH_correction_array)[selfaminoacidindex].covariance_matrix_self[j][k]+=Cposframe[j]*Cposframe[k];
                            }
                        }
                        (backboneH_correction_array)[selfaminoacidindex].count++;
                    }
                }
            }
        }
    }
}

void calculate_covariance_sidechain(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray, correction_structure *corrections){
    int i, j, k, l, found, aminoacidindex, my_corrected_atom_index, my_correcting_atom_index, corrected_atom_site, correcting_atom_site, my_correction_index;
    double_triple corrected_atom_pos, corrected_atom_pos_frame, correcting_atom_pos_ownframe, correcting_atom_pos_labframe;
    double pos_difference_av[3], correcting_atom_pos_corrected_frame[3];
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
        if(aminoacidindex>=0){
            for(l=0;l<amino_acid_list[aminoacidindex].N_corrections;l++){
                my_correction_index=amino_acid_list[aminoacidindex].correction_index[l];
                my_corrected_atom_index=amino_acid_list[aminoacidindex].corrected_atom_index[l];
                my_correcting_atom_index=amino_acid_list[aminoacidindex].correcting_atom_index[l];
                corrected_atom_site=amino_acid_list[aminoacidindex].sitemap[my_corrected_atom_index];
                correcting_atom_site=amino_acid_list[aminoacidindex].sitemap[my_correcting_atom_index];
                if((residuearray[i].atomfound[my_corrected_atom_index]==1)&&(residuearray[i].atomfound[my_correcting_atom_index]==1)){
                    if((cgresiduearray[i].sitedefined[corrected_atom_site]==1)&&(cgresiduearray[i].sitedefined[correcting_atom_site]==1)){
                        if(corrected_atom_site==correcting_atom_site){
                            printf("corrected and correcting atoms are in the same site (%i)!\n", corrected_atom_site);
                            printf("corrected atom %s and correcting atom %s in %s\n", amino_acid_list[aminoacidindex].atomnames[my_corrected_atom_index], amino_acid_list[aminoacidindex].atomnames[my_correcting_atom_index], amino_acid_list[aminoacidindex].resname);
                            exit(1);
                        }
                        corrected_atom_pos=subtract_double_triple(residuearray[i].atomposition[my_corrected_atom_index], cgresiduearray[i].coord[corrected_atom_site].pos);
                        
                        corrected_atom_pos_frame.x=dot_product(corrected_atom_pos, cgresiduearray[i].coord[corrected_atom_site].ex);
                        corrected_atom_pos_frame.y=dot_product(corrected_atom_pos, cgresiduearray[i].coord[corrected_atom_site].ey);
                        corrected_atom_pos_frame.z=dot_product(corrected_atom_pos, cgresiduearray[i].coord[corrected_atom_site].ez);
                        pos_difference_av[0]=corrected_atom_pos_frame.x-amino_acid_list[aminoacidindex].avpos[my_corrected_atom_index].x;
                        pos_difference_av[1]=corrected_atom_pos_frame.y-amino_acid_list[aminoacidindex].avpos[my_corrected_atom_index].y;
                        pos_difference_av[2]=corrected_atom_pos_frame.z-amino_acid_list[aminoacidindex].avpos[my_corrected_atom_index].z;
                        
                        correcting_atom_pos_ownframe=amino_acid_list[aminoacidindex].avpos[my_correcting_atom_index];
                        correcting_atom_pos_labframe=cgresiduearray[i].coord[correcting_atom_site].pos;
                        correcting_atom_pos_labframe=add_double_triple(correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ex, correcting_atom_pos_ownframe.x));
                        correcting_atom_pos_labframe=add_double_triple(correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ey, correcting_atom_pos_ownframe.y));
                        correcting_atom_pos_labframe=add_double_triple(correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ez, correcting_atom_pos_ownframe.z));
                        correcting_atom_pos_labframe=subtract_double_triple(correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].pos);
                        //recenter_double_triple(&correcting_atom_pos_labframe, pbc);
                        correcting_atom_pos_corrected_frame[0]=dot_product(correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ex);
                        correcting_atom_pos_corrected_frame[1]=dot_product(correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ey);
                        correcting_atom_pos_corrected_frame[2]=dot_product(correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ez);
                        
                        for(j=0;j<3;j++){
                            corrections[my_correction_index].predicted_pos_av[j]+=correcting_atom_pos_corrected_frame[j];
                            corrections[my_correction_index].predicted_pos_var[j]+=pow(correcting_atom_pos_corrected_frame[j], 2);
                            corrections[my_correction_index].pos_difference_av[j]+=pos_difference_av[j];
                            for(k=0;k<3;k++){
                                corrections[my_correction_index].covariance_matrix_cross[j][k]+=correcting_atom_pos_corrected_frame[k]*pos_difference_av[j];
                                corrections[my_correction_index].covariance_matrix_self[j][k]+=correcting_atom_pos_corrected_frame[j]*correcting_atom_pos_corrected_frame[k];
                            }
                        }
                        corrections[my_correction_index].count++;
                    }
                }
            }
        }
    }
}

void calculate_correction_from_covariance(correction_structure *pcorrection, char *corrections_file){
    
    //  Calculate linear correction from covariance; this accounts for cis-trans dihedral rotation
    
    FILE *outp;
    outp=fopen(corrections_file, "w");
    int i, j, k;
    declare_matrix(double, inverse_covariance_matrix, 3, 3);
    for(i=0;i<3;i++){
        (*pcorrection).predicted_pos_av[i]/=(1.*(*pcorrection).count);
        (*pcorrection).predicted_pos_var[i]/=(1.*(*pcorrection).count);
        (*pcorrection).pos_difference_av[i]/=(1.*(*pcorrection).count);
        (*pcorrection).intercept[i]=(*pcorrection).pos_difference_av[i];
    }
	fprintf(outp, "Correcting atom position in corrected atom frame:\n");
    fprintf(outp, "Average:\t\t%.16f\t%.16f\t%.16f\n", (*pcorrection).predicted_pos_av[0], (*pcorrection).predicted_pos_av[1], (*pcorrection).predicted_pos_av[2]);
    fprintf(outp, "Standard deviation:\t%.16f\t%.16f\t%.16f\n", sqrt((*pcorrection).predicted_pos_var[0]-pow((*pcorrection).predicted_pos_av[0], 2)), sqrt((*pcorrection).predicted_pos_var[1]-pow((*pcorrection).predicted_pos_av[1], 2)), sqrt((*pcorrection).predicted_pos_var[2]-pow((*pcorrection).predicted_pos_av[2], 2)));
	fprintf(outp, "\nChange in corrected atom position in corrected atom frame (after applying correction):\n");
    fprintf(outp, "Average:\t\t%.16f\t%.16f\t%.16f\n", (*pcorrection).pos_difference_av[0], (*pcorrection).pos_difference_av[1], (*pcorrection).pos_difference_av[2]);
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            (*pcorrection).covariance_matrix_cross[i][j]=(*pcorrection).covariance_matrix_cross[i][j]/(1.*(*pcorrection).count)-(*pcorrection).pos_difference_av[i]*(*pcorrection).predicted_pos_av[j];
            (*pcorrection).covariance_matrix_self[i][j]=(*pcorrection).covariance_matrix_self[i][j]/(1.*(*pcorrection).count)-(*pcorrection).predicted_pos_av[i]*(*pcorrection).predicted_pos_av[j];
        }
    }
	fprintf(outp, "\nCross correlation between corrected atom position and correcting atom position:\n");
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_cross[0][0], (*pcorrection).covariance_matrix_cross[0][1], (*pcorrection).covariance_matrix_cross[0][2]);
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_cross[1][0], (*pcorrection).covariance_matrix_cross[1][1], (*pcorrection).covariance_matrix_cross[1][2]);
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_cross[2][0], (*pcorrection).covariance_matrix_cross[2][1], (*pcorrection).covariance_matrix_cross[2][2]);
	fprintf(outp, "\nSelf cross correlation for correcting atom position:\n");
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_self[0][0], (*pcorrection).covariance_matrix_self[0][1], (*pcorrection).covariance_matrix_self[0][2]);
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_self[1][0], (*pcorrection).covariance_matrix_self[1][1], (*pcorrection).covariance_matrix_self[1][2]);
    fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).covariance_matrix_self[2][0], (*pcorrection).covariance_matrix_self[2][1], (*pcorrection).covariance_matrix_self[2][2]);
    
    invert_3x3_matrix((*pcorrection).covariance_matrix_self, inverse_covariance_matrix);
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            (*pcorrection).slope[i][j]=0;
            for(k=0;k<3;k++){
                (*pcorrection).slope[i][j]+=inverse_covariance_matrix[j][k]*(*pcorrection).covariance_matrix_cross[i][k];
            }
        }
    }
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            (*pcorrection).intercept[i]-=(*pcorrection).slope[i][j]*(*pcorrection).predicted_pos_av[j];
        }
    }    
	fprintf(outp, "\nCorrected atom position += B + M x (Correcting atom position)\n");
	fprintf(outp, "\nB:\n%.16f\t%.16f\t%.16f\n", (*pcorrection).intercept[0], (*pcorrection).intercept[1], (*pcorrection).intercept[2]);
	fprintf(outp, "\nM:\n%.16f\t%.16f\t%.16f\n", (*pcorrection).slope[0][0], (*pcorrection).slope[0][1], (*pcorrection).slope[0][2]);
	fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).slope[1][0], (*pcorrection).slope[1][1], (*pcorrection).slope[1][2]);
	fprintf(outp, "%.16f\t%.16f\t%.16f\n", (*pcorrection).slope[2][0], (*pcorrection).slope[2][1], (*pcorrection).slope[2][2]);
    free_matrix(inverse_covariance_matrix, 3);
    fclose(outp);
}

void calculate_rmsd_separate_corrections(amino_acid_struct *amino_acid_list, residuedata *residuearray, cgresidue *cgresiduearray, int Nresidues, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array, correction_structure *sidechain_corrections, char *pdb_file, char *large_disagreement_file_specific, char *entry_average_file){
    
    //  Calculate root mean square displacement between original all-atom configuration and all-atom configuration generated from coarse-grained representation
    
    int i, j, k, found, site, aminoacidindex, nextaminoacidindex, previousaminoacidindex, Oindex, Nindex, Hindex, Cindex, dormsd, my_corrected_atom_index, my_correcting_atom_index, corrected_atom_site, correcting_atom_site, myindex, entry_average_heavy_count_residuespecific=0, entry_average_H_count_residuespecific=0, Cterminus, index, correctionindexwithinaa;
    double distancesquared, msd_entry_average_heavy_residuespecific=0, fourthmoment_entry_average_heavy_residuespecific=0, msd_entry_average_H_residuespecific=0, fourthmoment_entry_average_H_residuespecific=0;
    double_triple position, position_frame, predicted_position, difference, predicted_correcting_atom_pos_ownframe, predicted_correcting_atom_pos_labframe, predicted_correcting_atom_pos_thisframe;
    
    declare_array(double, intercept_O, 3);
    declare_array(double, intercept_H, 3);
    declare_matrix(double, slope_O, 3, 3);
    declare_matrix(double, slope_H, 3, 3);
    
    FILE *outp;
    for(i=0;i<Nresidues;i++){
        aminoacidindex=residuearray[i].residnumber;
		if(aminoacidindex>=0){
            for(j=0;j<3;j++){
                intercept_O[j]=carbonyl_correction_array[aminoacidindex].intercept[j];
                intercept_H[j]=backboneH_correction_array[aminoacidindex].intercept[j];
                for(k=0;k<3;k++){
                    slope_O[j][k]=carbonyl_correction_array[aminoacidindex].slope[j][k];
                    slope_H[j][k]=backboneH_correction_array[aminoacidindex].slope[j][k];
                }
            }
			Oindex=amino_acid_list[aminoacidindex].Oindex;
			Hindex=amino_acid_list[aminoacidindex].Hindex;
			for(j=0;j<residuearray[i].Natoms;j++){
				if(residuearray[i].atomfound[j]==1){
                    if((cgresiduearray[i].Cterminussite>-1)&&(amino_acid_list[aminoacidindex].Cterminusindex[j]>-1)){
                        site=cgresiduearray[i].Cterminussite;
                        Cterminus=1;
                    }
                    else{
                        site=amino_acid_list[aminoacidindex].sitemap[j];
                        Cterminus=0;
                    }
					if(site>=0){
						if(cgresiduearray[i].sitedefined[site]==1){
							position=subtract_double_triple(residuearray[i].atomposition[j], cgresiduearray[i].coord[site].pos);
							position_frame.x=dot_product(position, cgresiduearray[i].coord[site].ex);
							position_frame.y=dot_product(position, cgresiduearray[i].coord[site].ey);
							position_frame.z=dot_product(position, cgresiduearray[i].coord[site].ez);
                            
                            if(Cterminus==1) index=amino_acid_list[aminoacidindex].Natoms+amino_acid_list[aminoacidindex].Cterminusindex[j];
                            else index=j;
                            predicted_position=amino_acid_list[aminoacidindex].avpos[index];
                            
                            if((Cterminus==0)&&(j==Oindex)){
                                if(cgresiduearray[i].Cterminussite==-1){
                                    dormsd=0;
                                    if(residuearray[i].atomfound[Oindex]==1){
                                        if((i<Nresidues-1)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                                            if(cgresiduearray[i+1].sitedefined[0]==1){
                                                dormsd=1;
                                                nextaminoacidindex=residuearray[i+1].residnumber;
                                                Nindex=amino_acid_list[nextaminoacidindex].Nindex;
                                                predicted_correcting_atom_pos_ownframe=amino_acid_list[nextaminoacidindex].avpos[Nindex];
                                                predicted_correcting_atom_pos_labframe.x=0;
                                                predicted_correcting_atom_pos_labframe.y=0;
                                                predicted_correcting_atom_pos_labframe.z=0;
                                                predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ex, predicted_correcting_atom_pos_ownframe.x));
                                                predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ey, predicted_correcting_atom_pos_ownframe.y));
                                                predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ez, predicted_correcting_atom_pos_ownframe.z));
                                                predicted_correcting_atom_pos_labframe=subtract_double_triple(add_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i+1].coord[0].pos), cgresiduearray[i].coord[0].pos);
                                                predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ex);
                                                predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ey);
                                                predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ez);
                                                predicted_position.x+=intercept_O[0];
                                                predicted_position.x+=(slope_O[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                                predicted_position.x+=(slope_O[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                                predicted_position.x+=(slope_O[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                                predicted_position.y+=intercept_O[1];
                                                predicted_position.y+=(slope_O[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                                predicted_position.y+=(slope_O[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                                predicted_position.y+=(slope_O[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                                predicted_position.z+=intercept_O[2];
                                                predicted_position.z+=(slope_O[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                                predicted_position.z+=(slope_O[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                                predicted_position.z+=(slope_O[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                            }
                                        }
                                    }
                                }
                            }
                            else if(j==Hindex){
                                dormsd=0;
                                if(residuearray[i].atomfound[Hindex]==1){
                                    if((i>0)&&(residuearray[i].chainid==residuearray[i-1].chainid)){
                                        if(cgresiduearray[i-1].sitedefined[0]==1){
                                            dormsd=1;
                                            previousaminoacidindex=residuearray[i-1].residnumber;
                                            Cindex=amino_acid_list[previousaminoacidindex].Cindex;
                                            predicted_correcting_atom_pos_ownframe=amino_acid_list[previousaminoacidindex].avpos[Cindex];
                                            predicted_correcting_atom_pos_labframe.x=0;
                                            predicted_correcting_atom_pos_labframe.y=0;
                                            predicted_correcting_atom_pos_labframe.z=0;
                                            predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ex, predicted_correcting_atom_pos_ownframe.x));
                                            predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ey, predicted_correcting_atom_pos_ownframe.y));
                                            predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ez, predicted_correcting_atom_pos_ownframe.z));
                                            predicted_correcting_atom_pos_labframe=subtract_double_triple(add_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i-1].coord[0].pos), cgresiduearray[i].coord[0].pos);
                                            predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ex);
                                            predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ey);
                                            predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ez);
                                            predicted_position.x+=intercept_H[0];
                                            predicted_position.x+=(slope_H[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                            predicted_position.x+=(slope_H[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                            predicted_position.x+=(slope_H[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                            predicted_position.y+=intercept_H[1];
                                            predicted_position.y+=(slope_H[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                            predicted_position.y+=(slope_H[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                            predicted_position.y+=(slope_H[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                            predicted_position.z+=intercept_H[2];
                                            predicted_position.z+=(slope_H[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                            predicted_position.z+=(slope_H[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                            predicted_position.z+=(slope_H[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                        }
                                    }
                                }
                            }
                            else if(amino_acid_list[aminoacidindex].atom_correction_index_withinaa[j]>=0){
                                correctionindexwithinaa=amino_acid_list[aminoacidindex].atom_correction_index_withinaa[j];
                                dormsd=0;
                                my_corrected_atom_index=amino_acid_list[aminoacidindex].corrected_atom_index[correctionindexwithinaa];
                                my_correcting_atom_index=amino_acid_list[aminoacidindex].correcting_atom_index[correctionindexwithinaa];
                                corrected_atom_site=amino_acid_list[aminoacidindex].sitemap[my_corrected_atom_index];
                                correcting_atom_site=amino_acid_list[aminoacidindex].sitemap[my_correcting_atom_index];
                                if((residuearray[i].atomfound[my_corrected_atom_index]==1)&&(residuearray[i].atomfound[my_correcting_atom_index]==1)){
                                    if((cgresiduearray[i].sitedefined[corrected_atom_site]==1)&&(cgresiduearray[i].sitedefined[correcting_atom_site]==1)){
                                        dormsd=1;
                                        predicted_correcting_atom_pos_ownframe=amino_acid_list[aminoacidindex].avpos[my_correcting_atom_index];
                                        predicted_correcting_atom_pos_labframe=cgresiduearray[i].coord[correcting_atom_site].pos;
                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ex, predicted_correcting_atom_pos_ownframe.x));
                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ey, predicted_correcting_atom_pos_ownframe.y));
                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ez, predicted_correcting_atom_pos_ownframe.z));
                                        predicted_correcting_atom_pos_labframe=subtract_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].pos);
                                        predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ex);
                                        predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ey);
                                        predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ez);
                                        myindex=amino_acid_list[aminoacidindex].atom_correction_index[j];
                                        predicted_position.x+=sidechain_corrections[myindex].intercept[0];
                                        predicted_position.x+=(sidechain_corrections[myindex].slope[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                        predicted_position.x+=(sidechain_corrections[myindex].slope[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                        predicted_position.x+=(sidechain_corrections[myindex].slope[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                        predicted_position.y+=sidechain_corrections[myindex].intercept[1];
                                        predicted_position.y+=(sidechain_corrections[myindex].slope[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                        predicted_position.y+=(sidechain_corrections[myindex].slope[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                        predicted_position.y+=(sidechain_corrections[myindex].slope[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                        predicted_position.z+=sidechain_corrections[myindex].intercept[2];
                                        predicted_position.z+=(sidechain_corrections[myindex].slope[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                        predicted_position.z+=(sidechain_corrections[myindex].slope[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                        predicted_position.z+=(sidechain_corrections[myindex].slope[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                    }
                                }
                            }
                            else dormsd=1;
                            if(dormsd==1){
                                difference=subtract_double_triple(position_frame, predicted_position);
                                difference.x=difference.x*difference.x;
                                difference.y=difference.y*difference.y;
                                difference.z=difference.z*difference.z;
                                amino_acid_list[aminoacidindex].msd_residuespecific[index]=add_double_triple(amino_acid_list[aminoacidindex].msd_residuespecific[index], difference);
                                distancesquared=difference.x+difference.y+difference.z;
                                if(distancesquared>distancesquared_error_threshold){
                                    outp=fopen(large_disagreement_file_specific, "a");
                                    fprintf(outp, "%s\t%s\t%i\t%i\t%s\t%i\t%s\t%g\t%f\t%f\t%f\t%f\t%f\t%f\n", pdb_file, residuearray[i].chainidchar, residuearray[i].resseq, i, amino_acid_list[aminoacidindex].resname, j, amino_acid_list[aminoacidindex].atomnames[j], sqrt(distancesquared), position_frame.x, position_frame.y, position_frame.z, predicted_position.x, predicted_position.y, predicted_position.z);
                                    fclose(outp);
                                }
                                amino_acid_list[aminoacidindex].fourthmoment_residuespecific[index]+=pow(distancesquared, 2);
                                if(amino_acid_list[aminoacidindex].element[index]==1){
                                    msd_entry_average_H_residuespecific+=distancesquared;
                                    fourthmoment_entry_average_H_residuespecific+=pow(distancesquared, 2);
                                    entry_average_H_count_residuespecific++;
                                }
                                else{
                                    msd_entry_average_heavy_residuespecific+=distancesquared;
                                    fourthmoment_entry_average_heavy_residuespecific+=pow(distancesquared, 2);
                                    entry_average_heavy_count_residuespecific++;
                                }
                                amino_acid_list[aminoacidindex].msdcount[index]++;
                            }
						}
					}
				}
			}
		}
	}
    msd_entry_average_heavy_residuespecific=msd_entry_average_heavy_residuespecific/(1.*entry_average_heavy_count_residuespecific);
    fourthmoment_entry_average_heavy_residuespecific=3./5.*fourthmoment_entry_average_heavy_residuespecific/(1.*entry_average_heavy_count_residuespecific)/pow(msd_entry_average_heavy_residuespecific, 2);
    msd_entry_average_heavy_residuespecific=sqrt(msd_entry_average_heavy_residuespecific);
    msd_entry_average_H_residuespecific=msd_entry_average_H_residuespecific/(1.*entry_average_H_count_residuespecific);
    fourthmoment_entry_average_H_residuespecific=3./5.*fourthmoment_entry_average_H_residuespecific/(1.*entry_average_H_count_residuespecific)/pow(msd_entry_average_H_residuespecific, 2);
    msd_entry_average_H_residuespecific=sqrt(msd_entry_average_H_residuespecific);
    outp=fopen(entry_average_file, "a");
    fprintf(outp, "%s", pdb_file);
    if(entry_average_heavy_count_residuespecific>0) fprintf(outp, "\t%.8f\t%.8f\t%i", msd_entry_average_heavy_residuespecific, fourthmoment_entry_average_heavy_residuespecific, entry_average_heavy_count_residuespecific, entry_average_H_count_residuespecific);
    else fprintf(outp, "\tN/A\t\tN/A\t\t%i", entry_average_heavy_count_residuespecific);
    if(entry_average_H_count_residuespecific>0) fprintf(outp, "\t%.8f\t%.8f\t%i\n", msd_entry_average_H_residuespecific, fourthmoment_entry_average_H_residuespecific, entry_average_H_count_residuespecific);
    else fprintf(outp, "\tN/A\t\tN/A\t\t%i\n", entry_average_H_count_residuespecific);
    fclose(outp);
    
    free(intercept_O);
    free(intercept_H);
    free_matrix(slope_O, 3);
    free_matrix(slope_H, 3);
}

void output_rmsd_onlyspecific(char *filename, amino_acid_struct *amino_acid_list, int nrecords, char **backbone_atom_types, int N_backbone_atom_types, char **Cterminusatomnames){
    int i, j, k, found, totalcount_heavy=0, totalcount_H=0, backboneaveragecount_heavy=0, backboneaveragecount_H=0, inbackbone, sitecountsum, sidesitecountsum, potentialatomsum, atomsperresidue, backbone_atom_type_O;
    double mag, totalaverage_heavy=0, totalaverage_H=0, backboneaveragemsd_heavy=0, backboneaveragemsd_H=0;
    declare_array(double, residueaveragemsd_H, Naminoacids);
    declare_array(double, residueaveragemsd_H_sidechain, Naminoacids);
    declare_array(double, residueaveragemsd_heavy, Naminoacids);
    declare_array(double, residueaveragemsd_heavy_sidechain, Naminoacids);
    declare_matrix(double, msdmag, Naminoacids, maxatomsperresidue);
    declare_array(int, residueaveragecount_H, Naminoacids);
    declare_array(int, residueaveragecount_H_sidechain, Naminoacids);
    declare_array(int, residueaveragecount_heavy, Naminoacids);
    declare_array(int, residueaveragecount_heavy_sidechain, Naminoacids);
    declare_array_nozero(double_triple, atomaveragemsd, N_backbone_atom_types);
    for(i=0;i<N_backbone_atom_types;i++){
        atomaveragemsd[i].x=0;
        atomaveragemsd[i].y=0;
        atomaveragemsd[i].z=0;
        if(strcmp(backbone_atom_types[i], " O  ")==0) backbone_atom_type_O=i;
    }
    declare_array(int, atomaveragecount, N_backbone_atom_types);
    declare_array(int, atomaveragesitecount, N_backbone_atom_types);
    declare_array_nozero(double_triple, Cterminus_atomaveragemsd, amino_acid_list[0].NCterminusatoms);
    for(i=0;i<amino_acid_list[0].NCterminusatoms;i++){
        Cterminus_atomaveragemsd[i].x=0;
        Cterminus_atomaveragemsd[i].y=0;
        Cterminus_atomaveragemsd[i].z=0;
    }
    declare_array(int, Cterminus_atomaveragecount, amino_acid_list[0].NCterminusatoms);
    declare_matrix(int, atomspersite, Naminoacids, maxsitesperresidue);
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].Natoms;j++){
            found=0;
            if(strcmp(amino_acid_list[i].atomnames[j], " OXT")==0) found=1;
            else if(strcmp(amino_acid_list[i].atomnames[j], " HXT")==0) found=1;
            else if(strcmp(amino_acid_list[i].atomnames[j], " H2 ")==0) found=1;
            else if(strcmp(amino_acid_list[i].atomnames[j], " H1 ")==0) found=1;
            else if(strcmp(amino_acid_list[i].atomnames[j], " H3 ")==0) found=1;
            if(found==0){
                atomspersite[i][amino_acid_list[i].sitemap[j]]++;
            }
        }
    }
    declare_array(int, backboneH, N_backbone_atom_types);
    for(i=0;i<N_backbone_atom_types;i++) backboneH[i]=-1;
    FILE *outp;
    outp=fopen(filename, "w");
    for(i=0;i<Naminoacids;i++){
        for(j=0;j<amino_acid_list[i].Natoms+amino_acid_list[i].NCterminusatoms;j++){
            if(amino_acid_list[i].msdcount[j]>0){
                
                //  Sum over atoms in residue before divide by N and sqrt
                
                if(amino_acid_list[i].element[j]==1){
                    residueaveragemsd_H[i]+=(amino_acid_list[i].msd_residuespecific[j].x+amino_acid_list[i].msd_residuespecific[j].y+amino_acid_list[i].msd_residuespecific[j].z);
                    residueaveragecount_H[i]+=amino_acid_list[i].msdcount[j];
                }
                else{
                    residueaveragemsd_heavy[i]+=(amino_acid_list[i].msd_residuespecific[j].x+amino_acid_list[i].msd_residuespecific[j].y+amino_acid_list[i].msd_residuespecific[j].z);
                    residueaveragecount_heavy[i]+=amino_acid_list[i].msdcount[j];
                }
                
                //  Sum over residues for atoms in backbone
                
                inbackbone=0;
                if((j<amino_acid_list[i].Natoms)&&(amino_acid_list[i].sitemap[j]==0)){
                    inbackbone=1;
                }
                
                if(inbackbone==1){
                    
                    //  Calculate average over amino acids for backbone sites
                    
                    found=0;
                    for(k=0;k<N_backbone_atom_types;k++){
                        if(strcmp(amino_acid_list[i].atomnames[j], backbone_atom_types[k])==0){
                            found=1;
                            break;
                        }
                    }
                    if(found==0) my_exit("Backbone atom type not found in output_rmsd!");
                    if(backboneH[k]==0){
                        if(amino_acid_list[i].element[j]<2){
                            my_exit("backbone H/heavy not consistent!\n");
                        }
                    }
                    else if(backboneH[k]==1){
                        if(amino_acid_list[i].element[j]!=1){
                            my_exit("backbone H/heavy not consistent!\n");
                        }
                    }
                    else{
                        if(amino_acid_list[i].element[j]==1) backboneH[k]=1;
                        else backboneH[k]=0;
                    }
                    atomaveragemsd[k]=add_double_triple(atomaveragemsd[k], amino_acid_list[i].msd_residuespecific[j]);
                    atomaveragecount[k]+=amino_acid_list[i].msdcount[j];
                    atomaveragesitecount[k]+=amino_acid_list[i].sitecount[0];
                }
                else if(j>=amino_acid_list[i].Natoms){
                    
                    //  C-terminus atoms
                    
                    k=j-amino_acid_list[i].Natoms;
                    Cterminus_atomaveragemsd[k]=add_double_triple(Cterminus_atomaveragemsd[k], amino_acid_list[i].msd_residuespecific[j]);
                    Cterminus_atomaveragecount[k]+=amino_acid_list[i].msdcount[j];
                }
                else{
                    if(amino_acid_list[i].element[j]==1){
                        residueaveragemsd_H_sidechain[i]+=(amino_acid_list[i].msd_residuespecific[j].x+amino_acid_list[i].msd_residuespecific[j].y+amino_acid_list[i].msd_residuespecific[j].z);
                        residueaveragecount_H_sidechain[i]+=amino_acid_list[i].msdcount[j];
                    }
                    else{
                        residueaveragemsd_heavy_sidechain[i]+=(amino_acid_list[i].msd_residuespecific[j].x+amino_acid_list[i].msd_residuespecific[j].y+amino_acid_list[i].msd_residuespecific[j].z);
                        residueaveragecount_heavy_sidechain[i]+=amino_acid_list[i].msdcount[j];
                    }
                }
                
                //  Now divide by N and sqrt
                
                amino_acid_list[i].msd_residuespecific[j]=scalar_multiply_double_triple(amino_acid_list[i].msd_residuespecific[j], 1./(1.*amino_acid_list[i].msdcount[j]));
                msdmag[i][j]=amino_acid_list[i].msd_residuespecific[j].x+amino_acid_list[i].msd_residuespecific[j].y+amino_acid_list[i].msd_residuespecific[j].z;
                if(amino_acid_list[i].fourthmoment_residuespecific[j]!=0){
                    amino_acid_list[i].fourthmoment_residuespecific[j]=3./5.*amino_acid_list[i].fourthmoment_residuespecific[j]/(1.*amino_acid_list[i].msdcount[j])/pow(msdmag[i][j], 2);
                }
                msdmag[i][j]=sqrt(msdmag[i][j]);
                amino_acid_list[i].msd_residuespecific[j].x=sqrt(amino_acid_list[i].msd_residuespecific[j].x);
                amino_acid_list[i].msd_residuespecific[j].y=sqrt(amino_acid_list[i].msd_residuespecific[j].y);
                amino_acid_list[i].msd_residuespecific[j].z=sqrt(amino_acid_list[i].msd_residuespecific[j].z);
            }
        }
        sitecountsum=0;
        potentialatomsum=0;
        atomsperresidue=0;
        sidesitecountsum=0;
        for(j=0;j<amino_acid_list[i].Nsites;j++){
            sitecountsum+=(amino_acid_list[i].sitecount[j]);
            potentialatomsum+=(amino_acid_list[i].sitecount[j]*atomspersite[i][j]);
            atomsperresidue+=atomspersite[i][j];
            if(j>0) sidesitecountsum+=(amino_acid_list[i].sitecount[j]);
        }
        for(j=0;j<amino_acid_list[i].Natoms;j++){
            if(amino_acid_list[i].msdcount[j]>0){
                fprintf(outp, "%s\tRes-spec\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%i\t%i\t%.8f\n", amino_acid_list[i].resname, amino_acid_list[i].atomnames[j], amino_acid_list[i].avpos[j].x, amino_acid_list[i].avpos[j].y, amino_acid_list[i].avpos[j].z, amino_acid_list[i].msd_residuespecific[j].x, amino_acid_list[i].msd_residuespecific[j].y, amino_acid_list[i].msd_residuespecific[j].z, msdmag[i][j], amino_acid_list[i].msdcount[j], amino_acid_list[i].sitecount[amino_acid_list[i].sitemap[j]], amino_acid_list[i].fourthmoment_residuespecific[j]);
            }
        }
        for(j=amino_acid_list[i].Natoms;j<amino_acid_list[i].Natoms+amino_acid_list[i].NCterminusatoms;j++){
            if(amino_acid_list[i].msdcount[j]>0){
                k=j-amino_acid_list[i].Natoms;
                fprintf(outp, "%s\tRes-spec\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%i\t?\t%.8f\n", amino_acid_list[i].resname, Cterminusatomnames[k], amino_acid_list[i].avpos[j].x, amino_acid_list[i].avpos[j].y, amino_acid_list[i].avpos[j].z, amino_acid_list[i].msd_residuespecific[j].x, amino_acid_list[i].msd_residuespecific[j].y, amino_acid_list[i].msd_residuespecific[j].z, msdmag[i][j], amino_acid_list[i].msdcount[j], amino_acid_list[i].fourthmoment_residuespecific[j]);
            }
        }
        if(residueaveragecount_heavy[i]>0) fprintf(outp, "%s\tRes-spec\tHeavy\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", amino_acid_list[i].resname, sqrt(residueaveragemsd_heavy[i]/residueaveragecount_heavy[i]), residueaveragecount_heavy[i], sitecountsum);
        else fprintf(outp, "%s\tRes-spec\tHeavy\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%i\t%i\n", amino_acid_list[i].resname, residueaveragecount_heavy[i], sitecountsum);
        if(residueaveragecount_H[i]>0) fprintf(outp, "%s\tRes-spec\tH\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", amino_acid_list[i].resname, sqrt(residueaveragemsd_H[i]/residueaveragecount_H[i]), residueaveragecount_H[i], sitecountsum);
        else fprintf(outp, "%s\tRes-spec\tH\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%i\t%i\n", amino_acid_list[i].resname, residueaveragecount_H[i], sitecountsum);
        if((residueaveragecount_heavy[i]+residueaveragecount_H[i])>0) fprintf(outp, "%s\tRes-spec\tOverall\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\tN/A\t\t%i\t%i\t%i\n\n", amino_acid_list[i].resname, sqrt((residueaveragemsd_heavy[i]+residueaveragemsd_H[i])/(residueaveragecount_heavy[i]+residueaveragecount_H[i])), residueaveragecount_heavy[i]+residueaveragecount_H[i], sitecountsum, potentialatomsum, atomsperresidue, amino_acid_list[i].Nsites);
        else fprintf(outp, "%s\tRes-spec\tOverall\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%i\t%i\tN/A\t\t%i\t%i\t%i\n\n", amino_acid_list[i].resname, residueaveragecount_heavy[i]+residueaveragecount_H[i], sitecountsum, potentialatomsum, atomsperresidue, amino_acid_list[i].Nsites);
        totalaverage_heavy+=residueaveragemsd_heavy[i];
        totalaverage_H+=residueaveragemsd_H[i];
        totalcount_heavy+=residueaveragecount_heavy[i];
        totalcount_H+=residueaveragecount_H[i];
        if(residueaveragecount_heavy_sidechain[i]>0) fprintf(outp, "%s\tSidechain\tHeavy\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", amino_acid_list[i].resname, sqrt(residueaveragemsd_heavy_sidechain[i]/residueaveragecount_heavy_sidechain[i]), residueaveragecount_heavy_sidechain[i], sidesitecountsum);
        if(residueaveragecount_H_sidechain[i]>0) fprintf(outp, "%s\tSidechain\tH\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", amino_acid_list[i].resname, sqrt(residueaveragemsd_H_sidechain[i]/residueaveragecount_H_sidechain[i]), residueaveragecount_H_sidechain[i], sidesitecountsum);
        if(residueaveragecount_heavy_sidechain[i]+residueaveragecount_H_sidechain[i]>0) fprintf(outp, "%s\tSidechain\tOverall\t\tN/A\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", amino_acid_list[i].resname, sqrt((residueaveragemsd_heavy_sidechain[i]+residueaveragemsd_H_sidechain[i])/(residueaveragecount_heavy_sidechain[i]+residueaveragecount_H_sidechain[i])), residueaveragecount_heavy_sidechain[i]+residueaveragecount_H_sidechain[i], sidesitecountsum);
        fprintf(outp, "\n");
    }
    sitecountsum=0;
    potentialatomsum=0;
    sidesitecountsum=0;
    for(i=0;i<Naminoacids;i++){
        sitecountsum+=(amino_acid_list[i].sitecount[0]);
        for(j=0;j<amino_acid_list[i].Nsites;j++) potentialatomsum+=(amino_acid_list[i].sitecount[j]*atomspersite[i][j]);
        for(j=1;j<amino_acid_list[i].Nsites;j++) sidesitecountsum+=(amino_acid_list[i].sitecount[j]);
    }
    for(i=0;i<N_backbone_atom_types;i++){
        if(atomaveragecount[i]>0){
            if(backboneH[i]==1){
                backboneaveragemsd_H+=(atomaveragemsd[i].x+atomaveragemsd[i].y+atomaveragemsd[i].z);
                backboneaveragecount_H+=atomaveragecount[i];
            }
            else{
                backboneaveragemsd_heavy+=(atomaveragemsd[i].x+atomaveragemsd[i].y+atomaveragemsd[i].z);
                backboneaveragecount_heavy+=atomaveragecount[i];
            }
            atomaveragemsd[i].x=atomaveragemsd[i].x/(1.*atomaveragecount[i]);
            atomaveragemsd[i].y=atomaveragemsd[i].y/(1.*atomaveragecount[i]);
            atomaveragemsd[i].z=atomaveragemsd[i].z/(1.*atomaveragecount[i]);
            mag=atomaveragemsd[i].x+atomaveragemsd[i].y+atomaveragemsd[i].z;
            atomaveragemsd[i].x=sqrt(atomaveragemsd[i].x);
            atomaveragemsd[i].y=sqrt(atomaveragemsd[i].y);
            atomaveragemsd[i].z=sqrt(atomaveragemsd[i].z);
            mag=sqrt(mag);
            fprintf(outp, "Average\tRes-spec\t%s\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%.8f\t%.8f\t%.8f\t%i\t%i\n", backbone_atom_types[i], atomaveragemsd[i].x, atomaveragemsd[i].y, atomaveragemsd[i].z, mag, atomaveragecount[i], atomaveragesitecount[i]);
        }
    }
    for(i=0;i<amino_acid_list[0].NCterminusatoms;i++){
        if(Cterminus_atomaveragecount[i]>0){
            if(amino_acid_list[0].element[amino_acid_list[0].Natoms+i]==1){
                backboneaveragemsd_H+=(Cterminus_atomaveragemsd[i].x+Cterminus_atomaveragemsd[i].y+Cterminus_atomaveragemsd[i].z);
                backboneaveragecount_H+=Cterminus_atomaveragecount[i];
            }
            else{
                backboneaveragemsd_heavy+=(Cterminus_atomaveragemsd[i].x+Cterminus_atomaveragemsd[i].y+Cterminus_atomaveragemsd[i].z);
                backboneaveragecount_heavy+=Cterminus_atomaveragecount[i];
            }
            Cterminus_atomaveragemsd[i].x=Cterminus_atomaveragemsd[i].x/(1.*Cterminus_atomaveragecount[i]);
            Cterminus_atomaveragemsd[i].y=Cterminus_atomaveragemsd[i].y/(1.*Cterminus_atomaveragecount[i]);
            Cterminus_atomaveragemsd[i].z=Cterminus_atomaveragemsd[i].z/(1.*Cterminus_atomaveragecount[i]);
            mag=Cterminus_atomaveragemsd[i].x+Cterminus_atomaveragemsd[i].y+Cterminus_atomaveragemsd[i].z;
            Cterminus_atomaveragemsd[i].x=sqrt(Cterminus_atomaveragemsd[i].x);
            Cterminus_atomaveragemsd[i].y=sqrt(Cterminus_atomaveragemsd[i].y);
            Cterminus_atomaveragemsd[i].z=sqrt(Cterminus_atomaveragemsd[i].z);
            mag=sqrt(mag);
            fprintf(outp, "Average\tRes-spec\t%s\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%.8f\t%.8f\t%.8f\t%i\n", Cterminusatomnames[i], Cterminus_atomaveragemsd[i].x, Cterminus_atomaveragemsd[i].y, Cterminus_atomaveragemsd[i].z, mag, Cterminus_atomaveragecount[i]);
        }
    }
    fprintf(outp, "Average\tRes-spec\tBB hvy\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt(backboneaveragemsd_heavy/backboneaveragecount_heavy), backboneaveragecount_heavy, sitecountsum);
    fprintf(outp, "Average\tRes-spec\tBB H\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt(backboneaveragemsd_H/backboneaveragecount_H), backboneaveragecount_H, sitecountsum);
    fprintf(outp, "Average\tRes-spec\tBB\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt((backboneaveragemsd_heavy+backboneaveragemsd_H)/(backboneaveragecount_heavy+backboneaveragecount_H)), backboneaveragecount_heavy+backboneaveragecount_H, sitecountsum);
    fprintf(outp, "Average\tRes-spec\tSC hvy\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt((totalaverage_heavy-backboneaveragemsd_heavy)/(totalcount_heavy-backboneaveragecount_heavy)), totalcount_heavy-backboneaveragecount_heavy, sidesitecountsum);
    fprintf(outp, "Average\tRes-spec\tSC H\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt((totalaverage_H-backboneaveragemsd_H)/(totalcount_H-backboneaveragecount_H)), totalcount_H-backboneaveragecount_H, sidesitecountsum);
    fprintf(outp, "Average\tRes-spec\tSC\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt((totalaverage_heavy+totalaverage_H-backboneaveragemsd_heavy-backboneaveragemsd_H)/(totalcount_heavy+totalcount_H-backboneaveragecount_heavy-backboneaveragecount_H)), totalcount_heavy+totalcount_H-backboneaveragecount_heavy-backboneaveragecount_H, sidesitecountsum);
    fprintf(outp, "Average\tRes-spec\tHeavy\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt(totalaverage_heavy/totalcount_heavy), totalcount_heavy, sitecountsum+sidesitecountsum);
    fprintf(outp, "Average\tRes-spec\tH\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\n", sqrt(totalaverage_H/totalcount_H), totalcount_H, sitecountsum+sidesitecountsum);
    fprintf(outp, "Average\tRes-spec\tOverall\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\tN/A\t\t%.8f\t%i\t%i\tN/A\t\t%i\n\n", sqrt((totalaverage_heavy+totalaverage_H)/(totalcount_heavy+totalcount_H)), totalcount_heavy+totalcount_H, sitecountsum+sidesitecountsum, potentialatomsum);
    fclose(outp);
    free(residueaveragemsd_heavy);
    free(residueaveragemsd_heavy_sidechain);
    free(residueaveragemsd_H);
    free(residueaveragemsd_H_sidechain);
    free(msdmag);
    free(residueaveragecount_heavy);
    free(residueaveragecount_heavy_sidechain);
    free(residueaveragecount_H);
    free(residueaveragecount_H_sidechain);
    free(atomaveragemsd);
    free(atomaveragecount);
    free(atomaveragesitecount);
    free_matrix(atomspersite, Naminoacids);
    free(backboneH);
}

void output_file_type_statistics(int *Ntype, char *filename){
    FILE *outp;
    outp=fopen(filename, "w");
    if(Ntype[0]>0) fprintf(outp, "Single model of protein consisting only of natural amino acids\t%i\n", Ntype[0]);
    if(Ntype[1]>0) fprintf(outp, "Ensemble of models\t\t\t\t\t\t%i\n", Ntype[1]);
    if(Ntype[2]>0) fprintf(outp, "Single model of molecule with non-amino-acid residues\t\t%i\n", Ntype[2]);
    if(Ntype[3]>0) fprintf(outp, "Unknown atom found\t\t\t\t\t\t%i\n", Ntype[3]);
	if(Ntype[4]>0) fprintf(outp, "No DBREF entry\t\t\t\t\t\t\t%i\n", Ntype[4]);
	if(Ntype[5]>0) fprintf(outp, "No ATOM or DBREF entry\t\t\t\t\t\t%i\n", Ntype[5]);
	if(Ntype[6]>0) fprintf(outp, "Residue indices out of order within a single chain\t\t%i\n", Ntype[6]);
	if(Ntype[7]>0) fprintf(outp, "Found terminal atom types in non-terminal residues\t\t%i\n", Ntype[7]);
	if(Ntype[8]>0) fprintf(outp, "File not found\t\t\t\t\t\t\t%i\n", Ntype[8]);
    fclose(outp);
}

void extract_piece_of_string(char *instring, char *outstring, int start, int end){
    int i;
    for(i=0;i<=end-start;i++){
        outstring[i]=instring[start-1+i];
    }
    outstring[end-start+1]='\0';
}

void print_error(int type){
    if(type==1) printf("Eensemble!");
    else if(type==2) printf("Non-amino-acid!");
    else if(type==3) printf("Unknown atom!");
    else if(type==4) printf("No DBREF!");
    else if(type==5) printf("No DBREF!");
    else if(type==6) printf("Residue indicies out of order!");
    else if(type==7) printf("Terminal atom in non terminal residue!");
    else if(type==8) printf("File not found!");
    else printf("Bad type in output_error_filename");
}

void input_cg_coords(char *filename, amino_acid_struct *amino_acid_list, char **Cterminusatomnames){
    int i=0, j=0, k, found, read, Cterminusindex;
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, firstword, 7);
    declare_array_nozero(char, resname, 3);
    declare_array_nozero(char, atomname, 4);
    FILE *inp;
    inp=fopen(filename, "r");
    read=mygetline(linechar, maxlinelength, inp);
    while(i<Naminoacids){
        while(j<amino_acid_list[i].Natoms){
            //  Assuming that there were no zero counts in generating filename
            
            extract_piece_of_string(linechar, resname, 1, 3);
            extract_piece_of_string(linechar, atomname, 5, 8);
            if(strcmp(resname, amino_acid_list[i].resname)!=0){
                i++;
                j=0;
                break;
            }
            if(strcmp(atomname, amino_acid_list[i].atomnames[j])!=0){
                j++;
            }
            else{
                sscanf(linechar, "%*s\t%*s\t%lf\t%lf\t%lf", &(amino_acid_list[i].avpos[j].x), &(amino_acid_list[i].avpos[j].y), &(amino_acid_list[i].avpos[j].z));
                j++;
                read=mygetline(linechar, maxlinelength, inp);
            }
        }
        
        //  C terminus atoms
        
        while(j<amino_acid_list[i].Natoms+amino_acid_list[i].NCterminusatoms){
            Cterminusindex=j-amino_acid_list[i].Natoms;
            extract_piece_of_string(linechar, resname, 1, 3);
            extract_piece_of_string(linechar, atomname, 5, 8);
            if(strcmp(resname, amino_acid_list[i].resname)!=0){
                i++;
                j=0;
                break;
            }
            if(strcmp(atomname, Cterminusatomnames[Cterminusindex])!=0){
                j++;
            }
            else{
                sscanf(linechar, "%*s\t%*s\t%lf\t%lf\t%lf", &(amino_acid_list[i].avpos[j].x), &(amino_acid_list[i].avpos[j].y), &(amino_acid_list[i].avpos[j].z));
                j++;
                read=mygetline(linechar, maxlinelength, inp);
            }
        }
        if(j==amino_acid_list[i].Natoms+amino_acid_list[i].NCterminusatoms){
            j=0;
            i++;
        }
    }
    
    fclose(inp);
    free(linechar);
    free(firstword);
    free(resname);
    free(atomname);
}

void input_correction(correction_structure *pcorrection, char *corrections_file){
    int i=0, j=0, read;
    declare_array_nozero(char, linechar, maxlinelength);
    FILE *inp;
    inp=fopen(corrections_file, "r");
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    for(i=0;i<3;i++){
        read=mygetline(linechar, maxlinelength, inp);
    }
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    for(i=0;i<3;i++){
        read=mygetline(linechar, maxlinelength, inp);
    }
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    sscanf(linechar, "%lf\t%lf\t%lf", &((*pcorrection).intercept[0]), &((*pcorrection).intercept[1]), &((*pcorrection).intercept[2]));
    read=mygetline(linechar, maxlinelength, inp);
    read=mygetline(linechar, maxlinelength, inp);
    for(i=0;i<3;i++){
        read=mygetline(linechar, maxlinelength, inp);
        sscanf(linechar, "%lf\t%lf\t%lf", &((*pcorrection).slope[i][0]), &((*pcorrection).slope[i][1]), &((*pcorrection).slope[i][2]));
    }
    
    fclose(inp);
    free(linechar);
}

void create_molecular_files_separate_corrections(char *pdbfile, char *directory, char *filebase, amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray, cgresidue *cgresiduearray, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array, correction_structure *sidechain_corrections, int do_amino_acid_examples){
    int i, j, k, foundatom=0, read, aminoacidindex, Oindex, Hindex, Nindex, Cterminus, index, correctionindexwithinaa, my_corrected_atom_index, corrected_atom_site, correcting_atom_site, found, myindex, my_correcting_atom_index, site, doprint, previousaminoacidindex, nextaminoacidindex, Cindex, foundallatomsinresidue, foundNterminus=0, foundCterminus=0;
    double_triple predicted_position, predicted_correcting_atom_pos_ownframe, predicted_correcting_atom_pos_labframe, predicted_correcting_atom_pos_thisframe, predicted_position_labframe, *predicted_position_labframe_array, dummydouble;
    FILE *inp, *outp;
    declare_array_nozero(char, newpdbfile, maxlinelength);
    declare_array_nozero(char, amino_base, maxlinelength);
    sprintf(newpdbfile, "%s/%s.pdb", directory, filebase);
    declare_array_nozero(char, linechar, maxlinelength);
    declare_array_nozero(char, label, 6);
    declare_array_nozero(char, resseqchar, 4);
    declare_matrix_nozero(char, positionchar, 3, 8);
    declare_array(int, foundamino, Naminoacids);
    if(do_amino_acid_examples==1){
        allocate_array(double_triple, predicted_position_labframe_array, maxatomsperresidue);
    }
    inp=fopen(pdbfile, "r");
    outp=fopen(newpdbfile, "w");
	if(inp==NULL){
		printf("%s not found!\n", pdbfile);
		exit(1);
	}
    
    declare_array(double, intercept_O, 3);
    declare_array(double, intercept_H, 3);
    declare_matrix(double, slope_O, 3, 3);
    declare_matrix(double, slope_H, 3, 3);
    
    read=mygetline(linechar, maxlinelength, inp);
    while(read>0){
        extract_piece_of_string(linechar, label, 1, 6);
        if(strcmp(label, "ATOM  ")==0){
            if(foundatom==0){
                
                //  Print model-predicted atom coordinates
                
                for(i=0;i<Nresidues;i++){
                    aminoacidindex=residuearray[i].residnumber;
                    if(aminoacidindex>=0){
                        for(j=0;j<3;j++){
                            intercept_O[j]=carbonyl_correction_array[aminoacidindex].intercept[j];
                            intercept_H[j]=backboneH_correction_array[aminoacidindex].intercept[j];
                            for(k=0;k<3;k++){
                                slope_O[j][k]=carbonyl_correction_array[aminoacidindex].slope[j][k];
                                slope_H[j][k]=backboneH_correction_array[aminoacidindex].slope[j][k];
                            }
                        }
                        Oindex=amino_acid_list[aminoacidindex].Oindex;
                        Hindex=amino_acid_list[aminoacidindex].Hindex;
                        foundallatomsinresidue=1;
                        for(j=0;j<residuearray[i].Natoms;j++){
                            if(residuearray[i].atomfound[j]==1){
                                
                                if((cgresiduearray[i].Cterminussite>-1)&&(amino_acid_list[aminoacidindex].Cterminusindex[j]>-1)){
                                    site=cgresiduearray[i].Cterminussite;
                                    Cterminus=1;
                                }
                                else{
                                    site=amino_acid_list[aminoacidindex].sitemap[j];
                                    Cterminus=0;
                                }
                                if(site>=0){
                                    if(cgresiduearray[i].sitedefined[site]==1){
                                        if(Cterminus==1) index=amino_acid_list[aminoacidindex].Natoms+amino_acid_list[aminoacidindex].Cterminusindex[j];
                                        else index=j;
                                        
                                        predicted_position=amino_acid_list[aminoacidindex].avpos[index];
                                        
                                        if((Cterminus==0)&&(j==Oindex)){
                                            doprint=0;
                                            if(residuearray[i].atomfound[j]==1){
                                                if((i<Nresidues-1)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                                                    if(cgresiduearray[i+1].sitedefined[0]==1){
                                                        doprint=1;
                                                        nextaminoacidindex=residuearray[i+1].residnumber;
                                                        Nindex=amino_acid_list[nextaminoacidindex].Nindex;
                                                        predicted_correcting_atom_pos_ownframe=amino_acid_list[nextaminoacidindex].avpos[Nindex];
                                                        predicted_correcting_atom_pos_labframe.x=0;
                                                        predicted_correcting_atom_pos_labframe.y=0;
                                                        predicted_correcting_atom_pos_labframe.z=0;
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ex, predicted_correcting_atom_pos_ownframe.x));
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ey, predicted_correcting_atom_pos_ownframe.y));
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[0].ez, predicted_correcting_atom_pos_ownframe.z));
                                                        predicted_correcting_atom_pos_labframe=subtract_double_triple(add_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i+1].coord[0].pos), cgresiduearray[i].coord[0].pos);
                                                        predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ex);
                                                        predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ey);
                                                        predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ez);
                                                        predicted_position.x+=intercept_O[0];
                                                        predicted_position.x+=(slope_O[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.x+=(slope_O[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.x+=(slope_O[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                                        predicted_position.y+=intercept_O[1];
                                                        predicted_position.y+=(slope_O[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.y+=(slope_O[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.y+=(slope_O[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                                        predicted_position.z+=intercept_O[2];
                                                        predicted_position.z+=(slope_O[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.z+=(slope_O[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.z+=(slope_O[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                                    }
                                                }
                                            }
                                        }
                                        else if(j==Hindex){
                                            doprint=0;
                                            if(residuearray[i].atomfound[Hindex]==1){
                                                if((i>0)&&(residuearray[i].chainid==residuearray[i-1].chainid)){
                                                    if(cgresiduearray[i-1].sitedefined[0]==1){
                                                        doprint=1;
                                                        previousaminoacidindex=residuearray[i-1].residnumber;
                                                        Cindex=amino_acid_list[previousaminoacidindex].Cindex;
                                                        predicted_correcting_atom_pos_ownframe=amino_acid_list[previousaminoacidindex].avpos[Cindex];
                                                        predicted_correcting_atom_pos_labframe.x=0;
                                                        predicted_correcting_atom_pos_labframe.y=0;
                                                        predicted_correcting_atom_pos_labframe.z=0;
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ex, predicted_correcting_atom_pos_ownframe.x));
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ey, predicted_correcting_atom_pos_ownframe.y));
                                                        predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[0].ez, predicted_correcting_atom_pos_ownframe.z));
                                                        predicted_correcting_atom_pos_labframe=subtract_double_triple(add_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i-1].coord[0].pos), cgresiduearray[i].coord[0].pos);
                                                        predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ex);
                                                        predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ey);
                                                        predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[0].ez);
                                                        predicted_position.x+=intercept_H[0];
                                                        predicted_position.x+=(slope_H[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.x+=(slope_H[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.x+=(slope_H[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                                        predicted_position.y+=intercept_H[1];
                                                        predicted_position.y+=(slope_H[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.y+=(slope_H[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.y+=(slope_H[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                                        predicted_position.z+=intercept_H[2];
                                                        predicted_position.z+=(slope_H[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                                        predicted_position.z+=(slope_H[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                                        predicted_position.z+=(slope_H[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                                    }
                                                }
                                            }
                                        }
                                        else if(amino_acid_list[aminoacidindex].atom_correction_index_withinaa[j]>=0){
                                            correctionindexwithinaa=amino_acid_list[aminoacidindex].atom_correction_index_withinaa[j];
                                            doprint=0;
                                            my_corrected_atom_index=amino_acid_list[aminoacidindex].corrected_atom_index[correctionindexwithinaa];
                                            my_correcting_atom_index=amino_acid_list[aminoacidindex].correcting_atom_index[correctionindexwithinaa];
                                            corrected_atom_site=amino_acid_list[aminoacidindex].sitemap[my_corrected_atom_index];
                                            correcting_atom_site=amino_acid_list[aminoacidindex].sitemap[my_correcting_atom_index];
                                            if((residuearray[i].atomfound[my_corrected_atom_index]==1)&&(residuearray[i].atomfound[my_correcting_atom_index]==1)){
                                                if((cgresiduearray[i].sitedefined[corrected_atom_site]==1)&&(cgresiduearray[i].sitedefined[correcting_atom_site]==1)){
                                                    doprint=1;
                                                    predicted_correcting_atom_pos_ownframe=amino_acid_list[aminoacidindex].avpos[my_correcting_atom_index];
                                                    predicted_correcting_atom_pos_labframe=cgresiduearray[i].coord[correcting_atom_site].pos;
                                                    predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ex, predicted_correcting_atom_pos_ownframe.x));
                                                    predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ey, predicted_correcting_atom_pos_ownframe.y));
                                                    predicted_correcting_atom_pos_labframe=add_double_triple(predicted_correcting_atom_pos_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[correcting_atom_site].ez, predicted_correcting_atom_pos_ownframe.z));
                                                    predicted_correcting_atom_pos_labframe=subtract_double_triple(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].pos);
                                                    predicted_correcting_atom_pos_thisframe.x=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ex);
                                                    predicted_correcting_atom_pos_thisframe.y=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ey);
                                                    predicted_correcting_atom_pos_thisframe.z=dot_product(predicted_correcting_atom_pos_labframe, cgresiduearray[i].coord[corrected_atom_site].ez);
                                                    myindex=amino_acid_list[aminoacidindex].atom_correction_index[j];
                                                    predicted_position.x+=sidechain_corrections[myindex].intercept[0];
                                                    predicted_position.x+=(sidechain_corrections[myindex].slope[0][0]*predicted_correcting_atom_pos_thisframe.x);
                                                    predicted_position.x+=(sidechain_corrections[myindex].slope[0][1]*predicted_correcting_atom_pos_thisframe.y);
                                                    predicted_position.x+=(sidechain_corrections[myindex].slope[0][2]*predicted_correcting_atom_pos_thisframe.z);
                                                    predicted_position.y+=sidechain_corrections[myindex].intercept[1];
                                                    predicted_position.y+=(sidechain_corrections[myindex].slope[1][0]*predicted_correcting_atom_pos_thisframe.x);
                                                    predicted_position.y+=(sidechain_corrections[myindex].slope[1][1]*predicted_correcting_atom_pos_thisframe.y);
                                                    predicted_position.y+=(sidechain_corrections[myindex].slope[1][2]*predicted_correcting_atom_pos_thisframe.z);
                                                    predicted_position.z+=sidechain_corrections[myindex].intercept[2];
                                                    predicted_position.z+=(sidechain_corrections[myindex].slope[2][0]*predicted_correcting_atom_pos_thisframe.x);
                                                    predicted_position.z+=(sidechain_corrections[myindex].slope[2][1]*predicted_correcting_atom_pos_thisframe.y);
                                                    predicted_position.z+=(sidechain_corrections[myindex].slope[2][2]*predicted_correcting_atom_pos_thisframe.z);
                                                }
                                            }
                                        }
                                        else doprint=1;
                                        if(doprint==1){
                                            
                                            if(residuearray[i].resseq<10) sprintf(resseqchar, "   %i", (residuearray[i].resseq));
                                            else if(residuearray[i].resseq<100) sprintf(resseqchar, "  %i", (residuearray[i].resseq));
                                            else if(residuearray[i].resseq<1000) sprintf(resseqchar, " %i", (residuearray[i].resseq));
                                            else sprintf(resseqchar, "%i", (residuearray[i].resseq));
                                            
                                            predicted_position_labframe=cgresiduearray[i].coord[site].pos;
                                            predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[site].ex, predicted_position.x));
                                            predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[site].ey, predicted_position.y));
                                            predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i].coord[site].ez, predicted_position.z));
                                            if(do_amino_acid_examples==1){
                                                predicted_position_labframe_array[j]=predicted_position_labframe;
                                            }
                                            
                                            if(predicted_position_labframe.x<0){
                                                if(predicted_position_labframe.x>-10){
                                                    sprintf(positionchar[0], "  %.3f", (predicted_position_labframe.x));
                                                }
                                                else if(predicted_position_labframe.x>-100){
                                                    sprintf(positionchar[0], " %.3f", (predicted_position_labframe.x));
                                                }
                                                else if(predicted_position_labframe.x>-1000){
                                                    sprintf(positionchar[0], "%.3f", (predicted_position_labframe.x));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.x);
                                                    exit(1);
                                                }
                                            }
                                            else{
                                                if(predicted_position_labframe.x<10){
                                                    sprintf(positionchar[0], "   %.3f", (predicted_position_labframe.x));
                                                }
                                                else if(predicted_position_labframe.x<100){
                                                    sprintf(positionchar[0], "  %.3f", (predicted_position_labframe.x));
                                                }
                                                else if(predicted_position_labframe.x<1000){
                                                    sprintf(positionchar[0], " %.3f", (predicted_position_labframe.x));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.x);
                                                    exit(1);
                                                }
                                            }
                                            if(predicted_position_labframe.y<0){
                                                if(predicted_position_labframe.y>-10){
                                                    sprintf(positionchar[1], "  %.3f", (predicted_position_labframe.y));
                                                }
                                                else if(predicted_position_labframe.y>-100){
                                                    sprintf(positionchar[1], " %.3f", (predicted_position_labframe.y));
                                                }
                                                else if(predicted_position_labframe.y>-1000){
                                                    sprintf(positionchar[1], "%.3f", (predicted_position_labframe.y));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.y);
                                                    exit(1);
                                                }
                                            }
                                            else{
                                                if(predicted_position_labframe.y<10){
                                                    sprintf(positionchar[1], "   %.3f", (predicted_position_labframe.y));
                                                }
                                                else if(predicted_position_labframe.y<100){
                                                    sprintf(positionchar[1], "  %.3f", (predicted_position_labframe.y));
                                                }
                                                else if(predicted_position_labframe.y<1000){
                                                    sprintf(positionchar[1], " %.3f", (predicted_position_labframe.y));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.y);
                                                    exit(1);
                                                }
                                            }
                                            if(predicted_position_labframe.z<0){
                                                if(predicted_position_labframe.z>-10){
                                                    sprintf(positionchar[2], "  %.3f", (predicted_position_labframe.z));
                                                }
                                                else if(predicted_position_labframe.z>-100){
                                                    sprintf(positionchar[2], " %.3f", (predicted_position_labframe.z));
                                                }
                                                else if(predicted_position_labframe.z>-1000){
                                                    sprintf(positionchar[2], "%.3f", (predicted_position_labframe.z));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.z);
                                                    exit(1);
                                                }
                                            }
                                            else{
                                                if(predicted_position_labframe.z<10){
                                                    sprintf(positionchar[2], "   %.3f", (predicted_position_labframe.z));
                                                }
                                                else if(predicted_position_labframe.z<100){
                                                    sprintf(positionchar[2], "  %.3f", (predicted_position_labframe.z));
                                                }
                                                else if(predicted_position_labframe.z<1000){
                                                    sprintf(positionchar[2], " %.3f", (predicted_position_labframe.z));
                                                }
                                                else{
                                                    printf("position=%f!\n", predicted_position_labframe.z);
                                                    exit(1);
                                                }
                                            }
                                            fprintf(outp, "%s%s %s %s %s%s%s   %s%s%s%s%s          %s%s\n", label, residuearray[i].serialchar[j], amino_acid_list[aminoacidindex].input_atomnames[j], amino_acid_list[aminoacidindex].resname, residuearray[i].chainidchar, resseqchar, residuearray[i].icodechar[j], positionchar[0], positionchar[1], positionchar[2], residuearray[i].occupancychar[j], residuearray[i].tempfactorchar[j], residuearray[i].elementchar[j], residuearray[i].chargechar[j]);
                                            
                                        }
                                        else{
                                            foundallatomsinresidue=0;
                                        }
                                    }
                                }
                            }
                            else if(residuearray[i].atomfound[j]==0){
                                
                                //  Incomplete if non-terminal atom type not found, but okay if terminal atom type not found
                                
                                if(j<amino_acid_list[aminoacidindex].Natoms-number_amino_hydrogen_atom_types){
                                    if((strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H  ")!=0)&&(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " H2 ")!=0)&&(strcmp(amino_acid_list[aminoacidindex].atomnames[j], " HXT")!=0)){
                                        foundallatomsinresidue=0;
                                    }
                                }
                            }
                        }
                        if(foundallatomsinresidue==1){
                            if(do_amino_acid_examples==1){
                                if(foundamino[aminoacidindex]==0){
                                    if((i>0)&&(i<Nresidues)&&(residuearray[i].chainid==residuearray[i-1].chainid)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                                        if(residuearray[i].Natoms+2>maxatomsperresidue){
                                            my_exit("Too many atoms in create_molecular_files");
                                        }
                                        previousaminoacidindex=residuearray[i-1].residnumber;
                                        if(previousaminoacidindex>-1){
                                            index=amino_acid_list[previousaminoacidindex].Cindex;
                                            if(residuearray[i-1].atomfound[index]==1){
                                                site=amino_acid_list[previousaminoacidindex].sitemap[index];
                                                if(cgresiduearray[i-1].sitedefined[site]==1){
                                                    predicted_position=amino_acid_list[previousaminoacidindex].avpos[index];
                                                    predicted_position_labframe=cgresiduearray[i-1].coord[site].pos;
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ex, predicted_position.x));
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ey, predicted_position.y));
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ez, predicted_position.z));
                                                    predicted_position_labframe_array[residuearray[i].Natoms]=predicted_position_labframe;
                                                }
                                            }
                                            else foundallatomsinresidue=0;
                                        }
                                        else foundallatomsinresidue=0;
                                        nextaminoacidindex=residuearray[i+1].residnumber;
                                        if(nextaminoacidindex>-1){
                                            index=amino_acid_list[nextaminoacidindex].Nindex;
                                            if(residuearray[i+1].atomfound[amino_acid_list[residuearray[i+1].residnumber].Nindex]==1){
                                                site=amino_acid_list[nextaminoacidindex].sitemap[index];
                                                if(cgresiduearray[i+1].sitedefined[site]==1){
                                                    predicted_position=amino_acid_list[nextaminoacidindex].avpos[index];
                                                    predicted_position_labframe=cgresiduearray[i+1].coord[site].pos;
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ex, predicted_position.x));
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ey, predicted_position.y));
                                                    predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ez, predicted_position.z));
                                                    predicted_position_labframe_array[residuearray[i].Natoms+1]=predicted_position_labframe;
                                                }
                                            }
                                            else foundallatomsinresidue=0;
                                        }
                                        else foundallatomsinresidue=0;
                                        if(foundallatomsinresidue==1){
                                            sprintf(amino_base, "%s.%s", filebase, amino_acid_list[aminoacidindex].resname);
                                            output_residue_molecular_files(directory, amino_base, amino_acid_list[aminoacidindex], residuearray[i], cgresiduearray[i], predicted_position_labframe_array, residuearray[i-1].atomposition[amino_acid_list[residuearray[i-1].residnumber].Cindex], residuearray[i+1].atomposition[amino_acid_list[residuearray[i+1].residnumber].Nindex], 1, 1, 0);
                                            foundamino[aminoacidindex]=1;
                                        }
                                    }
                                }
                                if(foundNterminus==0){
                                    if((i==0)||(residuearray[i].chainid!=residuearray[i-1].chainid)){
                                        if((i<Nresidues)&&(residuearray[i].chainid==residuearray[i+1].chainid)){
                                            foundallatomsinresidue=1;
                                            nextaminoacidindex=residuearray[i+1].residnumber;
                                            if(nextaminoacidindex>-1){
                                                index=amino_acid_list[nextaminoacidindex].Nindex;
                                                if(residuearray[i+1].atomfound[amino_acid_list[residuearray[i+1].residnumber].Nindex]==1){
                                                    site=amino_acid_list[nextaminoacidindex].sitemap[index];
                                                    if(cgresiduearray[i+1].sitedefined[site]==1){
                                                        predicted_position=amino_acid_list[nextaminoacidindex].avpos[index];
                                                        predicted_position_labframe=cgresiduearray[i+1].coord[site].pos;
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ex, predicted_position.x));
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ey, predicted_position.y));
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i+1].coord[site].ez, predicted_position.z));
                                                        predicted_position_labframe_array[residuearray[i].Natoms+1]=predicted_position_labframe;
                                                    }
                                                }
                                                else{
                                                    foundallatomsinresidue=0;
                                                }
                                            }
                                            else{
                                                foundallatomsinresidue=0;
                                            }
                                            if(foundallatomsinresidue==1){
                                                sprintf(amino_base, "%s.Nterm.%s", filebase, amino_acid_list[aminoacidindex].resname);
                                                output_residue_molecular_files(directory, amino_base, amino_acid_list[aminoacidindex], residuearray[i], cgresiduearray[i], predicted_position_labframe_array, dummydouble, residuearray[i+1].atomposition[amino_acid_list[residuearray[i+1].residnumber].Nindex], 0, 1, 0);
                                                foundNterminus=1;
                                            }
                                        }
                                    }
                                }
                                if(foundCterminus==0){
                                    if((i==Nresidues-1)||(residuearray[i].chainid!=residuearray[i+1].chainid)){
                                        if((i>0)&&(residuearray[i].chainid==residuearray[i-1].chainid)){
                                            foundallatomsinresidue=1;
                                            previousaminoacidindex=residuearray[i-1].residnumber;
                                            if(previousaminoacidindex>-1){
                                                index=amino_acid_list[previousaminoacidindex].Cindex;
                                                if(residuearray[i-1].atomfound[index]==1){
                                                    site=amino_acid_list[previousaminoacidindex].sitemap[index];
                                                    if(cgresiduearray[i-1].sitedefined[site]==1){
                                                        predicted_position=amino_acid_list[previousaminoacidindex].avpos[index];
                                                        predicted_position_labframe=cgresiduearray[i-1].coord[site].pos;
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ex, predicted_position.x));
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ey, predicted_position.y));
                                                        predicted_position_labframe=add_double_triple(predicted_position_labframe, scalar_multiply_double_triple(cgresiduearray[i-1].coord[site].ez, predicted_position.z));
                                                        predicted_position_labframe_array[residuearray[i].Natoms]=predicted_position_labframe;
                                                    }
                                                }
                                                else foundallatomsinresidue=0;
                                            }
                                            else foundallatomsinresidue=0;
                                            if(foundallatomsinresidue==1){
                                                sprintf(amino_base, "%s.Cterm.%s", filebase, amino_acid_list[aminoacidindex].resname);
                                                output_residue_molecular_files(directory, amino_base, amino_acid_list[aminoacidindex], residuearray[i], cgresiduearray[i], predicted_position_labframe_array, residuearray[i-1].atomposition[amino_acid_list[residuearray[i-1].residnumber].Cindex], dummydouble, 1, 0, 0);
                                                foundCterminus=1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
            }
            foundatom=1;
        }
        else{
            fprintf(outp, "%s", linechar);
        }
        read=mygetline(linechar, maxlinelength, inp);
    }
    int total=0;
    for(i=0;i<Naminoacids;i++){
        if(foundamino[i]==1) total++;
    }
    printf("%i amino acids found\n", total);
    if(foundNterminus==1) printf("N terminus found\n");
    if(foundCterminus==1) printf("C terminus found\n");
    fclose(inp);
    fclose(outp);
    free(linechar);
    free(label);
    free(resseqchar);
    free_matrix(positionchar, 3);
    free(newpdbfile);
    free(foundamino);
    free(amino_base);
    if(do_amino_acid_examples==1){
        free(predicted_position_labframe_array);
    }
    free(intercept_O);
    free(intercept_H);
    free_matrix(slope_O, 3);
    free_matrix(slope_H, 3);
}

void output_residue_molecular_files(char *directory, char *filebase, amino_acid_struct myaminoacid, residuedata myresidue, cgresidue mycgresidue, double_triple *predictedatoms, double_triple lastCpos, double_triple nextNpos, int doleft, int doright, int doCterminus){
    int i, j, k, heavycode, Natoms=0, Hcount=0;
    double Hrmsd=0, heavyrmsd=0;
    declare_3d_tensor(int, foundmatrix, maxsitesperresidue, 3, 2);
    declare_array(int, correctiontype, myresidue.Natoms);
    declare_array(int, atomcounter, myresidue.Natoms);
    declare_array_nozero(char, element, 1);
    declare_array_nozero(char, originalxyzfilename, maxstringlength);
    declare_array_nozero(char, originalxyzfile, maxstringlength);
    declare_array_nozero(char, roundtripxyzfilename, maxstringlength);
    declare_array_nozero(char, roundtripxyzfile, maxstringlength);
    declare_array_nozero(char, original_and_sites_tclfile, maxstringlength);
    declare_array_nozero(char, roundtrip_and_sites_tclfile, maxstringlength);
    declare_array_nozero(char, original_tclfile, maxstringlength);
    declare_array_nozero(char, roundtrip_tclfile, maxstringlength);
    sprintf(originalxyzfilename, "%s.original.xyz", filebase);
    sprintf(originalxyzfile, "%s/%s.original.xyz", directory, filebase);
    sprintf(roundtripxyzfilename, "%s.roundtrip.xyz", filebase);
    sprintf(roundtripxyzfile, "%s/%s.roundtrip.xyz", directory, filebase);
    sprintf(original_and_sites_tclfile, "%s/%s.original.sites.tcl", directory, filebase);
    sprintf(roundtrip_and_sites_tclfile, "%s/%s.roundtrip.sites.tcl", directory, filebase);
    sprintf(original_tclfile, "%s/%s.original.tcl", directory, filebase);
    sprintf(roundtrip_tclfile, "%s/%s.roundtrip.tcl", directory, filebase);
    FILE *originalxyzoutp, *roundtripxyzoutp, *originaltcloutp, *roundtriptcloutp, *originalandsitestcloutp, *roundtripandsitestcloutp;
    
    correctiontype[myaminoacid.Hindex]=1;
    correctiontype[myaminoacid.Cindex]=-1;
    correctiontype[myaminoacid.Oindex]=1;
    correctiontype[myaminoacid.Nindex]=-1;
    for(i=0;i<myaminoacid.N_corrections;i++){
        correctiontype[myaminoacid.corrected_atom_index[i]]=1;
        correctiontype[myaminoacid.correcting_atom_index[i]]=-1;
    }
    for(i=0;i<myresidue.Natoms;i++){
        if(myresidue.atomfound[i]==1){
            atomcounter[i]=Natoms;
            Natoms++;
            if(myaminoacid.element[i]==1){
                Hcount++;
                Hrmsd+=norm(subtract_double_triple(myresidue.atomposition[i], predictedatoms[i]));
            }
            else heavyrmsd+=norm(subtract_double_triple(myresidue.atomposition[i], predictedatoms[i]));
            if(myaminoacid.element[i]==1) heavycode=0;
            else heavycode=1;
            foundmatrix[myaminoacid.sitemap[i]][correctiontype[i]+1][heavycode]=1;
        }
    }
    
    originaltcloutp=fopen(original_tclfile, "w");
    roundtriptcloutp=fopen(roundtrip_tclfile, "w");
    originalandsitestcloutp=fopen(original_and_sites_tclfile, "w");
    roundtripandsitestcloutp=fopen(roundtrip_and_sites_tclfile, "w");
    
    fprintf(originaltcloutp, "mol new %s\n", originalxyzfilename);
    fprintf(roundtriptcloutp, "mol new %s\n", roundtripxyzfilename);
    fprintf(originalandsitestcloutp, "mol new %s\n", originalxyzfilename);
    fprintf(roundtripandsitestcloutp, "mol new %s\n", roundtripxyzfilename);
    
    setup_tcl(originaltcloutp, foundmatrix, doleft, doright, doCterminus);
    setup_tcl(roundtriptcloutp, foundmatrix, doleft, doright, doCterminus);
    setup_tcl(originalandsitestcloutp, foundmatrix, doleft, doright, doCterminus);
    setup_tcl(roundtripandsitestcloutp, foundmatrix, doleft, doright, doCterminus);
    
    for(i=0;i<myresidue.Natoms;i++){
        if(myresidue.atomfound[i]==1){
            for(j=0;j<myaminoacid.nconect[i];j++){
                if(myaminoacid.conect[i][j]<i){
                    if(myresidue.atomfound[myaminoacid.conect[i][j]]==1){
                        print_bond(originaltcloutp, atomcounter[i], atomcounter[myaminoacid.conect[i][j]]);
                        print_bond(roundtriptcloutp, atomcounter[i], atomcounter[myaminoacid.conect[i][j]]);
                        print_bond(originalandsitestcloutp, atomcounter[i], atomcounter[myaminoacid.conect[i][j]]);
                        print_bond(roundtripandsitestcloutp, atomcounter[i], atomcounter[myaminoacid.conect[i][j]]);
                    }
                }
            }
            if(doleft==1){
                if(i==myaminoacid.Nindex){
                    
                    //  N-lastC bond
                    
                    print_bond(originaltcloutp, atomcounter[i], Natoms);
                    print_bond(roundtriptcloutp, atomcounter[i], Natoms);
                    print_bond(originalandsitestcloutp, atomcounter[i], Natoms);
                    print_bond(roundtripandsitestcloutp, atomcounter[i], Natoms);
                }
                if(doright==1){
                    if(i==myaminoacid.Cindex){
                        
                        //  C-nextN bond
                        
                        print_bond(originaltcloutp, atomcounter[i], Natoms+1);
                        print_bond(roundtriptcloutp, atomcounter[i], Natoms+1);
                        print_bond(originalandsitestcloutp, atomcounter[i], Natoms+1);
                        print_bond(roundtripandsitestcloutp, atomcounter[i], Natoms+1);
                    }
                }
            }
            else if(doright==1){
                if(i==myaminoacid.Cindex){
                    //  C-nextN bond
                    print_bond(originaltcloutp, atomcounter[i], Natoms);
                    print_bond(roundtriptcloutp, atomcounter[i], Natoms);
                    print_bond(originalandsitestcloutp, atomcounter[i], Natoms);
                    print_bond(roundtripandsitestcloutp, atomcounter[i], Natoms);
                }
            }
        }
    }
    
    print_sites(originalandsitestcloutp, mycgresidue);
    print_sites(roundtripandsitestcloutp, mycgresidue);
    
    fclose(originaltcloutp);
    fclose(roundtriptcloutp);
    fclose(originalandsitestcloutp);
    fclose(roundtripandsitestcloutp);
    
    originalxyzoutp=fopen(originalxyzfile, "w");
    roundtripxyzoutp=fopen(roundtripxyzfile, "w");
    
    fprintf(originalxyzoutp, "%i\n%s\t%i\t%s\t%.8f\t%.8f\n", Natoms+doleft+doright, myresidue.chainidchar, myresidue.resseq, myaminoacid.resname, sqrt(heavyrmsd/(Natoms-Hcount)), sqrt(Hrmsd/Hcount));
    fprintf(roundtripxyzoutp, "%i\n%s\t%i\t%s\t%.8f\t%.8f\n", Natoms+doleft+doright, myresidue.chainidchar, myresidue.resseq, myaminoacid.resname, sqrt(heavyrmsd/(Natoms-Hcount)), sqrt(Hrmsd/Hcount));
    for(i=0;i<myresidue.Natoms;i++){
        if(myresidue.atomfound[i]==1){
            
            if((mycgresidue.Cterminussite>-1)&&(myaminoacid.Cterminusindex[i]>-1)){
                print_position(originalxyzoutp, "X", myresidue.atomposition[i]);
                print_position(roundtripxyzoutp, "X", predictedatoms[i]);
            }
            else if(myaminoacid.sitemap[i]==0){
                if(correctiontype[i]==1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "A", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "A", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "B", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "B", predictedatoms[i]);
                    }
                }
                else if(correctiontype[i]==-1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "D", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "D", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "E", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "E", predictedatoms[i]);
                    }
                }
                else{
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "F", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "F", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "G", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "G", predictedatoms[i]);
                    }
                }
            }
            else if(myaminoacid.sitemap[i]==1){
                if(correctiontype[i]==1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "I", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "I", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "J", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "J", predictedatoms[i]);
                    }
                }
                else if(correctiontype[i]==-1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "K", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "K", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "L", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "L", predictedatoms[i]);
                    }
                }
                else{
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "M", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "M", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "P", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "P", predictedatoms[i]);
                    }
                }
            }
            else if(myaminoacid.sitemap[i]==2){
                if(correctiontype[i]==1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "Q", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "Q", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "R", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "R", predictedatoms[i]);
                    }
                }
                else if(correctiontype[i]==-1){
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "T", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "T", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "U", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "U", predictedatoms[i]);
                    }
                }
                else{
                    if(myaminoacid.element[i]==1){
                        print_position(originalxyzoutp, "V", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "V", predictedatoms[i]);
                    }
                    else{
                        print_position(originalxyzoutp, "W", myresidue.atomposition[i]);
                        print_position(roundtripxyzoutp, "W", predictedatoms[i]);
                    }
                }
            }
            else{
                my_exit("Too many sites in output_residue_molecular_files\n");
            }
        }
    }
    if(doleft==1){
        print_position(originalxyzoutp, "X", lastCpos);
        print_position(roundtripxyzoutp, "X", predictedatoms[myresidue.Natoms]);
    }
    if(doright==1){
        print_position(originalxyzoutp, "X", nextNpos);
        print_position(roundtripxyzoutp, "X", predictedatoms[myresidue.Natoms+1]);
    }
    fclose(originalxyzoutp);
    fclose(roundtripxyzoutp);
    
    free(originalxyzfilename);
    free(originalxyzfile);
    free(roundtripxyzfilename);
    free(roundtripxyzfile);
    free(original_and_sites_tclfile);
    free(roundtrip_and_sites_tclfile);
    free(original_tclfile);
    free(roundtrip_tclfile);
    free(correctiontype);
    free(atomcounter);
    free_3d_tensor(foundmatrix, maxsitesperresidue, 3);
}

void print_position(FILE *outp, char *element, double_triple position){
    fprintf(outp, "%s\t%.8f\t%.8f\t%.8f\n", element, position.x, position.y, position.z);
}

void print_bond(FILE *outp, int a, int b){
    fprintf(outp, "topo addbond %i %i\n", a, b);
}

void setup_tcl(FILE *outp, int ***foundmatrix, int doleft, int doright, int doCterminus){
    fprintf(outp, "topo clearbonds\n");
    fprintf(outp, "set sel [atomselect top \" name A\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name B\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name D\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name E\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name F\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name G\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name I\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name J\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name K\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name L\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name M\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name P\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name Q\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name R\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name T\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name U\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name V\"]\n");
    fprintf(outp, "$sel set radius %f\n", Hradius);
    fprintf(outp, "set sel [atomselect top \" name W\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    fprintf(outp, "set sel [atomselect top \" name X\"]\n");
    fprintf(outp, "$sel set radius %f\n", heavyradius);
    if(foundmatrix[0][2][0]==1) fprintf(outp, "color Name A %s\n", color0ed);
    if(foundmatrix[0][2][1]==1) fprintf(outp, "color Name B %s\n", color0ed);
    if(foundmatrix[0][0][0]==1) fprintf(outp, "color Name D %s\n", color0ing);
    if(foundmatrix[0][0][1]==1) fprintf(outp, "color Name E %s\n", color0ing);
    if(foundmatrix[0][1][0]==1) fprintf(outp, "color Name F %s\n", color0);
    if(foundmatrix[0][1][1]==1) fprintf(outp, "color Name G %s\n", color0);
    if(foundmatrix[1][2][0]==1) fprintf(outp, "color Name I %s\n", color1ed);
    if(foundmatrix[1][2][1]==1) fprintf(outp, "color Name J %s\n", color1ed);
    if(foundmatrix[1][0][0]==1) fprintf(outp, "color Name K %s\n", color1ing);
    if(foundmatrix[1][0][1]==1) fprintf(outp, "color Name L %s\n", color1ing);
    if(foundmatrix[1][1][0]==1) fprintf(outp, "color Name M %s\n", color1);
    if(foundmatrix[1][1][1]==1) fprintf(outp, "color Name P %s\n", color1);
    if(foundmatrix[2][2][0]==1) fprintf(outp, "color Name Q %s\n", color2ed);
    if(foundmatrix[2][2][1]==1) fprintf(outp, "color Name R %s\n", color2ed);
    if(foundmatrix[2][0][0]==1) fprintf(outp, "color Name T %s\n", color2ing);
    if(foundmatrix[2][0][1]==1) fprintf(outp, "color Name U %s\n", color2ing);
    if(foundmatrix[2][1][0]==1) fprintf(outp, "color Name V %s\n", color2);
    if(foundmatrix[2][1][1]==1) fprintf(outp, "color Name W %s\n", color2);
    if((doleft==1)||(doright==1)){
        fprintf(outp, "color Name X %s\n", colorneighbor);
    }
    fprintf(outp, "proc vmd_draw_arrow {mol start end color} {\n");
    fprintf(outp, "    # an arrow is made of a cylinder and a cone\n");
    fprintf(outp, "    set middle [vecadd $start [vecscale %f [vecsub $end $start]]]\n", stemfraction);
    fprintf(outp, "    graphics $mol color $color\n");
    fprintf(outp, "    graphics $mol cylinder $start $middle radius %f\n", cylinderradius);
    fprintf(outp, "    graphics $mol cone $middle $end radius %f\n", coneradius);
    fprintf(outp, "}\n");
}

void print_sites(FILE *outp, cgresidue mycgresidue){
    double_triple start, end;
    int i;
    declare_array_nozero(char, color, maxstringlength);
    for(i=0;i<mycgresidue.Nsites;i++){
        if(i==mycgresidue.Cterminussite) strcpy(color, colorneighbor);
        else{
            if(i==0) strcpy(color, color0);
            else if(i==1) strcpy(color, color1);
            else if(i==2) strcpy(color, color2);
            else my_exit("Too many sites in print_sites");
        }
        start=mycgresidue.coord[i].pos;
        end=add_double_triple(mycgresidue.coord[i].pos, scalar_multiply_double_triple(mycgresidue.coord[i].ex, arrowlength));
        fprintf(outp, "draw arrow {%f %f %f} {%f %f %f} %s\n", start.x, start.y, start.z, end.x, end.y, end.z, color);
        start=mycgresidue.coord[i].pos;
        end=add_double_triple(mycgresidue.coord[i].pos, scalar_multiply_double_triple(mycgresidue.coord[i].ey, arrowlength));
        fprintf(outp, "draw arrow {%f %f %f} {%f %f %f} %s\n", start.x, start.y, start.z, end.x, end.y, end.z, color);
        start=mycgresidue.coord[i].pos;
        end=add_double_triple(mycgresidue.coord[i].pos, scalar_multiply_double_triple(mycgresidue.coord[i].ez, arrowlength));
        fprintf(outp, "draw arrow {%f %f %f} {%f %f %f} %s\n", start.x, start.y, start.z, end.x, end.y, end.z, color);
    }
}
