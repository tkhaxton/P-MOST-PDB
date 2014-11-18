
# define maxlinelength 400
# define maxstringlength 200
# define Naminoacids 20
# define maxatomsperresidue 50
# define maxsiteatoms 6
# define maxNchains 100
# define residuebuffer 1000
# define max_sidechain_corrections 60
# define max_tetrahedra_per_residue 5
# define max_overlap_pairs_per_residue 5
# define minimum_atom_separation 0.1
# define maxNCdistance2 9.0
# define maxbondfactor 1.5
# define minbondfactor 0.6
# define maxelement 17
# define maxatomtypes 100
# define maxbackboneatomtypes 30
# define maxsitesperresidue 3
# define extra_backbone_atom_types 3
# define max_rings_per_residue 2
# define number_amino_hydrogen_atom_types 15
# define HN2offset 1
# define HM1offset 2
# define HM2offset 3
# define HM3offset 4
# define GLYoffset 5
# define PROoffset 10
# define max_angles_per_residue 60
# define maxangledifference (M_PI/6.)
# define O_C_N_angle 122.5
# define CA_C_N_angle 116.5
# define C_N_H_angle 123.0
# define minimum_tetrahedral_angle (M_PI/6.)
# define Hradius 1.2
# define heavyradius 1.7
# define color0 "red"
# define color0ed "red2"
# define color0ing "red3"
# define color1 "blue"
# define color1ed "blue2"
# define color1ing "blue3"
# define color2 "green"
# define color2ed "green2"
# define color2ing "green3"
# define colorneighbor "black"
# define cylinderradius 0.1
# define coneradius 0.2
# define stemfraction 0.8
# define arrowlength 1.0
# define max_aromatic_hydrogens 6
# define max_overlapping_atoms 1000
# define max_atoms_in_radial_set 5

# define distancesquared_error_threshold 1

typedef struct chaindata_tag{
    int seqbegin;
    int seqend;
	int dbseqend;
} chaindata;

typedef struct residuedata_tag{
    int residnumber;
    int chainid;
    int Natoms;             //  Total number of possible atoms for the residue
    int *atomfound;         //  0: not found; 1: found; 2: not found because it is a terminal atom type and the residue is not terminal
    double_triple *atomposition;
    char *chainidchar;
    int resseq;
    int Nterminusflag;
    char **serialchar;
    char **icodechar;
    char **occupancychar;
    char **tempfactorchar;
    char **elementchar;
    char **chargechar;
} residuedata;

typedef struct amino_acid_struct_tag{
	char *resname;
 	int Natoms;
	char **atomnames;
	char **input_atomnames;
    int Nsites;
    int *sitecode;          //  Code for how the site is defined; e.g. sitecode[i]=0 denotes triad based on siteatomindex[i][0], siteatomindex[i][1], siteatomindex[i][2]
    int *sitecount;
    int **siteatomindex;
    int *sitemap;           //  Denotes which site atom belongs to
    int *Cterminusindex;    //  Maps to -1 (default) or index among C terminus site atoms
    int NCterminusatoms;
    int *Cterminussiteatomindex;
    double_triple *avpos;
    double_triple *varpos;
    double_triple *msd_residuespecific;
    double *fourthmoment_residuespecific;
    int *msdcount;
    int *count;
    int Oindex;
    int Nindex;
    int Cindex;
    int CAindex;
    int Hindex;
    int aminoHindex;
    int OXTindex;
    //int HXTindex;
    int N_corrections;
    int *correction_index;  //  Maps from index within amino acid to overall index used for array of correction structures; -1 denotes no correction
    int *corrected_atom_index;
    int *correcting_atom_index;
    int *atom_correction_index; //  Maps from atom to overall index used for array of correction structures
    int *atom_correction_index_withinaa;
    int *nconect;
    int **conect;
    int *element;
    int N_indistinguishable_sets;
    int *indistinguishable_set_size;
    int **indistinguishable_atoms;
    double *first_dividing_angle;
    int N_aromatic_hydrogens;
    int *aromatic_hydrogen_list;
    int N_radial_sets;
    int *radial_set_size;
    int **radial_atoms;
    int *radial_count;
    double_double *radial_average;
} amino_acid_struct;

typedef struct cgcoord_tag{
    double_triple pos;
    double_triple ex;
    double_triple ey;
    double_triple ez;
} cgcoord;

typedef struct cgresidue_tag{
    int Nsites;
    int *sitedefined;
    cgcoord *coord;
    int Cterminussite;      //  -1 by default; otherwise maps to C-terminus site (OXT, O, HXT)
} cgresidue;

typedef struct correction_structure_tag{
    double *predicted_pos_av;
    double *predicted_pos_var;
    double *pos_difference_av;
    double **covariance_matrix_cross;
    double **covariance_matrix_self;
    double **slope;
    double *intercept;
    int count;
} correction_structure;

typedef struct atomlookup_tag{
    char *filename;
    int residue;
    int atom;
} atomlookup;

void create_amino_acid_struct(char *directory, char *site_director_atoms_file, char *site_atoms_file, char *corrections_parameter_file, char *elements_file, char *Cterminus_atoms_file, char *indistinguishable_atoms_file, char *radial_atoms_file, amino_acid_struct **pamino_acid_list, int *pN_sidechain_corrections, char ***psidechain_correction_identifiers, char ***pbackbone_atom_types, int *pN_backbone_atom_types, int *pbackbone_Oindex, int *pbackbone_Nindex, int *pbackbone_Hindex, int *pbackbone_Cindex, char ***pCterminusatomnames, char *aromatic_hydrogens_file);
void change_atom_names_in_amino_acid_struct(amino_acid_struct *amino_acid_list, char *change_atom_names_file, int *pN_backbone_atom_types, char **backbone_atom_types);
void initialize_correction_structure(correction_structure *pmystructure);
int read_pdb_to_residuearray(char *filename, amino_acid_struct *amino_acid_list, residuedata **presiduearray, int *pNresidues, int record_all_data, int dbref);
void calculate_average_radial_distances(amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray);
void output_error_filename(char *file, char *base, int type, int *pfirst);
void reassign_Nterminus_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite, int record_all_data);
int assign_atomfound_terminii(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list);
void discard_overlapping_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite, atomlookup *overlapping_atoms_list, int *pN_overlapping_atoms);
void setup_bond_length_structure(char *covalent_bond_lengths_file, int ***pfoundbond, double ***pmaxbondlength, double ***pminbondlength);
void check_covalent_bond_lengths(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int **foundbond, double **maxbondlength, double **minbondlength, char *output_file, char *pdb_file, int dowrite);
void setup_bond_angles_structure(char *bond_angles_file, amino_acid_struct *amino_acid_list, int **pNangles, double ***pangles, int ****pangle_index);
void check_bond_angles(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int *Nangles, double **angles, int ***angle_index, char *output_file, char *pdb_file, int dowrite);
void discard_bad_tetrahedra(residuedata *residuearray, int Nresidues, int *Ntetrahedra, int ***tetrahedraindices, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite);
void setup_rings_structure(char *rings_file, amino_acid_struct *amino_acid_list, int **pNrings, int ***pringsize, int ***pnumber_dependent_atoms, int ****patomindex);
void check_rings(residuedata *residuearray, int Nresidues, int *Nrings, int **ringsize, int **number_dependent_atoms, int ***atomindex, amino_acid_struct *amino_acid_list, char *output_file, char *pdb_file, int dowrite);
void print_indistinguishable_xyz_files(amino_acid_struct *amino_acid_list, residuedata *residuearray, int Nresidues, char *base, int *pfirst);
void resolve_indistinguishable_atoms(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int dowrite, char *output_file, char *pdb_file, int record_all_data);
void resolve_aromatic_hydrogens(residuedata *residuearray, int Nresidues, amino_acid_struct *amino_acid_list, int dowrite, char *output_file, char *pdb_file, int record_all_data);
void map_to_cgmodel(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue **pcgresiduearray, int docountsites);
void record_moments_relative_to_oriented_cg_sites(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray);
void calculate_Cterminus_bond_lengths_and_angles(amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray, double **Cterminus_avs, int *pCterminus_count);
void output_Cterminus_bond_lengths_and_angles(char *filename, double **Cterminus_avs, int Cterminus_count);
void free_cgresiduearray(cgresidue **pcgresiduearray, int Nresidues);
void free_residuearray(residuedata **presiduearray, int Nresidues, int record_all_data);
void add_header_to_xyz_files(char *base, amino_acid_struct *amino_acid_list);
void average_and_output_radial_positions(char *filename, amino_acid_struct *amino_acid_list);
void average_and_output_moments_relative_to_oriented_cg_sites(char *filename, amino_acid_struct *amino_acid_list, char **backbone_atom_types, int N_backbone_atom_types, char **Cterminusatomnames);
void calculate_covariance_separate(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array);
void calculate_covariance_sidechain(int Nresidues, residuedata *residuearray, amino_acid_struct *amino_acid_list, cgresidue *cgresiduearray, correction_structure *corrections);
void calculate_correction_from_covariance(correction_structure *pcarbonyl_correction, char *corrections_file);
void calculate_rmsd_separate_corrections(amino_acid_struct *amino_acid_list, residuedata *residuearray, cgresidue *cgresiduearray, int Nresidues, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array, correction_structure *sidechain_corrections, char *pdb_file, char *large_disagreement_file_specific, char *entry_average_file);
void output_rmsd_onlyspecific(char *filename, amino_acid_struct *amino_acid_list, int nrecords, char **backbone_atom_types, int N_backbone_atom_types, char **Cterminusatomnames);
void output_file_type_statistics(int *Ntype, char *filename);

void extract_piece_of_string(char *instring, char *outstring, int start, int end);
void extract_piece_of_string_lowercase(char *instring, char *outstring, int start, int end);
void print_error(int type);

void input_cg_coords(char *filename, amino_acid_struct *amino_acid_list, char **Cterminusatomnames);
void input_correction(correction_structure *pcorrection, char *corrections_file);
void create_molecular_files_separate_corrections(char *pdbfile, char *directory, char *filebase, amino_acid_struct *amino_acid_list, int Nresidues, residuedata *residuearray, cgresidue *cgresiduearray, correction_structure *carbonyl_correction_array, correction_structure *backboneH_correction_array, correction_structure *sidechain_corrections, int do_amino_acid_examples);

void output_residue_molecular_files(char *directory, char *filebase, amino_acid_struct myaminoacid, residuedata myresidue, cgresidue mycgresidue, double_triple *predictedatoms, double_triple lastCpos, double_triple nextNpos, int doleft, int doright, int doCterminus);
void print_position(FILE *outp, char *element, double_triple position);
void print_bond(FILE *outp, int a, int b);
void setup_tcl(FILE *outp, int ***foundmatrix, int doleft, int doright, int doCterminus);
void print_sites(FILE *outp, cgresidue mycgresidue);
int convert_entry_to_integer(char *entry_id);




