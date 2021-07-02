#include <KBaseFeatureValues.spec>
#include <KBaseAssembly.spec>

/*
This module is for storing microbial growth phenotype data, e.g., from
FEBA project or other ENIGMA experiments
*/

module KBaseRBTnSeq {

    /*
    A handle id from the Handle Service for a shock node.
    @id handle
    */
    typedef string handle_id;
    
    /*
    @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;


    /*
    @id ws KBaseRBTnSeq.RBTS_InputGenesTable
    */
    typedef string genes_table_ref;
    
    /*
    @id ws KBaseFile.SingleEndLibrary 
    */
    typedef string fastq_ref;
   
    /*
    A list of fastq_refs
    */
    typedef list<fastq_ref> fastqs;

    /*
    A header for a  column
    */
    typedef string col_header;
    
    
    /*
    The list of column headers
    */
    typedef list<col_header> col_list;

   

    /*
    file_type - KBaseRBTnSeq.RBTS_InputGenesTable, the name of the file type.
    input_genes_table - handle that allows to download file, and get info re. shock node, shock url,
    handle_type - the type of the handle. This should always be ‘shock’.
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    column_header_list - a list of the headers of the columns, the length of this 
        list should be the num of columns in the file. Currently: 
        <"locusId", "sysName", "type", "scaffoldId", "begin", "end", "strand", "name", "desc", "GC", "nTA">
        making a total of 11 columns.
    column_headers_str - a string; comma-separated column headers for the file
    num_lines - the number of lines in the file - keeps track of the general size
    related_genome_ref -  the genome which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
   
    @optional related_genome_ref 
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws num_lines as num_lines
    @metadata ws column_headers_str as column_headers_str
    @metadata ws related_genome_ref as related_genome_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    */

    typedef structure {
        string file_type;
        handle_id input_genes_table;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        col_list column_header_list;
        string column_headers_str;
        string num_lines;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
    } RBTS_InputGenesTable;


    /*
    file_type - KBaseRBTnSeq.RBTS_Model, the name of the file type.
    utc_created - the Coordinated Universal Time of creation
    standard_model_name - e.g. 'Sc_Tn5', or 'pKMW3_universal' 
    model_string - The model string: e.g.  nnnnnnGATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGACGGCCGGCCAGACCGGGGACTTATCAGCCAACCTGT
    past_end_string - The Past end string: e.g. TATGTGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTAATTCTTGAAGA
   
    @metadata ws utc_created as utc_created
    @metadata ws standard_model_name as standard_model_name
    @metadata ws model_string as model_string
    @metadata ws past_end_string as past_end_string
    */

    typedef structure {
        string file_type;
        string utc_created;
        string standard_model_name;
        string model_string;
        string past_end_string;
    
    } RBTS_Model;




    /*
    file_type - KBaseRBTnSeq.RBTS_PoolFile, the name of the file type.
    poolfile - handle that allows to download file, and get info re. shock node, shock url,
    handle_type - the type of the handle. This should always be ‘shock’.
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    column_header_list - a list of the headers of the columns, the length of this 
        list should be the num of columns in the file. Currently: 
        <"barcode", "rcbarcode", "nTot", "n", "scaffold", "strand", 
        "pos", "n2", "scaffold2", "strand2", "pos2", "nPastEnd">
        making a total of 12 columns.
    column_headers_str - a string; comma-separated column headers for the file
    num_lines - the number of lines in the file - keeps track of the general size
    model_used - The string of the model used to create this poolfile
    related_genes_table_ref - the genes table which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    fastqs_used - the fastqs which were used to create the poolfile
    fastqs_used_str - comma separated string with refs of fastqs used to create file. 
    description - A description given by the uploader as to what the
        pool file means.
    
    @optional fastqs_used
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws num_lines as num_lines
    @metadata ws model_used as model_used
    @metadata ws column_headers_str as column_headers_str
    @metadata ws related_genes_table_ref as related_genes_table_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws fastqs_used_str as fastqs_used_str
    @metadata ws description
    */
    typedef structure {
        string file_type;
        handle_id poolfile;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        col_list column_header_list;
        string column_headers_str;
        string num_lines;
        genes_table_ref related_genes_table_ref;
        string related_organism_scientific_name;
        string model_used;
        fastqs fastqs_used; 
        string fastqs_used_str; 
        string description;
    
    } RBTS_PoolFile;


  






    /*
    file_type KBaseRBTnSeq.RBTS_PoolCount
    handle_id will be poolcount file handle
    handle_type - the type of the handle. This should always be ‘shock’.
    column_header_list will be
        barcode, rcbarcode, scaffold, strand, pos, and an unknown number of columns
    column_headers_str - a string; comma-separated column headers for the file
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    set_name - the name of the set
    num_lines - the number of lines in the file - keeps track of the general size
    related_genes_table_ref - the genes table which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    fastqs_used - the fastqs used to create the poolcount file
    poolfile_ref - the ref for the poolfile used to create the poolcount file
    description - A description given by the uploader as to what the
        pool file means.
    
    @optional poolfile_ref fastqs_used
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws set_name
    @metadata ws column_headers_str as column_headers_str
    @metadata ws related_genes_table_ref as related_genes_table_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws fastqs_used_str as fastqs_used_str
    @metadata ws description
    @metadata ws num_lines
    */
    typedef structure {

        string file_type;
        handle_id poolcount;
        string handle_type;
        col_list column_header_list;
        string column_headers_str;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        string set_name;
        string num_lines;
        genes_table_ref related_genes_table_ref;
        string related_organism_scientific_name;
        fastqs fastqs_used; 
        string fastqs_used_str; 
        string poolfile_ref;
        string description;
    
    } RBTS_PoolCount;


    /*
    file_type KBaseRBTnSeq.RBTS_ExperimentsTable
    exps_file will be experiments file handle
    handle_type - the type of the handle. This should always be ‘shock’.
    column_header_list will have required parts:
        SetName, Index, Description, Date_pool_expt_started
        Not required headers have the following sometimes:
            Drop, Person, Mutant Library,
            gDNA plate, gDNA well, Sequenced at, Media, Growth Method,
            Group, Temperature, pH, Liquid v. solid, Aerobic_v_Anaerobic, Shaking,
            Condition_1, Concentration_1, Units_1, Condition_2, Concentration_2,
            Units_2, Timecourse, Timecourse Sample, Growth Plate ID, 
            Growth Plate wells, StartOK, EndOD, Total Generations
    column_headers_str - a string; comma-separated column headers for the file
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    num_lines - the number of lines in the file - keeps track of the general size
    related_genes_table_ref -  the genes_table which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    description - A description given by the uploader as to what the
        pool file means.

    @optional poolfile_ref
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws column_headers_str as column_headers_str
    @metadata ws related_genes_table_ref as related_genes_table_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws description
    @metadata ws num_lines
    */
    typedef structure {

        string file_type;
        handle_id expsfile;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        col_list column_header_list;
        string column_headers_str;
        string num_lines;
        genes_table_ref related_genes_table_ref;
        string related_organism_scientific_name;
        string poolfile_ref;
        string description;

    } RBTS_ExperimentsTable;












    typedef int bool;
        
    /*
    @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;
    
    /*
    @id ws KBaseGenomes.Contig
    */
    typedef string contig_ref;
    
    /*
    Replace with KBaseFile object?
    @id ws KBaseAssembly.SingleEndLibrary
    */
    typedef string reads_ref;
    
    /*
    TnSeq barcode model, showing how the barcodes / transposon
    were designed
    */
    typedef tuple <string read_template, string flanking_sequence> tnseq_model;
    /*
    A MappedReads object stores the mapping of reads to a genome.
    Unique and non-unique read positions are stored in arrays indexed using
    the contig index.  The last set of reads in each of these arrays
    corresponds to "past end" reads.  Handle should be replaced by
    a KBaseFile.Handle, but this is not registered yet.
    */
    typedef structure {
        genome_ref genome;
        reads_ref reads;
        tnseq_model model;
        KBaseAssembly.Handle mapped_reads_file;
        list<list<int>> unique_insert_pos_by_contig;
        list<list<int>> nonunique_insert_pos_by_contig;
    } MappedReads;
    
    /*
    @id ws KBaseRBTnSeq.MappedReads
    */
    typedef string mapped_reads_ref;
    
    /*
    @id ws KBaseBiochem.Media
    */
    typedef string media_ref;
    
    /*
    Reference to a Feature object of a genome in the workspace
    @id subws KBaseGenomes.Genome.features.[*].id
    */
    typedef string feature_ref;
    
    /*
    enum: insertion, deletion, substitution
    The latter is not strictly necessary, but convenient to avoid storing
    two separate events.
    */
    typedef string change_type;
    
    /*
      A Delta is a description of a single change to a strain.  A series
      of Deltas defines the transition from one Strain to another.  For
      sequenced insertions or substitutions, give the 0-indexed position
      on the contig where the insertion/substitution begins, in the +
      direction.  For sequenced deletions and substitutions, give the
      position and length.  The position of all Deltas should be
      calculated relative to the parent strain (derived_from_strain), so
      that the Deltas could be applied in any order.
    @optional change feature contig_index sequence position length
    */
    typedef structure {
        string description;
        change_type change;
        feature_ref feature;
        int contig_index;
        string sequence;
        int position;
        int length;
    } Delta;
    /*
    @id ws KBaseRBTnSeq.Strain
    */
    typedef string strain_ref;
    
    /*
      A Strain is a particular genetic variant of an organism.  Optionally,
      it may be:
        * derived from another Strain (e.g., as an engineered mutant)
        * sequenced
        * a wild-type example of a Genome
      If a strain is "wild type" it should have a non-null genome_ref and a
      null derived_from_strain.  If not wild type (or otherwise not
      characterized), genome_ref should be set to null.  (If genome_ref
      were not null, a large pool of mutant strains would have too many references
      to the genome in our current data model.)
    @optional description genome derived_from_strain deltas
    */
    typedef structure {
        string name;
        string description;
        genome_ref genome;
        strain_ref derived_from_strain;
        list<Delta> deltas;
    } Strain;
    /*
    @id ws KBaseRBTnSeq.Pool
    */
    typedef string pool_ref;
    
    /*
      A Pool is a collection of barcoded strains.  Barcodes, tags, etc should
      be stored as Deltas in each strain.
    @optional pool_hit_file pool_unhit_file pool_surprise_file
    */
    typedef structure {
        genome_ref genome;
        mapped_reads_ref mapped_reads;
        KBaseAssembly.Handle pool_file;
        KBaseAssembly.Handle pool_hit_file;
        KBaseAssembly.Handle pool_unhit_file;
        KBaseAssembly.Handle pool_surprise_file;
        list<Strain> strains;
        list<int> counts;
    } Pool;
    /*
     Reference to a compound object in a biochemistry
     @id subws KBaseBiochem.Biochemistry.compounds.[*].id
    */
    typedef string compound_ref;
    
    /*
      A Condition is something that is added to particular aliquots in
      a growth experiment, in addition to the media.  e.g., it may be a stress
      condition, or a nutrient.  Compound is needed if the condition is
      addition of a chemical in the KBase Biochemistry database.
    @optional concentration units compound
    */
    typedef structure {
        string name;
        float concentration;
        string units;
        compound_ref compound;
    } Condition;
    /*
    @id ws KBaseRBTnSeq.Condition
    */
    typedef string condition_ref;
    
    /*
      GrowthParameters describes all the conditions a particular aliquot
      was subjected to in an experiment
    @optional description gDNA_plate gDNA_well index media growth_method group temperature pH isLiquid isAerobic shaking growth_plate_id growth_plate_wells startOD endOD total_generations
    */
    typedef structure {
        string description;
        string gDNA_plate;
        string gDNA_well;
        string index;
        media_ref media;
        string growth_method;
        string group;
        float temperature;
        float pH;
        bool isLiquid;
        bool isAerobic;
        string shaking;
        string growth_plate_id;
        string growth_plate_wells;
        float startOD;
        float endOD;
        float total_generations;
    } GrowthParameters;
    /*
    @id ws KBaseRBTnSeq.GrowthParameters
    */
    typedef string growth_parameters_ref;
    
    /*
      A BarSeqExperiment is an experiment in which a pool is grown in
      several parallel aliquots (e.g., wells or tubes), each potentially
      treated with a different set of conditions (but usually it is just one
      Condition)
      BarSeqExperiment object represents only one such condition
    @optional person mutant_lib_name tnseq_pool
    */
    typedef structure {
        string name;
        string person;
        string mutant_lib_name;
        string start_date;
        string sequenced_at;
        GrowthParameters growth_parameters;
        list<Condition> conditions; /* usually one or two conditions */
        pool_ref tnseq_pool;
    } BarSeqExperiment;
    
    /*
    @id ws KBaseRBTnSeq.BarSeqExperiment
    */
    typedef string barseq_experiment_ref;
    
    /*
      Number of times a barcode (i.e. a strain) was detected by sequencing a pool at beginning (refernce state)
      and at the end of GrowthParameters, and a calculated log ratio of strain abundance relative to a starting
      condition.
      feature_index - Genome.features[index] optional, if strain_index is NA and log ratio corresponds to a gene, not a strain
      strain_index - index of a strain in Pool.strains list
      count_begin - number of instances of the strain identified from sequencing at the beginning of experiment
      count_end - at the end of experiment
      norm_log_ratio - normalized log ratio between count_end and count_begin
    */
    typedef tuple<int feature_index,int strain_index,int count_begin,int count_end,float norm_log_ratio> bar_seq_result;
    
    /*
     A single (i.e. one condition) BarSeq experiment
     experiment - describes the BarSeq experiment
     bar_seq_result - list of counts. Can be per gene or per strain.
    */
    typedef tuple<BarSeqExperiment experiment, list< bar_seq_result > results> bar_seq_exp;
    
    /*
      BarSeqExperimentResults stores the log ratios calculated from
      a BarSeqExperiment.  There is one log ratio per strain or per gene per
      GrowthParameters.
      May contain a number of BarSeqExperiment done on same species (i.e. 
      usually it would be dozens and even more than a hundred barseq experiments 
      performed under various conditions starting from a the same library).
      The raw data, including start/end mutant counts, is stored in:
         list<bar_seq_exp> experiments
      It does not require to have values for all genes/loci and the number of features does not have 
      to be the same for all experiments.
    
      features_by_experiments - is a 2D matrix that contains only normalized log ratios from experiments[*]['results' = 1][*]['norm_log_ratio' = 5 ]
      Therefore, there is a redundancy, we store log ratios twice. 'experiments' is used to store all raw data and allows quickly retrieve a single condition. 
      'features_by_experiments' is used to support visualization widgets and other methods that work with KBaseFeatureValues.FloatMatrix2D. 
      It also allows to quickly retrieve log ratios for all conditions per gene.
    */
    typedef structure {
        genome_ref genome;
        mapping<int feature_index, string feature_id> feature_index_to_id;
        list<bar_seq_exp> experiments;
        KBaseFeatureValues.FloatMatrix2D features_by_experiments;    
        mapping<string, int> col_to_index; /* column name of 'features_by_experiments' to index within 'experiments' mapping */
        mapping<string, int> row_to_index; /* row names (usually genes) to feature_index, i.e. reverse of 'feature_index_to_id' */
    } BarSeqExperimentResults;
    /*
    @id ws KBaseRBTnSeq.BarSeqExperimentResults
    */
    typedef string barseq_experiment_results_ref;
    /*
      Computes essential genes from a TnSeq pool
      Input: pool_ref - reference to a Pool
      Output: list of genes 
    */
    funcdef essential_genes(pool_ref) returns (list<feature_ref>) authentication required;
    /*
      Computes essential genes from a BarSeq experiment
      Input: barseq_experiment_results_ref - reference to a BarSeqExperimentResults object
      Output: list of genes 
    */
    funcdef essential_genes(barseq_experiment_results_ref) returns (list<feature_ref>) authentication required;
    /*
      Computes gene fitness within a TnSeq pool
      Input: tnseq_library_ref - reference to a TnSeqLibrary
      Output: list of genes with their fitness
    */
    funcdef gene_fitness(pool_ref) returns (list<tuple<feature_ref gene, float fitness>>) authentication required;
    /*
      Computes gene fitness from a BarSeq experiment
      Input: barseq_experiment_results_ref - reference to a BarSeqExperimentResults object
      Output: list of genes with their fitness
    */
    funcdef gene_fitness(barseq_experiment_results_ref) returns (list<tuple<feature_ref gene, float fitness>>) authentication required;

};


