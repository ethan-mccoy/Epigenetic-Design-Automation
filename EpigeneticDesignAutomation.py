from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import GEOparse
from dataclasses import dataclass
from typing import List
import copy


@dataclass(frozen=True)
class GEOData:
    """
    TODO: Use this class to adapt GEOProcessor for function/model pattern and simplify the code
    """
    is_dataset = bool # True if GDS, False if GSM
    id = str # GSM or GDS ID
    table = pd.DataFrame # Table of data
    platform = str # Platform ID
    standardized_table = pd.DataFrame # Table of data with platform probe IDs converted to gene IDs
    filtered_by = tuple(str, str) # (column, value) tuple of column and value to filter by

class GEOProcessor:
    """
    Functionality: A helper class for TargetFinder that processes GEO data
    inputs: 
    - data_directory: str (directory to store GEO data)
    - gds_id: str (GEO DataSet ID)
    - indicator_col: str (column to filter by)
    - indicator_value: str (value to filter by)
    - gsm_id: str (GEO Sample ID)
    - genes_of_interest: list of genes of interest
    - gds: GEOparse object (GEO DataSet)
    """

    def __init__(self, data_directory="./geo_data"):
        self.data_directory = data_directory

    def get_filtered_dataset(self, gds_id, indicator_col, indicator_value):
        """
        Retrieves and filters a GEO DataSet (GDS). Commonly used to filter for control samples.
        Output: pandas dataframe (GDS dataset with only the filtered samples)
        """
        gds = GEOparse.get_GEO(geo = gds_id, destdir=self.data_directory, silent = True)
        std_gds_df = self.geo_with_gene_id(gds) 
        filtered_gds_df = self.filter_GDS_by_indicator(gds, indicator_col, indicator_value)
        filtered_gds_df = filtered_gds_df.merge(std_gds_df, on=['ID_REF', 'IDENTIFIER'], how='left')
        return filtered_gds_df
    
    def get_sample_dataset(self, gsm_id):
        """
        Retrieves a GEO sample dataset.
        Output: pandas dataframe (GSM dataset with gene IDs and expression values)
        """
        gsm = GEOparse.get_GEO(geo = gsm_id, destdir=self.data_directory, silent = True)
        gsm_df = self.geo_with_gene_id(gsm)
        return gsm_df

    def geo_with_gene_id(self, geo):
        """
        Parses epigenetic data for standardization; converts probe IDs to gene IDs.
        Output: pandas dataframe of GEO table with gene IDs and expression values
        """
        #TODO: Make it smarter, remove all the if elses. Find a way to identify the title of the gene name column in each platform
        # Maybe possible with GPT3.5 API but could overcomplicate things and cause bad outputs
        platform_column = 'platform_id' if 'platform_id' in geo.metadata else 'platform'
        platform_df = GEOparse.get_GEO(geo = geo.metadata[platform_column][0], destdir="./geo_data", silent = True).table
        geo_df = geo.table
        if ('Gene Symbol' in platform_df.columns):
            geo_df['IDENTIFIER'] = geo_df.merge(platform_df, left_on='ID_REF', right_on='ID')['Gene Symbol']
        elif ('GENE_SYMBOL' in platform_df.columns):
            geo_df['IDENTIFIER'] = geo_df.merge(platform_df, left_on='ID_REF', right_on='ID')['GENE_SYMBOL']
        elif ('ILMN_Gene' in platform_df.columns):
            geo_df['IDENTIFIER'] = geo_df.merge(platform_df, left_on='ID_REF', right_on='ID')['ILMN_Gene']
        geo_df['GB_ACC'] = geo_df.merge(platform_df, left_on='ID_REF', right_on='ID')['GB_ACC']
        return geo_df

    def filter_GDS_by_indicator(self, gds, indicator_col, indicator_value):
        """
        Filters a GDS dataset based on a specific indicator.
        Output: pandas DataFrame of filtered data
        """
        indicator = gds.columns
        epigenes = gds.table
        gsm_to_indicator = dict(zip(indicator.index, indicator[indicator_col]))
        indicated_cols = ['ID_REF', 'IDENTIFIER']
        for col, label in gsm_to_indicator.items():
            if label == indicator_value:
                indicated_cols.append(col)
        indicated_df = epigenes[indicated_cols]
        return indicated_df
    



class TargetFinder:
    """
    A class that helps find genes with different methylation levels between populations
    Input: 
    - target_gsm_df: pandas dataframe of target GSM
    - control_gds_df: pandas dataframe of control GDS
    - genes_of_interest: list of genes of interest
    Output: dictionary of target genes to methylation status
    """

    def __init__(self):
        pass

    def run(self, target_gsm_df, control_gds_df, genes_of_interest):
        target_gene_to_methylation = {}
        for gene_id in genes_of_interest:
            target_gene_to_methylation[gene_id] = self.assess_gene_methylation(target_gsm_df, control_gds_df, gene_id)

        return target_gene_to_methylation
    
    def assess_gene_methylation(self, target_gsm_df, control_gds_df, gene_id):
        """
        Assesses methylation status of a single gene. 
        Input: gene_id: str (gene ID)
        Output: str (methylation status of the gene, 'hypo' or 'hyper')
        """
        target_gene = target_gsm_df[target_gsm_df['IDENTIFIER'] == gene_id]
        control_gene = control_gds_df[control_gds_df['IDENTIFIER'] == gene_id]

        for i in range(len(target_gene)):
            control_mean = control_gene.iloc[i, 2:].mean()
            target_value = target_gene['VALUE'].iloc[i]
            percent_difference = abs((target_value - control_mean) / control_mean) * 100
            if percent_difference > 10:
                return 'hypo' if target_value < control_mean else 'hyper'
        return None
    

@dataclass(frozen=True)
class Gene:
    """TODO: Use this class to adapt TargetFinder for function/model pattern"""
    refseq_id: str
    gene_id: str
    seq = str
    promoter_seq = str
    pam_to_protospacer = dict

class SequenceProcessor:
    """A class that helps process DNA sequences"""

    def __init__(self):
        pass

    def get_promoter_sequence(self, refseq_id):
        """
        Attempts to find the promoter sequence of a gene from NCBI given its refseq id
        Input: str (refseq id of a gene)
        Output: str (promoter sequence of the gene)
        """
        try:
            handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            promoter_end_pos = 300 # Default if no CDS listed
            for feature in record.features:
                if feature.type == "CDS":
                    promoter_end_pos = feature.location.start
            seq = record.seq
            promoter_seq = seq[0 : promoter_end_pos]
            return promoter_seq
        except Exception as e:
            print(f"Error fetching sequence for {refseq_id}: {str(e)}")
            return None

    def find_protospacers(self, promoter_seq):
        """
        Searches the forward strand of the promoter sequence for PAM sites with protospacers that have G at position 1 and A/T at position 17
        If none are found, finds protospacers with G at position 20.
        Input: 
        - promoter_seq: str (promoter sequence of a gene)
        Returns: dictionary of PAM sites and their protospacers
        """
        pam_to_protospacer = {}
        for i in range(20, len(promoter_seq) - 2):
            if promoter_seq[i+2:i+3] == 'G' and promoter_seq[i+1:i+2] == 'G':
                protospacer = promoter_seq[i-20:i]
                protospacer = str(protospacer)  # Removes seq object
                if protospacer[0] == 'G' and (protospacer[16] == 'A' or protospacer[16] == 'T'):
                    pam_to_protospacer[i] = protospacer.lower()

        if not pam_to_protospacer:
            for i in range(20, len(promoter_seq) - 2):
                if promoter_seq[i+2:i+3] == 'G' and promoter_seq[i+1:i+2] == 'G':
                    protospacer = promoter_seq[i-20:i]
                    protospacer = str(protospacer)  # Removes seq object
                    if protospacer[18] == 'G': 
                        pam_to_protospacer[i] = protospacer.lower()

        return pam_to_protospacer


@dataclass()
class EpigeneticConstruct:
    # Plasmid backbone info
    plasmid_name : str 
    plasmid_methylation_type : str
    contains_scaffold: bool
    downstream_U6_sites: list        # Downstream of U6/T7 promoters / gRNA scaffold if present, but upstream of dCas9
    upstream_promoter_sites: list     # Site for extra promoter control of Cas9
    # Target gene info
    target_gene : str 
    target_gene_forward_oligo : str
    target_gene_reverse_oligo : str
    
class OligoDesigner:
    """
    Functionality: Designs oligos for epigenetic constructs
    Inputs:
    - genes_to_methylation: dictionary of genes to methylation status
    - target_gsm_df: dataframe of target GSM
    Outputs: dictionary of genes to list of constructs
    """
    def __init__(self):
        self.sequence_processor = SequenceProcessor()

        Tet1 = EpigeneticConstruct('pdCas9-Tet1-CD', 'hypo', False, ['Acc65I', 'KpnI'], ['AarI', 'AgeI'], '', '', '')
        p300 = EpigeneticConstruct('pcDNA-dCas9-p300 Core', 'hypo', False, ['ClaI', 'BspDI'], ['SacII'], '', '', '')
        DNMT3A = EpigeneticConstruct('pdCas9-DNMT3A-EGFP', 'hyper', True, ['XbaI', 'Acc65I', 'KpnI'], ['AarI', 'AgeI'], '', '', '')
        MQ1 = EpigeneticConstruct('pcDNA3.1-dCas9-MQ1(Q147L)-EGFP', 'hyper', False, ['EcoRV', 'NotI', 'XbaI'], [], '', '', '')
        self.constructs = [Tet1, p300, DNMT3A, MQ1]

        #TODO: find online database of common restriction sites, add as txt file, read in as dict
        self.restriction_sites = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'NotI': 'GCGGCCGC',
            'XhoI': 'CTCGAG',
            'EcoRV': 'GATATC',
            'XbaI': 'TCTAGA',
            'Acc65I': 'GGTACC',
            'KpnI': 'GGTACC',
            'AarI': 'CACCTGC',
            'AgeI': 'ACCGGT',
            'ClaI': 'ATCGAT',
            'BspDI': 'ATCGAT',
            'SacII': 'CCGCGG',
        }

    def run(self, genes_to_methylation, target_gsm_df):
        """
        Design oligos based on methylation status
        Inputs:
        - genes_to_methylation: dictionary of genes to methylation status
        - target_gsm_df: dataframe of target GSM
        Outputs: dictionary of genes to list of constructs
        """
        target_genes = genes_to_methylation.keys()
        self.target_gsm_df = target_gsm_df
        id_to_seq = self.find_promoter_sequences(target_genes)
        id_to_protospacer = self.find_protospacers(id_to_seq)
        gene_to_constructs = self.design_oligos(self.constructs, id_to_protospacer, genes_to_methylation)
        return gene_to_constructs

    def find_protospacers(self, id_to_seq, chopchop = False):
        """
        Finds potential protospacers in sequences
        Inputs: 
        - id_to_seq: dictionary of gene IDs to promoter sequences
        - chopchop: boolean indicating whether to use chopchop to find protospacers
        Outputs: dictionary of gene IDs to protospacers
        """
        # Protospacers with G at 0 and A/T at 16
        if chopchop == False:
            id_to_protospacers = {}
            for id, seq in id_to_seq.items():
                protospacers = self.sequence_processor.find_protospacers(seq)
                id_to_protospacers[id] = protospacers
            return id_to_protospacers
        elif chopchop == True:
            id_to_protospacers = {}
            #for id, seq in id_to_seq.items():
                #TODO: Set up chopchop cmd line tool, many dependences

    def find_promoter_sequences(self, target_genes):
        """
        Finds promoter sequences for target genes
        Inputs: 
        - target_genes: list of target genes
        Outputs: dictionary of gene IDs to promoter sequences
        """
        id_to_seq = {}
        for gene_id in target_genes:
            gene_df = self.target_gsm_df[self.target_gsm_df['IDENTIFIER'] == gene_id]
            if gene_df.empty:
                print(f"Error: Gene ID {gene_id} not found in the dataframe.")
                continue

            if ('DETECTION P-VALUE' in gene_df.columns):
                present_transcript_id = gene_df.sort_values(by='DETECTION P-VALUE').iloc[0:1]['GB_ACC'].iloc[0]
            else:
                present_transcript_id = gene_df.sort_values(by='Detection Pval').iloc[0:1]['GB_ACC'].iloc[0]

            promoter_seq = self.sequence_processor.get_promoter_sequence(present_transcript_id)
            id_to_seq[gene_id] = promoter_seq
        return id_to_seq

    def design_oligos(self, constructs, id_to_protospacer, genes_to_methylation): 
        """
        Designs oligos for constructs
        Inputs: 
        - constructs: list of EpigeneticConstruct objects
        - id_to_protospacer: dictionary of gene IDs to protospacers
        - genes_to_methylation: dictionary of genes to methylation status
        Outputs: dictionary of genes to list of constructs
        """
        gene_to_constructs = {}
        for gene, methylation in genes_to_methylation.items():
            for construct in constructs:
                if construct.plasmid_methylation_type == methylation:
                    gRNA_scaffold = ''
                    if construct.contains_scaffold == False:
                        gRNA_scaffold = 'gtttaagagctatgctggaaacagcatagcaagtttaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgcttttttt' 
                    gene_construct = copy.deepcopy(construct)
                    gene_construct.target_gene = gene

                    if gene in id_to_protospacer:
                        protospacer = id_to_protospacer[gene]
                    else:
                        print(f"Warning: Gene {gene} not found in id_to_protospacer. Skipping.")
                        continue

                    restriction_site = gene_construct.downstream_U6_sites[0]
                    restriction_seq = self.restriction_sites.get(restriction_site)
                    restriction_seq_overhang = 'ataat' + restriction_seq

                    gene_construct.target_gene_forward_oligo = str(restriction_seq_overhang) + str(protospacer) + str(gRNA_scaffold)
                    gene_construct.target_gene_reverse_oligo = restriction_seq_overhang
                    if gene not in gene_to_constructs:
                        gene_to_constructs[gene] = []
                    gene_to_constructs[gene].append(gene_construct)
        return gene_to_constructs
    

class EpigeneticDesignAutomation:
    """
    Functionality: Automates the design of epigenetic constructs
    Inputs:
    - email: str (email for Entrez)
    - target_gsm_id: str (GSM ID of the target sample)
    - reference_geo_id: str (GSM or GDS ID of the reference sample)
    - genes_of_interest: list or dictionary of genes of interest
    - indicator_col: str (column to filter reference data by)
    - indicator_value: str (value to filter reference data by)
    Outputs: dictionary of genes to list of constructs
    """

    def __init__(self, email):
        Entrez.email = email
        self.geo_processor = GEOProcessor()
        self.target_finder = TargetFinder()
        self.oligo_designer = OligoDesigner()

    def run(self, target_gsm_id, reference_geo_id, genes_of_interest, indicator_col = None, indicator_value = None):
        target_gsm = self.get_target_gsm(target_gsm_id)
        ref_values = self.get_reference_methylation(reference_geo_id, indicator_col, indicator_value)
        genes_to_methylation = self.get_target_genes(target_gsm, ref_values, genes_of_interest)
        gene_to_constructs = self.oligo_designer.run(genes_to_methylation, target_gsm)
        
        return gene_to_constructs

    def get_target_gsm(self, target_gsm_id):
        """
        Retrieves the dataset for a specified target GSM ID.
        Output: pandas dataframe of the target GSM
        """
        try:
            target_gsm = self.geo_processor.get_sample_dataset(target_gsm_id)
            return target_gsm
        except Exception as e:
            raise ValueError(f"Not a valid GSM id: {target_gsm_id}. Error: {e}")
        
    def get_reference_methylation(self, reference_geo_id, indicator_col = None, indicator_value = None):
        """
        Retrieves and processes reference methylation data.
        Output: pandas dataframe of the reference data
        """
        if reference_geo_id.startswith("GDS"):
            return self.handle_gds_reference(reference_geo_id, indicator_col = indicator_col, indicator_value = indicator_value)
        elif reference_geo_id.startswith("GSM"):
            return self.handle_gsm_reference(reference_geo_id)
        else:
            raise ValueError("Invalid GEO ID. Please enter a valid GSM or GDS ID.")
        
    def handle_gds_reference(self, gds_id, indicator_col, indicator_value):
        """
        Handles retrieval and processing of GDS data.
        Output: pandas dataframe of the GDS data
        """
        filtered_gds = self.geo_processor.get_filtered_dataset(gds_id, indicator_col, indicator_value)
        ref_values = filtered_gds[['ID_REF', 'IDENTIFIER']].copy()
        ref_values['VALUE'] = filtered_gds.filter(regex='^GSM').mean(axis=1)
        return ref_values
        
    def handle_gsm_reference(self, reference_gsm_id):
        """
        Handles retrieval and processing of GSM data.
        Output: pandas dataframe of the GSM data
        """
        reference_gsm = self.geo_processor.get_sample_dataset(reference_gsm_id)
        reference_gsm = reference_gsm[reference_gsm['VALUE'].notna()]
        ref_values = reference_gsm[['ID_REF', 'IDENTIFIER', 'VALUE']].copy()
        return ref_values

    def get_target_genes(self, target_gsm, ref_values, genes_of_interest):
        """
        Creates a dictionary linking genes of interest to their methylation status (hypo/hyper)
        Input: 
        - ref_values: pandas dataframe of the reference data
        Output: Dictionary of target genes to methylation status
        """
        if isinstance(genes_of_interest, dict):
            genes_to_methylation = genes_of_interest
        elif isinstance(genes_of_interest, list):
            genes_to_methylation = self.target_finder.run(target_gsm, ref_values, genes_of_interest)
        else:
            raise ValueError("genes_of_interest must be a list or a dictionary.")
        return genes_to_methylation