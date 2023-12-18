from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import GEOparse
from dataclasses import dataclass
import copy

class GEOProcessor:
    """A helper class for TargetFinder that processes GEO data"""

    def __init__(self, data_directory="./geo_data"):
        self.data_directory = data_directory

    def get_filtered_dataset(self, gds_id, indicator_col, indicator_value):
        """
        Retrieves and filters a GEO DataSet (GDS).
        Input: GDS ID, column to filter, value to filter
        Returns: a pandas dataframe of the GDS dataset with only the filtered/control samples
        """
        gds = GEOparse.get_GEO(geo = gds_id, destdir=self.data_directory, silent = True)
        std_gds_df = self.geo_with_gene_id(gds) 
        filtered_gds_df = self.filter_GDS_by_indicator(gds, indicator_col, indicator_value)
        filtered_gds_df = filtered_gds_df.merge(std_gds_df, on=['ID_REF', 'IDENTIFIER'], how='left')
        return filtered_gds_df
    
    def get_sample_dataset(self, gsm_id):
        """
        Retrieves a GEO sample dataset.
        Input: GSM (Sample) ID
        Returns: a pandas dataframe of the GSM dataset with gene IDs and expression values
        """
        gsm = GEOparse.get_GEO(geo = gsm_id, destdir=self.data_directory, silent = True)
        gsm_df = self.geo_with_gene_id(gsm)
        return gsm_df

    def geo_with_gene_id(self, geo):
        """
        Parses epigenetic data for standardization; converts probe IDs to gene IDs.
        Input: a GEOparse object (Sample or Dataset)
        Returns: a pandas dataframe of the GEO object with gene IDs and expression values
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
        Input: GDS object, indicator column, indicator value
        Returns: a pandas DataFrame with filtered data
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
    """ A class that helps find genes with different methylation levels between populations"""

    def __init__(self):
        pass

    def run(self, target_gsm_df, control_gds_df, genes_of_interest):
        """
        Identifies target genes based on methylation levels.
        Input: DataFrames for target sample, control dataset, and list of genes of interest
        Returns: Dictionary of target genes with methylation status
        """
        target_gene_to_methylation = {}
        for gene_id in genes_of_interest:
            target_gene_to_methylation[gene_id] = self.assess_gene_methylation(target_gsm_df, control_gds_df, gene_id)

        return target_gene_to_methylation
    
    def assess_gene_methylation(self, target_sample, healthy_dataset, gene_id):
        """
        Assesses methylation status of a single gene.
        Input: DataFrames for target sample and control dataset, and a single gene ID
        Returns: Methylation status ('hypo' or 'hyper') of the gene
        """
        sample_gene = target_sample[target_sample['IDENTIFIER'] == gene_id]
        healthy_gene = healthy_dataset[healthy_dataset['IDENTIFIER'] == gene_id]

        for i in range(len(sample_gene)):
            healthy_mean = healthy_gene.iloc[i, 2:].mean()
            sample_value = sample_gene['VALUE'].iloc[i]
            percent_difference = abs((sample_value - healthy_mean) / healthy_mean) * 100

            if percent_difference > 10:
                return 'hypo' if sample_value < healthy_mean else 'hyper'

        return None
    
class SequenceProcessor:
    """A class that helps process DNA sequences"""

    def __init__(self):
        pass

    def get_promoter_sequence(self, refseq_id):
        """Attempts to find the promoter sequence of a gene from NCBI given its refseq id"""
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
        Input: promoter sequence of a gene
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
        """
        id_to_seq = {}
        for gene_id in target_genes:
            gene_df = self.target_gsm_df[self.target_gsm_df['IDENTIFIER'] == gene_id]

            if gene_df.empty:
                # Log an error and continue to the next gene_id
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
        Design oligos for constructs
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
    def __init__(self, email):
        Entrez.email = email
        self.geo_processor = GEOProcessor()
        self.target_finder = TargetFinder()
        self.oligo_designer = OligoDesigner()

    def run(self, target_gsm_id, reference_geo_id, genes_of_interest, indicator_col = None, indicator_value = None):
        target_gsm = self.get_target_gsm(target_gsm_id)
        ref_values = self.get_reference_methylation(reference_geo_id, indicator_col, indicator_value)
        genes_to_methylation = self.get_target_genes(target_gsm, ref_values, genes_of_interest)
        gene_to_constructs = self.get_constructs(genes_to_methylation, target_gsm)

        return gene_to_constructs

    def get_target_gsm(self, gsm_id):
        """
        Input: Target GSM ID
        Output: pandas dataframe of the target GSM
        """
        try:
            target_gsm = self.geo_processor.get_sample_dataset(gsm_id)
            return target_gsm
        except Exception as e:
            raise ValueError(f"Not a valid GSM id: {gsm_id}. Error: {e}")
        
    def get_reference_methylation(self, geo_id, indicator_col = None, indicator_value = None):
        if geo_id.startswith("GDS"):
            return self.handle_gds_reference(geo_id, indicator_col = indicator_col, indicator_value = indicator_value)
        elif geo_id.startswith("GSM"):
            return self.handle_gsm_reference(geo_id)
        else:
            raise ValueError("Invalid GEO ID. Please enter a valid GSM or GDS ID.")
        
    def handle_gds_reference(self, gds_id, indicator_col, indicator_value):
        filtered_gds = self.geo_processor.get_filtered_dataset(gds_id, indicator_col, indicator_value)
        ref_values = filtered_gds[['ID_REF', 'IDENTIFIER']].copy()
        ref_values['VALUE'] = filtered_gds.filter(regex='^GSM').mean(axis=1)
        return ref_values
        
    def handle_gsm_reference(self, reference_gsm_id):
        reference_gsm = self.geo_processor.get_sample_dataset(reference_gsm_id)
        reference_gsm = reference_gsm[reference_gsm['VALUE'].notna()]
        ref_values = reference_gsm[['ID_REF', 'IDENTIFIER', 'VALUE']].copy()
        return ref_values

    def get_target_genes(self, target_gsm, ref_values, genes_of_interest):
        """
        Input: Target GSM, reference values, and genes of interest (either as a list or a dict)
        Output: Dictionary of target genes to methylation status
        """
        if isinstance(genes_of_interest, dict):
            genes_to_methylation = genes_of_interest
        elif isinstance(genes_of_interest, list):
            genes_to_methylation = self.target_finder.run(target_gsm, ref_values, genes_of_interest)
        else:
            raise ValueError("genes_of_interest must be a list or a dictionary.")
        return genes_to_methylation
    
    def get_constructs(self, genes_to_methylation, target_gsm):
        return self.oligo_designer.run(genes_to_methylation, target_gsm)


class CLI():
    """Command line interface for EpigeneticDesignAutomation"""
    def __init__(self):
        self.eda = EpigeneticDesignAutomation(email="ethanmccoy@example.com")

    def run(self):
        target_gsm = self.get_target_gsm()
        ref_values = self.get_reference_methylation()
        genes_to_methylation = self.get_target_genes(target_gsm, ref_values)
        constructs = self.get_constructs(genes_to_methylation, target_gsm)

        return constructs

    def get_target_gsm(self):
        """
        Gets a target GSM from user input
        Returns a pandas dataframe of the target GSM
        """
        while True:
            target_gsm_id = input("\nEnter a target GSM id like GSM1324896: ")
            try:
                return self.eda.get_target_gsm(target_gsm_id)
            except ValueError as e:
                print(e)

    def get_reference_methylation(self):
        """
        Entering a GDS ID will filter the GDS for control samples to create an average methylation profile
        Entering a GSM ID will use its methylation profile
        """
        geo_id = input("Enter a GEO ID (GSM or GDS) like GSM1324896 or GDS5047: ")
        if geo_id.startswith("GDS"):
            return self.get_gds_reference_methylation(geo_id)
        elif geo_id.startswith("GSM"):
            return self.eda.handle_gsm_reference(geo_id)
        else:
            print("Invalid GEO ID. Please enter a valid GSM or GDS ID.")

    def get_gds_reference_methylation(self, gds_id):
        indicator_col = input("Enter the indicator column name: ")
        indicator_value = input("Enter the indicator value: ")
        return self.eda.get_reference_methylation(gds_id, indicator_col, indicator_value)
    
    def get_target_genes(self, target_gsm, ref_values):
        """
        Interacts with the user to get target genes.
        Returns a dictionary of target genes to methylation status.
        """
        print("Enter a list or a dict w/ known methylation status"
                "\nDict Example: {'ID4': 'hyper', 'CDKN2B': 'hypo', etc.}"
                "\nList Example: ['ID4', 'CDKN2B', 'CDKN1A', 'CDKN2A']")
        genes_of_interest = eval(input())
        return self.eda.get_target_genes(target_gsm, ref_values, genes_of_interest)
    
    def get_constructs(self, genes_to_methylation, target_gsm):
        return self.eda.get_constructs(genes_to_methylation, target_gsm)

# TODO: Add pickling and file management
def main():
    #cli = CLI()
    #constructs = cli.run()
    #print(constructs)
    print('')


if __name__ == '__main__':
    main()

