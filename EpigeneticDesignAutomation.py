from subprocess import call

call('pip install importlib_metadata', shell=True)
from importlib_metadata import version, PackageNotFoundError

def install_packages():
    REQUIRED_PACKAGES = [
        'Bio',
        'gdown',
        'pandas',
        'GEOparse'
    ]

    for package in REQUIRED_PACKAGES:
        try:
            dist_version = version(package)
            print('{} ({}) is installed'.format(package, dist_version))
        except PackageNotFoundError:
            print('{} is NOT installed'.format(package))
            call("pip install " + package, shell=True)



class AutoEdit:
    def __init__(self):
        print('')

    def find_oligos(self, gds_id, target_gsm_id, indicator_col, indicator_value, genes_of_interest):
        print('Begin Running AutoEdit')

        geo_processor = GEOProcessor()
        target_finder = TargetFinder()
        oligo_designer = OligoDesigner()

        print('Downloading and filtering GDS')
        filtered_gds = geo_processor.get_filtered_dataset(gds_id, indicator_col, indicator_value)
        target_gsm = geo_processor.get_sample_dataset(target_gsm_id)

        target_genes = target_finder.run(target_gsm, filtered_gds, genes_of_interest)
        print('Found target genes:', target_genes)

        oligos = oligo_designer.run(target_genes, target_gsm)
        print('Found oligos:', oligos)

        return oligos



class GEOProcessor:
    def get_filtered_dataset(self, gds_id, indicator_col, indicator_value):
        gds = GEOparse.get_GEO(geo = gds_id, destdir="./", silent = True)
        filtered_gds_df = self.filter_GDS_by_indicator(gds, indicator_col, indicator_value) # Keep only the indicated control samples
        return filtered_gds_df

    def get_sample_dataset(self, gsm_id):
        gsm = GEOparse.get_GEO(geo = gsm_id, destdir="./", silent = True)
        gsm_df = self.geo_with_gene_id(gsm)
        return gsm_df

    def geo_with_gene_id(self, geo):
        platform_column = 'platform_id' if 'platform_id' in geo.metadata else 'platform'
        platform_df = GEOparse.get_GEO(geo = geo.metadata[platform_column][0], destdir="./", silent = True).table
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
    def __init__(self):
        print('')

    def run(self, target_gsm_df, control_gds_df, genes_of_interest):
        targets = self.find_targets(target_gsm_df, control_gds_df, genes_of_interest)
        return targets

    def find_targets(self, target_sample, healthy_dataset, genes_of_interest):
        target_gene_to_methylation = {}
        for gene_id in genes_of_interest:
            sample_gene = target_sample[target_sample['IDENTIFIER'] == gene_id]
            healthy_gene = healthy_dataset[healthy_dataset['IDENTIFIER'] == gene_id]

            for i in range(len(sample_gene)):
                healthy_mean = healthy_gene.iloc[i, 2:].mean().mean()
                sample_value = sample_gene['VALUE'].iloc[i]

                # If greater than 10 % difference add to list
                percent_difference = abs((sample_value - healthy_mean) / healthy_mean) * 100
                if percent_difference > 10:
                    if sample_value > healthy_mean:
                        target_gene_to_methylation[gene_id] = 'hypo'
                    else:
                        target_gene_to_methylation[gene_id] = 'hyper'

        return target_gene_to_methylation

    def find_potential_genes_of_interest(self):
        print('TODO')



class OligoDesigner:
    def __init__(self):
        print('')

    def run(self, target_genes, target_gsm_df):
        sp = SequenceProcessor()

        # Find the promoter sequence of each target gene
        id_to_seq = {}
        for gene_id in target_genes:
            gene_df = target_gsm_df[target_gsm_df['IDENTIFIER'] == gene_id]

            if ('DETECTION P-VALUE' in gene_df.columns):
                present_transcript_id = gene_df.sort_values(by='DETECTION P-VALUE').iloc[0:1]['GB_ACC'].iloc[0]
            else:
                present_transcript_id = gene_df.sort_values(by='Detection Pval').iloc[0:1]['GB_ACC'].iloc[0]

            promoter_seq = sp.get_promoter_sequence(present_transcript_id)
            id_to_seq[gene_id] = promoter_seq

        # Find possible protospacers for each target promoter
        id_to_protospacers = {}
        for id, seq in id_to_seq.items():
              protospacers = sp.find_protospacers(seq)
              id_to_protospacers[id] = protospacers

        # Design sgRNA oligos to address hyper/hypomethylation in each target gene
        id_to_oligos = {}
        for gene_id, methylation in target_genes.items():
            protospacers = id_to_protospacers.get(gene_id)

            oligos = []
            for protospacer in protospacers.values():
                if methylation == 'hyper':
                    oligos.append(self.design_Tet1_oligos(protospacer))
                elif methylation == 'hypo':
                    oligos.append(self.design_DNMT3A_oligos(protospacer))

            id_to_oligos[gene_id] = oligos

        return id_to_oligos

    # U6 promoter and gRNA scaffold w/ MS2 binding sites is already on the plasmid
    def design_Tet1_oligos(self, protospacer):
        Acc65IPlus5 = 'ataatGGTACC'
        forward_oligo = Acc65IPlus5 + protospacer
        reverse_oligo = Acc65IPlus5
        return (forward_oligo, reverse_oligo)

    def design_DNMT3A_oligos(self, protospacer):
        Acc65IPlus5 = 'ataatGGTACC'
        forward_oligo = Acc65IPlus5 + protospacer
        reverse_oligo = Acc65IPlus5
        return (forward_oligo, reverse_oligo)


class SequenceProcessor:
    def get_promoter_sequence(self, refseq_id):
        try:
            # Fetch the sequence data from NCBI
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
        # Search for PAM sites on forward strand of promoter
        pam_sites = {}
        for i in range(20, len(promoter_seq) - 2):
            if promoter_seq[i+2:i+3] == 'G' and promoter_seq[i+1:i+2] == 'G':
                protospacer = promoter_seq[i-20:i]
                protospacer = str(protospacer) # Removes seq object
                if protospacer[0] == 'G' and (protospacer[16] == 'A' or protospacer[16] == 'T'):
                    pam_sites[i] = protospacer.lower()

        return pam_sites
    

def psoriasis_test(): 
    gds_id = 'GDS2518' # Plaque psoriasis dataset
    target_gsm_id = 'GSM154769' # Lesional psoriasis sample
    indicator_col = 'disease state'
    indicator_value = 'uninvolved'
    genes_of_interest = ['ID4', 'CDKN2B', 'CDKN1A', 'CDKN2A']

    ae = AutoEdit()
    oligos = ae.find_oligos(gds_id, target_gsm_id, indicator_col, indicator_value, genes_of_interest)

    return oligos

def cocaine_test(): 
    gds_id = 'GDS5047' # cocaine brains dataset
    target_gsm_id = 'GSM1324896' # cocaine abuse sample
    indicator_col = 'agent'
    indicator_value = 'control'
    genes_of_interest = ['CDKN1A', 'CCL2', 'PVALB', 'SLC6A3']

    geo_processor = GEOProcessor()
    target_finder = TargetFinder()
    oligo_designer = OligoDesigner()

    ae = AutoEdit()
    oligos = ae.find_oligos(gds_id, target_gsm_id, indicator_col, indicator_value, genes_of_interest)

    return oligos
    
def main():
    install_packages()
    from Bio import Entrez
    from Bio import SeqIO
    import gdown
    import pandas as pd
    import GEOparse

    Entrez.email = "ethanmccoy@example.com"  # Replace with your email
    openAI_auth = ''

    print('Hi, whats your name?')

    name = input()

    print('Hello', name)

    psoriasis_oligos = psoriasis_test()
    #cocaine_oligos = cocaine_test()

    print(psoriasis_oligos)

if __name__ == '__main__':
    main()