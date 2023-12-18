from dataclasses import dataclass
from EpigeneticDesignAutomation import *
from typing import List

@dataclass(frozen=True)
class Test:
    name : str
    gds_id: str
    target_gsm_id: str
    indicator_col: str
    indicator_value: str
    genes_of_interest: List[str]

class Tests():
    def __init__(self):
        self.eda = EpigeneticDesignAutomation(email="example@email.com")
        self.tests = [
            Test(name = 'Psoriasis test', gds_id='GDS2518', target_gsm_id='GSM154769', 
                 indicator_col='disease state', indicator_value='uninvolved', 
                 genes_of_interest=['ID4', 'CDKN2B', 'CDKN1A', 'CDKN2A']),
            Test(name = 'Cocaine abuse test', gds_id='GDS5047', target_gsm_id='GSM1324896', 
                 indicator_col='agent', indicator_value='control', 
                 genes_of_interest=['CDKN1A', 'CCL2', 'PVALB', 'SLC6A3'])
        ]

    def run(self):
        for test in self.tests:
            self.run_test(test)

    def run_test(self, test: Test):
        gene_to_constructs = self.eda.run(test.target_gsm_id, test.gds_id, test.genes_of_interest, 
                                  indicator_col=test.indicator_col, indicator_value=test.indicator_value)
        print(f"\nTest results for {test.name}:")
        for gene, constructs in gene_to_constructs.items():
            print('\n' + gene)
            for construct in constructs:
                print('\nplasmid:', construct.plasmid_name)
                print('restriction site:', construct.downstream_U6_sites[0])
                print('forward oligo:', construct.target_gene_forward_oligo)
                print('reverse oligo:', construct.target_gene_reverse_oligo)

def main():
    tests = Tests()
    tests.run()

if __name__ == '__main__':
    main()