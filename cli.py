from EpigeneticDesignAutomation import *

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
    
def main():
    cli = CLI()
    constructs = cli.run()
    print(constructs)

if __name__ == '__main__':
    main()