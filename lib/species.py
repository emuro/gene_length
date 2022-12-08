# python3
# ################################################################## #
# main_species.py (C) March-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: The Class Species
#
# ################################################################## #
import pandas as pd
from lib import constants as c


class Species:
    """A single species
    """
    def __init__(self, name="", division="", assembly="", db="", db_v=0):
        self.name = name
        self.division = division
        self.assembly = assembly
        self.db = db
        self.db_version = str(db_v)

        self.base_path = ""
        self.base_annotation_path = ""
        self.root_annotation_path = ""
        self.trunk_annotation_path = ""
        self.annotation_file = ""
        self.db_annotation_file = ""

        self.base_biotype_distribution_path = ""
        self.root_biotype_distribution_path = ""
        self.trunk_biotype_distribution_path = ""

        self.base_genes_path = ""
        self.root_genes_path = ""
        self.trunk_genes_path = ""

        self.min_number_of_genes_per_biotype = ""
        self.number_of_biotypes = 0
        self.biotype_distribution = pd.DataFrame(columns=['biotype', 'number'])

    def init_annotation_os(self, base_path="", root_annot_path="", trunk_annot_path="",
                           annot_file="", db_annot_file=""):
        self.base_path = base_path
        self.base_annotation_path = self.base_path
        self.root_annotation_path = root_annot_path
        self.trunk_annotation_path = trunk_annot_path
        self.annotation_file = annot_file
        self.db_annotation_file = db_annot_file

    def init_biotype_distribution_os(self, base_path="", root_biotype_distr_path="",
                                     trunk_biotype_distribution_path=""):
        self.base_biotype_distribution_path = base_path
        self.root_biotype_distribution_path = root_biotype_distr_path
        if not trunk_biotype_distribution_path:
            self.trunk_biotype_distribution_path = \
                self.trunk_annotation_path.replace("/gtf/", '/' + c.BIOTYPE_DISTR_PATH_NAME + '/')
        else:
            self.trunk_biotype_distribution_path = trunk_biotype_distribution_path

    def init_biotype_distribution_from_df(self, dist, minimum=0):
        self.min_number_of_genes_per_biotype = minimum
        self.biotype_distribution = dist.copy()
        self.number_of_biotypes = len(dist.biotype)

    def init_biotype_distribution_from_string(self, dist_string, minimum=0):
        self.min_number_of_genes_per_biotype = minimum

        biotypes_dist = dist_string.split(",")
        df_dist = pd.DataFrame([], columns=['biotype','number'])
        for s in biotypes_dist:
            df_dist.loc[len(df_dist.index)] = s.split(":")
        self.biotype_distribution = df_dist.copy()
        self.number_of_biotypes = len(df_dist.biotype)
        del df_dist

    def init_genes(self, base_path, root_path, trunk_path):
        self.base_genes_path = base_path
        self.root_genes_path = root_path
        self.trunk_genes_path = trunk_path

    def show_info(self, bool_show_number_of_biotypes=0):
        print(f'Species,\n\tname: {self.name}\n\tdivision: {self.division}\n\tassembly: {self.assembly}'
              f'\n\tdb: {self.db}\n\tdb_version: {self.db_version}'
              f'\n\tbase_annotation_path: {self.base_annotation_path}'
              f'\n\troot_annotation_path: {self.root_annotation_path}'
              f'\n\ttrunk_annotation_path: {self.trunk_annotation_path}'
              f'\n\tannotation_file: {self.annotation_file}'
              f'\n\tdb_annotation_file: {self.db_annotation_file}'
              #
              f'\n\tbase_biotype_distribution_path: {self.base_biotype_distribution_path}'
              f'\n\troot_biotype_distribution_path: {self.root_biotype_distribution_path}'
              f'\n\ttrunk_biotype_distribution_path: {self.trunk_biotype_distribution_path}'
              #
              f'\n\tbase_genes_path: {self.base_genes_path}'
              f'\n\troot_genes_path: {self.root_genes_path}'
              f'\n\ttrunk_genes_path: {self.trunk_genes_path}'
              )
        if bool_show_number_of_biotypes:
            print(f'\t#\n\tmin_number_of_genes_per_biotype: {self.min_number_of_genes_per_biotype}'
                  f'\n\tnumber_of_biotypes: {self.number_of_biotypes}'
                  f'\n\tbiotype_distribution:\n\t{self.biotype_distribution}')

    def string_basic_info(self, header=0):
        if header:
            return "species\tassembly\tdivision\tdb\tdb_version"
        else:
            return(f'{self.name}\t{self.assembly}\t{self.division}\t{self.db}\t{self.db_version}')

    def string_biotype_distribution_paths(self, header=0):
        if header:
            return "root_biotype_distribution_path\ttrunk_biotype_distribution_path"
        else:
            return(f'{self.root_biotype_distribution_path}\t{self.trunk_biotype_distribution_path}')

    def string_with_the_biotype_distribution(self, header=0):
        if header:
            return "min_number_of_genes_per_biotype\tnumber_of_biotypes\tsorted_biotype_distribution"
        else:
            str_out = f'{self.min_number_of_genes_per_biotype}\t{self.number_of_biotypes}\t'
            for index, row in self.biotype_distribution.iterrows():
                str_out += row['biotype'] + ":" + str(row['number']) + ","
            str_out = str_out[:-1]
            return str_out

    def get_annotation_file_with_full_path(self):
        return self.base_path + self.root_annotation_path + self.trunk_annotation_path + self.annotation_file

    def get_biotype_distribution_full_path(self):
        return self.base_biotype_distribution_path + self.root_biotype_distribution_path + \
               self.trunk_biotype_distribution_path

    def get_genes__full_path(self):
        aux_trunk_path = self.trunk_genes_path
        return self.base_genes_path + self.root_genes_path + \
                aux_trunk_path.replace(self.name+"/", "")

    def get_genes_of_species__full_path(self):
        return self.base_genes_path + self.root_genes_path + \
               self.trunk_genes_path

    def get_file_name_of_biotype_distribution_histogram(self):
        return "biotype_distribution_of_species_" + self.name + "_" + str(self.division) + "_" + \
               self.db + "_" + str(self.db_version) + ".png"

    def get_file_name_of_genes_of_species_and_biotype(self, biotype=""):
        return biotype + ".genes." + self.name + "." + str(self.division) + "." + \
               self.db + "." + str(self.db_version) + ".tsv"

    def get_species_file_name_of_genes_count_foreach_contig(self):
        return "count_genes_in_contigs." + self.name + "." + str(self.division) + "." + \
               self.db + "." + str(self.db_version) + ".tsv"

    def get_biotype_file_name_of_genes_count_foreach_contig(self, biotype=""):
        return "count_genes_in_contigs." + biotype + "." + self.name + "." + str(self.division) + "." + \
               self.db + "." + str(self.db_version) + ".tsv"

    def get_file_name_of_genes_of_species_and_biotype__full_path(self,biotype=""):
        return self.get_genes_of_species__full_path() + self.get_file_name_of_genes_of_species_and_biotype(biotype)
