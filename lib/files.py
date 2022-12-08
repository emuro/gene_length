# python3
# ################################################################## #
# files.py (C) March-2021 Mainz.
# Author: Enrique M. Muro
# ################################################################## #
#
# ------------------------------------------------------------------------
# Project: geneLength
#
# Purpose: module to deal with files
#
# ################################################################## #
import sys
sys.path.append('./lib/')

import pandas as pd


def creates_an_empty_file(input_file=""):
    ''' inputs:
            -input_file,
                The file (with path)
        output:
            No ouput
        Note: creates an empty file
     '''

    try:
        f = open(input_file, "a")
        f.close()  # creates an empty file
    except Exception as e:
        print (__name__, "\t", creates_an_empty_file.__name__, ",", sep="")
        print("\t", e)
        sys.exit()
    return

def add_line_to_existing_file(input_file="", new_line=""):
    ''' inputs:
            -input_file,
                The file (with path)
        Note: The file must be already created
     '''

    try:
        f = open(input_file, "a")
        f.write(new_line + "\n")
        f.close()  # creates an empty file
    except Exception as e:
        print (__name__, "\t", add_line_to_existing_file.__name__, ",", sep="")
        print("\t", e)
        sys.exit()
    return


def get_df_from_tsv_file(input_file=""):
    ''' inputs:
            -input_file,
                The file (with path) where the data is saved
        output:
            -df_data
                A pandas data frame with the file data
        Note: It test if the file exists
     '''

    try:
        df_data = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError as e:
        print(e)
        sys.exit()
    except Exception as e:
        print(e)
        sys.exit()
    else:
        if 0:
            pd.options.display.max_columns=None
            print(df_data.head(10))
        if 0:
            print(df_data.columns)
        if 0:
            print(df_data.shape)
    finally:
        pass

    return df_data


def get_df_from_tsv_file_or_create_file(input_file=""):
    ''' inputs:
            -input_file,
                The file (with path) where the data was presumably saved
        output:
            -df_data
                A pandas data frame with the file's data

        Note: If the file does not exists, just create it
     '''

    try:
        df_data = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError as e:
        print(get_df_from_tsv_file_or_create_file.__name__, ",", sep="")
        print("\t", e)
        creates_an_empty_file(input_file) # creates an empty file
        df_data = pd.DataFrame()
    except Exception as e:
        print(get_df_from_tsv_file_or_create_file.__name__,",",sep="")
        print("\t",e, ". File: (", input_file, ")", sep="")
        df_data=pd.DataFrame()
    if 0:
        pd.options.display.max_columns=None
        print(df_data.head(10))
    if 0:
        print(df_data.columns)
    if 0:
        print(df_data.shape)

    return df_data


def save_columns_givenDf(df, cols=[], BOOL_SAVE_TO_FILES=0, output_path="", output_file_stem=""):
    ''' given the next inputs:
            - df Dataframe with the information (for instance each row can be a gene)
              NOTE: the lengths should be sorted numerically first
            -columns: the columns I wand to display or save
            -BOOL_SAVE_TO_FILE:    1 (yes, save it); 0 (no, just display it)
            -output_file: the file to save the data
    '''
    ## !! HACERLO ENTERO
    stem__out_file = output_path + output_file_stem
    if 0: # just some info about the df
        print("--->output file:", output_path + output_file_stem)
        print("\n\t...save_length_info_givenDf")
        print("\t", len(df), "number of genes")
        print("\t", df.columns)
        pd.set_option('display.max_rows', None)
        print(df.head(2))
        print(df.tail(2))
    if 1 and BOOL_SAVE_TO_FILES: # contigs
        s = df["contig"].value_counts().reset_index()
        s.columns = ['contig','counts']
        s.to_csv(stem__out_file+".contig"+".tsv", sep="\t", index = False)

    if 1 and BOOL_SAVE_TO_FILES: # contigs
        f = open(stem__out_file+".length"+".tsv", 'w') # lengths
        line = "" # header
        for c in cols:
            line = line+c+"\t"
        line = line.rstrip() + "\n"
        for index, row in df.iterrows():  # content
            for c in cols:
                line = line+str(row[c])+"\t"
            line = line.rstrip() + "\n"
        f.writelines(line)
        f.close()


def get_list_from_column_from_tsv_file(input_file="", column=""):
    """
        From a *.tsv file with headers, get a list from an specified column

        Parameters
        ----------
        input_file: str (*.tsv),
            The file (with path) where the data was presumably saved
        column: str
            The name of the column to get the list
        Output:
        ------
        list
            Data of the column from the data of the input_file
     """

    try:
        df_data = pd.read_csv(input_file, sep="\t")
    except FileNotFoundError as e:
        print(get_list_from_column_from_tsv_file.__name__, ",", sep="")
        print("\t", e)
    except Exception as e:
        print(get_list_from_column_from_tsv_file.__name__, ",", sep="")
        print("\t", e, ". File: (", input_file, ")", sep="")

    return list(df_data[column])