import pandas as pd



def save_columns_givenDf(df, cols=[], BOOL_SAVE_TO_FILES=0, output_path="", output_file_stem=""):
    ''' given the next inputs:
            - df Dataframe with the information (for instance each row can be a gene)
              NOTE: the lengths should be sorted numerically first
            -columns: the columns I wand to display or save
            -BOOL_SAVE_TO_FILE:    1 (yes, save it); 0 (no, just display it)
            -output_file: the file to save the data

        Changes: 16.2.21. I change the code and now saves the files faster
        but seems there is some new space at the end of lines. Do not how why (?)
    '''
    stem__out_file = output_path + output_file_stem
    if 0: # just some info about the df
        print("--->output file:", output_path + output_file_stem)
        print("\n\t...save_length_info_givenDf")
        print("\t", len(df), "number of species")
        print("\t", df.columns)
        pd.set_option('display.max_rows', None)
        print(df.head(2))
        print(df.tail(2))
    if 1 and BOOL_SAVE_TO_FILES: # contigs
        s = df["contig"].value_counts().reset_index()
        s.columns = ['contig','counts']
        s.to_csv(stem__out_file+".contig"+".tsv", sep="\t", index = False)

    # 16.2.2019 this has been replaced by the next lines
    # it was extremely slow
    # if 0 and BOOL_SAVE_TO_FILES: # contigs
    #     f = open(stem__out_file+".length"+".tsv", 'w') # lengths
    #     line = "" # header
    #     for c in cols:
    #         line = line+c+"\t"
    #     line = line.rstrip() + "\n"
    #     for index, row in df.iterrows():  # content
    #         for c in cols:
    #             line = line+str(row[c])+"\t"
    #         line = line.rstrip() + "\n"
    #     f.writelines(line)
    #     f.close()

    if 1 and BOOL_SAVE_TO_FILES:  # contigs
        df[cols].to_csv(stem__out_file+".length"+".tsv", sep='\t', index=False)
