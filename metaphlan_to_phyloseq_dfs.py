import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import re


def clean_merged_table(df,merged_file_name):
    '''Cleans the merged table by removing file extenstion suffix which gets added to sample columns during merge.'''
    df = pd.read_table(df,sep='\t',engine='python')
    df.columns = df.columns.str.replace(merged_file_name, '')
    return df

def filter_rows_by_taxa(df,rank):
    '''Filters all metaphlan rows by ID to return only selected taxanomic rank i.e. specify f for family, s for species'''
    df=df[df['ID'].str.contains(r'\|'+rank+'__[^|]*$')|(df['ID']=='UNKNOWN')]
    df.reset_index(drop=True,inplace=True)
    return df

def get_taxa_columns(df,rank):
    '''Splits ID into taxanomic ranks to make taxa table''' 
    df_taxa = df['ID'].str.split('|',expand=True)
    taxa_cols = ["Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain"]
    taxa_dict = {'Kingdom':1,"Phylum":2,"Class":3,"Order":4,"Family":5,"Genus":6,"Species":7,"Strain":8}
    value = taxa_dict.get(rank)
    taxa_cols=taxa_cols[0:value]
    df_taxa.columns=taxa_cols
    for col in df_taxa.columns:
        df_taxa[col]=df_taxa[col].apply(trim_taxa_names)    
    #df_taxa=df_taxa.drop_duplicates(subset=df_taxa.columns, keep='last').reset_index(drop=True)
    otu_index = []
    for i in range(0, len(df)):
        otu_index.append("Otu"+str(i))
    df_taxa['Otu']=otu_index 
    return df_taxa

def trim_taxa_names(x):
    '''Removes leading characters before taxa ID e.g. s__ '''
    match = re.sub(r'^[kpcofgs]__',"",str(x))
    return match

def get_sample_cols(df):
    '''Finds and returns sample columns in dataframe, presuming sample names contain a number'''
    r = re.compile(r'^.*[0-9].*$') #match column names that contain a number anywhere
    sample_cols=[]
    for col in df:
        if(r.match(col)):
            sample_cols.append(col)
    return sample_cols

def create_sample_df(abun_matrix):
    sample_cols=get_sample_cols(abun_matrix)
    sample_df=pd.DataFrame({'Sample':sample_cols})
    sample_df['Behaviour'] = sample_df['Sample'].apply(get_behaviour)
    # Add extra column to sample df so phyloseq ordination plots behave
    sample_df['Type'] = 'murine'
    return sample_df

def add_otu_primary_key(df):
    '''Adds otu primary key column to dataframe'''
    otu_index = []
    for i in range(0, len(df)):
        otu_index.append("Otu"+str(i))
    df['Otu']=otu_index 
    return df

def get_behaviour(x):
    '''Creates sample_data columns based on sample name starting with letter'''
    if x.startswith('E'):
        return "Extinction"
    elif x.startswith('R'):
        return "Resilient"
    elif x.startswith('S'):
        return "Susceptible"
    else:
        return ''


if __name__ == "__main__":
    # Create abundance matrix
    df_k_1=clean_merged_table("data/pool1_profile_known.txt","_known_profiled_metagenome")
    df_k_2=clean_merged_table("data/pool2_profile_known.txt","_known_profiled_metagenome")
    df_known = pd.merge(df_k_1,df_k_2,on='ID',how="outer")
    abun_matrix=add_otu_primary_key(df_known)
    
    # Create taxa dataframe 
    rank_matrix=filter_rows_by_taxa(abun_matrix,"s")
    taxa_cols=get_taxa_columns(rank_matrix,"Species")
    
    # Remove ID column from abundance matrix after creating taxa cols
    abun_matrix.drop("ID",axis=1,inplace=True)
    # Create sample df
    sample_df=create_sample_df(abun_matrix)
    
    print(abun_matrix.head())
    print(taxa_cols.head())
    print(sample_df.head())
    
    # Export to csv to import into R
    abun_matrix.to_csv("species_known_abundance.csv",index=False)
    taxa_cols.to_csv("species_known_taxa.csv",index=False)
    sample_df.to_csv("sample_df.csv",index=False)
 
    
    