#!/usr/bin/env python
# coding: utf-8

# In[1]:


#ignoring warnings
import warnings
warnings.filterwarnings('ignore')


# In[2]:


#importing Library
import os

#printing the current directory
datapath = os.getcwd()


# In[3]:


#Joing the required file to the datapath
file1_path = datapath + "\\B_subtilis_168.ptt"

#print
file1_path

#Importing pandas package as pd
import pandas as pd

#Loading the ppt file as pandas dataframe using read_csv function
df = pd.read_csv(file1_path, sep='\t', names = ['a','b','c','d','e','f','g','h','i'])


# In[4]:


#defining a function called "process_dataframe" that does necessary changes to dataframe
def process_dataframe(df):
    
    #Subsetting the dataframe and including the rows from second and so on
    df = df.iloc[2:]
    
    #Assigning row 1 values as column headers
    df.columns = df.iloc[0]
   
    # Droping the first row from the dataset
    df = df.drop(index=[2])
    
    #Reseting the index
    df = df.reset_index(drop=True)
    
    #Spliting the data from location column into two different columns
    df[['start', 'stop']] = df['Location'].str.split('.', 1, expand=True)
    
    # Removing the period'.' from the 'stop' column
    df['stop'] = df['stop'].str.replace('.', '')
    
    #Converting the stop column values from string to intergers
    df['stop'] = pd.to_numeric(df['stop'])
    
    #Converting the start column values from string to intergers
    df['start'] = pd.to_numeric(df['start'])
    
    #Sorting the dataframe based on start column
    df = df.sort_values(by='start')
    
    return df


# In[5]:


#Processing the dataframe using the predefined function
df = process_dataframe(df)

#Dropping the columns that are not required for further analysis
df.drop(['Location', 'Code', 'COG', 'Product', 'PID', 'Synonym', 'Length'], axis=1, inplace=True)


# In[6]:


#defining a function to calculate the number of operons and the operons
def get_operons(df, column_name):
    
    #Initializing a dictionary with keys as positive and negative strand
    operons = {'+': [[]], '-': [[]]}
    
    #creating a variable to hold the gene value to later use it in calculating the distance
    prev_gene = None

    # iterating over each row and ignoring the index value
    for _, gene in df.iterrows():
        
        #assigning the strand value from each row to the variable
        strand = gene['Strand']
        
        #subsetting an operon using the strand value
        current_operon = operons[strand][-1]

        #if loop to clasify if a gene belongs to the same operon by calculating distance and checking the strand
        if prev_gene is None or (gene['start'] - prev_gene['stop']) > 50 or gene['Strand'] != prev_gene['Strand']:
            
            #appending the gene if the distance is less than 50 and is on the same strand as prev gene
            operons[strand].append([gene[column_name]])
            
        else:
            
            #else creating a new list
            current_operon.append(gene[column_name])

        # updating the prev_gene variable
        prev_gene = gene

    #for loop to print the output in desired manner
    for strand, operon_list in operons.items():
        
        #calculating the total number of operons belonging to +/- strand
        operon_count = len(operon_list) - 1
        
        #print statement which prints the total number of operons on a strand along with the operons
        print(f"Operons for {strand} strand ({operon_count}):")
        for i, operon in enumerate(operon_list[1:], start=1):
            print(f"operon{i}: {operon}")


# In[7]:


get_operons(df, 'Gene')


# In[8]:


#Joing the required file to the datapath
file2_path = datapath + "\\E_coli_K12_MG1655.ptt"

#print
file2_path

#Loading the ppt file as pandas dataframe using read_csv function
df = pd.read_csv(file2_path, sep='\t', names = ['a','b','c','d','e','f','g','h','i'])

#Processing the dataframe using the predefined function
df = process_dataframe(df)

#Dropping the columns that are not required for further analysis
df.drop(['Location', 'Code', 'COG', 'Product', 'PID', 'Synonym', 'Length'], axis=1, inplace=True)


# In[9]:


get_operons(df, 'Gene')


# In[10]:


#Joing the required file to the datapath
file3_path = datapath + "\\Halobacterium_NRC1.ptt"

#print
file3_path

#Loading the ppt file as pandas dataframe using read_csv function
df = pd.read_csv(file3_path, sep='\t', names = ['a','b','c','d','e','f','g','h','i'])

##Processing the dataframe using the predefined function
df = process_dataframe(df)

#Dropping the columns that are not required for further analysis
df.drop(['Location', 'Code', 'COG', 'Product', 'PID', 'Gene', "Length"], axis=1, inplace=True)


# In[11]:


###Used synonym instead of gene names as few of the gene names were missing####
get_operons(df, 'Synonym')


# In[12]:


#Joing the required file to the datapath
file4_path = datapath + "\\Synechocystis_PCC6803_uid159873.ptt"

#print
file4_path

#Loading the ppt file as pandas dataframe using read_csv function
df = pd.read_csv(file4_path, sep='\t', names = ['a','b','c','d','e','f','g','h','i'])

#Processing the dataframe using the predefined function
df = process_dataframe(df)

#Dropping columns that are not required for further analysis
df.drop(['Location', 'Code', 'COG', 'Product', 'PID', 'Synonym', 'Length'], axis=1, inplace=True)


# In[13]:


get_operons(df, 'Gene')


# In[14]:


#Joing the required file to the datapath
gff_file = datapath + "\\2088090036.gff"

#print
gff_file

#Loading the gff file as pandas dataframe using read_csv function
df = pd.read_csv(gff_file, sep='\t', names = ['a','b','c','d','e','f','g','h','i'])

#Renaming the required columns
df = df.rename(columns = {"a" : "contig" , "d" : "start" , "e" : "stop" , "g" : "Strand"})

#Subsetting the dataframe
df = df.iloc[1:-1]

#Resetting the index
df = df.reset_index(drop=True)

#Dropping columns that are not needed for further analysis
df.drop(['b', 'c', 'f', 'h', 'i'], axis=1, inplace=True)

# sort the DataFrame by the 'Start' column
df = df.sort_values(by='start')


# In[15]:


get_operons(df, 'contig')

