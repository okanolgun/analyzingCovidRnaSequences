#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 


metadata = pd.read_csv("ncbi_datasets.csv")


# In[2]:


metadata.head()


# In the first row, you can see the "RefSeq", the reference sequence. 
# The first thing we can say that, it has 29 000 nucleotides. It shows us that how a virus can be so complex for nucleotides  

# In[3]:


metadata.shape


# In[ ]:





# In[4]:


metadata["Collection Data"] = pd.to_datetime(metadata["Collection Date"])

#pass the collection date into that function
#it ll convert it from a string to a date time 


# In[5]:


metadata.dtypes

#checked "collection date" is like how we want


# In[6]:


metadata.columns = [c.lower().replace(" ","_",) for c in metadata.columns]

#uppercases and spaces gone. 
#metadata.columns gives you all of the columns names


# In[7]:


metadata


# In[ ]:





# In[8]:


metadata["continent"] = metadata["geo_location"].str.replace(";.+","",regex=True) 

#started to data analysis 
#created a continent column. because it is together with continent and city/country
#";.+","",regex=True means that replace the semicolon and ANYTHING after it with NOTHING 
#and we tell pandas yo use a regular expression


# In[9]:


metadata

#we can see "continent" columns is created


# In[ ]:





# In[10]:


metadata.groupby("continent").apply(lambda x: x.sort_values("collection_date").iloc[0]) 

#we looked first coded sequence that sequenced on each continent 
#we group the data by continent, sorted them by time and applied a function to only take first row 
#for example, we can see in Africa row, first sqeuence was collected in 2020 february 


# In[11]:


metadata["continent"].value_counts()

#we looked how many sequences collected in each continent


# We can say that, in North America, there were so many nucleotides that collected. Almost twice that of its closest competitor, Europe.
# 
# This can give us an idea of how much the virus mutates in the US.

# In[12]:


metadata["nucleotide_length"].describe()

#to figure out the longest nucleotide length is ad what the shortest is


# We can observe from this table what was the variation in the virus.  In general, the covid virus has around 29000 nucleotides along the covid rna sequence.  

# In[13]:


metadata[metadata["nucleotide_length"] == metadata["nucleotide_length"].min()]

#this means that, search metadata for any row where the nucleotide length, equals the minimum lucleotide length


# Now, if we look at this row we can see that the "completeness" column says it is a complete genome. So this is probably just a typo and they meant to out an extra zero or some other digit in here. 
# 

# In[ ]:





# In[14]:


sample_month = pd.Series(pd.to_datetime(metadata["collection_date"]).values.astype("<M8[M]"))
 
#we can look at how many samples were collected by month. 
#we basically changes the type of the column.

#"<M8[M]", actually this code isn't a thing that comes up very often 
#It's converting it to a specific kind of date value. 
#And this date value is basically just year and month.

#The reason why we are doing this is it's going to make our sample of plotting how many samples are collected by month


# In[15]:


sample_month.value_counts().sort_index().plot()

#to just summarize how many samples collected in each month. 


# We can see the activity increased in January 2021, and then it increased in December 2022 again. Perhaps the relaxation of quarantine rules around the world caused these increases.
# 
# We could do a lot more exploratory data analysis if we wanted to. 

# In[ ]:





# In[16]:


metadata[metadata["sequence_type"] == "RefSeq"] 


# We pulled out the reference sequence, "RefSeq" here, which is the first covid genome. 

# In[17]:


metadata[metadata["isolate_name"].str.contains("Delta").fillna(False)] 


# We looked specifically at "Delta" variant. In "isolate_name" column contains information about which variant it is. We searched the column and find if there is any "Delta" variant.
# 
# Some of the isolate names are missing and we want to fill "na" with false. 
# 
# We can see the sequences of the delta variant after we run that code line. 

# In[18]:


metadata[metadata["isolate_name"].str.contains("Omicron").fillna(False)] 


# We did the same thing for "Omicron" variant and we observed that there is only one "Omicron" sequence as well. 

# In[ ]:





# In[19]:


sequences = ["NC_045512.2","OL467832.1","OM061695.1","OM095411.1"]
human_names = ["reference","base","delta","omicron"] 

#"NC_045512.2": we already looked that sequence. it is the first covid genome. 
#"OL467832.1": this is the first sequence from north america, just for we have another base sequence to look at. Just to see how that compares to the reference sequence.  
#"OM061695.1": for delta, we grabbed this sequence because this is tagged "Delta-1"
#"OM095411.1": for omicron, we grabbed this sequence because this is tagged "Omicron-1"

#we created a list to be reference. the first is a reference sequence, then we have the base sequence, delta and omicron. 
#these are 4 sequences that we are going to work with 


# In[20]:


selected_sequences = metadata[metadata["nucleotide_accession"].isin(sequences)]

#we can just look at the metadata for thoso sequences. the id is called "nucleotide_session".
#and we'll basically look for any rows in metadata where the nucleotide accession is in this list(sequences)


# In[21]:


selected_sequences 

#with this line we can see the sequences we pulled out.


# 

# In[ ]:





# Now we are going to use a library called "biopyhton" which has a really handy method which helps you look up RNA sequences and DNA sequences by accession number and then download them. 
# 
# If you didn't install this package yet, you can run the code "!pip install biopyhton" 
# 

# In[22]:


from Bio import Entrez 
Entrez.email = "okanolgun980@gmail.com"

#We need to set an email so I set my email but you can change it. The reason why we are doing this, whenever you download a sequence,
#It sends your email to the National Institutes of Health (NIH), the organization which we downloaded our dataset. 

#It is just a security thing to do to let them know who you are and who is trying to get the datas from them. 


# In[ ]:





# Now let them contact with you, if you are doing anything weird or unsafe with how you're accessing the site. 
# 

# In[23]:


def download_sequence(id_code) : 
    handle = Entrez.esearch(db="nucleotide", term=id_code, retmax="1") 
    #this esearch func. will search the nih's database for a given sequence  
    #we're going to tell it to search the nucleotide database  
    #and search it using our id code which is this nucleotide accession number 
    #retmex="1" is for returning only a single value. 
    
    record = Entrez.read(handle)
    #we need to read from that search and will see result of search
    
    handle = Entrez.efetch(db="nucleotide", id=record["IdList"][0], rettype="fasta", retmode="text") 
    #we're doing another search for nucleotide database 
    #we'll pass in our id.  
    #record["IdList"][0], this is just internals to the library
    #rettype="fasta", we want to get the sequence in "fasta" format 
    #retmode="text", we want the results in text format. 
    
    return handle.read() 
    #gives us back the sequence. 


    # ** fasta is just a specific format for RNA and DNA sequences. 


# Actually, download_sequence(id_code) is a really high level method. It took a while for me to understand. 
# 
# This function will take in an id number for a nucleotide sequence and give us back the sequence. 

# In[ ]:





# In[24]:


sequence_data = {} 
for sequence in sequences : 
    sequence_data[sequence] = {"fasta": download_sequence(sequence)} 
    
#It creates a loop to loop though all of our sequences and for each one.
#it's going to download the sequences from the nih's database and store it into this dictionary. 

#It might take a little bit to run because it has to talk to an external site and do some downloading 


# In[25]:


# sequence_data


# If you run the last comment line, it's going to be really long and you can see the id number and the actual sequence for the covid virus. It will be on the top. 
# 
# Also, we can see the sequences for delta and omicton variants as well. So we can analyze them in just a second. 

# In[ ]:





# You may noticed that the data is in "fasta" format. We actually want to parse set that and just turn it into something we can read more easily. With the next for loop, we'll be able to do that. 

# In[26]:


#we're going to use the "fasta parser" from bioPython
from Bio import SeqIO 
import io 

#We gooing to loop through the sequence data dictionary. (k)eys are the ids and the (v)alues the actual sequences 
for k, v in sequence_data.items() :
    f = io.StringIO(v["fasta"])
    #basically cerating a file object using the "fasta" sequence
    #this is basically says parse has to turn it into a file. Because bioPyhton only works with files. 
    
    sequence_data[k]["parsed"] = list(SeqIO.parse(f,"fasta"))[0] 
    #and we going to parse that using the bioPyhton parser
    #created a new key in our dictionary, called "parsed". This will be the sequence that is parsed 
    #this parsing will create several different objects and we only want the first one. [0] was for that. 
    
    


# StringIO gives you file-like access to strings, so you can use an existing module that deals with a file and change almost nothing and make it work with strings.

# In[27]:


# sequence_data


# If you run the comment line, you can see that we have a parsed field now. I am not running it now because it'll be a very long output. 

# In[ ]:





# In[28]:


sequence_data["NC_045512.2"]["parsed"]


# This is the actual parsed sequence record for RefSeq. For example we can align as one of the different methods of analyzing data. 

# ![image.png](attachment:image.png) 

# Actually, there is a process behind the alignment like I showed you the picture above. It's a little bit tricky. We need to figure out where the sequences align and where they dont. We need to understand which part of these sequences aligns to the other part and which part doesn't align. 
# 
# Our goal is to figure out which portions of the sequence align because if we can understand that, we can figure out how different the sequences are. We'll calculate a score for how different the sequences are. 

# In[ ]:





# In[29]:


from Bio import Align
aligner = Align.PairwiseAligner() 

#we created this pairwise aligner class that will basically help us align multiple sequences.


# In[30]:


aligner.algorithm 

#our aligner is using an algorithm to actually align the sequences. 
#This algorithm is basically how the computer can look at two 'Rna' sequences and figure out which parts overlap.


# In the output we're seeing the 'Needleman-Wunsch'. 
# The Needleman-Wunsch algorithm is an algorithm used in bioinformatics to align protein or nucleotide sequences. 

# In[31]:


score = aligner.score(sequence_data["NC_045512.2"]["parsed"].seq,sequence_data["OM061695.1"]["parsed"].seq)

#we can use the aligner to find a score. We'll find the alignment between reference sequence and one other sequence. 

#["NC_045512.2"]["parsed"], for reference sequence 

#.seq, gives us the actual sequence from the parsed data structure. 
# It will be look something like that : seq=Seq('ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGT...AAA' 

#["OM061695.1"]["parsed"], for delta-1 sequence 


# In[32]:


score


# We can see our alignment score here. This score essentially tells you how aligned the sequences are in relation to the number of nucleotides in the sequence. 

# In[33]:


len(sequence_data["NC_045512.2"]["parsed"]) 

#the length of our reference sequence. 


# In[34]:


29818.0 / 29903 


# We can see that these two sequences are pretty fairly aligned. 99% aligned so it's pretty close. 
# 
# We can say that the virus didn't mutate/change a lot between first being sequenced in China and then later being sequenced a couple months later in North America. 
# 
# We observed that how a virus mutated over time which it's really cool because you can predict that how fast a virus can change.

# In[ ]:





# In[35]:


import numpy as np 

comparisons = np.zeros((4,4))

#created a 4x4 empty matrix 


# In[36]:


comparisons


# In[37]:


for i in range(0,4):
    for j in range(0,i+1):
        score = aligner.score(sequence_data[sequences[i]]["parsed"].seq,sequence_data[sequences[j]]["parsed"].seq)
        comparisons[i,j] = score
        
#this inner for loops will compare each sequence with another sequence 
# and it'll assign the score values to '0' places one by one


# In[38]:


comparisons_df = pd.DataFrame(comparisons, columns=human_names, index=human_names)


# In[39]:


comparisons_df 


# In[40]:


comparisons_df = pd.DataFrame(comparisons, columns=human_names, index=human_names)


# In[41]:


comparisons_df


# In[42]:


comparisons_df.iloc[:,0] / 29903


# We got a the kind of percentage overlap between each of the variants and the reference sequence. 
# 
# So we can see, "base" was the sequence from North America two months after reference sequence. The overlap is very high for them. 
# 
# "delta" sequence overlap is still very high but a little bit lower than "base". 
# 
# It shows us, as time goes on, the virus is mutuating more and more. It becomes more different than the reference sequence. We call it "mutuation". 

# In[ ]:





# In[43]:


seq1 = sequence_data["NC_045512.2"]["parsed"].seq
seq2 = sequence_data["OM061695.1"]["parsed"].seq

delta_alignments = aligner.align(seq1, seq2) 


# In[44]:


delta_alignments


# In[45]:


delta_alignment = delta_alignments[0]


# In[46]:


delta_alignment.aligned 


# We got these result for checking the mutation points for our sequences. 
# 
# Actually, this is a tuple of tuples. Or you can just think that it is a list of lists. 
# 
# It shows you the overlaps between the sequences. 
# 
# We can through this list to figure out where the insertions and deletions.

# In[47]:


seq1_end = None 
seq2_end = None

for alignments in zip(delta_alignment.aligned[0], delta_alignment.aligned[1]) :
    if seq1_end and seq2_end : 
        seq1_mismatch = seq1[seq1_end: alignments[0][0]]
        seq2_mismatch = seq2[seq2_end: alignments[1][0]] 
        print("1: {}".format(seq1_mismatch))
        print("2: {}".format(seq2_mismatch)) 
    
    seq1_end = alignments[0][1] 
    seq2_end = alignments[1][1]
        

#the whole for loop is looping through all of the alignments. and saying "if it's the first alignment, then skip if statement". 
# and figure out where the end of the alignment is and assign it. 

#anything from the end of the previous to the start of the next alignment we're going to print out
# because that's a mismatch and it's the more important one for us. 


# Now, it shows us all of the mismatches. There is not a lot but it doesn't show us where  the mismatches are happening 
# 
# Next thing we can do is smarten this up a little bit for showing us to actual mismatches. 
# 
# We'll color it for making them more understandable. 

# In[48]:


from IPython.display import HTML 

#IPython.display can let us display data in html format.


# In[49]:


def color_print(s, color="black") : 
    return "<span style='color:{}'>{}</span".format(color, s)  


#This function just displays data in color for us. 
#we'll style this with some HTML. 
#It'll return HTML and it will display that html in JupyterNotebook. 


# In[ ]:





# Now we need to check/improve this is existing sequence of code. 
# 
# We'll show the type of mismatch and in order to do that, we gonna check the lengths of the mismatches. 
# 
# for the sequence 1, if the mismatch is zero and for sequence two, the mismatch is several characters then it's an insertion. 
# Other way around it is a deletion. 
# If they match in length, then it is a substitution. 

# In[50]:


#copied from before 

seq1_end = None 
seq2_end = None
display_seq = [] #we need to print them separately. so created an empty list for 
                 #keep track of the all of the things we need. 

for alignments in zip(delta_alignment.aligned[0], delta_alignment.aligned[1]) :
    if seq1_end and seq2_end : 
        seq1_mismatch = seq1[seq1_end: alignments[0][0]]
        seq2_mismatch = seq2[seq2_end: alignments[1][0]] 
        if len(seq2_mismatch) == 0:
            display_seq.append(color_print(seq1_mismatch, "red")) #"red" for deletion
        elif len(seq1_mismatch) == 1:
            display_seq.append(color_print(seq2_mismatch, "green")) #"green" for insertion 
        else : 
            display_seq.append(color_print(seq2_mismatch, "blue")) #blue for substitution
        
    
    display_seq.append(seq1[alignments[0][0]:alignments[0][1]])
    
    
    seq1_end = alignments[0][1] 
    seq2_end = alignments[1][1]


# In[51]:


display_seq = [str(i) for i in display_seq] 


# In[52]:


display(HTML('<br>'.join(display_seq))) 

#we joined our display sequence


# With this output, we can see the full sequence. 
# 
# Black is where the sequences align, 
# 
# Red is the deletion, means that, something from the original sequence was deleted. In this case the reference sequence was deleted by the delta variant. 
# 
# Green is where something was inserted by the delta variant, 
# 
# Blue is a substitution.  
# 
# 
# You can see exactly where those occurred all.
