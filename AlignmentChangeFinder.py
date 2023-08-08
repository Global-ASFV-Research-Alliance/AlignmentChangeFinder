import pandas as pd
import argparse

#python AlignmentChangeFinder.py -c ASFV_Example.csv -m ..\\..\\muscle3.8.31_i86win32.exe -s "Georgia"

argParser = argparse.ArgumentParser()

argParser.add_argument("-c", "--csv", type= str, help="Input .csv Location")
argParser.add_argument("-m", "--muscle", type= str, help="path to MUSCLE installation")
argParser.add_argument("-o", "--oldgenomes", type=str, help=".txt file of the names of old (reference) genomes",default='OldGenomes.txt')
argParser.add_argument("-n", "--newgenomes", type=str, help=".txt file of the names of new genomes", default='NewGenomes.txt')
argParser.add_argument("-s", "--special", type=str, help="Name of special reference", default=None)

### Optional Inputs
argParser.add_argument("-f", "--folder", type= str, help="Destination Folder", default=".\\Output")
argParser.add_argument("-a", "--aminoacid", type=str, help = "How to display Amino Acids (options = 'Ref', 'Majority')", default='Default')
argParser.add_argument("-p", "--position", type=str, help = "How to display position (options = 'Ref', 'Self', 'Align')", default='Default')
argParser.add_argument("-d", "--dropna", default=False)

args = argParser.parse_args()

my_file = open(args.oldgenomes, "r")
data = my_file.read()
OldList = data.split("\n")
my_file.close()

my_file = open(args.newgenomes, "r")
data = my_file.read()
NewList = data.split("\n")
my_file.close()

if args.position == 'Default':
    if args.special != None:
        args.position = 'Ref'
    else:
        args.position = 'Self'
        
if args.aminoacid == 'Default':
    if args.special != None:
        args.aminoacid = 'Ref'
    else:
        args.aminoacid = 'Majority'

#From FastaTransformer Package
def FastaToAlignment(input_file = 'output_endcount.fasta', output_file = 'output_aligned.fasta', musclepath = 'muscle3.8.31_i86win32.exe', gapextpenalty = -1.0, gapopenpenalty = -10.0):
    
    """
    Takes a .fasta file and uses the Muscle program via MuscleCommandLine to output an alignment file. Inputs are filepath locations, except for gapextpenalty and gapopenpenalty which are aspects of Muscle. It can be applied to an entire folder of .fasta files using :func:`GeneBankToAlignmentBank`.
    
    Documentation for .fasta files can be found at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ and examples provided at https://en.wikipedia.org/wiki/FASTA_format#Description_line.

    Parameters
    ---
    input_file : String of a path to an input .fasta file.

    output_file : String of a path to an output .fasta file.
    
    musclepath : String of a path to the Muscle installation, a .exe file.

    gapextpenalty : Negative Float. Default value of -1.0.

    gapopenpenalty : Negative Float. Default value of -10.0.

    Returns
    ---
    An alignment .fasta file.

    Examples
    ---
    >>> input_file
    >>> .\\FastaBank\\Gene01.fasta
    >Genome01
    AAAABC
    >Genome02
    AAABC

    >>> output_file
    >>> .\\AlignmentBank\\Gene01.fasta
    >Genome01
    AAAABC
    >Genome02
    AAA-BC

    """

    from Bio.Align.Applications import MuscleCommandline
    import os
    with open(output_file, 'w') as fp:
        pass
    
    muscle_exe =  musclepath
    in_file = input_file
    out_file = output_file
    muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file, gapextend = gapextpenalty, gapopen = gapopenpenalty, quiet = True) 
    os.system(str(muscle_cline))

def FolderPathFixer(FolderPath):
    """
    An internal function used to fix folder paths.
    
    Parameters
    ---
    A path to a folder.
    
    Returns
    ---
    The folder path, if not made, will be made. Additionally, a "\\" will be added at the end to allow more items to the end.
    """
    import os
    os.makedirs(FolderPath, exist_ok=True)
    return FolderPath + "\\"

def MakeGeneBank(input_csv, FastaBank, Seq = "Sequence", ID = "Genome", Grouping = "Gene", DropDuplicates = "Partial", DeleteStar = True):

    """
    Take an input csv with Genomes, their genes, and their sequences and split them into a series of .fasta files with the title of the gene and each entry described as a geneome. The output of this function can be fed into :func:`GeneBankToAlignmentBank`.

    Documentation for .fasta files can be found at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ and examples provided at https://en.wikipedia.org/wiki/FASTA_format#Description_line.
    
    Parameters
    ---
    input_csv : A CSV or pd.DataFrame containing columns for a sequence, genome, and genes to be converted into Fasta files.
    
    FastaBank : String of a path to a folder to place the .fasta files in.
    >>> Example
    "Project_AD-Unique Protein Sequences in new Genomes\\FastaBank\\"
    
    Seq : A string indicating the column of the CSV used to indicate the sequence. Default "Sequence".
    
    ID : A stirng indicating the column of the CSV used to indicate the Accession, Genome, or Identifier of the Sequence. Default "Genome".
    
    Grouping : A string indicating the column of the CSV used to split the sequences into different .fasta files. Default "Gene".
    
    DropDuplicates : Indicates if Duplicate Sequences in a Grouping (Such as Gene) should be dropped. Increases computational time. Default Partial.
        All : Remove all identical Seq in a Grouping
        Partial : Remove idential entries only if ID, Seq, and Grouping are the Same
    
    DeleteStar : Indicates if values after an * should be deleted from each sequence, to only keep the portion before the sequences end. Default True.

    Returns
    ---
    Fills the selected FastaBank folder by splitting the rows of the input_csv.

    It also returns a dataframe of all duplicate genes in each genome to allow bugfixing.

    Examples
    ---
    >>> MakeGeneBank('FixSamTSV Code\\p72 Genotyping Website\\Data\\currated_ASFV_db_v02.csv','FixSamTSV Code\\p72 Genotyping Website\\Data\\FastaBank')

    """

    import pandas as pd
    if type(input_csv) == type("ABC"):
        input_df = pd.read_csv(input_csv)
    else:
        input_df = input_csv
    try:
        input_df[Grouping] = input_df[Grouping].str.replace('\*', "", regex=True)
    except:
        pass
    if DropDuplicates == "All":
        input_df = input_df.drop_duplicates([Grouping, Seq])
    if DropDuplicates == "Partial":
        input_df = input_df.drop_duplicates([Grouping, Seq, ID])

    Groups = input_df.groupby(Grouping)

    keys = Groups.groups.keys()

    FastaBank = FolderPathFixer(FastaBank)

    BugBank = pd.DataFrame(columns=[Seq, Grouping, ID])
    for name, group in Groups:
        
        try:
            if len(group[group.duplicated([ID], keep=False)]) > 0:
                print(group[group.duplicated([ID], keep=False)].drop(['Unnamed: 0'], axis=1))
                BugBank = pd.concat([BugBank,(group[group.duplicated([ID], keep=False)]).drop(['Unnamed: 0'], axis=1)])
        except:
            try:
                if len(group[group.duplicated([ID], keep=False)]) > 0:
                    BugBank = pd.concat([BugBank,group[group.duplicated([ID], keep = False)]])
                    group = group.drop_duplicates([ID], keep=False)
            except:
                print("Failed to run, likely duplicates with " + name + " , skipping")
                continue
            
        for index, row in group.iterrows():        
            outputname_fasta = FastaBank + str(name).replace('/', '-') + ".fasta"
            try:
                if DeleteStar == True:
                    row[Seq] = row[Seq].split('*', 1)[0]

                with open(outputname_fasta, "a") as f:
                        print(">" + str(row[ID]) + "\n" + str(row[Seq]).replace("[", "").replace("]", "").replace("'", ""), file = f)
            except:
                try: 
                    print("Failed to Format")
                    print(row)
                except:
                    print("Failed to Format")
                    print(name)
                continue
    return BugBank

def GeneBankToAlignmentBank(FastaBank, AlignmentBank, musclepath = 'muscle3.8.31_i86win32.exe', gapextpenalty = -1.0, gapopenpenalty = -10.0):
    
    """
    Take a folder of .fasta files seperated by gene, and create a matching folder of alignment files using the :func:`FastaToAlignment` function.

    Documentation for .fasta files can be found at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ and examples provided at https://en.wikipedia.org/wiki/FASTA_format#Description_line.

    Parameters
    ---
    FastaBank : String of a path to a folder containing the .fasta files. A FastaBank can be created using the :func:`MakeGeneBank` function.
    >>> Example
    "Project_AD-Unique Protein Sequences in new Genomes\\FastaBank"

    AlignmentBank : String of a path to a folder to place the alignments corresponding with those in the FastaBank.
    >>> Example
    "Project_AD-Unique Protein Sequences in new Genomes\\AlignmentBank"

    musclepath : String of a path to the Muscle installation to be used in the wrapped function :func:`FastaToAlignment`.

    gapextpenalty : Negative Integer. Default value of -1.0.

    gapopenpenalty : Negative Integer. Default value of -10.0.

    Returns
    ---
    Fills the selected AlignmentBank folder with .fasta files corresponding to those in FastaBank.

    Examples
    ---
    >>> GeneBankToAlignmentBank('WebSite_LSDV\\FastaBank','WebSite_LSDV\\AlignmentBank')

    """
    import glob
    
    FastaBank = FolderPathFixer(FastaBank)
    AlignmentBank = FolderPathFixer(AlignmentBank)
    
    my_files = glob.glob(FastaBank + '*')
    for fasta in my_files:
        ouputfasta = fasta.replace(FastaBank, AlignmentBank)
        FastaToAlignment(fasta, ouputfasta, musclepath, gapextpenalty = gapextpenalty,gapopenpenalty = gapopenpenalty)

def InputToList(Object, Keyword = None):

    """Take a list, a pd.Series, or a pd.DataFrame with a specific keyword and return a list of unique strings in that list. Used in this package to prevent failures caused by different uplaod types. Used as a starting function in other functions to help put parameters into a correct format.
    
    Parameters
    ---
    Object : a list, pd.Series, or pd.Dataframe with Keyword set to a value other than None.

    Keyword : Used to designate the column to be used in the case that a pd.DataFrame is submitted.
    
    Returns
    ---
    A list of unique strings.
    """

    import pandas as pd
    if type(Object) is list:
            Output = list(set(Object))
    elif type(Object) is pd.Series:
        Output = list(Object.unique())
    elif type(Object) is pd.DataFrame:
        Output = list(Object[Keyword].unique())
    return Output

def AlignmentChangeFinder(AlignmentBank, NewGenomeNames, OldGenomeNames, SpecialReference = None):
    """
    A function intended to find novel (unique) changes in genes between a set of new genomes and reference ('old') genomes using a set of .fasta files. The function returns a specially formatted dataframe.

    Documentation for .fasta files can be found at https://www.ncbi.nlm.nih.gov/genbank/fastaformat/ and examples provided at https://en.wikipedia.org/wiki/FASTA_format#Description_line.

    Parameters
    ---
    AlignmentBank : String of the path to a folder containing .fasta files with the file named after a gene and the sequences named after the genomes. An AlignmentBank can be generated using the :func:`GeneBankToAlignmentBank` function in this package. 
    >>> Example
    "Project_AD-Unique Protein Sequences in new Genomes\\AlignmentBank"
    
    NewGenomeNames : A list of strings containing the names of genomes to be tested. Must match the genome names used in the AlignmentBank. Can also be submitted as a pd.Series (such as by using df['GenomeColumnName'] for the df that generated the Fasta\\AlignmentBanks).
    >>> Example
    ['Cameroon_2016_C1', 'Cameroon_2016_C5', 'Cameroon_2017_C-A2', 'Cameroon_2018_C02', 'Cameroon_2018_C-F3']
    
    OldGenomeNames : A list of strings containing the names of genomes to be references against. Must match the genome names used in the AlignmentBank. Can also be submitted as a pd.Series (such as by using df['GenomeColumnName'] for the df that generated the Fasta\\AlignmentBanks).
    >>> Example
    ['ASFV-G (Georgia-2007)', 'Malawi Lil-20/1 (Malawi: Chalaswa-1983)', 'L60 (Portugal-1960)', 'BA71V (Spain-1971)', 'Benin 97/1 (Benin-1997)', 'E75 (Spain-1975)', 'OURT 88/3 (Portugal-1988)', 'Warmbaths (South Africa: Warmbaths-1987)', 'Warthog (Namibia-1980)', 'Ken05/Tk1 (Kenya-2005)', 'Ken06.Bus (Kenya-2006)', 'Kenya 1950 (Kenya-1950)', 'Mkuzi 1979 (South Africa: Mkuzi Game Reserve-1979)', 'Tengani 62 (Malawi: Tengani-1962)', 'Pretorisuskop/96/4 (South Africa: Kruger National Park-1996)', 'NHV (Portugal-1968)']

    SpecialReference : If measurements are desired against a specific old/reference geneome, set this equal to the name of a string in OldGenomes. Default value is None.
    >>> Example
    'Benin 97/1 (Benin-1997)'

    Returns
    ---
    The function returns a dataframe that contains all data. This dataframe can be fed into :func:`AlignmentChangeFinderSelector` and/or :func:`AlignmentChangeFinderCleanup` to clean up and select wanted results.
    >>> df
                Gene1   Gene2
    NewGeneome1 list11  list12
    NewGeneome2 list21  list22

    Each list contains elements indicating unique differences between the 'New' Genomes and the Old/reference Genomes. Each element is a formatted string with the following parts, with each part name in parenthesis followed by its data:

    RefAA : A letter indicating the ammino acid changed from in the SpecialReference version of the gene. Only present if a SpecialReference is provided.

    RefPos : A number indicating the position in the SpecialReference version of the gene where the change took place (this number is different from AlignPos because it does not count gaps). Only present if a SpecialReference is provided. Positions begin at 1, not 0.

    Old : The Ammino Acids present in all the reference strains, by number. Each Ammino acid is followed by #, then the number of ammino acids at the Alignposition of that type, then . and the next ammino acid present.
    
    AlignPos : A number indicating the position in the Alignment version of the gene where the change took place (counting gaps). Only present if a SpecialReference is provided. Positions begin at 1, not 0.

    SelfPos : A number indicating the position in the New Genome version of the gene where the change took place (this number is different from AlignPos because it does not count gaps). Positions begin at 1, not 0.

    SelfAA : A letter indicating the ammino acid changed to in the new Genome.

    >>> Example
    ['(RefAA)V (RefPos)104 (Old)V#16.-#12.S#5 (AlignPos)105 (SelfPos)104 (SelfAA)I']

    """
    import Bio.SeqIO as SeqIO
    import pandas as pd
    import glob

    NewGenomeNames = InputToList(NewGenomeNames, 'Genomes')
    OldGenomeNames = InputToList(OldGenomeNames, 'Genomes')
    AlignmentBank = FolderPathFixer(AlignmentBank)
    
    if SpecialReference != None:
        if SpecialReference not in OldGenomeNames:
            print("SpecialReference Must be in OldGenomeNames, To Avoid Error SpecialReference set to None")

    my_files = glob.glob(AlignmentBank + '*')
    FinalFrame = pd.DataFrame(NewGenomeNames, columns=['id']).set_index('id')
    NoHistoricGenes = []
    NoNewGenes = []
    for fasta in my_files:

        #fasta = AlignmentBank+'B962L.fasta'
        gene = fasta.replace(AlignmentBank, "").replace(".fasta", "").replace(".fa", "")

        df = pd.DataFrame(columns=['id', 'sequence'])
        for Genome in SeqIO.parse(fasta,format='fasta'):
            df.loc[len(df.index)] = [Genome.description, str(Genome.seq)]
        
        Old_df = df[df['id'].isin(OldGenomeNames)].set_index('id')
        if len(Old_df) == 0:
            NoHistoricGenes.append(gene)
            continue
        
        if SpecialReference != None:
            try:
                Spec_Seq = df[df['id'] == SpecialReference].set_index('id').iloc[0]['sequence']
            except:
                Spec_Seq = None
        else:
            Spec_Seq = None

        New_df = df[df['id'].isin(NewGenomeNames)].set_index('id')
        if len(New_df) == 0:
            NoNewGenes.append(gene)
        
        #Duplicate Entires Check
        if len(df[df.duplicated(['id'], keep=False)]) > 0:
            print("Gene " + gene + " contains duplicate entries, gene skipped to avoid crash.")
            print(df[df.duplicated(['id'], keep=False)])
            continue

        AllSequences = []
        for sequence in list(map(''.join, zip(*Old_df['sequence']))):
            AllSequences.append('.'.join([i + "#" + str(sequence.count(i)) for i in set(sequence)]))
        
        def FindFlaw(Sequence, CompareSequences, Spec_Seq = None):
            Flaw = []
            try:
                for i in range(len(Sequence)):
                    if Sequence[i] not in CompareSequences[i]:
                        
                        if Spec_Seq == None:
                            SSinfo = ""
                        else:
                            SSinfo =  "(RefAA)"+Spec_Seq[i] + " " + "(RefPos)"+str(len(Spec_Seq[:i+1].replace("-",""))) + " "

                        Flaw.append(SSinfo + "(Old)" + CompareSequences[i]+ " " + "(AlignPos)" + str(i+1) + " " + "(SelfPos)" + str(len(Sequence[:i+1].replace("-",""))) + " " + "(SelfAA)" + Sequence[i])
                if len(Sequence) == 0:
                    return "NA"
                if len(Flaw) == 0:
                    return ""
                return Flaw
            except:
                print("Failed At FindFlaw on " + fasta)
                return "ERROR"
        
        try:
            results = pd.DataFrame()
            results[str(gene)] = New_df['sequence'].map(lambda x: FindFlaw(x, AllSequences, Spec_Seq))
            FinalFrame = pd.concat([FinalFrame,results], axis=1)
        except:
            print("Failed at Final Frame Concatonation at " + fasta)
            continue
    
    if len(NoHistoricGenes) > 0:
        print("The following genes did not have a historic entry and were skipped: ")
        print(NoHistoricGenes)
    if len(NoNewGenes) > 0:
        print("The following genes did not have a new entry, they were not skipped but their rows will have na. These can be removed in AlignmentChangeFinderCleanup by setting DropNA = True: ")
        print(NoNewGenes)
    
    return FinalFrame

def Rangemaker(info, FirstAA, FirstPositions):
        
        """
        A function for cleaning up a list object within a pd.DataFrame produced by :func:`AlignmentChangeFinder` into a more interpretable form. This function primarily exists for development purposes. To convert an entire pd.DataFrame, use its wrapper, :func:`AlignmentChangeFinderCleanup`.

        Parameters
        ---
        info : A list entry in a pd.DataFrame constructed by :func:`AlignmentChangeFinder`.

        FirstAA : Indicates which Ammino Acid should be shown for the result. Options include:

            'Ref' : Use the Ammino Acid indicated by (RefAA). Default option.

            'Majority' : Uses the most common Ammino Acid in the (Old) substring.

        FirstPositions : Indicates which positional marker should be shown for the result. Options include:

            'Ref' :  Use the position from the (RefPos) substring. Does not include gaps. Default option.

            'Self' : Use the data in the (SelfPos) substring. Does not include gaps.

            'Align' : Use the data in the (AlignPos) substring.

        Returns
        ---
        A cleaned up string in place of the list. The string is the chosen ammino acid, chosen position, and final ammino acid (followed by a comma and any other aread). In the case there are continous positions in the list, the elements are combined into a single entry. If either the chosen ammino acid or the final position are the same indicator, they are converted to one object (for example, AAAA173TQKY would become A173TQKY).

        >>> An example list before cleanup.
        ['(RefAA)C (RefPos)352 (Old)Y#3.C#10.-#3 (AlignPos)386 (SelfPos)334 (SelfAA)G', '(RefAA)P (RefPos)353 (Old)S#2.P#11.-#3 (AlignPos)387 (SelfPos)335 (SelfAA)R', '(RefAA)K (RefPos)356 (Old)K#8.E#3.-#5 (AlignPos)390 (SelfPos)338 (SelfAA)H', '(RefAA)C (RefPos)358 (Old)Y#1.H#1.C#9.-#5 (AlignPos)392 (SelfPos)340 (SelfAA)G', '(RefAA)S (RefPos)366 (Old)S#8.-#8 (AlignPos)400 (SelfPos)348 (SelfAA)P', '(RefAA)E (RefPos)368 (Old)E#7.-#9 (AlignPos)402 (SelfPos)350 (SelfAA)K', '(RefAA)S (RefPos)369 (Old)S#6.T#1.-#9 (AlignPos)403 (SelfPos)351 (SelfAA)P', '(RefAA)Y (RefPos)370 (Old)Y#8.-#8 (AlignPos)404 (SelfPos)352 (SelfAA)C', '(RefAA)S (RefPos)371 (Old)S#8.-#8 (AlignPos)405 (SelfPos)353 (SelfAA)P']
        >>> Changed Example
        'CP352-353GR, K356H, C358G, S366P, ESYS368-371KPCP'

        Warnings
        ---
        If there are no references in the info list, then the output ['NoReference'] will be given instead.

        """

        import pandas as pd
        from itertools import groupby
        from operator import itemgetter
        import pandas as pd

        if FirstAA or FirstPositions == 'Ref':
            if '(RefPos)' not in info[0]:
                return ['NoReference']

        Positions = []
        AAFirst = []
        AALast = []       
        
        for string in info:
            splitstring = string.split(" ")
            AALast.append(list(filter(lambda x: '(SelfAA)' in x, splitstring))[0].replace('(SelfAA)',""))
            if FirstPositions == 'Self':
                Positions.append(int(list(filter(lambda x: '(SelfPos)' in x, splitstring))[0].replace('(SelfPos)',"")))
            if FirstPositions == 'Ref':
                Positions.append(int(list(filter(lambda x: '(RefPos)' in x, splitstring))[0].replace('(RefPos)',"")))
            if FirstPositions == 'Align':
                Positions.append(int(list(filter(lambda x: '(AlignPos)' in x, splitstring))[0].replace('(AlignPos)',"")))
            if FirstAA == 'Ref':
                AAFirst.append(list(filter(lambda x: '(RefAA)' in x, splitstring))[0].replace('(RefAA)',""))
            if FirstAA == 'Majority':
                Old = list(filter(lambda x: '(Old)' in x, splitstring))[0].replace('(Old)',"").split('.')
                OldAA = []
                OldCount = []
                for OldEntry in Old:
                    OldAA.append(OldEntry.split('#')[0])
                    OldCount.append(int(OldEntry.split('#')[1]))
                AAFirst.append(OldAA[OldCount.index(max(OldCount))])
        
        #Identify Ranges
        
        def RangeMakerInside(AAFirst, Positions, AALast):
            listdf = pd.DataFrame(list(zip(AAFirst, Positions, AALast)), columns=['AAFirst','Positions','AALast']).sort_values('Positions', ascending=True)    
            data = sorted(set(Positions))
            Ranges = []
            for k, g in groupby(enumerate(data), lambda ix : ix[0] - ix[1]):
                Rangedf = listdf.loc[listdf['Positions'].isin(list(map(itemgetter(1), g)))]

                AAFirstRange = ''.join(list(Rangedf['AAFirst']))
                if len(set(AAFirstRange)) == 1:
                    AAFirstRange = ''.join(set(AAFirstRange))
                AALastRange = ''.join(list(Rangedf['AALast']))
                if len(set(AALastRange)) == 1:
                    AALastRange = ''.join(set(AALastRange))
                if AAFirstRange == "-":
                    AAFirstRange = "Ins"
                if AALastRange == "-":
                    AALastRange = "Del"
                
                if len(Rangedf) > 1:
                    Ranges.append(AAFirstRange + str(Rangedf['Positions'].iloc[0]) + "-" + str(Rangedf['Positions'].iloc[-1]) + AALastRange)
                else:
                    Ranges.append(AAFirstRange + str(Rangedf['Positions'].iloc[0]) + AALastRange)
            return Ranges

        return RangeMakerInside(AAFirst, Positions, AALast)

def AlignmentChangeFinderCleanup(df, FirstAA = 'Ref', FirstPositions = 'Ref', GeneList = None, DropNA = True):
    
    """
    A function for cleaning up a pd.DataFrame produced by :func:`AlignmentChangeFinder`. It is a wrapper for the :func:`Rangemaker` function, which should be viewed for further explanation.

    Parameters
    ---
    df : A pd.DataFrame constructed by :func:`AlignmentChangeFinder`.
    
    FirstAA : Indicates which Ammino Acid should be shown for the result. Options include:

        'Ref' : Use the Ammino Acid indicated by (RefAA). Default option.

        'Majority' : Uses the most common Ammino Acid in the (Old) substring.

    FirstPositions : Indicates which positional marker should be shown for the result. Options include:

        'Ref' : Use the position from the (RefPos) substring. Does not include gaps. Default option.

        'Self' : Use the data in the (SelfPos) substring. Does not include gaps.

        'Align' : Use the data in the (AlignPos) substring.
        
    Genelist : A list, or pd.Series object with a list of Genes to be kept. All other genes will be removed.

    DropNa : Drop all genes only containing NaN entries. Does not include Genes that are empty because there were no genes missing in the Genomes and differences were not found.

    """

    import numpy as np
    import pandas as pd
    
    if GeneList != None:
        GeneList = InputToList(GeneList, 'Gene')
        df = df[GeneList]
    if DropNA == True:
        df = df.dropna(how='all', axis=1)

    def list2Str(lst, FirstAA, FirstPositions):
        if type(lst) is list: # apply conversion to list columns
            return ", ".join(map(str, Rangemaker(lst, FirstAA, FirstPositions)))
        else:
            return lst
    
    df_clean = df.apply(lambda x: [list2Str(i, FirstAA, FirstPositions) for i in x]).replace(np.nan,"na") 

    return df_clean

#Remaining Code
print("Making Gene Bank")
MakeGeneBank(input_csv = pd.read_csv(args.csv), FastaBank=args.folder + "\\FastaBank", ID="Genome", Seq="Sequence", Grouping="Gene")
print("Making Alignment Bank")
GeneBankToAlignmentBank(FastaBank=args.folder + "\\FastaBank",AlignmentBank= args.folder + "\\AlignmentBank", musclepath=args.muscle)
print("Searching for Novel Changes")
df1 = AlignmentChangeFinder(args.folder + "\\AlignmentBank", NewGenomeNames = NewList, OldGenomeNames = OldList, SpecialReference = args.special)
df1.transpose().to_csv(args.folder + "\\rawoutput.csv")
print("Cleaning Novel Changes")
df_sparkles = AlignmentChangeFinderCleanup(df1, FirstAA = args.aminoacid, FirstPositions = args.position, DropNA=args.dropna)
df_sparkles.transpose().to_csv(args.folder + "\\output.csv")
print("Run complete. For more details on internal functions, see the github page for the FastaTransformer package at https://github.com/Global-ASFV-Research-Alliance/FastaTransformer")