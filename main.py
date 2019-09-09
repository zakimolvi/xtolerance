""""
ALLELIC CROSS TOLERANCE SCORING

This script uses neoepitope seq input and HLA typing to:
1. find matches in an unmutated database according to the scoring criteria
2. Calculate NetMHCpan binding preds, %rank for neo and unmut epitopes
3. Output, for each neoepitope, matching unmut at each allele and the % rank value.

Modalities of attaining a matching score:
1 1 1 1 1
1 1 0.25 1 1
1 1 1 0.25 1
1 1 0.5 0.5 1
"""
import pandas as pd
from Bio import Align
from mhctools import NetMHCpan4

def predictBinding(hla, seqs):
    """
    Use NetMHCpan4.0 backend to make predictions for unmut and neo
    :param hla: patient hla typing
    :param seqs: peptide sequences
    :return: df containing binding predictions
    """
    predictor = NetMHCpan4(alleles=hla)
    predictor._check_hla_alleles(alleles=hla)

    bindPreds = predictor.predict_peptides_dataframe(seqs)
    return bindPreds

def prepPtHLA(fname='allHLAtypings.csv'):
    """
    We must create a dictionary of all the patient HLA typings
    return: HLA typing as query-able dictionary
    ex. pt_hla[PID] = [HLA typing as list]
    """
    import csv

    with open(fname, encoding='utf-8-sig') as csv_data:
        reader = csv.DictReader(csv_data)
        pt_hla = {}

        for row in reader:
            key = row['PID']
            val = row['HLA'].split(',')

            pt_hla[key] = val

    return pt_hla

def neoPreds():
    """Submit all neoepitopes for NetMHCpan4 predictions.
    -combinedInput.csv is long-form summary from all cohorts
    """

    #get all PIDs in a list
    df = pd.read_csv('combinedInput.csv')
    samples = df.Sample.drop_duplicates().tolist()
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665') #pt CR3665 excluded in automation bc A33:51 bug with local installation of netmhcpan
    for ID in dictHLA.keys():
        print('Current pt:'+ ID)
        hla = dictHLA[ID]
        seqs = df[df.Sample == ID]['MT.Peptide'].tolist()
        output = predictBinding(hla=hla, seqs=seqs)
        output.to_csv("neo_preds/"+ID+"_mtpeptidePreds.csv")
#neoPreds()


def filterNeo():
    """
    Open neo predictions and filter <2.5% rank
    -save to csv
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')

    for ID in dictHLA.keys():
        print(ID)
        df = pd.read_csv("neo_preds/"+ID+"_mtpeptidePreds.csv")
        df = df.drop('Unnamed: 0', axis=1)
        df = df[df.percentile_rank <= 2.50]
        df.to_csv("neo_preds_filtered/"+ID+"_mtpeptidePredsFiltered.csv")

#filterNeo()

def subMatrix(fname='matrix.txt'):
    """
    Takes tab-delimited text file and creates dictionary substitution matrix
    :param fname: substitution/scoring matrix, 'matrix.txt'
    :return: scoring matrix dict
    """
    fl = [] #array that will hold lines of txt input
    with open(fname) as f:
        for line in f:
            fl.append(line.split())

    matrix = {}
    for j in range(1, len(fl[:][0])+1):
        for i in range(0,len(fl[0])):
            #1st AA
            key = (fl[j][0], fl[0][i])
            val = float(fl[j][i+1])

            matrix[key] = val
    return matrix
    ## Write matrix to file for examination
    #with open('matrix_test.txt', 'w') as f:
    #     for key in matrix.keys():
    #         f.write("%s,%s\n"%(key,matrix[key]))


def matchAndScore(neo, db):
    """
    Find peptides in unmutated DB that match at P4,5,8 using regex
    then finds those with score >=1.25 at P6+7 using scoring matrix.

    :param neo: neoepitope string input
    :param db: unmutated db as pandas DF
    :return: Matching peptides with neoepitope as pandas DF
    """
    import re
    pattern = '^([A-Z]){3}'+neo[3]+neo[4]+'([A-Z]){2}'+neo[7]
    #get pos 4, 5, 8 matches
    pos_matches = db[db.str.contains(pattern, regex=True)].tolist()

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = subMatrix()
    unmut_matches = []
    for match in pos_matches:
        score6 = aligner.score(neo[5], match[5])
        score7 = aligner.score(neo[6], match[6])
        score67 = aligner.score(neo[5:7], match[5:7]) #compare P6 and P7

        if (score67 >= 1.0 and score6 >=0.25 and score7 >=0.25):
            unmut_matches.append(match)

    scoredDF = pd.DataFrame()
    scoredDF['Peptide'] = unmut_matches
    scoredDF['Neoepitope'] = neo
    scoredDF['P4,5,8'] = neo[3]+neo[4]+neo[7]
    return scoredDF

def neoMatchAndScore():
    """This part takes a few hours to run.

    For each patient, open their binding neoepitope DF and find unmutated matches
    by means of matchAndScore, then save them in a DF in neo_wt_matches/
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    unmut = pd.read_excel("unmutated database.xlsx", squeeze=True, index_col=None, header=None)
    for ID in dictHLA.keys():
        print("Current pt: "+ID)
        df = pd.read_csv('neo_preds_filtered/'+ID+'_mtpeptidePredsFiltered.csv')
        score_all = pd.DataFrame() #holds all scored matches for pt
        for seq in df['peptide'].drop_duplicates():
            scored = matchAndScore(seq, unmut)
            score_all = score_all.append(scored)
        score_all.to_csv("neo_wt_matches/"+ID+"_scored_matches.csv")
#neoMatchAndScore()

def wtPreds():
    """"
    Take wt matches for each pt and get binding preds by netmhcpan4
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    for ID in dictHLA.keys():
        hla = dictHLA[ID]
        print("Current pt: "+ID)
        df = pd.read_csv('neo_wt_matches/'+ID+'_scored_matches.csv')
        if df.empty is False:
            seqs = df['Peptide'].drop_duplicates().tolist()
            preds = predictBinding(hla=hla, seqs=seqs)
            preds.to_csv('wt_preds/'+ID+'_wt_preds.csv')

#wtPreds()

def pivotWtMt():
    """
    Takes NetMHC predictions for wt and neoepitopes and pivots them to wide format
    -filters wt preds for <=2.5% rank
    -save to (wt, neo)_preds_filtered_pivoted/
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    for ID in dictHLA.keys():
        hla = dictHLA[ID]
        print("Current pt: " + ID)
        if (ID == '16') or (ID == '42') or (ID=='71'): #manually exclude these patients with no wt matches
            continue
        #open wt preds, filter, and pivot
        df = pd.read_csv('wt_preds/'+ID+'_wt_preds.csv')
        df = df.drop('Unnamed: 0', axis=1)
        df = df[df.percentile_rank <= 2.50]
        df = df.pivot(index='peptide', columns='allele', values = ['affinity', 'percentile_rank'])
        #df.to_csv('wt_preds_filtered_pivoted/'+ID+'_wt_preds_filtered_pivoted.csv')
        df.to_pickle('wt_preds_filtered_pivoted/'+ID+'_wt_preds_filtered_pivoted.pkl')

        #open filtered neo preds, and pivot
        df = pd.read_csv('neo_preds_filtered/'+ID+'_mtpeptidePredsFiltered.csv')
        df = df.drop('Unnamed: 0', axis=1)
        df = df.pivot(index='peptide', columns='allele', values=['affinity', 'percentile_rank'])
        #df.to_csv('neo_preds_filtered_pivoted/'+ID+'_neo_preds_filtered_pivoted.csv')
        df.to_pickle('neo_preds_filtered_pivoted/' + ID + '_neo_preds_filtered_pivoted.pkl')
#pivotWtMt()

def neoWtPreds():
    """
    Open each neo preds DF, search for its matching wt peptides, get wt preds and append
    them to the neo DF.
    -pt 16+42+71 in Mellman had no wt matches for neoepitopes so we exclude them
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    dictHLA.pop('16')
    dictHLA.pop('42')
    dictHLA.pop('71')
    for ID in dictHLA.keys():
        hla = dictHLA[ID]
        print("Current pt: " + ID)

        db = pd.read_csv('neo_wt_matches/'+ID+'_scored_matches.csv')
        db = db.drop('Unnamed: 0', axis=1)
        neo = pd.read_pickle('neo_preds_filtered_pivoted/'+ID+'_neo_preds_filtered_pivoted.pkl')
        neo['quality'] = 'neo'
        neo['Neoepitope'] = neo.index
        wt = pd.read_pickle('wt_preds_filtered_pivoted/'+ID+'_wt_preds_filtered_pivoted.pkl')
        wt['quality'] = 'unmut'

        for seq in neo.index:
            xdb = db[db.Neoepitope == seq] #slice of wt:neo DB that has matching wt peptides for neo
            if xdb.empty == False:
                wtslice = wt[wt.index.isin(xdb.Peptide)] #take slice of wt_preds DF that has neo-wt matches of interest
                wtslice['Neoepitope'] = seq
                print(wtslice)
                neo = neo.append(wtslice) #append matching wt preds to neo DF
        neo = neo.sort_values(by=['Neoepitope', 'quality']) #unmut matches are under neo
        neo.to_csv('neo_wt_preds/' + ID + '_joint.csv')
        neo.to_pickle('neo_wt_preds/'+ID+'_joint.pkl')
#neoWtPreds()


def fixedNeoBools():
    """"
    Add TRUE or FALSE to neoepitopes with matching unmutated peptides.
    """
    import numpy as np
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    dictHLA.pop('16')
    dictHLA.pop('42')
    dictHLA.pop('71')
    for ID in dictHLA.keys():
        print("Current pt: " + ID)

        df = pd.read_pickle('neo_wt_preds/' + ID + '_joint.pkl')
        neo = df[df[('quality', '')] == "neo"]  # neo df
        unmut = df[df[('quality', '')] == "unmut"]  # unmutdf matching neo

        for seq in neo.index:
            print("Current Neo: "+seq)
            #seq_matches = []
            wtmatch = unmut[unmut[('Neoepitope', '')] == seq]
            neoSP = df[df[('Neoepitope', '')] == seq] #current neoepitope

            for allele in wtmatch['percentile_rank'].columns.tolist():
                #iterating over all the alleles in the neoepitope's predictions
                print(allele)
                col = allele+' match' #new column name for allele matches
                neo_NA = neoSP[('percentile_rank', allele)].isna()
                if neo_NA.iloc[0] == True: #if it's true that neoepitope has NaN %rank at this point, it has no matches bc >2.5% rank
                    df.set_value(seq, col, False)
                    #seq_matches.append(False)
                else:
                    rank = neoSP[neoSP.index == seq][('percentile_rank', allele)].iloc[0]
                    score_bool = (wtmatch[('percentile_rank')] <=5*rank) & (wtmatch[('percentile_rank')] >=0.2*rank)
                    #if score_bool has any True values, then the neo is a match at the current allele
                    if score_bool.any(axis=None): #looking at ALL alleles to find unmut binders that could tolerize
                        df.set_value(seq, col, True)
                        #seq_matches.append(True)
                    else:
                        df.set_value(seq, col, False)
                        #seq_matches.append(False)
        df.to_csv('finalMatches/'+ID+'.csv')
        df.to_pickle('finalMatches/'+ID+'.pkl')

#fixedNeoBools()

def sameAlleleBools():
    """
    Test without allelic xtolerance to find matches that are only at the
    same allele as the one the neoepitope is predicted to bind
    -run for %rank as before
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    dictHLA.pop('16')
    dictHLA.pop('42')
    dictHLA.pop('71')
    for ID in dictHLA.keys():
        print("Current pt: " + ID)

        df = pd.read_pickle('neo_wt_preds/' + ID + '_joint.pkl')
        neo = df[df[('quality', '')] == "neo"]  # neo df
        unmut = df[df[('quality', '')] == "unmut"]  # unmutdf matching neo

        for seq in neo.index:
            print("Current Neo: " + seq)
            # seq_matches = []
            wtmatch = unmut[unmut[('Neoepitope', '')] == seq]
            neoSP = df[df[('Neoepitope', '')] == seq]  # current neoepitope

            for allele in wtmatch['percentile_rank'].columns.tolist():
                print(allele)
                col = allele + ' match'
                neo_NA = neoSP[('percentile_rank', allele)].isna()
                if neo_NA.iloc[0] == True:  # if it's true that neoepitope has NaN %rank at this point, it has no matches bc >2.5% rank
                    df.set_value(seq, col, False)

                else: #if the neoepitope does bind at the given allele
                    rank = neoSP[neoSP.index == seq][('percentile_rank', allele)].iloc[0]
                    score_bool = (wtmatch[('percentile_rank')] <= 5 * rank) & (
                                wtmatch[('percentile_rank')] >= 0.2 * rank)
                    # if score_bool has any True values, then the neo is a match at the current allele
                    if score_bool[allele].any(axis=None): #looking at SAME allele to find tolerizing unmut binders
                        df.set_value(seq, col, True)
                    else:
                        df.set_value(seq, col, False)

        df.to_csv('sameAlleleMatches/' + ID + '.csv')
        df.to_pickle('sameAlleleMatches/' + ID + '.pkl')

#sameAlleleBools()

def sameAlleleBoolsByAffinity():
    """
    Test without allelic xtolerance to find matches that are ONLY at the
    same allele as the one the neoepitope is predicted to bind
    -run for absolute affinity
    """
    dictHLA = prepPtHLA()
    dictHLA.pop('CR3665')
    dictHLA.pop('16')
    dictHLA.pop('42')
    dictHLA.pop('71')
    for ID in dictHLA.keys():
        print("Current pt: " + ID)

        df = pd.read_pickle('neo_wt_preds/' + ID + '_joint.pkl')
        neo = df[df[('quality', '')] == "neo"]  # neo df
        unmut = df[df[('quality', '')] == "unmut"]  # unmutdf matching neo

        for seq in neo.index:
            print("Current Neo: " + seq)
            # seq_matches = []
            wtmatch = unmut[unmut[('Neoepitope', '')] == seq]
            neoSP = df[df[('Neoepitope', '')] == seq]  # current neoepitope

            for allele in wtmatch['affinity'].columns.tolist():
                print(allele)
                col = allele + ' match'
                neo_NA = neoSP[('affinity', allele)].isna()
                if neo_NA.iloc[0] == True:  # if it's true that neoepitope has NaN %rank at this point, it has no matches bc >2.5% rank
                    df.set_value(seq, col, False)

                else: #if the neoepitope does bind at the given allele
                    rank = neoSP[neoSP.index == seq][('affinity', allele)].iloc[0]
                    score_bool = (wtmatch[('affinity')] <= 5 * rank) & (
                                wtmatch[('affinity')] >= 0.2 * rank)
                    # if score_bool has any True values, then the neo is a match at the current allele
                    if score_bool[allele].any(axis=None): #looking at SAME allele to find tolerizing unmut binders
                        df.set_value(seq, col, True)
                    else:
                        df.set_value(seq, col, False)

        df.to_csv('sameAlleleMatchesByAffinity/' + ID + '.csv')
        df.to_pickle('sameAlleleMatchesByAffinity/' + ID + '.pkl')

#sameAlleleBoolsByAffinity()