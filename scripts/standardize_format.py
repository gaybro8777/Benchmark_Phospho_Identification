import re
import pandas as pd

#standardized maxquant
def standardized_maxquant_msms(data):
    df_maxquant = data[data['Phospho (STY)']!=0]
    df_maxquant=df_maxquant[['Scan number','Charge','Sequence','Length','Proteins','Modified sequence',
                             'Phospho (STY) Probabilities','Phospho (STY)','Retention time']]
    pep=df_maxquant['Modified sequence'].tolist()
    phospep=[]
    for i in pep:
        p1 = re.compile(r'p', re.S)
        p2 = re.compile(r'[(](ox)[)]', re.S)
        p3 = re.compile(r'[(](ac)[)]', re.S)
        p4 = re.compile(r'[(](ph)[)]', re.S)
        tar_index1 = [m.start() for m in re.finditer(p1, i)]
        tar_index2 = [m.start() for m in re.finditer(p2, i)]
        tar_index3 = [m.start() for m in re.finditer(p4, i)]
        li=list(i)
        for index in tar_index1:
            li[index+1]=li[index+1].lower()
        for index in tar_index2:
            li[index-1]=li[index-1].lower()
        for index in tar_index3:
            li[index-1]=li[index-1].lower()
        li="".join(li)       
        no_c=re.sub(p4,'',li)
        no_c=re.sub(p2,'',no_c)
        no_c=re.sub(p3,'',no_c)
        no_c=re.sub(p1,'',no_c)
        phospep.append(no_c[1:-1])
    df_maxquant['phosphopep']=pd.Series(phospep).values  
    phosprob=df_maxquant['Phospho (STY) Probabilities'].tolist()
    phosnum=df_maxquant['Phospho (STY)'].tolist()
    accept=[]
    for i in range(len(phosprob)):
        p = re.compile(r'[(](.*?)[)]', re.S)
        tar_index1 = [m.start() for m in re.finditer(p, phosprob[i])]
        all_num=re.findall(p,phosprob[i])
        hi_index=[x for x in range(len(all_num)) if float(all_num[x])>0.75]
        if len(hi_index)>=phosnum[i]:
            accept.append('y')
        else:
            accept.append('n')
    df_maxquant['Acceptance']=pd.Series(accept).values  
    df_maxquant.columns=['scan','charge','sequence','length','protein','modified sequence',
                         'phospho_prob','phospho_num','retention time','phosphopep','Acceptance']
    return(df_maxquant)


#standardized Michi
def standardized_Michi_psm(data):
    data.index=data.Spectrum
    Michi = data[['Peptide','Phospho Site Localization','Entry Name','Charge','Modified Peptide','Retention','Assigned Modifications','PeptideProphet Probability']]
    # Michi_all.index=Michi_all.Spectrum
    # Michi = Michi_all[['Peptide','Phospho Site Localization','Entry Name','Charge','Modified Peptide','Retention','Assigned Modifications','PeptideProphet Probability']]
    idx=[]
    for i in range(len(Michi)):
        if str(Michi['Phospho Site Localization'].iloc[i])=='nan':
            continue
        else:
            idx.append(i)
    Michi=Michi.iloc[idx]
    exp =[]
    scannum=[]
    for i in Michi.index:
        exp.append(str.split(i,'.')[0])
        scannum.append(str.split(i,'.')[1])
    scannum=[int(i) for i in scannum]
    row1=pd.Series(exp,list(Michi.index))
    row2=pd.Series(scannum,list(Michi.index))
    Michi['experiment'] = pd.Series(exp).values
    Michi['scan'] = pd.Series(scannum).values
    Michi.index=Michi.experiment

    phoprob_list=Michi["Phospho Site Localization"].tolist()
    pho_list=Michi["Modified Peptide"].tolist()
    Mich_phosphopep = []
    for i in range(len(pho_list)):
        p = re.compile(r'\[(.*?)\]', re.S)
        tar_index = [m.start() for m in re.finditer(p, pho_list[i])]
        if len(tar_index)>0:
            li=list(pho_list[i])
            for index in tar_index:
                    li[int(index)-1]=li[int(index)-1].lower()
            li="".join(li)       
            no_c=re.sub(p,'',li)
            no_c=list(no_c)
            no_c=''.join(no_c[1:])
            Mich_phosphopep.append(no_c)
    accept=[]
    for i in range(len(phoprob_list)):
        p = re.compile(r'[(](.*?)[)]', re.S)
        tar_index = [m.start() for m in re.finditer(p, phoprob_list[i])]
        all_num=re.findall(p,phoprob_list[i])
        hi_index=[x for x in range(len(all_num)) if float(all_num[x])>0.75]
        num = Mich_phosphopep[i].count('s')+Mich_phosphopep[i].count('t')+Mich_phosphopep[i].count('y')
        if len(hi_index)>=num:
            accept.append('y')
        else:
            accept.append('n')
    Michi['phosphopep']=pd.Series(Mich_phosphopep).values
    Michi['Acceptance']=pd.Series(accept).values
    Michi.columns=['sequence','phospho_prob','protein','charge',
                 'modified sequence','retention time','assigned modifications','pep_prob','experiment','scan','phosphopep','Acceptance']
    return(Michi)

#standardized PNNL
def standardized_PNNL_report(data):
    PNNL=data[["Charge",'site','maxAScore','Dataset','RefSeq.D','Scan','cleanSeq','isDecoy','pepSeq','RetentionTime','MSGFScore','exp']]
#     PNNL_exp=[]
#     index=[]
#     for i in range(len(PNNL_report)):
#         item = PNNL_report.Dataset[i].split('_')
#         if item[5] == 'B1S1':
#             PNNL_exp.append('01CPTAC_UCEC_P_PNNL_20170922_B1S1_'+item[6])
#             index.append(i)
#     PNNL=PNNL_report.loc[index]
#     PNNL['exp'] = pd.Series(PNNL_exp).values
    seq_PNNL=[]
    pho_list=list(PNNL['cleanSeq'])
    cleanseq=[]
    for i in range(len(pho_list)):
        p = re.compile(r'\*', re.S)
        tar_index = [m.start() for m in re.finditer(p, pho_list[i])]
        li=list(pho_list[i])
        if len(tar_index)>0:
            for index in tar_index:
                li[index-1]=li[index-1].lower()
            li="".join(li)       
            no_c=re.sub(p,'',li)
            cleanseq.append(pho_list[i])
            seq_PNNL.append(no_c)
    PNNL['phosphopep'] = pd.Series(seq_PNNL).values
    PNNL.drop_duplicates(keep=False,inplace=True) 
    PNNL.columns=["charge",'site','ascore','dataset','protein','scan','cleanSeq','isDecoy','sequence','retention time','MSGFScore','exp','phosphopep']
    PNNL.index=PNNL.exp
    return(PNNL)

#standardized CDAP
def standardized_CDAP_psm(data):
    CDAP=data[data.FullyLocalized.isin(['Y','N'])]
    CDAP=CDAP[['ScanNum','PhosphoRSPeptide','OriginalCharge','nPhospho','FullyLocalized',
               'PeptideSequence','MSGFScore','RTAtPrecursorHalfElution','Protein']]
    exp =[]
    for i in CDAP.index:
        exp.append(str.split(i,'.')[0])
    phoprob_list=CDAP["PhosphoRSPeptide"].tolist()
    seq=[]
    for i in phoprob_list:
        p = re.compile(r'\[(.*?)\]', re.S)
        no_c=re.sub(p,'',i)
        seq.append(no_c)
    CDAP['seq']=pd.Series(seq).values
    CDAP['exp']=pd.Series(exp).values
    pho_list=CDAP["PeptideSequence"].tolist()
    for i in range(len(pho_list)):
        p0 = re.compile(r'M\+15', re.S)
        p1 = re.compile(r'S\+79', re.S)
        p2 = re.compile(r'T\+79', re.S)
        p3 = re.compile(r'Y\+79', re.S)
        p0_ = re.compile(r'M', re.S)
        p1_= re.compile(r'S', re.S)
        p2_= re.compile(r'T', re.S)
        p3_= re.compile(r'Y', re.S)
        tar_index_m = [m.start() for m in re.finditer(p0, pho_list[i])]
        tar_index_M = [m.start() for m in re.finditer(p0_,pho_list[i])]
        CDAP_index_M=[m.start() for m in re.finditer(p0_, seq[i])]
        tar_index_s = [m.start() for m in re.finditer(p1, pho_list[i])]
        tar_index_S = [m.start() for m in re.finditer(p1_,pho_list[i])]
        CDAP_index_S=[m.start() for m in re.finditer(p1_, seq[i])]
        tar_index_t = [m.start() for m in re.finditer(p2, pho_list[i])]
        tar_index_T = [m.start() for m in re.finditer(p2_,pho_list[i])]
        CDAP_index_T=[m.start() for m in re.finditer(p2_, seq[i])]
        tar_index_y = [m.start() for m in re.finditer(p3, pho_list[i])]
        tar_index_Y = [m.start() for m in re.finditer(p3_,pho_list[i])]
        CDAP_index_Y=[m.start() for m in re.finditer(p3_, seq[i])]
        for t in range(len(tar_index_M)):
            if tar_index_M[t] in tar_index_m:
                seq[i]=list(seq[i])
                seq[i][int(CDAP_index_M[t])]=seq[i][int(CDAP_index_M[t])].lower()
                seq[i]=''.join(seq[i])
        for t in range(len(tar_index_S)):
            if tar_index_S[t] in tar_index_s:
                seq[i]=list(seq[i])
                seq[i][int(CDAP_index_S[t])]=seq[i][int(CDAP_index_S[t])].lower()
                seq[i]=''.join(seq[i])
        for t in range(len(tar_index_T)):
            if tar_index_T[t] in tar_index_t:
                seq[i]=list(seq[i])
                seq[i][int(CDAP_index_T[t])]=seq[i][int(CDAP_index_T[t])].lower()
                seq[i]=''.join(seq[i])
        for t in range(len(tar_index_Y)):
            if tar_index_Y[t] in tar_index_y:
                seq[i]=list(seq[i])
                seq[i][int(CDAP_index_Y[t])]=seq[i][int(CDAP_index_Y[t])].lower()
                seq[i]=''.join(seq[i])
    CDAP['phosphopep']=pd.Series(seq).values
    p1 = re.compile(r'\(.*\)', re.S)
    Pro=[]
    for i in range(len(CDAP)):
        pro=CDAP.Protein.iloc[i]
        Pro.append(';'.join([re.sub(p1,'',j) for j in pro.split(';')]))
    CDAP['protein']=pd.Series(Pro).values
    CDAP.columns=['scan','phospho_prob','charge','phospho_num','FullyLocalized','modified sequence','MSGFScore','retention time','Proteins','sequence','exp','phosphopep','protein']
    CDAP.index=CDAP.exp
    return(CDAP)

#Convert overlap, Michi, PNNL, df_maxquant_T to psm.txt
#Carbamidomethyl[C]
#TMT6plex[AnyN-term]
#TMT6plex[K]
def convert_to_psmtxt(data,name):
    data_pdeep3=data[['scan','charge','sequence','phosphopep']]
    data_pdeep3=data_pdeep3.drop_duplicates()
    pep_list=data_pdeep3['phosphopep'].tolist()
    modinf=[]
    for pep in pep_list:
        modinf_each=[]
        mod='0,TMT6plex[AnyN-term];'
        modinf_each.append(mod)
        if 'C' in pep:
            p_s = re.compile(r'C', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('C')+1)+','+'Carbamidomethyl[C];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'Carbamidomethyl[C];')
                    modinf_each.append(mod)
        if 'K' in pep:
            p_s = re.compile(r'K', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('K')+1)+','+'TMT6plex[K];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'TMT6plex[K];')
                    modinf_each.append(mod)
        if 's' in pep:
            p_s = re.compile(r's', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('s')+1)+','+'Phospho[S];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'Phospho[S];')
                    modinf_each.append(mod)
        if 'm' in pep:
            p_s = re.compile(r'm', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('m')+1)+','+'Oxidation[M];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'Oxidation[M];')
                    modinf_each.append(mod)
        if 'y' in pep:
            p_s = re.compile(r'y', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('y')+1)+','+'Phospho[Y];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'Phospho[Y];')
                    modinf_each.append(mod)       
        if 't' in pep:
            p_s = re.compile(r't', re.S)
            tar_index_p = [m.start() for m in re.finditer(p_s, pep)]
            if len(tar_index_p)==1:
                mod=str(str(pep.index('t')+1)+','+'Phospho[T];')
                modinf_each.append(mod)
            else:
                for inde in tar_index_p:
                    mod=str(str(inde+1)+','+'Phospho[T];')
                    modinf_each.append(mod)
        modinf_each="".join(modinf_each)
        modinf.append(modinf_each)
    data_pdeep3['modinfo']=pd.Series(modinf).values
    #outputStr = ''
    #for z in range(len(data_pdeep3)): 
    #    if (pd.isna(data_pdeep3.modinfo.iloc[z])==True):
    #        outputStr += str(data_pdeep3.index[z]) + "\t" + str(data_pdeep3['scan'].iloc[z]) + "\t"+  str(data_pdeep3['sequence'].iloc[z])+"\t"+"\t"+str(data_pdeep3.charge.iloc[z]) +'\n'
    #    else:
    #        outputStr += str(data_pdeep3.index[z]) + "\t" + str(data_pdeep3['scan'].iloc[z]) + "\t"+  str(data_pdeep3['sequence'].iloc[z])+"\t"+str(data_pdeep3.modinfo.iloc[z])+"\t"+str(data_pdeep3.charge.iloc[z]) +'\n'
    #file=''.join(['../data/pDeep3_data/psmtxt/8_5_2021/UCEC/UCEC_'+name+'_exp1_psm.txt'])
    #with open(file, "w") as outputter:
    #    outputter.write('raw_name'+"\t"+'scan'+"\t"+"peptide"+"\t"+'modinfo'+'\t'+'charge'+'\n')
    #    outputter.write(outputStr)

    #    outputStr = ''
    #for z in range(len(data_pdeep3)): 
    #    if (pd.isna(data_pdeep3.modinfo.iloc[z])==True):
    #        outputStr += str(data_pdeep3.index[z]) + "\t" + str(data_pdeep3['scan'].iloc[z]) + "\t"+  str(data_pdeep3['sequence'].iloc[z])+"\t"+"\t"+str(data_pdeep3.charge.iloc[z]) +'\n'
    #    else:
    #        outputStr += str(data_pdeep3.index[z]) + "\t" + str(data_pdeep3['scan'].iloc[z]) + "\t"+  str(data_pdeep3['sequence'].iloc[z])+"\t"+str(data_pdeep3.modinfo.iloc[z])+"\t"+str(data_pdeep3.charge.iloc[z]) +'\n'
    #file=''.join(['../data/PXD015284/pDeep2_data/psmtxt/'+name+'_exp1_psm.txt'])
    #file=''.join(['../data/PXD007145/pDeep2_data/psmtxt/'+name+'_exp1_psm.txt'])
    #with open(file, "w") as outputter:
    #    outputter.write('raw_name'+"\t"+'scan'+"\t"+"peptide"+"\t"+'modinfo'+'\t'+'charge'+'\n')
    #    outputter.write(outputStr)  



def prepare_phi_vali_AutoRT(data,path):
    rt_michi_only=[]
    phosphopep_michi_only=[]
    idx_michi_only=[]
    for i in range(len(data)):
        if [(data.index[i],int(data.scan.iloc[i]),
             data.phosphopep.iloc[i])][0] in Michi_only_PXD015284:
            seq=data['phosphopep'].iloc[i].replace('m','M')
            seq=seq.replace('s','S')
            seq=seq.replace('t','T')
            seq=seq.replace('y','Y')
            if seq in seq_upper:
                continue
            if len(list(data['phosphopep'].iloc[i]))<=48:
                rt_michi_only.append(data['retention time'].iloc[i]/60)
                phosphopep_michi_only.append(data['phosphopep'].iloc[i])
                idx_michi_only.append(i)
                print(i)
    michi_Autort=data.iloc[idx_michi_only]
    michi_Autort.to_pickle(path+'/michi_Autort.pkl')
    cp_michi_only1=list(zip(phosphopep_michi_only,rt_michi_only))
    cp_michi_only1=list(set(cp_michi_only1))
    phosphopep_michi=[]
    rt_michi=[]
    for i in range(len(cp_michi_only1)):
        pep=cp_michi_only1[i][0].replace('m','1')
        pep=pep.replace('s','2')
        pep=pep.replace('t','3')
        pep=pep.replace('y','4')
        if list(pep)[0]=='n':
            pep=''.join(pep[1:])
            phosphopep_michi.append(pep)
        else:
            phosphopep_michi.append(pep)
        rt_michi.append(cp_michi_only1[i][1])
    outputStr = ''
    for z in range(len(phosphopep_michi)): 
        outputStr += str(phosphopep_michi[z]) + "\t" + str(rt_michi[z]) +'\n'
    file=''.join([path+'/Michi_prediction_2.tsv'])
    with open(file, "w") as outputter:
        outputter.write('x'+"\t"+'y'+'\n')
        outputter.write(outputStr)

def prepare_MQ_vali_AutoRT(data,path):
    #MaxQuant prediction, peptide < 48 AA
    rt_maxquant_only=[]
    phosphopep_maxquant_only=[]
    idx_maxquant_only=[]
    for i in range(len(data)):
        if [(data.index[i],data['scan'].iloc[i],
            data['phosphopep'].iloc[i])][0] in maxquant_only_PXD015284:
            seq=data['phosphopep'].iloc[i].replace('m','M')
            seq=seq.replace('s','S')
            seq=seq.replace('t','T')
            seq=seq.replace('y','Y')
            if seq in seq_upper:
                continue
            if len(list(data['phosphopep'].iloc[i]))<=48:
                rt_maxquant_only.append(data['retention time'].iloc[i])
                phosphopep_maxquant_only.append(data['phosphopep'].iloc[i])
                idx_maxquant_only.append(i)
                print(i)
    maxquant_Autort=data.iloc[idx_maxquant_only]
    maxquant_Autort.to_pickle(path+'/maxquant_Autort.pkl')
    cp_maxquant_only1=list(zip(phosphopep_maxquant_only,rt_maxquant_only))
    cp_maxquant_only1=list(set(cp_maxquant_only1))
    phosphopep_MaxQuant=[]
    rt_MaxQuant=[]
    for i in range(len(cp_maxquant_only1)):
        pep=cp_maxquant_only1[i][0].replace('m','1')
        pep=pep.replace('s','2')
        pep=pep.replace('t','3')
        pep=pep.replace('y','4')
        if list(pep)[0]=='n':
            pep=''.join(pep[1:])
            phosphopep_MaxQuant.append(pep)
        else:
            phosphopep_MaxQuant.append(pep)
        rt_MaxQuant.append(cp_maxquant_only1[i][1])
    outputStr = ''
    for z in range(len(phosphopep_MaxQuant)): 
        outputStr += str(phosphopep_MaxQuant[z]) + "\t" + str(rt_MaxQuant[z]) +'\n'
    file=''.join([path+'/MaxQuant_prediction.tsv'])
    with open(file, "w") as outputter:
        outputter.write('x'+"\t"+'y'+'\n')
        outputter.write(outputStr)       
