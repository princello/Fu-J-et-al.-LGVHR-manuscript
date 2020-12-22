import numpy as np
import pandas as pd

def treatcr(tcr_file):
    tcr_data = pd.read_csv(tcr_file)
    tcr_data = tcr_data[(tcr_data['chain'] == 'TRB') & (tcr_data['v_gene'] != 'None') & (tcr_data['j_gene'] != 'None') & (tcr_data['raw_consensus_id'] != 'None')]
    tcr_data = tcr_data.drop_duplicates(subset='barcode',keep=False)
    tcr_data['#clones'] = tcr_data['raw_clonotype_id'].map(tcr_data['raw_clonotype_id'].value_counts())
    tcr_data['V1'] = tcr_data['v_gene'].str.split('-').str[0].str[4:].str.zfill(2)
    tcr_data['V2'] = tcr_data['v_gene'].str.split('-').str[1].str.zfill(2).fillna('01')
    tcr_data['J1'] = tcr_data['j_gene'].str.split('-').str[0].str[4:].str.zfill(2)
    tcr_data['J2'] = tcr_data['j_gene'].str.split('-').str[1].str.zfill(2).fillna('01')
    return tcr_data[['barcode','raw_clonotype_id', 'cdr3', 'cdr3_nt','#clones','V1','V2','J1','J2']]

def treathvg(bulk_file, sample_name, hvg_file):
    bulk_data = pd.read_csv(bulk_file, sep='\t')
    hvg = [i[:-1] for i in open(hvg_file,'r').readlines()]
    #sample_data = bulk_data[bulk_data['sample_name'] == sample_name]
    #hvg_data = sample_data[(sample_data['rearrangement'].isin(hvg)) & (sample_data['v_resolved'].str[:2] != 'un') & (sample_data['j_resolved'].str[:2] != 'un')]
    sample_data = bulk_data[bulk_data['sample_name'].isin(sample_name)]
    
    hvg_data = sample_data[(sample_data['rearrangement'].isin(hvg)) & (sample_data['v_resolved'].str[:2] != 'un') & (sample_data['j_resolved'].str[:2] != 'un')]
    hvg_data = hvg_data.drop_duplicates(subset='rearrangement',keep='first')
    hvg_data.insert(1, 'V1' , hvg_data['v_resolved'].str.split('-').str[0].str[5:])
    hvg_data.insert(1, 'V2' , hvg_data['v_resolved'].str.split('-').str[1].str[:2].fillna('01'))
    hvg_data.insert(1, 'J1' , hvg_data['j_resolved'].str.split('-').str[0].str[5:])
    hvg_data.insert(1, 'J2' , hvg_data['j_resolved'].str.split('-').str[1].str[:2].fillna('01'))
    #new add
    return hvg_data[['rearrangement','amino_acid','V1','V2','J1','J2','productive_frequency']]

def treathvg_all(bulk_file):#,bulk_file_nr):
    bulk_data = pd.read_csv(bulk_file, sep='\t')
    #bulk = pd.read_csv(bulk_file_nr, sep='\t')
    #hvg = list(bulk[bulk['Present In']>=1]['Rearrangement'])
    hvg_data = bulk_data[(bulk_data['v_resolved'].str[:2] != 'un') & (bulk_data['j_resolved'].str[:2] != 'un')]  #(bulk_data['rearrangement'].isin(hvg)) & 
    hvg_data = hvg_data.drop_duplicates(subset='rearrangement',keep='first')
    hvg_data.insert(1, 'V1' , hvg_data['v_resolved'].str.split('-').str[0].str[5:])
    hvg_data.insert(1, 'V2' , hvg_data['v_resolved'].str.split('-').str[1].str[:2].fillna('01'))
    hvg_data.insert(1, 'J1' , hvg_data['j_resolved'].str.split('-').str[0].str[5:])
    hvg_data.insert(1, 'J2' , hvg_data['j_resolved'].str.split('-').str[1].str[:2].fillna('01'))
    return hvg_data[['rearrangement','amino_acid','V1','V2','J1','J2','productive_frequency']]

def link(tcr_data,hvg_data):
    indexs = []
    ans = []
    #hvg_sel = pd.DataFrame(columns=hvg_data.columns) 
    for index, row in tcr_data.iterrows():
        record = hvg_data[hvg_data['rearrangement'].str.count(row['cdr3_nt']) == 1]
        if len(record) != 0:
            ans_temp = []
            for ind, rec in record.iterrows():
                if rec['V1']  == row['V1'] and rec['V2'] == row['V2'] and rec['J1']  == row['J1'] and rec['J2']  == row['J2']:
                    indexs.append(index)
                    #hvg_sel = hvg_sel.append(rec)
                    ans_temp.append(rec['rearrangement'])
            if ans_temp != []:
                ans.append(ans_temp)
    ans_list = ['; '.join(temp) for temp  in ans]
    indexs = list(set(indexs))
    indexs.sort()
    tcr_sel = tcr_data.loc[indexs]
    tcr_sel['Adaptive NT sequences'] = ans_list
    return tcr_sel

def addsample(pt_output, bulk_file_nr):
    bulk = pd.read_csv(bulk_file_nr,sep = '\t')
    rear_list = list(pt_output['Adaptive NT sequences'])
    samples = []
    for rear in rear_list:
        record = []
        for rea in rear.split('; '):
            record.extend(list(bulk.columns[list(np.where(bulk[bulk['Rearrangement'] == rea].values[0][1:]>0)[0]+1)][2:]))
        samples.append(list(set(record)))
    ans_list = ['; '.join(temp) for temp  in samples]
    pt_output['samples'] = ans_list
    return pt_output


def classcell(bulk_file,sample_name,hvg_head,hvg_file,tcr_file):
    tcr_data = treatcr(tcr_file)
    output = []
    j = [2,3,0,1,6,7,4,5]
    for i in range(len(hvg_file)):
        hvg_data = treathvg(bulk_file, [sample_name[i],sample_name[j[i]]], hvg_head+hvg_file[i])
        tcr_sel = link(tcr_data,hvg_data)
        output.append(tcr_sel)
        print (len(tcr_sel))
    return output

def classcell(bulk_file,sample_name,hvg_head,hvg_file,tcr_file):
    tcr_data = treatcr(tcr_file)
    output = []
    j = [2,3,0,1,6,7,4,5]
    for i in range(len(hvg_file)):
        hvg_data = treathvg(bulk_file, sample_name, hvg_head+hvg_file[i])
        tcr_sel = link(tcr_data,hvg_data)
        output.append(tcr_sel)
        print (len(tcr_sel))
    return output


def classcell_all(bulk_file,bulk_file_nr, tcr_file):
    tcr_data = treatcr(tcr_file)
    hvg_data = treathvg_all(bulk_file)#,bulk_file_nr)
    tcr_sel = link(tcr_data,hvg_data)
    return tcr_sel


labels = ['CD4_HvG', "CD4_H'vG",'CD4_nonHvG', "CD4_nonH'vG", 'CD8_HvG', "CD8_H'vG", 'CD8_nonHvG', "CD8_nonH'vG"]
#prepost = ['CD4_HvG','CD4_nonHvG', 'CD8_HvG','CD8_nonHvG', "preTx_Unmappable" ,"CD4_H'vG","CD4_nonH'vG","CD8_H'vG", "CD8_nonH'vG",'postTx_Unmappable']
#cell_file = 'mj006_cells.csv'
#labels = ['CD4_HvG', "CD4_GvH",'CD4_nonHvG', "CD4_nonGvH", 'CD8_HvG', "CD8_GvH", 'CD8_nonHvG', "CD8_nonGvH"]

def makelabel(cell_file,output,labels,tcr_file,post=True):
    cells = pd.read_csv(cell_file)
    cells['Unnamed: 0'] = cells['Unnamed: 0'].str.split('-1').str[0] #alt
    tcr_data = treatcr(tcr_file)
    barcode = list(tcr_data['barcode'].str.split('-').str[0])
    clonotype = list(tcr_data['raw_clonotype_id'])
    dct = dict(zip(barcode,clonotype))
    cells['clonotype'] = ['None'] * len(cells)
    for barc in barcode:
        cells.loc[cells['Unnamed: 0'] == barc, 'clonotype'] = dct[barc]
    cells['pre_label'] = ['Un'] * len(cells)
    if post == False:
        for i in range(8):
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'pre_label'] = labels[i]
    else:
        cells['post_label'] = ['Un'] * len(cells)
        cells['pre_post_label'] = ['Un'] * len(cells)
        for i in [6,4,2,0]:
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'pre_label'] = labels[i]
        for i in [7,5,3,1]:
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'post_label'] = labels[i]
        cells['pre_post_label'] = cells['pre_label'] +'; ' +cells['post_label']
        cells.loc[~cells['pre_post_label'].isin(sel_labels), 'pre_post_label'] = 'Others'
    cells['Unnamed: 0'] = [i+'-1' for i in list(cells['Unnamed: 0'])]
    return cells

def makelabel_cells(cells,output,labels,tcr_file,post=True):
    tcr_data = treatcr(tcr_file)
    barcode = list(tcr_data['barcode'].str.split('-').str[0])
    clonotype = list(tcr_data['raw_clonotype_id'])
    dct = dict(zip(barcode,clonotype))
    cells['clonotype'] = ['None'] * len(cells)
    for barc in barcode:
        cells.loc[cells['Unnamed: 0'] == barc, 'clonotype'] = dct[barc]
    cells['pre_label'] = ['Un'] * len(cells)
    if post == False:
        for i in range(4):
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'pre_label'] = labels[i]
    else:
        cells['post_label'] = ['Un'] * len(cells)
        cells['pre_post_label'] = ['Un'] * len(cells)
        for i in [6,4,2,0]:
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'pre_label'] = labels[i]
        for i in [7,5,3,1]:
            cells.loc[cells['Unnamed: 0'].isin(list(output[i].barcode.str.split('-').str[0])), 'post_label'] = labels[i]
        cells['pre_post_label'] = cells['pre_label'] +'; ' +cells['post_label']
        cells.loc[~cells['pre_post_label'].isin(sel_labels), 'pre_post_label'] = 'Others'
    return cells



#pre_post label
sel_labels = ["CD4_HvG; CD4_H'vG", 
"Un; CD4_H'vG",
"CD4_HvG; Un",
"CD4_HvG; CD4_nonH'vG",
"CD4_nonHvG; CD4_nonH'vG",

"CD8_HvG; CD8_H'vG",
"Un; CD8_H'vG",
"CD8_HvG; Un",
"CD8_HvG; CD8_nonH'vG",
"CD8_nonHvG; CD8_nonH'vG",

"Un; Un"]

#pre_post separate
labels = ['CD4_HvG', "CD4_H'vG",'CD4_nonHvG', "CD4_nonH'vG", 'CD8_HvG', "CD8_H'vG", 'CD8_nonHvG', "CD8_nonH'vG"]

def table2(output,labels):
    table = pd.DataFrame(columns=output[0].columns)
    table.insert(10, 'Pre', [])
    table.insert(11, 'Post',[])
    for i in [7,6,5,4,3,2,1,0]:
        lab = labels[i]
        output_data = output[i]
        for ind, out in output_data.iterrows():
            out = list(out)
            out[0] = out[0].split('-')[0]
            if lab.count("'") == 0:
                if out[0] not in list(table['barcode']):
                    out.append(lab)
                    out.append('Un')
                    table.loc[-1] = out
                    table.index = table.index + 1
                    table = table.sort_index()
                else:
                    table.loc[table['barcode'] == out[0], 'Pre'] = lab
                    #table.loc[table['barcode'] == out[0], '#clones'] = out[4]
                    #table.loc[table['barcode'] == out[0], 'Adaptive NT sequences'] = out[9]
            else:
                if out[0] not in list(table['barcode']):
                    out.append('Un')
                    out.append(lab)
                    table.loc[-1] = out
                    table.index = table.index + 1
                    table = table.sort_index()
                else:
                    table.loc[table['barcode'] == out[0], 'Post'] = lab
                    #table.loc[table['barcode'] == out[0], '#clones'] = out[4]
                    #table.loc[table['barcode'] == out[0], 'Adaptive NT sequences'] = out[9]
    return table
    
def count_clone(cl,gvh=False):
    for i in range(0,8,2):
        la = labels[i]
        print (la, len(cl[cl['pre_label'] == la]), len(set(cl[cl['pre_label'] == la]['clonotype'])))
    print ('pre_un', len(cl[(cl['pre_label'] == 'Un') & (cl['clonotype'] !='None')]), len(set(cl[(cl['pre_label'] == 'Un') & (cl['clonotype'] !='None')]['clonotype'])))
    if gvh == False:
        for i in range(1,9,2):
            la = labels[i]
            print (la, len(cl[cl['post_label'] == la]), len(set(cl[cl['post_label'] == la]['clonotype'])))
        print ('post_un', len(cl[(cl['post_label'] == 'Un') & (cl['clonotype'] !='None')]), len(set(cl[(cl['post_label'] == 'Un') & (cl['clonotype'] !='None')]['clonotype'])))
    else:
        for i in range(1,9,2):
            la = labels[i]
            print (la, len(cl[cl['pre_label'] == la]), len(set(cl[cl['pre_label'] == la]['clonotype'])))
        print ('pre_un', len(cl[(cl['pre_label'] == 'Un') & (cl['clonotype'] !='None')]), len(set(cl[(cl['pre_label'] == 'Un') & (cl['clonotype'] !='None')]['clonotype'])))


def table2(output,labels):
    table = pd.DataFrame(columns=output[0].columns)
    table.insert(10, 'HvG', [])
    table.insert(11, 'GvH',[])
    for i in [7,6,5,4,3,2,1,0]:
        lab = labels[i]
        output_data = output[i]
        for ind, out in output_data.iterrows():
            out = list(out)
            out[0] = out[0].split('-')[0]
            if lab.count("GvH") == 0:
                if out[0] not in list(table['barcode']):
                    out.append(lab)
                    out.append('Un')
                    table.loc[-1] = out
                    table.index = table.index + 1
                    table = table.sort_index()
                else:
                    table.loc[table['barcode'] == out[0], 'HvG'] = lab
                    #table.loc[table['barcode'] == out[0], '#clones'] = out[4]
                    #table.loc[table['barcode'] == out[0], 'Adaptive NT sequences'] = out[9]
            else:
                if out[0] not in list(table['barcode']):
                    out.append('Un')
                    out.append(lab)
                    table.loc[-1] = out
                    table.index = table.index + 1
                    table = table.sort_index()
                else:
                    table.loc[table['barcode'] == out[0], 'GvH'] = lab
                    #table.loc[table['barcode'] == out[0], '#clones'] = out[4]
                    #table.loc[table['barcode'] == out[0], 'Adaptive NT sequences'] = out[9]
    return table


cell_file = 'hvg/mj001_cell_post.csv'
cl1 = makelabel(cell_file,pt15_post_output_2,labels,pt15_tcr_file, True)
count_clone(cl1)
cl1.to_csv('hvg/mj001_post.csv')

cell_file = 'hvg/mj002_cell_post.csv'
cl2 = makelabel(cell_file,pt13_iel_post_output_2,labels,pt13_iel_tcr_file,True)
count_clone(cl2)
cl2.to_csv('hvg/mj002_post.csv')

cell_file = 'hvg/mj003_cell_post.csv'
cl3 = makelabel(cell_file,pt13_lpl_post_output_2,labels,pt13_lpl_tcr_file,True)
count_clone(cl3)
cl3.to_csv('hvg/mj003_post.csv')

cell_file = 'hvg/mj005_cell_post.csv'
cl5 = makelabel(cell_file,pt14_post_output,labels,pt14_tcr_file,True)
count_clone(cl5)
cl5.to_csv('hvg/mj005_post.csv')

cell_file = 'hvg/mj006_cell_post.csv'
cl6 = makelabel(cell_file,pt21_post_output_2,labels,pt21_tcr_file,True)
count_clone(cl6)
cl6.to_csv('hvg/mj006_post.csv')

cell_file = 'hvg/mj008_cell_post.csv'
cl8 = makelabel(cell_file,pt4_iel_post_output,labels,pt4_iel_tcr_file,True)
count_clone(cl8)
cl8.to_csv('hvg/mj008_post.csv')

cell_file = 'hvg/mj009_cell_post.csv'
cl9 = makelabel(cell_file,pt4_lpl_post_output,labels,pt4_lpl_tcr_file,True)
count_clone(cl9)
cl9.to_csv('hvg/mj009_post.csv')


cell_file = 'int/cell/mj008_cells.csv'
cl8 = makelabel(cell_file,pt4_iel_post_output,labels,pt4_iel_tcr_file,True)
count_clone(cl8)
cl8.to_csv('int/cell/mj008_post.csv')

cell_file = 'int/cell/mj009_cells.csv'
cl9 = makelabel(cell_file,pt4_lpl_post_output,labels,pt4_lpl_tcr_file,True)
count_clone(cl9)
cl9.to_csv('int/cell/mj009_post.csv')



labels = ['CD4_HvG', "CD4_GvH",'CD4_nonHvG', "CD4_nonGvH", 'CD8_HvG', "CD8_GvH", 'CD8_nonHvG', "CD8_nonGvH"]
cell_file = 'hvg/mj014_cell_post.csv'
cl14 = makelabel(cell_file,pt16_post_output,labels,pt16_tcr_file,False)
count_clone(cl14,True)
cl14.to_csv('hvg/mj014_post.csv')



labels = ['CD4_HvG', "CD4_GvH",'CD4_nonHvG', "CD4_nonGvH", 'CD8_HvG', "CD8_GvH", 'CD8_nonHvG', "CD8_nonGvH"]
cell_file = 'hvg/mj013_cell_post.csv'
cl13 = makelabel(cell_file,pt19_post_output_2,labels,pt19_tcr_file,False)
count_clone(cl13,True)
cl13.to_csv('hvg/mj013_post.csv')

labels = ['CD4_HvG', "CD4_GvH",'CD4_nonHvG', "CD4_nonGvH", 'CD8_HvG', "CD8_GvH", 'CD8_nonHvG', "CD8_nonGvH"]
cell_file = 'hvg/mj012_cell_post.csv'
cl12 = makelabel(cell_file,pt18_post_output_2,labels,pt18_tcr_file,False)
count_clone(cl12,True)
cl12.to_csv('hvg/mj012_post.csv')


labels = ['CD4_HvG', "CD4_GvH",'CD4_nonHvG', "CD4_nonGvH", 'CD8_HvG', "CD8_GvH", 'CD8_nonHvG', "CD8_nonGvH"]
cell_file = 'hvg/mj010_cell_post.csv'
cl10 = makelabel(cell_file,pt18_post_output_10,labels,pt18_tcr_file,False)
count_clone(cl10,True)
cl10.to_csv('hvg/mj010_post.csv')


table = table2(pt15_post_output_2,labels)
table.to_csv('mj001_table2.tsv',sep='\t')

table = table2(pt13_iel_post_output_2,labels)
table.to_csv('mj002_table2.tsv',sep='\t') 
table = table2(pt13_lpl_post_output_2,labels)
table.to_csv('mj003_table2.tsv',sep='\t') 

table = table2(pt14_post_output,labels)
table.to_csv('mj005_table2.tsv',sep='\t') 

table = table2(pt21_post_output_2,labels)
table.to_csv('mj006_table2.tsv',sep='\t')

table = table2(pt4_iel_post_output,labels)
table.to_csv('mj008_table2.tsv',sep='\t') 
table = table2(pt4_lpl_post_output,labels)
table.to_csv('mj009_table2.tsv',sep='\t') 

table = table2(pt18_post_output,labels)
table.to_csv('mj010_table2.tsv',sep='\t')

table = table2(pt19_post_output,labels)
table.to_csv('mj011_table2.tsv',sep='\t')



cell = pd.read_csv('/Users/zw2595/Desktop/cells.csv') 


cells_mj001 = cell[cell['x'] == 'Pt15_POD1194']

cells_mj002 = cell[cell['x'] == 'Pt13_POD1032_IEL'] 

cells_mj003 = cell[cell['x'] == 'Pt13_POD1032_LPL'] 

cells_mj005 = cell[cell['x'] == 'Pt14_POD1764']
cells_mj006 = cell[cell['x'] == 'Pt21_POD626'] 
cells_mj007 = cell[cell['x'] == 'D251']


cells_mj001['Unnamed: 0'] = cells_mj001['Unnamed: 0'].str.split('_').str[0]
cmj1 = makelabel_cells(cells_mj001, pt15_post_output_2, labels, pt15_tcr_file) 

cells_mj002['Unnamed: 0'] = cells_mj002['Unnamed: 0'].str.split('_').str[0] 
cmj2 = makelabel_cells(cells_mj002, pt13_iel_post_output_2, labels, pt13_iel_tcr_file) 

cells_mj003['Unnamed: 0'] = cells_mj003['Unnamed: 0'].str.split('_').str[0]  
cmj3 = makelabel_cells(cells_mj003, pt13_lpl_post_output_2, labels, pt13_lpl_tcr_file)

cells_mj005['Unnamed: 0'] = cells_mj005['Unnamed: 0'].str.split('_').str[0]
cmj5 = makelabel_cells(cells_mj005, pt14_post_output, labels, pt14_tcr_file)

cells_mj006['Unnamed: 0'] = cells_mj006['Unnamed: 0'].str.split('_').str[0]  
cmj6 = makelabel_cells(cells_mj006, pt21_post_output_2, labels, pt21_tcr_file)

cells_mj007 = cell[cell['x'] == 'D251'] 
cells_mj007['Unnamed: 0'] = cells_mj007['Unnamed: 0'].str.split('-').str[0] 




clono = list(cmj1.clonotype)+list(cmj5.clonotype)+list(cmj2.clonotype)+list(cmj3.clonotype)+list(cmj6.clonotype)+['None']*len(cells_mj007)
pre = list(cmj1.pre_label)+list(cmj5.pre_label)+list(cmj2.pre_label)+list(cmj3.pre_label)+list(cmj6.pre_label)+['Un']*len(cells_mj007)
post = list(cmj1.post_label)+list(cmj5.post_label)+list(cmj2.post_label)+list(cmj3.post_label)+list(cmj6.post_label)+['Un']*len(cells_mj007)

pre_post = list(cmj1.pre_post_label)+list(cmj5.pre_post_label)+list(cmj2.pre_post_label)+list(cmj3.pre_post_label)+list(cmj6.pre_post_label)+['Un; Un']*len(cells_mj007)


cell['clonotype'] = clono
cell['pre_label'] = pre 
cell['post_label'] = post 
cell['pre_post_label'] = pre_post 




cell_file = 'int/cell/mj008_cells.csv'
cl8 = makelabel(cell_file,pt4_iel_post_output,labels,pt4_iel_tcr_file,True)
count_clone(cl8)
cl8.to_csv('int/cell/mj008.csv')

cell_file = 'int/cell/mj009_cells.csv'
cl9 = makelabel(cell_file,pt4_lpl_post_output,labels,pt4_lpl_tcr_file,True)
count_clone(cl9)
cl9.to_csv('int/cell/mj009.csv')


