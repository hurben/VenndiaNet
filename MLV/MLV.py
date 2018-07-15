import sys, os
import json
import numpy as np
import pandas as pd
from itertools import chain, combinations, compress
import FL_MLV_Methods

def powerset(iterable):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    xs = list(iterable)
    # note we return an iterator rather than a list
    return chain.from_iterable(combinations(xs,n) for n in range(len(xs)+1))



def load_state(fp):
  with open('work/%s/work.json' % fp,'r') as f:
    return json.load(f)

def save_state(j):
  fp = 'work/%s/work.json' % j['id']
  with open(fp,'w') as f:
    return json.save(j, f)

# input: directory (work.json, files)
# output: set information (set name - [gene name]), gene information (gene name - gene info)
def generate_venn(path):
  s = load_state(path)
  dfs = []
  for fp in s['files']:
    df_genes = pd.read_csv('work/%s/f/%s' % (s['id'], fp), sep='\t')
    # change column name of last one
    new_col_name = df_genes.columns.tolist()
    if ('TS_ID' not in new_col_name):
      new_col_name[0] = 'TS_ID'
    if ('GENESYMBOL' not in new_col_name):
      new_col_name[1] = 'GENESYMBOL'
    new_col_name[-1] = fp
    df_genes.columns = new_col_name
    dfs.append(df_genes)
  # join genes and p-value to single dataframe
  df_integrate = dfs[0]
  for df in dfs[1:]:
    df_integrate = pd.merge(df_integrate, df, on=['TS_ID', 'GENESYMBOL'], how='outer')
  # save integrated groups
  df_integrate.to_csv('work/%s/genes.txt' % path, sep='\t')
  # gather group information
  groups = {}     # exclude intersection
  groups_in = {}  # include intersection
  for x in powerset(s['files']):
    if (len(x) == 0):
      continue
    groups_in['__'.join(sorted(x))] = []
  for _idx,row in df_integrate.iterrows():
    # groups dict
    groups_set = sorted(row[row.notna()].index.tolist()[2:])
    groupname = '__'.join(groups_set)
    if (groupname not in groups):
      groups[groupname] = []
    groups[groupname].append(str(row['TS_ID']))
    # groups_in dict
    for gns in powerset(groups_set):
      if (len(gns) == 0):
        continue
      nm = str(row['TS_ID'])
      groups_in['__'.join(sorted(gns))].append(nm)
  # save group information (gene - group dataframe)
  gene_cond_data = []
  genenames = []
  for cond,genes in groups.items():
    for g in genes:
      genenames.append(g)
      gene_cond_data.append( (0,cond,'',0) )
  df = pd.DataFrame(gene_cond_data, index=genenames, columns=['id','cond','name','pvalue'])
  df.to_csv('work/%s/pvalue.txt'%s['id'])
  # save group information (each gene contains multiple group(cond))
  with open('work/%s/group.txt'%s['id'],'w') as f:
    #for (df,fn) in zip(dfs,s['files']):
    #  f.write('%s\t%d\t%s\n' % (fn, len(genes), ','.join(genes)))
    for cond,genes in groups_in.items():
      if (len(genes) <= 0):
        continue
      f.write('%s\t%d\t%s\n' % (cond, len(genes), ','.join(genes)))


# input: directory (genes.txt), select_set, subtract_set
# output: updated selection info (sel.txt)
def generate_sel(path, sel_set, sub_set):
  s = load_state(path)
  groups = {}
  with open('work/%s/group.txt'%s['id'],'r') as f:
    for l in f:
      x = l.split();
      groups[x[0]] = x[2].split(',')
  genes = set()
  for set_ in sel_set:
    genes.update( groups[set_] )
  for set_ in sub_set:
    genes = genes - set(groups[set_])
  # write new selected gene files
  # (filter out previous pvalue.txt)
  df = pd.read_csv('work/%s/pvalue.txt'%s['id'], index_col=0)
  print df
  df_sel = df[list(genes)]
  df_sel.to_csv('work/%s/pvalue_sel.txt'%s['id'],sep='\t')




def REMOVE_FC(gene):
  gene = gene.split('_')[0]
  return gene

def CHECK_UNABLE_GENES(condition_dict, ppi_gene_list):
  #CHECK unable genes
  unable_gene_list = []
  for condition in condition_dict.keys():
    deg_genes = condition_dict[condition]

    for gene in deg_genes:
      gene = REMOVE_FC(gene) #NEED
      if gene not in ppi_gene_list:
        unable_gene_list.append(gene)

  unable_gene_list = list(set(unable_gene_list))
  return unable_gene_list

def do_adj_mat_RWR(params, cut_off=0.2):
  ppi_dict = params.ppi_dict
  unable_gene_list = params.unable_gene_list
  condition_dict = params.condition_dict

  DEG_topology_dict = FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part1(condition_dict, ppi_dict, unable_gene_list)
  #DEG_topology_dict does not contain unable_gene
  condition_specific_gene_entry_dict = FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part2(DEG_topology_dict,condition_dict, "work/%s/.BH"%params.name, unable_gene_list)
  #condition_specific_gene_entry_dict does not contain unable_gene
  # ---
  # now make p0 vector
  intersection_gene_dict = FL_MLV_Methods.Create_intersection_dict(condition_dict, unable_gene_list)
  #intersection_dict does not include unable_gene
  FL_MLV_Methods.Venndianet_RWR().Create_RWR_p0_vector(condition_specific_gene_entry_dict, intersection_gene_dict, "work/%s/.BH"%params.name)
  # ---
  # run RWR
  condition_list = condition_specific_gene_entry_dict.keys()
  FL_MLV_Methods.Venndianet_RWR().Run_RWR("work/%s/.BH"%params.name, condition_list)

def do_adj_mat_neighbor(params, condition_dict, cut_off):
  for condition in condition_dict.keys():
    if '__' not in condition:
      deg_genes = condition_dict[condition]

      for deg_gene in deg_genes:
        if method == "neighbor":
          result_dict = FL_MLV_Methods.Venndianet_Neighbor().Condition_Exclusive(deg_gene, condition, condition_dict,ppi_dict, unable_gene_list, result_dict)

  return result_dict



def fill_params(params):
  # ---
  # reload data, and fill dataframe.
  # (tabbed splitted file, index exists.)
  for cond in condition_list.keys():
    df = pd.read_csv("work/%s/.BH.%s.rwr"%(params.name,cond), index_col=0, sep='\t')
    for gn,row in df.iterrows():
      params.df[gn]['pvalue'] = row.iloc[0]
  # ---
  # extract graph data from adjacency matrix.
  params.edges = {}
  for cond in condition_list.keys():
    with open("work/%s/.BH.%s.adj.matrix"%(params.name,cond),'r') as f:
      gns = f.readline().strip().split('\t')
      count_remove_mark = 1   # to remove edge duplication
      for l in f:
        x = l.strip.split('\t')
        gn, gmark = x[0], [_==1 for _ in x[1:]]
        for i in range(count_remove_mark):
          gmark[i] = False
        targets = compress(gns, gmark)
        for target in targets:
          params.edges.append( { 'source':gn, 'target':target, 'weight':1 } )
        count_remove_mark += 1
  # ---
  # fill node data
  params.genes = []
  for gn, row in params.df.iterrows():
    params.genes.append( { 'name': gn, 'group':row['cond'], 'value':row['pvalue'] } )

# df_pvalue: row genename, value is p-value
class FL_MLV_params(object):
  def __init__(self, n, df_pvalue):
    self.name = n
    self.df = df_pvalue
    self.genes = []   # graph component, to be filled
    self.edges = []   # graph component, to be filled
    self.gene_list = self.df.index.tolist()
    # fill condition_dict
    self.condition_dict = {}  # does it really necessary?
    for gn,row in self.df.iterrows():
      cond = row['cond']
      if (cond not in self.condition_dict):
        self.condition_dict[cond] = []
      self.condition_dict[cond].append(gn)
    # fill unable gene list
    ppi_profile_file = open('data/STRING_PPI_Mus_musculus_Symbol.txt','r')
    ppi_profile_file_readlines = ppi_profile_file.readlines()
    self.ppi_dict, self.ppi_gene_list = FL_MLV_Methods.StringDB_to_dict()
    self.unable_gene_list = CHECK_UNABLE_GENES(self.condition_dict, self.ppi_gene_list)

# input: directory (sel info; sel.txt)
# output: generated graph file (RWR; random walk probability)
def generate_graph(path):
  s = load_state(path)
  df = pd.read_csv('work/%s/pvalue_sel.txt'%s['id'],sep='\t',index_col=0)
  # process random walk and fill pvalue
  FL_MLV_params = FL_MLV_params(s['id'], df)
  do_adj_mat_RWR(params)
  fill_params(params)
  # save graph and df_node
  df.to_csv('work/%s/pvalue_sel.txt'%s['id'],sep='\t')
  with open('work/%s/graph.json'%s['id'],'w') as f:
    json.dump({'nodes':params.nodes,'edges':params.edges}, f)



# read venn diagram, return as dict.
def read_venn(id_):
  groups = {}
  with open('work/%s/group.txt'%id_,'r') as f:
    for l in f:
      x = l.split();
      groups[x[0]] = x[2].split(',')
  r = []
  for g,gns in groups.items():
    grps = g.split('__')
    r.append({
      'sets': grps,
      'figure': len(gns),
      'label': ' and '.join(grps),
      'size': len(gns)
    })
  return r


##
# Manual commands
##
if (__name__=="__main__"):
  tasks = {
    'generate_venn' : generate_venn,
    'generate_sel':generate_sel,
    'generate_graph': generate_graph
  }

  import argparse
  parser = argparse.ArgumentParser(description='Process some integers.')
  parser.add_argument('task', help='task name to proceed (%s)' % ','.join(tasks.keys()))
  parser.add_argument('path', help='destination path including work.json and data')
  parser.add_argument('--sel_set')
  parser.add_argument('--sub_set')
  parser.add_argument('--method', help='argument of RWR method')
  args = parser.parse_args()
  if (args.task == 'generate_sel'):
    generate_sel(args.path, args.sel_set.split(','), args.sub_set.split(','))
  else:
    tasks[ args.task ](args.path)
