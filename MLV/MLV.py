import sys, os
import json
import numpy as np
import pandas as pd
from itertools import chain, combinations, compress
import FL_MLV_Methods_v6 as FL_MLV_Methods

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

		if ('ID' not in new_col_name):
			new_col_name[0] = 'ID'
		if ('GENESYMBOL' not in new_col_name):
			new_col_name[1] = 'GENESYMBOL'

		new_col_name[-1] = fp
		df_genes.columns = new_col_name
		dfs.append(df_genes)
	
	# join genes and p-value to single dataframe
	df_integrate = dfs[0]
	for df in dfs[1:]:
		df_integrate = pd.merge(df_integrate, df, on=['ID', 'GENESYMBOL'], how='outer')
	# save integrated groups
	df_integrate.to_csv('work/%s/genes.txt' % path, sep='\t')

	# gather group information
	groups = {}			# exclude intersection
	groups_in = {}	# include intersection
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
		groups[groupname].append(str(row['GENESYMBOL']))
		# groups_in dict
		for gns in powerset(groups_set):
			if (len(gns) == 0):
				continue
			nm = str(row['GENESYMBOL'])
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
		for cond,genes in groups.items():

#		for cond,genes in groups_in.items():
#19.03.20
#why? groups_in

			if (len(genes) <= 0):
				continue
			f.write('%s\t%d\t%s\n' % (cond, len(genes), ','.join(genes)))

	# save .mlv file for further process
	mlv_file_dir = 'work/' + str(s['id']) + '/' + str(s['id']) +'.mlv'
	input_list_dir = 'work/' + str(s['id']) + '/input.list'
	input_list_txt = open(input_list_dir,'w')
	
	for fp in s['files']:
		input_list_txt.write('work/' + str(s['id']) + '/f/' + str(fp) + '\n')
	input_list_txt.close()
	
	make_mlv_cmd = 'python //var/www/htdocs/MLV_WEB/MLV/source_code/MAKE_MLV_FILE.py -i %s -o %s' % (input_list_dir, mlv_file_dir)
	os.system(make_mlv_cmd)

#MAKE RWR SEEDS
def select_venn_group(path, selection):
	# liver.HFD.deg.sig.gene.matrix,liver.HFD.deg.sig.gene.matrix__muscle.HFD.deg.sig.gene.matrix,muscle.HFD.deg.sig.gene.matrix
	s = load_state(path)

	sel_list = selection
	sel_list = str(sel_list)
	sel_list = sel_list.split(',')

	#pvalue.txt changed to seeds.txt
	#create seed gene list from condition
	mlv_file = 'work/%s/%s.mlv'% (s['id'],s['id'])
	mlv_condition_dict, mlv_gene_list = FL_MLV_Methods.MLV_Profile_to_dict(mlv_file)

	rwr_seed_txt = open('work/%s/seeds.txt'%s['id'], 'w')
	rwr_seed_txt.write('#Selected Conditions:\n')
	for condition in sel_list:
		rwr_seed_txt.write('#%s\n' % condition)

	#This process is essential
	#Because the order of CONDITION1_CONDITION2_CONDITION3 is different with dongwon
	for condition in sel_list:

		if '__' not in condition:

			if condition in mlv_condition_dict.keys():
				gene_list_from_condition = mlv_condition_dict[condition]

				for gene in gene_list_from_condition:
					rwr_seed_txt.write('%s\n' % gene)

		if '__' in condition:

			if condition in mlv_condition_dict.keys():

				if condition in mlv_condition_dict.keys():
					gene_list_from_condition = mlv_condition_dict[condition]

					for gene in gene_list_from_condition:
						rwr_seed_txt.write('%s\n' % gene)

			if condition not in mlv_condition_dict.keys():
				condition_sep_list = condition.split('__')
				condition_sep_list = sorted(condition_sep_list)

				for mlv_condition in mlv_condition_dict.keys():
					mlv_condition_sep_list = mlv_condition.split('__')
					mlv_condition_sep_list = sorted(mlv_condition_sep_list)

					if condition_sep_list == mlv_condition_sep_list:
						gene_list_from_condition = mlv_condition_dict[mlv_condition]

						for gene in gene_list_from_condition:
							rwr_seed_txt.write('%s\n' % gene)


def Bootstrap(path, selection):

	s = load_state(path)

	mlv_file = 'work/%s/%s.mlv'% (s['id'],s['id'])
	adj_file = 'work/%s/%s.adj.matrix'% (s['id'],s['id'])
	seed_file = 'work/%s/seeds.txt'% s['id']

	mlv_condition_dict, mlv_gene_list = FL_MLV_Methods.MLV_Profile_to_dict(mlv_file)

	StringDB_dict, PPI_gene_list = FL_MLV_Methods.StringDB_to_dict()
 	unable_gene_list = FL_MLV_Methods.Check_Unable_Genes(mlv_condition_dict, StringDB_dict)
	DEGtopology_dict, common_gene_list = FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part1(mlv_condition_dict, StringDB_dict, unable_gene_list)
	FL_MLV_Methods.Venndianet_RWR().Create_RWR_Condition_AdjacencyMatrix_part2(DEGtopology_dict, common_gene_list, adj_file, unable_gene_list)

	seed_gene_list = []
	print len(DEGtopology_dict.keys())
	open_seed_file = open(seed_file,'r')
	seed_file_readlines = open_seed_file.readlines()
	for i in range(len(seed_file_readlines)):
		read = seed_file_readlines[i]
		if '#' not in read[0]:
			gene = read.replace('\n','')
			seed_gene_list.append(gene)
	print "(DEBUG) Number of Seed genes: ", len(seed_gene_list)

	FL_MLV_Methods.Venndianet_RWR().Create_RWR_p0_vector(DEGtopology_dict, seed_gene_list, s['id'])
	#FL_MLV_Methods.Graph_Methods().Create_Graph_Info(s['id'])




def Run_RWR(path):

	s = load_state(path)
	work_id = s['id']

	#path_to_rwr_script = '//var/www/htdocs/MLV_WEB/MLV/RWR.R'
	path_to_rwr_script = '//var/www/htdocs/MLV/analysis/ppi/protocol/RWR.R'
	
	p0_vector_dir = 'work/%s/%s.p0.vector'% (work_id,work_id)
	adj_matrix_dir = 'work/%s/%s.adj.matrix'% (work_id,work_id)

	result_file = 'work/%s/RWR.result'% (work_id)
	Rscript_RWR_cmd = 'Rscript %s %s %s %s'% (path_to_rwr_script, adj_matrix_dir, p0_vector_dir, result_file)
	os.system(Rscript_RWR_cmd)


def After_RWR(path):

	s = load_state(path)
	work_id = s['id']

	rwr_dir = 'work/%s/RWR.result'% work_id
	parameter_info_dir = 'work/%s/parameter.txt'% work_id

	cutoff, sorting_order = FL_MLV_Methods.Parameter_Handling().Get_Parameter(parameter_info_dir)
	cutoff_gene_list = []

	if cutoff == 'None':
		FL_MLV_Methods.Graph_Methods().Create_Graph_Info_DEFAULT(s['id'])

	if cutoff != 'None':
		cutoff_gene_list = FL_MLV_Methods.Parameter_Handling().Cutoff_RWR(rwr_dir, cutoff, sorting_order)
		FL_MLV_Methods.Graph_Methods().Create_Graph_Info_CUTOFF(s['id'], cutoff_gene_list)

	return cutoff, sorting_order, cutoff_gene_list


#//var/www/htdocs/MLV/analysis/ppi/protocol/summerize_global_usage.py
def Create_Result_Summary(path, cutoff, sort_order):

	s = load_state(path)
	work_id = s['id']

	rwr_dir = 'work/%s/RWR.result'% work_id
	mlv_dir = 'work/%s/%s.mlv'% (work_id,work_id)
	output_dir = 'work/%s/%s.result.summary'% (work_id, work_id)
	output_txt = open(output_dir,'w')

	condition_id_dict, gene_condition_dict, condition_gene_dict = FL_MLV_Methods.Result_Summary().MLV_to_Dict(mlv_dir)
	rwr_dict = FL_MLV_Methods.Result_Summary().RWR_Result_to_Dict(rwr_dir)
	gene_rank_dict = {}

	if sort_order != 'ascend':
		sorted_rwr_dict = sorted(rwr_dict.items(), key=lambda x: x[1], reverse=True)
	if sort_order == 'ascend':
		sorted_rwr_dict = sorted(rwr_dict.items(), key=lambda x: x[1])

	if cutoff != 'None':
		result_length = len(rwr_dict.keys())
		cutoff = float(cutoff)
		cutoff = int(result_length * cutoff)


	#Header
	for conditionID in condition_gene_dict.keys():
		output_txt.write('ConditionID\tConditionName\n')
		output_txt.write('%s\t%s\n'% (conditionID, condition_id_dict[conditionID]))
	output_txt.write('Rank\tGene\tProbablity\tCondition\n')


	#Main
	for i in range(len(sorted_rwr_dict)):
		gene = sorted_rwr_dict[i][0]

		prob = rwr_dict[gene]
		output_txt.write('%s\t%s\t%s'% (i, gene, prob))

#modified 19.03.20
#issue : number is zero-based
#soluction : change it to 1 -based
#		gene_rank_dict[gene] = i
		gene_rank_dict[gene] = i + 1

		condition = gene_condition_dict[gene]
		output_txt.write('\t%s\n'% condition)

	output_txt.close()

	return gene_rank_dict, gene_condition_dict, condition_id_dict

def Create_Rank_Gene_dict(gene_rank_dict):

	rank_gene_dict = {}

	for gene in gene_rank_dict.keys():
		rank = gene_rank_dict[gene]
		rank_gene_dict[rank] = gene

	return rank_gene_dict


def Create_Rank_Gene_dict_without_seeds(rank_gene_dict, seed_gene_list):

	rank_idx = 0
	rank_gene_without_seed_dict = {}

	for rank in rank_gene_dict.keys():
		gene = rank_gene_dict[rank]

		if gene not in seed_gene_list:
			rank_gene_without_seed_dict[rank_idx] = gene
			rank_idx += 1

	return rank_gene_without_seed_dict


def Gene_Color_Rank_dict(rank_gene_dict, sorting_order, work_id):


	seed_info_dir = 'work/%s/seeds.txt' % work_id
	seed_file = open(seed_info_dir,'r')
	seed_file_readlines = seed_file.readlines()

	seed_gene_list = []

	for i in range(len(seed_file_readlines)):
		read = seed_file_readlines[i]
		read = read.replace('\n','')
		if '#' not in read:
			seed_gene_list.append(read)


	gene_rank_color_dict = {}

	#hex_color = ['#BD0013','#BB192D','#BA3347','#B94D61','#B8677B','#B78096','#B69AB0','#B5B4CA','#B4CEE4','#B3E8FF']
	#hex_color = ['#BD0013','#C41C2D','#CB3847','#D35561','#DA717B','#E18D96','#E9AAB0','#F0C6CA','#F7E2E4','#FFFFFF']
	#hex_color = ['#BD0013','#C3192A','#CA3342','#D04C59','#D76671','#DE7F89','#E499A0','#EBB2B8','#F1CCCF','#F8E5E7','#FFFFFF']
	hex_color = ['#BD0013','#CA3342','#D76671','#E499A0','#F1CCCF','#FFFFFF']
	rgb_color = []

	for hex_code in hex_color:
		value = hex_code.replace('#','')
		lv = len(value)
		rgb = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))
		rgb = str(rgb)
		rgb_color.append('rgb%s' % rgb)


	if sorting_order == 'ascend':
		rgb_color = rgb_color[::-1]

	maximum_len_genes = len(rank_gene_dict.keys()) 
	#chunck_size = int(round(maximum_len_genes / 10))
	chunck_size = int(round(maximum_len_genes / 5))
	rank_idx = 0

	#for color_idx in range(11):
	for color_idx in range(6):
		left_chunck = chunck_size * color_idx
		right_chunck = chunck_size * (color_idx + 1)
		print left_chunck, right_chunck

		for rank_idx in range(left_chunck, right_chunck):
			try :
				gene = rank_gene_dict[rank_idx]
				if gene in seed_gene_list:
					gene_rank_color_dict[gene] = 'rgb(0, 0, 0,0)'
				else:
					gene_rank_color_dict[gene] = rgb_color[color_idx]
			except KeyError:
				break


	return gene_rank_color_dict, seed_gene_list



def Create_Gene_Color_dict(cutoff_gene_list, gene_condition_dict, condition_id_dict):
	#Would be better if i used dataframe..
	colorset_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
	gene_color_dict = {}

	for gene in gene_condition_dict.keys():
		condition_id = gene_condition_dict[gene]
		condition = condition_id_dict[condition_id]

		if gene in cutoff_gene_list:

			if '__' not in condition:
				try :
					color_id = int(condition_id) - 1
					color_code = colorset_list[color_id]
				except IndexError:
					color_code = '#262626'
			if '__' in condition:
				color_code = '#262626'
		else:
			color_code = '#e5e5e5'

		gene_color_dict[gene] = color_code

	return gene_color_dict	




# read venn diagram, return as list.
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

# read combination map, return as {combinations, conditions}.
def read_combination(id_):
	groups = []
	df = pd.read_csv('work/%s/pvalue.txt'%id_,index_col=0)
	conds = set()
	for name,group in df.groupby('cond'):
		groups.append( {'name':name, 'size':group.shape[0]} )
		conds.update( set(name.split('__')) )
	return {'rows':groups, 'cols':list(conds)}

def read_graph(id_):
	with open('work/%s/graph.json'%id_,'r') as f:
		return json.load(f)


##
# Manual commands
##
if (__name__=="__main__"):
	tasks = {
		'generate_venn' : generate_venn,
		'generate_sel':generate_sel,
		'generate_graph': generate_graph,
		'select_venn_group': select_venn_group,
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
		tasks[ args.task ](args.path, args.sel_set.split(','), args.sub_set.split(','))
	elif (args.task == 'select_venn_group'):
		tasks[ args.task ](args.path, args.sel_set.split(','))
	else:
		tasks[ args.task ](args.path)
