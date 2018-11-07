def MLV_Profile_to_dict(mlv_file):

	open_mlv_file = open(mlv_file,'r')
	mlv_profile_readlines = open_mlv_file.readlines()
	
	condition_dict = {}
	condition_id_dict = {}
	gene_list = []
	data_check_point = 0

	for i in range(1,len(mlv_profile_readlines)):
		read = mlv_profile_readlines[i]

		if data_check_point == 0:
			read = read.replace('\n','')
			if '#' not in read[0]:
				token = read.split('\t')
				condition_id = token[0]
				condition = token[1]

				condition_id_dict[condition_id] = condition

		if data_check_point == 1:
			read = read.replace('\n','')
			token = read.split('\t')
			condition = token[0]

			for j in range(1, len(token)):
				gene = token[j]

				try : condition_dict[condition].append(gene)
				except KeyError: condition_dict[condition] = [gene]

				if gene not in gene_list:
					gene = gene.split('_')[0]
					gene_list.append(gene)

		if mlv_profile_readlines[i][0:3] == '###':
			data_check_point = 1

	for condition_id in condition_id_dict.keys():

		condition = condition_id_dict[condition_id]
		condition_gene_list = condition_dict[condition_id]
		condition_dict[condition] = condition_gene_list
		
	for condition_id in condition_id_dict.keys():
		condition_dict.pop(condition_id,'None')

	return condition_dict, gene_list

def Check_Unable_Genes(condition_dict, ppi_gene_list):
	#CHECK unable genes
	unable_gene_list = []
	for condition in condition_dict.keys():
		deg_genes = condition_dict[condition]

		for gene in deg_genes:
			if gene not in ppi_gene_list:
				unable_gene_list.append(gene)

	unable_gene_list = list(set(unable_gene_list))
	return unable_gene_list

def StringDB_to_dict():
	print '[FL_MLV] START StringDB_to_dict'

	StringDB_dir = '//var/www/htdocs/MLV_WEB/data/MED_STRING_PPI_Mus_musculus_Symbol.txt'
	StringDB_open = open(StringDB_dir,'r')
	StringDB_readlines = StringDB_open.readlines()

	StringDB_dict = {}
	Gene_list = []

	for i in range(len(StringDB_readlines)):
		read = StringDB_readlines[i]
		read = read.replace('\n','')
		token = read.split('\t')
		node_1 = token[0]
		node_2 = token[1]
		#if node_1 not in StringDB_dict.keys():
		try :
			StringDB_dict[node_1].append(node_2)
		except KeyError:
			StringDB_dict[node_1] = [node_2]

		Gene_list.append(node_1)
		Gene_list.append(node_2)

	Gene_list = list(set(Gene_list))
	SortGene_list = sorted(Gene_list)

	return StringDB_dict, SortGene_list

class Venndianet_RWR():

	def Create_RWR_Condition_AdjacencyMatrix_part1(self, Condition_dict, StringDB_dict, UnableGene_list):
	#Express condition should be unified. 
	#ex) only UP DEG or only DOWN DEG
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : START'
		DEGtopology_dict = {} 

		common_gene_list = [] #shares common node
		
		for condition in Condition_dict.keys():

			gene_list = Condition_dict[condition]
			common_gene_list = list(set(common_gene_list + gene_list))

		common_gene_list = list(set(common_gene_list) - set(UnableGene_list))
		common_gene_list = sorted(common_gene_list)
		print '##################'
		print len(common_gene_list)
		print '##################'


		for node1 in StringDB_dict.keys():
			if node1 in common_gene_list:
			#[1] If PPI node 1 is DEG gene

				DEGtopology_dict[node1] = []
				node2_list = StringDB_dict[node1]

				for node2 in node2_list:
					if node2 in common_gene_list:
						#[2] And PPI node 2 is DEG gene, append
						DEGtopology_dict[node1].append(node2)
	
		##debug
		debug_txt = open('debug.txt','w')
		for key in DEGtopology_dict.keys():
			debug_txt.write('%s\n' % key)
			dict_gene_list = DEGtopology_dict[key]
			debug_txt.write('\n'.join(dict_gene_list))

					
		print "[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : END"

		return DEGtopology_dict, common_gene_list
		#DEGtopology_dict : contains every node-edges of every conditions.
		#common_gene_list : every genes of every conditions


	def Create_RWR_Condition_AdjacencyMatrix_part2(self, DEGtopology_dict, common_gene_list, adj_file, unable_gene_list):
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : START'
		adj_matrix_file_name = adj_file
		adj_matrix_txt = open(adj_matrix_file_name,'w')

		print '[FL_stringDB_RWR] Number of genes for Adjacency matrix : ', len(common_gene_list)
		for node1 in common_gene_list:
			adj_matrix_txt.write('\t' + str(node1))
		adj_matrix_txt.write('\n')

		for node1 in common_gene_list:
			adj_matrix_txt.write(str(node1))
			#[1] write row name

			for node2 in common_gene_list:
				if node1 == node2:
				#[2-1] check current position. is duped
					adj_matrix_txt.write('\t0')

				if node1 != node2:
				#[2-2] if not duped,
					node2_list = DEGtopology_dict[node1]

					if node2 in node2_list:
						adj_matrix_txt.write('\t1')
					if node2 not in node2_list:
						adj_matrix_txt.write('\t0')

			adj_matrix_txt.write('\n')
		adj_matrix_txt.close()

	
	def Create_RWR_p0_vector(self, DEG_topology_dict, seed_gene_list, work_id):

		p0_vector_dir = 'work/%s/%s.p0.vector'% (work_id,work_id)
		p0_vector_txt = open(p0_vector_dir, 'w')

		for gene in DEG_topology_dict.keys():
			p0_vector_txt.write(str(gene))

			#[2] if gene is in intersection, it is considered as seed = p0
			if gene in seed_gene_list:
				p0_vector_txt.write('\t1\n')

			#[2] if gene is not in intersection, it is not considered as seed = p0
			if gene not in seed_gene_list:
				p0_vector_txt.write('\t0\n')


	#NOT IN USE
	def Run_RWR(self, result_file_string):
		import os
		path_to_rwr_script = '//var/www/htdocs/MLV/analysis/ppi/protocol/RWR.R'

		p0_vector_file = str(result_file_string) + '.p0.vector'
		adj_matrix_file = str(result_file_string) + '.adj.matrix'

		rwr_result = str(result_file_string) + '.rwr'
		Rscript_RWR_cmd = 'Rscript ' + str(path_to_rwr_script) + ' ' + str(adj_matrix_file) + ' ' + str(p0_vector_file) + ' ' + str(rwr_result)
		print Rscript_RWR_cmd
		os.system(Rscript_RWR_cmd)

class Parameter_Handling():

	def Get_Parameter(self, parameter_info_dir):
		
		cutoff = 'None'

		parameter_info_open = open(parameter_info_dir,'r')
		parameter_info_readlines = parameter_info_open.readlines()

		parameter_string = parameter_info_readlines[0]
		parameter_string = parameter_string.replace('\n','')
		cutoff = parameter_string

		sort_order_string = parameter_info_readlines[1]
		sort_order = sort_order_string.replace('\n','')
	
		parameter_info_open.close()

		return cutoff, sort_order
	
	def Cutoff_RWR(self, rwr_dir, cutoff, sort_order):

		cutoff = float(cutoff)

		rwr_open = open(rwr_dir,'r')
		rwr_readlines = rwr_open.readlines()
		rwr_dict = {}
		gene_list = []

		for i in range(1, len(rwr_readlines)):
			read = rwr_readlines[i]
			read = read.replace('\n','')
			token = read.split()

			gene = token[0]
			prob = token[1]
			
			if gene not in rwr_dict:
				rwr_dict[gene] = prob
			else:
				print '(debug) CRITICAL ERROR in Cutoff_RWR or INPUT'

		if sort_order == 'descend':
			sorted_rwr_dict = sorted(rwr_dict.items(), key=lambda x: x[1], reverse=True)
		if sort_order == 'ascend':
			sorted_rwr_dict = sorted(rwr_dict.items(), key=lambda x: x[1])

		result_length = len(rwr_dict.keys())
		cutoff = int(result_length * cutoff)

		for i in range(cutoff):
			gene = sorted_rwr_dict[i][0]

			if gene not in gene_list:
				gene_list.append(gene)
			else:
				print '(debug) CRITICAL ERROR in Cutoff_RWR or INPUT'

		return gene_list

	



class Graph_Methods():

	def Create_Graph_Info_DEFAULT(self, work_id):

		mlv_file = 'work/%s/%s.mlv'% (work_id, work_id)
		graph_json_file = 'work/%s/graph.json'% (work_id)
		graph_json_open = open(graph_json_file,'w')

		json_string = ''

		condition_dict, gene_list = MLV_Profile_to_dict(mlv_file)
		StringDB_dict, SortGene_list = StringDB_to_dict()

		json_string = Graph_Methods().Create_Json_Template(json_string, 1)
		json_string = Graph_Methods().Create_Node_Json(condition_dict, json_string)
		json_string = Graph_Methods().Create_Json_Template(json_string, 2)
		json_string = Graph_Methods().Create_Edge_Json(gene_list, StringDB_dict, json_string)
		json_string = Graph_Methods().Create_Json_Template(json_string, 3)

		graph_json_open.write(json_string)
		graph_json_open.close()

	def Create_Graph_Info_CUTOFF(self, work_id, cutoff_gene_list):
		print cutoff_gene_list

		mlv_file = 'work/%s/%s.mlv'% (work_id, work_id)
		graph_json_file = 'work/%s/graph.json'% (work_id)
		graph_json_open = open(graph_json_file,'w')

		json_string = ''

		condition_dict, gene_list = MLV_Profile_to_dict(mlv_file)
		rep_condition_dict = {}

		#refinement
		for condition in condition_dict.keys():
			condition_gene_list = condition_dict[condition]
			rep_condition_gene_list = condition_gene_list[:]

			for gene in condition_gene_list:
				if gene not in cutoff_gene_list:
					rep_condition_gene_list.remove(gene)

			rep_condition_dict[condition] = rep_condition_gene_list


		StringDB_dict, SortGene_list = StringDB_to_dict()

		condition_dict = rep_condition_dict
		print rep_condition_dict

		json_string = Graph_Methods().Create_Json_Template(json_string, 1)
		json_string = Graph_Methods().Create_Node_Json(condition_dict, json_string)
		json_string = Graph_Methods().Create_Json_Template(json_string, 2)
		json_string = Graph_Methods().Create_Edge_Json(gene_list, StringDB_dict, json_string)
		json_string = Graph_Methods().Create_Json_Template(json_string, 3)

		graph_json_open.write(json_string)
		graph_json_open.close()

	def Create_Graph_Info_CUTOFF_V2(self, work_id, cutoff_gene_list):
		print cutoff_gene_list

		mlv_file = 'work/%s/%s.mlv'% (work_id, work_id)
		graph_json_file = 'work/%s/graph.json'% (work_id)
		graph_json_open = open(graph_json_file,'w')

		json_string = ''

		condition_dict, gene_list = MLV_Profile_to_dict(mlv_file)
		StringDB_dict, SortGene_list = StringDB_to_dict()

		json_string = Graph_Methods().Create_Json_Template(json_string, 1)
		json_string = Graph_Methods().Create_Node_Json_ver_CUTOFF(condition_dict, json_string, cutoff_gene_list)
		json_string = Graph_Methods().Create_Json_Template(json_string, 2)
		json_string = Graph_Methods().Create_Edge_Json(gene_list, StringDB_dict, json_string)
		json_string = Graph_Methods().Create_Json_Template(json_string, 3)

		graph_json_open.write(json_string)
		graph_json_open.close()


	def Create_Json_Template(self, json_string, stage):
		if stage == 1:
			json_string += '{"nodes": ['
		if stage == 2:
			json_string += '],"edges": ['
		if stage == 3:
			json_string += ']}'
		return json_string

	def Create_Node_Json(self, condition_dict, json_string):

		for condition in condition_dict.keys():
			gene_list = condition_dict[condition]
			if len(gene_list) != 0:
				for gene in gene_list:
					json_string += '{"group": "%s", "name": "%s", "value": %s},'% (condition, gene, 1)
			
		json_string = json_string[:len(json_string) -1] #removing unnecessaory colon(,)
		return json_string

	def Create_Node_Json_ver_CUTOFF(self, condition_dict, json_string, cutoff_gene_list):

		for condition in condition_dict.keys():
			gene_list = condition_dict[condition]
			if len(gene_list) != 0:
				for gene in gene_list:
					if gene not in cutoff_gene_list:
						value = 0 
					else:
						value = 1

					json_string += '{"group": "%s", "name": "%s", "value": %s},'% (condition, gene, value)
			
		json_string = json_string[:len(json_string) -1] #removing unnecessaory colon(,)
		return json_string



	def Create_Edge_Json(self, gene_list, StringDB_dict, json_string):
	
		for source_gene in gene_list:

			if source_gene in StringDB_dict.keys():
				target_gene_list = StringDB_dict[source_gene]

				for target_gene in target_gene_list:
					if target_gene in gene_list:
						json_string += '{"source": "%s", "target": "%s", "weight": %s},'% (source_gene, target_gene, 1)

		json_string = json_string[:len(json_string) -1] #removing unnecessaory colon(,)
		return json_string



class Result_Summary():

	def MLV_to_Dict(self, mlv_file):

		condition_id_dict = {}
		gene_condition_dict = {}
		condition_gene_dict = {}

		mlv_open = open(mlv_file,'r')
		mlv_readlines = mlv_open.readlines()

		check_point = 0

		for i in range(1, len(mlv_readlines)):

			if check_point == 1:

				read = mlv_readlines[i]
				read = read.replace('\n','')
				token = read.split('\t')
				conditionID = token[0]

				for j in range(1,len(token)):
					gene = token[j]
					gene = gene.split('_')[0]

					try : condition_gene_dict[conditionID].append(gene)
					except KeyError: condition_gene_dict[conditionID] = [gene]

					gene_condition_dict[gene] = conditionID

			if check_point == 0:

				if mlv_readlines[i][0:3] != '###':

					read = mlv_readlines[i]
					read = read.replace('\n','')
					token = read.split('\t')

					conditionID = token[0]
					condition_tag = token[1]
					condition_id_dict[conditionID] = condition_tag

				else:
					check_point = 1

		return condition_id_dict, gene_condition_dict, condition_gene_dict

	def RWR_Result_to_Dict(self, rwr_file):

		rwr_dict = {}

		rwr_open = open(rwr_file,'r')
		RWR_readlines = rwr_open.readlines()

		for i in range(1, len(RWR_readlines)):
			read = RWR_readlines[i]
			read = read.replace('\n','')
			token = read.split()

			gene = token[0]
			prob = token[1]

			rwr_dict[gene] = float(prob)

		return rwr_dict

