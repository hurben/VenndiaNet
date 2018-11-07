###################################
#Module library for venn_ppi_RWR.py
#
#
#MAKE adjacency matrix for R
#We will create matrix with same probablity

def WELCOME_MESSAGE():
	print '********************************'
	print 'Venndia-net                     '
	print 'v0.0.1'
	print '        > S T A R T <'
	print '\n'
	print '********************************'

def REMOVE_FC(gene):
    gene = gene.split('_')[0]
    return gene

def REMOVE_FC_FROM_LIST(gene_fc_list):

    gene_list = []
    for gene_fc in gene_fc_list:
        gene = REMOVE_FC(gene_fc)
        gene_list.append(gene)
    return gene_list


def MLV_Profile_to_dict(mlv_profile_readlines):

	condition_dict = {}
	gene_list = []

	for i in range(len(mlv_profile_readlines)):
		read = mlv_profile_readlines[i]
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

	return condition_dict, gene_list

def StringDB_to_dict():
	print '[FL_MLV] START StringDB_to_dict'

	StringDB_dir = '//var/www/htdocs/MLV/analysis/ppi/STRING_PPI_Mus_musculus_Symbol.txt'
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
	print "[FL_MLV] Total # of Genes: ", len(SortGene_list)
	print "[FL_MLV] Total # of StringDB_dict Keys: ", len(StringDB_dict.keys())
	print '[FL_MLV] END StringDB_to_dict'

	return StringDB_dict, SortGene_list

def Create_intersection_dict(condition_dict, unable_gene_list):
	#purpose : intersection_dict[condition] = [intersection genes] where '__' not in condition
	#So, it makes intersection gene list of non-intersection condition.
	#example ) intersection_dict[muscle] = [geneA__conditionB, geneC__conditionC]  of course, __conditionN will not be included.
	#note : pure intersection genes
	intersection_dict = {}

	#[1] : loop condition, if non-intersection == condition
	for condition_i in condition_dict.keys():
		if '__' not in condition_i:
			intersection_dict[condition_i] = []

			#[2] : loop condition, if intersection == condition
			for condition_j in condition_dict.keys():
				if condition_j != condition_i and '__' in condition_j:
					condition_A = condition_j.split('__')[0]
					condition_B = condition_j.split('__')[1]
					if condition_A == condition_i or condition_B == condition_i:

						gene_list = REMOVE_FC_FROM_LIST(condition_dict[condition_j])
						intersection_dict[condition_i].extend(gene_list)
	#[3] 
	unable_intersection_gene_count = 0
	for condition in intersection_dict.keys():
		max_length = len(intersection_dict[condition])
		for unable_gene in unable_gene_list:
			if unable_gene in intersection_dict[condition]:
				unable_intersection_gene_count += 1
				intersection_dict[condition].remove(unable_gene)
		print '<', condition, '> includes', len(intersection_dict[condition]), 'intersection genes while',  unable_intersection_gene_count, 'genes are not in PPI. Meaning that initially there were ', max_length, 'genes'

	#intersection_dict does not include unable_gene
	return intersection_dict

class Venndianet_RWR():

	def Create_RWR_AdjacencyMatrix_part1(self, Gene_list, StringDB_dict, UnableGene_list):
		#Based on selected gene list
		#make DEG whole network topology
		print '[FL_stringDB_RWR] Create_RWR_AdjacencyMatrix START'
		DEGtopology_dict = {} 

		for node1 in StringDB_dict.keys():
			if node1 in Gene_list:
				if node1 not in UnableGene_list: 
					#[1] If PPI node 1 is DEG gene
					#Initialize topology dict
					try : dummy = DEGtopology_dict[node1]
					except KeyError: DEGtopology_dict[node1] = []

					node2_list = StringDB_dict[node1]

					for node2 in node2_list:
						if node2 in Gene_list:
							if node2 not in UnableGene_list: 
								#[2] And PPI node 2 is DEG gene, append
								DEGtopology_dict[node1].append(node2)

		print "[FL_stringDB_RWR] DEGtopology contains ", len(DEGtopology_dict.keys()), " Nodes"

		for node1 in DEGtopology_dict.keys():
			node2_list = DEGtopology_dict[node1]
			if len(node2_list) == 0:
				print node1, " does not have any edges"

		print "[FL_stringDB_RWR] Create_RWR_AdjacencyMatrix END"
		return DEGtopology_dict

	def Create_RWR_AdjacencyMatrix_part2(self, DEGtopology_dict,result_file_string):
		
		import os.path

		adj_matrix_file_name = str(result_file_string) + '.adj.matrix.txt'

		if os.path.isfile(adj_matrix_file_name) == False:

			adj_matrix_txt = open(adj_matrix_file_name,'w')

			gene_list = DEGtopology_dict.keys()
			gene_list = sorted(gene_list)

			#write header
			adj_matrix_txt.write(str(gene_list[0]))
			for i in range(1, len(gene_list)):
				node1 = gene_list[i]
				adj_matrix_txt.write('\t' + str(node1))
			adj_matrix_txt.write('\n')

			#write contents
			for node1 in gene_list:
				#[1] write row name
				adj_matrix_txt.write(str(node1))

				for node2 in gene_list:

					#[2] check current position. is duped
					if node1 == node2:
						adj_matrix_txt.write('\t0')

					if node1 != node2:
					#[3] if not duped, 
						node2_list = DEGtopology_dict[node1]

						#[4] check whether node2 exists in node1 list
						if node2 in node2_list:
							adj_matrix_txt.write('\t1')
						if node2 not in node2_list:
							adj_matrix_txt.write('\t0')

				adj_matrix_txt.write('\n')
			adj_matrix_txt.close()

		else:
			print "[FL_stringDB_RWR][Create_RWR_AdjacencyMatrix_part2] NO need to do this process"

		return gene_list

	def Create_RWR_Condition_AdjacencyMatrix_part1(self, Condition_dict, StringDB_dict, UnableGene_list):
	#Based on selected gene list
	#Disected Topology
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : START'
		DEGtopology_dict = {} 

		for condition in Condition_dict.keys():
			if '__' not in condition:
				Gene_list = Condition_dict[condition]
				Gene_list = REMOVE_FC_FROM_LIST(Gene_list)

				for node1 in StringDB_dict.keys():
					if node1 in Gene_list:
					#[1] If PPI node 1 is DEG gene

						DEGtopology_dict[condition, node1] = []
						node2_list = StringDB_dict[node1]

						for node2 in node2_list:
							if node2 in Gene_list:
								#[2] And PPI node 2 is DEG gene, append
								DEGtopology_dict[condition, node1].append(node2)
								

		print "[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : END"

		#node2 != unable_gene
		return DEGtopology_dict

	def Create_RWR_Condition_AdjacencyMatrix_part2(self, DEGtopology_dict, Condition_dict, result_file_string, unable_gene_list):
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : START'
		condition_specific_gene_entry_dict = {}

		for condition in Condition_dict.keys():
			if '__' not in condition:
				Gene_list = Condition_dict[condition]
				Gene_list = REMOVE_FC_FROM_LIST(Gene_list)
				Gene_list = list(set(Gene_list))

				#ESSENTIAL
				for gene in unable_gene_list:
					if gene in Gene_list:
						Gene_list.remove(gene)
				#ESSENTIAL

				Gene_list = list(sorted(Gene_list))

				adj_matrix_file_name = str(result_file_string) + '.' + str(condition) + '.adj.matrix'
				adj_matrix_txt = open(adj_matrix_file_name,'w')
				adj_matrix_txt.write('\t' + str(Gene_list[0]))
				
				condition_specific_gene_entry_dict[condition] = Gene_list

				print '[FL_stringDB_RWR] ', condition, '> Number of genes for Adjacency matrix : ', len(Gene_list)
				for i in range(1, len(Gene_list)):
					node1 = Gene_list[i]
					adj_matrix_txt.write('\t' + str(node1))
				adj_matrix_txt.write('\n')

				for node1 in Gene_list:
					adj_matrix_txt.write(str(node1))
					#[1] write row name

					for node2 in Gene_list:
						if node1 == node2:
						#[2-1] check current position. is duped
							adj_matrix_txt.write('\t0')

						if node1 != node2:
						#[2-2] if not duped, 
							node2_list = DEGtopology_dict[condition, node1]

							#[4] check whether node2 exists in node1 list
							if node2 in node2_list:
								adj_matrix_txt.write('\t1')
							if node2 not in node2_list:
								adj_matrix_txt.write('\t0')

					adj_matrix_txt.write('\n')
				adj_matrix_txt.close()

		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : END'
		#unable gene list excluded
		return condition_specific_gene_entry_dict




	def Create_RWR_Condition_AdjacencyMatrix_not_default_part1(self, Condition_dict, StringDB_dict, UnableGene_list):
	#Express condition should be unified. 
	#ex) only UP DEG or only DOWN DEG
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : START'
		DEGtopology_dict = {} 

		common_gene_list = [] #shares common topolgoy
		
		for condition in Condition_dict.keys():

			if '__' not in condition:

				gene_list = Condition_dict[condition]
				gene_list = REMOVE_FC_FROM_LIST(gene_list)
				#Does not matter. Because every FC status is UP or DOWN
				common_gene_list = list(set(common_gene_list + gene_list))

		common_gene_list = list(set(common_gene_list) - set(UnableGene_list))
		print len(common_gene_list)

		for node1 in StringDB_dict.keys():
			if node1 in common_gene_list:
			#[1] If PPI node 1 is DEG gene

				DEGtopology_dict[node1] = []
				node2_list = StringDB_dict[node1]

				for node2 in node2_list:
					if node2 in gene_list:
						#[2] And PPI node 2 is DEG gene, append
						DEGtopology_dict[node1].append(node2)
						
		print "[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 1 : END"

		return DEGtopology_dict, common_gene_list
		#DEGtopology_dict : contains every node-edges of every conditions.
		#common_gene_list : every genes of every conditions


	def Create_RWR_Condition_AdjacencyMatrix_not_default_part2(self, DEGtopology_dict, common_gene_list, result_file_string, unable_gene_list):
		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : START'
		adj_matrix_file_name = str(result_file_string) + '.adj.matrix'
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

		print '[FL_stringDB_RWR] Create_RWR_Condition_AdjacencyMatrix STEP 2 : END'


	def Create_RWR_p0_vector(self, condition_specific_gene_entry_dict, intersection_dict, result_file_string): #
		#[1] : loop condition of adjacency matrixs
		#condition = number of adjacency matrix, typically should not include '__'
		for condition in condition_specific_gene_entry_dict.keys():
			p0_vector_txt = open(str(result_file_string) + '.' + str(condition) + '.p0.vector', 'w')

			for gene in condition_specific_gene_entry_dict[condition]:
				p0_vector_txt.write(str(gene))

				#[2] if gene is in intersection, it is considered as seed = p0
				if gene in intersection_dict[condition]:
					p0_vector_txt.write('\t1\n')
				#[2] if gene is not in intersection, it is not considered as seed = p0
				if gene not in intersection_dict[condition]:
					p0_vector_txt.write('\t0\n')
	
	def Set_Seed_Genes(self, condition_dict, selected_condition, unable_gene_list):
		seed_gene_list = []

		selected_condition_list = selected_condition.split('__')

		print 'Total Condition selected :',len(selected_condition_list)

		for condition in selected_condition_list:
			print condition , '#########'
			temp_list = condition_dict[condition]
			set_b = set(temp_list)

			if len(seed_gene_list) != 0:
				
				set_a = set(seed_gene_list)
				seed_gene_list = list(set_a.intersection(set_b))

			if len(seed_gene_list) == 0:
				seed_gene_list = condition_dict[condition]

		print len(seed_gene_list)

#########HARD CODED. NEED TO INCREASE REUSEABLITY.
#########REUSEABLITY DOES NOT EXISTs
#		for condition in condition_dict.keys():
#			if '__' not in condition:
#				if condition not in selected_condition_list:
#					print condition,'~~~~~~~~~~'
#					set_a = set(seed_gene_list)
#					set_b = set(condition_dict[condition])
#					seed_gene_list = list(set_a - set_b)

#		print len(seed_gene_list)
#################################################################################
#		set_a = set(seed_gene_list)
#		set_b = set(unable_gene_list)
#		seed_gene_list = list(set_a - set_b)
#		print len(seed_gene_list)

		seed_gene_list = REMOVE_FC_FROM_LIST(seed_gene_list)

		return seed_gene_list
		#seed gene dict shares same seed genes

	def Create_RWR_p0_vector_not_default(self, DEG_topology_dict, seed_gene_list, result_file_string): #

		p0_vector_txt = open(str(result_file_string) + '.p0.vector', 'w')

		for gene in DEG_topology_dict.keys():
			p0_vector_txt.write(str(gene))

			#[2] if gene is in intersection, it is considered as seed = p0
			if gene in seed_gene_list:
				p0_vector_txt.write('\t1\n')
			#[2] if gene is not in intersection, it is not considered as seed = p0
			if gene not in seed_gene_list:
				p0_vector_txt.write('\t0\n')

	def Run_RWR(self, result_file_string, condition_list):

		import os
		path_to_rwr_script = '//var/www/htdocs/MLV/analysis/ppi/protocol/RWR.R'

		for condition in condition_list:
			p0_vector_file = str(result_file_string) + '.' + str(condition) + '.p0.vector'
			adj_matrix_file = str(result_file_string) + '.' + str(condition) + '.adj.matrix'

			rwr_result = str(result_file_string) + '.' + str(condition) + '.rwr'
			Rscript_RWR_cmd = 'Rscript ' + str(path_to_rwr_script) + ' ' + str(adj_matrix_file) + ' ' + str(p0_vector_file) + ' ' + str(rwr_result)
			os.system(Rscript_RWR_cmd)

	def Run_RWR_not_default(self, result_file_string):
		import os
		path_to_rwr_script = '//var/www/htdocs/MLV/analysis/ppi/protocol/RWR.R'

		p0_vector_file = str(result_file_string) + '.p0.vector'
		adj_matrix_file = str(result_file_string) + '.adj.matrix'

		rwr_result = str(result_file_string) + '.rwr'
		Rscript_RWR_cmd = 'Rscript ' + str(path_to_rwr_script) + ' ' + str(adj_matrix_file) + ' ' + str(p0_vector_file) + ' ' + str(rwr_result)
		os.system(Rscript_RWR_cmd)


	def Cut_RWR(self, result_file_string, condition_list, unable_gene_list, cut_off):

		import operator
		rwr_cut_dict = {}

		for condition in condition_list:
			rwr_file = str(result_file_string) + '.' + str(condition) + '.rwr'
			rwr_open = open(rwr_file,'r')
			rwr_readlines = rwr_open.readlines()

			temp_dict = {}
			rwr_cut_dict[condition] = []

			for i in range(1, len(rwr_readlines)):
				read = rwr_readlines[i]
				read = read.replace('\n','')
				token = read.split()

				gene_symbol = token[0]
				value = token[1]

				try : value = float(value)
				except ValueError: None

				node_score = value
				temp_dict[gene_symbol] = node_score

			sorted_temp_list = sorted(temp_dict.items(), key=operator.itemgetter(1))
			#sorted_temp_list[n][0] = keys

			MAX_LEN = len(sorted_temp_list)
			#CUTOFF = 0.2 #20% BASED ON RESULTS OF HFD, RWR results.
			CUTOFF = float(cut_off) #20% BASED ON RESULTS OF HFD, RWR results.
			LEN_CUTOFF = int(MAX_LEN) - int(MAX_LEN * CUTOFF) 
			print condition, 'cutoff :', LEN_CUTOFF, ' from ' , MAX_LEN

			for i in range(LEN_CUTOFF):
				gene_symbol = sorted_temp_list[i][0]
				rwr_cut_dict[condition].append(gene_symbol)

		return rwr_cut_dict

	def Summarize_RWR(self, condition_list, condition_dict, rwr_cut_dict, intersection_dict, result_file_string):

	
		for condition in condition_list:

			gene_list = rwr_cut_dict[condition]
			original_gene_list = condition_dict[condition]
			intersection_gene_list =  intersection_dict[condition]


			result_file_txt = str(result_file_string) + '.' + str(condition) + '.rwr.summary'
			result_file_txt = open(result_file_txt,'w')

			intersection_gene_count = 0
			
			for gene in gene_list:
				if gene not in intersection_gene_list:

					for original_gene_with_status in original_gene_list:

						original_gene_without_status = REMOVE_FC(original_gene_with_status)

						if original_gene_without_status == gene:	
							result_file_txt.write(str(original_gene_with_status) + '\n')

				if gene in intersection_gene_list:
					#print 'CRITICAL ERROR IN Summarize_RWR'
					intersection_gene_count += 1

			print condition, ': From Total', len(intersection_gene_list), ' intersection genes. Removed ' + str(intersection_gene_count) + ' intersection genes from gene list of ' + str(len(gene_list))
			result_file_txt.close()


	def Summarize_RWR_planB(self, condition_list, condition_dict, rwr_cut_dict, intersection_dict, result_file_string):

	
		for condition in condition_list:

			gene_list = rwr_cut_dict[condition]
			original_gene_list = condition_dict[condition]
			intersection_gene_list =  intersection_dict[condition]


			result_file_txt = str(result_file_string) + '.' + str(condition) + '.rwr.summary.planB'
			result_file_txt = open(result_file_txt,'w')

			intersection_gene_count = 0
			
			for gene in gene_list:

				for original_gene_with_status in original_gene_list:

					original_gene_without_status = REMOVE_FC(original_gene_with_status)

					if original_gene_without_status == gene:	
						result_file_txt.write(str(original_gene_with_status) + '\t')

				if gene in intersection_gene_list:
					result_file_txt.write('1\n')
				if gene not in intersection_gene_list:
					result_file_txt.write('0\n')

				if gene in intersection_gene_list:
					#print 'CRITICAL ERROR IN Summarize_RWR'
					intersection_gene_count += 1

			print condition, ': From Total', len(intersection_gene_list), ' intersection genes. Removed ' + str(intersection_gene_count) + ' intersection genes from gene list of ' + str(len(gene_list))
			result_file_txt.close()





class Venndianet_Neighbor():

	def Condition_Exclusive(self, deg_gene, condition, condition_dict, ppi_dict, unable_gene_list, result_dict):
		#condition_dict[condition] = gene_up
		deg_gene_strip = REMOVE_FC(deg_gene)
		#deg_gene = TP53_UP
		#deg_gene_strip = TP53

		if deg_gene_strip not in unable_gene_list:
		#[1] Current condition gene is not in other condition
			for condition_i in condition_dict.keys():
				if condition_i != condition and '__' in condition_i:
				#[2] if condition_loop is not current condition
				#[3] from intersections,

					condition_1 = condition_i.split('__')[0]
					condition_2 = condition_i.split('__')[1]

					if condition_1 == condition or condition_2 == condition:
						#[4] that contains current condition.

						intersection_deg_genes = condition_dict[condition_i]
						#intersection_deg_genes = [TP53_UP, PI3K_DOWN]

						if deg_gene not in intersection_deg_genes:
						#[5] make sure deg_gene is not in intersection
							check = 1

							for intersection_deg_gene in intersection_deg_genes:
							#[6] make sure this gene is not included in 'any' neighbor of intersection genes
								intersection_deg_gene_strip = REMOVE_FC(intersection_deg_gene)
								if check == 1:
									if intersection_deg_gene_strip not in unable_gene_list:

										node2_genes = ppi_dict[intersection_deg_gene_strip]
										#sould be no error in here, because we filtered out as unable gene.

										if deg_gene_strip not in node2_genes: #pointing out exclusive
											check = 1

										if deg_gene_strip in node2_genes:
											check = 0

								if check == 0:
									#this deg gene is not available
									break #disable this loop					

							#After process, [6]
							if check == 1:
							#We know now that this deg_gene 
							#is not intersection gene
							#is not neighbor of intersection gene
								term = '.exclusive'
								condition = str(condition) + str(term)
								try : result_dict[condition].append(deg_gene)
								except KeyError : result_dict[condition] = [deg_gene]

							if check == 0:
								None
		return result_dict


