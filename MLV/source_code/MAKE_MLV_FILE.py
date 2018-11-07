import sys
import argparse

#MAKE .mlv file from inputlist

def REMOVE_FC(gene):
	gene = gene.split('_')[0]
	return gene

def REMOVE_FC_FROM_LIST(gene_fc_list):

	gene_list = []
	for gene_fc in gene_fc_list:
		gene = REMOVE_FC(gene_fc)
		gene_list.append(gene)
	return gene_list

def text_to_dict(file_dir, condition_gene_dict, condition_list):

	if '/' in file_dir:
		condition_name = file_dir.split('/')
		slash_n = len(condition_name)
		condition_name = condition_name[slash_n - 1]

	condition_list.append(condition_name)

	open_file = open(file_dir,'r')
	file_readlines = open_file.readlines()
	gene_list = []

	for i in range(1, len(file_readlines)):
		read = file_readlines[i]
		read = read.replace('\n','')
		read = read.replace('\r','')

		token = read.split('\t')
		gene_id = token[0]
		gene = token[1]


		if gene in gene_list:
			print 'Duplicated DEG :', gene

		if gene not in gene_list:
		#in order to avoid duplicates
			try : condition_gene_dict[condition_name].append(gene)
			except KeyError: condition_gene_dict[condition_name] = [gene]

		total_gene_list.append(gene)
		gene_list.append(gene)


	return condition_gene_dict, condition_list, total_gene_list

def organize_conditions(condition_gene_dict, condition_list, total_gene_list):

	number_of_cases = power(2, len(condition_list)) - 1
	
	condition_id_dict, condition_id_dict_reverse = define_condition_classes(condition_list, number_of_cases)
	venn_sheet_dict = create_venn_sheet(condition_id_dict, condition_id_dict_reverse, condition_gene_dict, total_gene_list)

	return venn_sheet_dict, condition_id_dict


def define_condition_classes(condition_list, number_of_cases):

	condition_id_dict = {}
	condition_id_dict_reverse = {}

	for i in range(len(condition_list)):

		if i == 0:
			condition_id = 1

		else:
			condition_id = power(2, i)

		condition = condition_list[i]
		condition_id_dict[condition_id] = [condition]
		condition_id_dict_reverse[condition] = condition_id

	return condition_id_dict, condition_id_dict_reverse

		
def create_venn_sheet(condition_id_dict, condition_id_dict_reverse, condition_gene_dict, gene_list):
#	condition_id_dict[id] = ['condition']
#	condition_id_dict['condition'] = id

	venn_sheet_dict = {}


	for gene in gene_list:

		gene_case_id = 0
		gene_condition_list = []

		for condition in condition_gene_dict:

			condition_gene_list = condition_gene_dict[condition]
		
			if gene in condition_gene_list:
				gene_case_id += condition_id_dict_reverse[condition]
				gene_condition_list.append(condition)

####
		if gene_case_id not in condition_id_dict.keys():
			condition_id_dict[gene_case_id] = gene_condition_list
			
		gene_condition_list_string = '__'.join(map(str,gene_condition_list))

		if gene_condition_list_string not in condition_id_dict_reverse.keys():
			condition_id_dict_reverse[gene_condition_list_string] = gene_case_id
####

		try : 
			if gene not in venn_sheet_dict[gene_case_id]:
				venn_sheet_dict[gene_case_id].append(gene)
		except KeyError: 
			venn_sheet_dict[gene_case_id] = [gene]

	return venn_sheet_dict	


def make_mlv_file(venn_sheet_dict, condition_id_dict, result_file):


	mlv_txt = open(result_file,'w')

	mlv_txt.write('ConditionID\tConditionList\n')

	condition_id_list = list(map(int, condition_id_dict.keys()))
	sorted_condition_id_list = sorted(condition_id_list)
	
	CR_sorted_condition_id_list = []
	#CR : CameraReady

	for condition_id in sorted_condition_id_list:

		gene_condition_list = condition_id_dict[condition_id]
		gene_condition_list_string = '__'.join(map(str,gene_condition_list))

		if '__' not in gene_condition_list_string:

			number_of_genes = len(venn_sheet_dict[condition_id])
			mlv_txt.write(str(condition_id) + '\t' + str(gene_condition_list_string) + '\t' +str(number_of_genes) + '\n')
			CR_sorted_condition_id_list.append(condition_id)

	for condition_id in sorted_condition_id_list:

		gene_condition_list = condition_id_dict[condition_id]
		gene_condition_list_string = '__'.join(map(str,gene_condition_list))

		if '__' in gene_condition_list_string:
			number_of_genes = len(venn_sheet_dict[condition_id])
			mlv_txt.write(str(condition_id) + '\t' + str(gene_condition_list_string) + '\t' +str(number_of_genes) + '\n')
			CR_sorted_condition_id_list.append(condition_id)



	mlv_txt.write('###\n')

	for condition_id in CR_sorted_condition_id_list:

		gene_list = venn_sheet_dict[condition_id]
		mlv_txt.write(str(condition_id))

		for gene in gene_list:
			mlv_txt.write('\t' + str(gene))
		mlv_txt.write('\n')


if __name__ == '__main__':
	
	import numpy
	power = numpy.power

	parser = argparse.ArgumentParser()

	parser.add_argument('-i', '--input_list', dest = 'input_list')
	parser.add_argument('-o', '--output', dest = 'result_file')

	args = parser.parse_args()
	input_list_file = args.input_list
	result_file = args.result_file

	condition_gene_dict = {}
	condition_list = []
	total_gene_list = []

	open_input_dir = open(input_list_file,'r')
	input_list_readlines = open_input_dir.readlines()
	result_txt = open(result_file,'w')

	for i in range(len(input_list_readlines)):

		read = input_list_readlines[i]
		read = read.replace('\n','')
		condition_file = read
		condition_gene_dict, condition_list, total_gene_list = text_to_dict(condition_file, condition_gene_dict, condition_list)	

	venn_sheet_dict, condition_id_dict = organize_conditions(condition_gene_dict, condition_list, total_gene_list)
	make_mlv_file(venn_sheet_dict,condition_id_dict, result_file)

	
