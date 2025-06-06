from django.http import HttpResponse, Http404
from django.shortcuts import render, redirect
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
import random, string, sys, os, json
import MLV
import threading

def main(request):
	return render(request, 'main.html')

def generate(request):
	def id_generator(size=12, chars=string.ascii_uppercase + string.digits):
		return ''.join(random.choice(chars) for _ in range(size))
	id_ = request.GET.get('id')
	if (not id_):
		id_ = id_generator()
	# make directory & initialize state for working
	os.mkdir('work/%s' % id_, 0775)
	with open('work/%s/work.json' % id_,'w') as f:
		json.dump({
			'id': id_,
			'state': 'not_uploaded',
			'cond_count': 0,
			'files': [],
			'result': 0
		}, f)
	print 'New generation of passing by: %s' % id_
	# now redirect to step_upload
	return redirect('step_upload', id=id_)

def delete(request, id):
	pass


# validation checker
def validate(id_):
	# check directory exists
	# check state file (json) exists
	if (not os.path.exists('work/%s' % id_)):
		return False
	if (not os.path.exists('work/%s/work.json' % id_)):
		return False
	return True
def load_state(id_):
	with open('work/%s/work.json' % id_, 'r') as f:
		return json.load(f)
def save_state(j):
	id_ = j['id']
	with open('work/%s/work.json' % id_, 'w') as f:
		json.dump(j, f)

def step_upload(request, id):
	if (not validate(id)):
		raise Http404
	s = load_state(id)

	if (request.method == "POST"):
		# prepare to save file
		if (os.path.exists('work/%s/f' % id) == True):

			clean_files_cmd = 'rm work/%s/f/*' % id
			os.system(clean_files_cmd)
			clean_files_cmd = 'rm work/%s/*.mlv' % id
			os.system(clean_files_cmd)

		if (not os.path.exists('work/%s/f' % id)):
			os.mkdir('work/%s/f' % id, 0775)
		s['files'] = []
		# save files and update state

		for _name, f in request.FILES.iteritems():
			if (not f):
				continue
			fn = f.name
			s['files'].append(fn)
			default_storage.save('work/%s/f/%s' % (id,fn), ContentFile(f.read()))

		s['state'] = 'uploaded'
		s['cond_count'] = len(s['files'])
		save_state(s)
		# generate set info file to draw venn diagram
		MLV.generate_venn(id)
		return redirect('step_select', id=id)

	return render(request, 'upload.html', {'id': id, 'state': s})

def step_select(request, id):

	print '< STEP SELECT > from views.py'

	if (not validate(id)):
		raise Http404
	s = load_state(id)
	# refuse to display if not uploaded
	if (s['state'] == 'not_uploaded'):
		return render(request, 'invalid.html', {'id': id, 'state': s})

	if (request.method == "POST"):
		# prepare to save selected genes
		selection = request.POST['selection']

		#Creates Seed Gene list for RWR : seeds.txt
		MLV.select_venn_group(id, selection)
		MLV.Bootstrap(id, selection)

		s['state'] = 'selected'
		save_state(s)
		return redirect('step_params', id=id)

	# draw venn diagram in case of count <= 4
	# else, draw in heatmap method (miRTarVis)
	venn = MLV.read_venn(id)
	comb = MLV.read_combination(id)
	return render(request, 'select.html', {
		'id': id, 'state': s, 'venn': venn, 'combination': comb
		})

def step_params(request, id):

	if (not validate(id)):
		raise Http404
	s = load_state(id)

	seed_info_txt = 'work/%s/seeds.txt'% id
	parameter_info_txt =  'work/%s/parameter.txt'% id
	selected_condition = ''

	if os.path.exists(seed_info_txt) == True:
		seed_info_open = open(seed_info_txt,'r')
		seed_info_readlines = seed_info_open.readlines()

		for i in range(1, len(seed_info_readlines)):
			read = seed_info_readlines[i]
			read = read.replace('\n','')
			if '#' in read:
				read = read.replace('#','')
				selected_condition += '%s<br>'% read


	if os.path.exists(seed_info_txt) == False:
		selected_condition = 'None'

	# refuse to display if not uploaded
	if (s['state'] == 'not_uploaded'):
		return render(request, 'invalid.html', {'id': id, 'state': s})

	if (request.method == "POST"):
		# process it on another thread and go to result page
		#retrieve parameters from here
		cut_off =  request.POST['cutoff']
		method =  request.POST['sort_method']
		parameter_info_open = open(parameter_info_txt,'w')
		parameter_info_open.write('%s\n'% cut_off )
		parameter_info_open.write('%s\n'% method )
		parameter_info_open.close()

		s['result'] = 0
		save_state(s)
		return redirect('result', id=id)

	return render(request, 'settings.html', {'id': id, 'state': s, 'selected_condition': selected_condition})

def RWR_thread(id):

	MLV.Run_RWR(id)

def result(request, id):

	if (not validate(id)):
		raise Http404
	s = load_state(id)
	work_id = s['id']
	# refuse to display if not uploaded
	if (s['state'] == 'not_uploaded'):
		return render(request, 'invalid.html', {'id': id, 'state': s})

	t = threading.Thread(target=RWR_thread, args=(id,))
	t.start()

	if (os.path.exists('work/%s/RWR.result' % id) == True):
		cutoff, sorting_order, cutoff_gene_list = MLV.After_RWR(id)
		gene_rank_dict, gene_condition_dict, condition_id_dict = MLV.Create_Result_Summary(id, cutoff, sorting_order)
		gene_color_dict = MLV.Create_Gene_Color_dict(cutoff_gene_list, gene_condition_dict, condition_id_dict)
		rank_gene_dict = MLV.Create_Rank_Gene_dict(gene_rank_dict)

		gene_rank_color_dict, seed_gene_list = MLV.Gene_Color_Rank_dict(rank_gene_dict, sorting_order, work_id)
		rank_gene_without_seed_dict = MLV.Create_Rank_Gene_dict_without_seeds(rank_gene_dict,seed_gene_list)

		manage_result_file_cmd ='cp work/%s/%s.result.summary //var/www/html/MLV_PUBLIC/upload/%s'% (work_id, work_id, work_id)
		os.system(manage_result_file_cmd)

	else:
		return render(request, 'invalid.html', {'id': id, 'state': s})

	# if result is not yet,
	# do automatic refresh
	return render(request, 'result.html', {'id': id, 'state': s, 'graph':json.dumps(MLV.read_graph(id)), 'gene_rank' : gene_rank_dict, 'rank_gene' : rank_gene_dict, 'gene_color' : gene_color_dict, 'gene_rank_color': gene_rank_color_dict, 'rank_gene_without_seed_dict': rank_gene_without_seed_dict})
