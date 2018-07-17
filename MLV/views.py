from django.http import HttpResponse, Http404
from django.shortcuts import render, redirect
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
import random, string, sys, os, json
import MLV

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
  if (not validate(id)):
    raise Http404
  s = load_state(id)
  # refuse to display if not uploaded
  if (s['state'] == 'not_uploaded'):
    return render(request, 'invalid.html', {'id': id, 'state': s})

  if (request.method == "POST"):
    # prepare to save selected genes
    selection = request.POST['selection']
    MLV.select_venn_group(id, selection)
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
  # refuse to display if not uploaded
  if (s['state'] == 'not_uploaded'):
    return render(request, 'invalid.html', {'id': id, 'state': s})

  if (request.method == "POST"):
    # process it on another thread and go to result page
    s['result'] = 0
    save_state(s)
    return redirect('result', id=id)

  return render(request, 'settings.html', {'id': id, 'state': s})

def result(request, id):
  if (not validate(id)):
    raise Http404
  s = load_state(id)
  # refuse to display if not uploaded
  if (s['state'] == 'not_uploaded'):
    return render(request, 'invalid.html', {'id': id, 'state': s})

  # if result is not yet,
  # do automatic refresh
  return render(request, 'result.html', {'id': id, 'state': s, 'graph':MLV.read_graph(id)})
