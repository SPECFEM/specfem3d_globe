# Create your views here.

from django import forms
from django.shortcuts import render_to_response, get_object_or_404
from datetime import datetime
from django.http import HttpResponse, HttpResponseRedirect, Http404
from Specfem3DGlobe.web.Specfem3DGlobe.models import Mesh, Model, Simulation, UserInfo
from cig.web.seismo.events.models import Event

def getres_checkbox(request,key):
	if request.has_key(key) and request.POST[key] == 'on':
		return True
	else: 
		return False

def index(request):
	user = login_check(request);
	prev_simulations = Simulation.objects.filter(user=user)
	form = forms.FormWrapper(Simulation.AddManipulator(), {}, {})

	return render_to_response('Specfem3DGlobe/home.html',
							 {'prev_simulations': prev_simulations,
							  'form': form,
							  'userinfo': user})

def setparam(request):
	if not request.has_key('events'):
		errors = 'you must select an event.'
		form = forms.FormWrapper(Simulation.AddManipulator(), {}, {})
		return render_to_response('Specfem3DGlobe/home.html', 
								  {'errors': errors,
								  'form': form})

	mesh_type = request.POST['mesh_type']
	if mesh_type == 'regional':
		template = 'Specfem3DGlobe/simulation_form_regional.html'
	else:
		template = 'Specfem3DGlobe/simulation_form_global.html'

	return render_to_response(template, 
							 {'events': request.POST['events'],
							  'simulation_type': request.POST['simulation_type']})

def detail(request, sim_id):
	sim = get_object_or_404(Simulation,id=sim_id)
	return render_to_response('Specfem3DGlobe/detail.html', {'sim': sim}) 

def delete(request,sim_id):
	sim = get_object_or_404(Simulation,id=sim_id)
	if sim:
		if sim.mesh:
			sim.mesh.delete()
		if sim.model:
			sim.model.delete()
		sim.delete()
	return HttpResponseRedirect('/specfem3dglobe/')

def create_simulation(request):
	# 
	# first, do some parameter validity checking!!!
	#
	# Do some checking here!
	user = login_check(request);

	if user == None:
		return HttpResponseRedirect('/specfem3dglobe/')

	#
	# mesh information
	#
	_nchunks = request.POST['mesh_nchunks']
	_nproc_xi = request.POST['mesh_nproc_xi']
	_nproc_eta = request.POST['mesh_nproc_eta']
	_nex_xi = request.POST['mesh_nex_xi']
	_nex_eta = request.POST['mesh_nex_eta']
	_save_files = getres_checkbox(request,'mesh_save_files')
	_type = 1
	if request.POST['mesh_type'] == 'regional':
		_type = 2
	_angular_width_eta = 0.0
	_angular_width_xi = 0.0
	_center_latitude = 0.0
	_center_longitude = 0.0
	_gamma_rotation_azimuth = 0.0
	if request.POST['mesh_type'] == 'regional':
		_angular_width_eta = request.POST['mesh_angular_width_xi']
		_angular_width_xi = request.POST['mesh_angular_width_eta']
		_center_latitude = request.POST['mesh_center_latitude']
		_center_longitude = request.POST['mesh_center_longitude']
		_gamma_rotation_azimuth = request.POST['mesh_gamma_rotation_azimuth']
		
	mesh = Mesh(    
					nchunks                         = _nchunks, 
					nproc_xi                        = _nproc_xi,
					nproc_eta                       = _nproc_eta,
					nex_xi                          = _nex_xi,
					nex_eta                         = _nex_eta,
					save_files                      = _save_files,
					type                            = _type,
					angular_width_eta               = _angular_width_eta,
					angular_width_xi                = _angular_width_xi,
					center_latitude                 = _center_latitude,
					center_longitude                = _center_longitude,
					gamma_rotation_azimuth          = _gamma_rotation_azimuth 
				)
	mesh.save()

	# 
	# model information
	#
	_type = request.POST['model_type']
	_oceans = getres_checkbox(request,'model_oceans')
	_gravity = getres_checkbox(request,'model_gravity')
	_attenuation = getres_checkbox(request,'model_attenuation')
	_topography = getres_checkbox(request,'model_topography')
	_rotation = getres_checkbox(request,'model_rotation')
	_ellipticity = getres_checkbox(request,'model_ellipticity')
	model = Model(  
					type                            = _type,
					oceans                          = _oceans,
					gravity                         = _gravity,
					attenuation                     = _attenuation,
					topography                      = _topography,
					rotation                        = _rotation,
					ellipticity                     = _ellipticity 
				)
	model.save()

	# 
	# simulation information
	#
	_date = datetime.now()
	_status = 2
	_record_length = request.POST['simulation_record_length']
	_receivers_can_be_buried = getres_checkbox(request,'simulation_receivers_can_be_buried')
	_print_source_time_function = getres_checkbox(request,'simulation_print_source_time_function')
	_save_forward = False
	_movie_surface = getres_checkbox(request,'simulation_movie_surface')
	_movie_volume = getres_checkbox(request,'simulation_movie_volume')
	_absorbing_conditions = False
	if request.POST['mesh_type'] == 'regional':
		_absorbing_conditions = True
	_ntstep_between_frames = 100
	_simulation_type = request.POST['simulation_type']

	user = UserInfo.objects.get(userid=request.session['userid'])

	# 
	# event information
	#
	simulation = Simulation(
					user                            = user,
					date                            = _date,
					mesh                            = mesh,
					model                           = model,
					status                          = _status,
					record_length                   = _record_length,
					receivers_can_be_buried         = _receivers_can_be_buried,
					print_source_time_function      = _print_source_time_function,
					save_forward                    = _save_forward,
					movie_surface                   = _movie_surface,
					movie_volume                    = _movie_volume,
					absorbing_conditions            = _absorbing_conditions,
					ntstep_between_frames           = _ntstep_between_frames,
					simulation_type                 = _simulation_type
				)
	simulation.save()
	simulation.events.add(Event.objects.get(id=request.POST['events']))

	return HttpResponseRedirect('/specfem3dglobe/')

def info(request, info_str):
	template = None
	if info_str == 'mesh':
		template = 'Specfem3DGlobe/mesh_info.html'
	elif info_str == 'model':
		template = 'Specfem3DGlobe/model_info.html'
	elif info_str == 'output_format':
		template = 'Specfem3DGlobe/output_format_info.html'
	elif info_str == 'movie':
		template = 'Specfem3DGlobe/movie_info.html'

	if template == None:
		raise Http404

	return render_to_response(template)

def simulation_pml(request, sim_id):
	from django.template import loader, Context

	response = HttpResponse(mimetype='text/xml')
	#response['Content-Disposition'] = 'attachment; filename=parameters.xml'

	# Get data from the database here.
	simulation = get_object_or_404(Simulation, id=sim_id)

	t = loader.get_template('Specfem3DGlobe/simulation.pml')
	c = Context({
		'simulation': simulation,
	})
	response.write(t.render(c))
	return response

def events_txt(request, sim_id):
	from django.template import loader, Context

	response = HttpResponse(mimetype='text/plain')

	# Get data from the database here.
	simulation = get_object_or_404(Simulation, id=sim_id)

	t = loader.get_template('Specfem3DGlobe/events.txt')
	c = Context({
		'events': simulation.events.all(),
	})
	response.write(t.render(c))
	return response

def stations_txt(request, sim_id):
	from django.template import loader, Context

	response = HttpResponse(mimetype='text/plain')

	# Get data from the database here.
	simulation = get_object_or_404(Simulation, id=sim_id)

	t = loader.get_template('Specfem3DGlobe/stations.txt')
	c = Context({
		'stations': simulation.stations.all(),
	})
	response.write(t.render(c))
	return response

def login_check(request):
	user = None
	try:
		user = UserInfo.objects.get(userid=request.session['userid'])
	except UserInfo.DoesNotExist:
		pass
	except KeyError:
		pass

	return user

def login(request):
	user = None
	error = None
	try:
		user = UserInfo.objects.get(userid=request.POST['userid'])
	except UserInfo.DoesNotExist:
		error = "user does not exist"
	except KeyError:
		error = "post error"

	if user and user.password == request.POST['password']:
		request.session['userid'] = user.userid
		return HttpResponseRedirect('/specfem3dglobe/')
	else:
		return HttpResponseRedirect('/specfem3dglobe/')

def logout(request):
	try:
		del request.session['userid']
	except KeyError:
		pass
	return HttpResponseRedirect('/specfem3dglobe/')

