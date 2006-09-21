# Create your views here.

from django import forms
from django.conf import settings
from django.contrib.auth import logout
from django.contrib.auth.decorators import user_passes_test
from django.core import validators
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render_to_response, get_object_or_404
from django.template.context import RequestContext
from cig.web.seismo.events.models import Event
from forms import SimulationTypeManipulator, SimulationWizardManipulator, RegistrationAddManipulator, RegistrationChangeManipulator
from models import Simulation
import os, os.path


OUTPUT_ROOT = os.path.join(settings.MEDIA_ROOT, 'Specfem3DGlobe', 'output')


# Create our own version of 'login_required' which redirects to our login page.
login_required = user_passes_test(lambda u: not u.is_anonymous(), "/specfem3dglobe/login")


def simulation_start_form():
	manipulator = SimulationTypeManipulator()
	new_data = {'mesh__nchunks': '6', 'simulation_type': '1' }
	return forms.FormWrapper(manipulator, new_data, {})

def index(request):
	recent_simulations = Simulation.objects.filter(user=request.user)
	return render_to_response('Specfem3DGlobe/home.html',
				  {'recent_simulations': recent_simulations,
				   'form': simulation_start_form()},
				  RequestContext(request, {}))
index = login_required(index)

def simulation_index(request):
	from django.views.generic.list_detail import object_list
	return object_list(request, Simulation.objects.filter(user=request.user),
			   allow_empty=True,
			   extra_context={'form': simulation_start_form()})
simulation_index = login_required(simulation_index)

def create_simulation(request):

	nchunks = int(request.POST['mesh__nchunks'])
	absorbing_conditions = "True"
	if nchunks == 1:
		template = 'Specfem3DGlobe/simulation_base.html'
	elif nchunks == 2:
		template = 'Specfem3DGlobe/simulation_form_2chunks.html'
	elif nchunks == 3:
		template = 'Specfem3DGlobe/simulation_form_3chunks.html'
	else: # nchunks == 6
		absorbing_conditions = "False"
		template = 'Specfem3DGlobe/simulation_form_global.html'

	manipulator = SimulationWizardManipulator(nchunks)

	if request.POST.has_key('blank'):
		# Return a new, blank form of the requested type.
		new_data = manipulator.flatten_data()
		new_data.update({
			# hidden fields
			'mesh__user': request.user.id,
			'mesh__nchunks': request.POST['mesh__nchunks'],
			'model__user': request.user.id,
			'user': request.user.id,
			'mesh': "dummy",
			'model': "dummy",
			'simulation_type': request.POST['simulation_type'],
			'absorbing_conditions': absorbing_conditions,
			})
		form = forms.FormWrapper(manipulator, new_data, {})
		return render_to_response(template,
					  { 'form': form },
					  RequestContext(request, {}))

	# User is POSTing data.

	new_data = request.POST.copy()

	# 
	# First, do some parameter validity checking!!!
	#
	# Do some checking here!
	#
	# This only checks for simple errors.
	#
	
	errors = manipulator.get_validation_errors(new_data)
	if errors:
		form = forms.FormWrapper(manipulator, new_data, errors)
		return render_to_response(template,
					  { 'form': form },
					  RequestContext(request, {}))

	manipulator.do_html2python(new_data)
	
	manipulator.save(new_data)
	
	return HttpResponseRedirect('/specfem3dglobe/')

create_simulation = login_required(create_simulation)

def delete(request, sim_id):
	from django.views.generic.create_update import delete_object
	post_delete_redirect = '/specfem3dglobe/'
 	if request.method == 'POST':
 		sim = get_object_or_404(Simulation, id=sim_id)
 		sim.mesh.delete()
 		sim.model.delete()
		sim.delete()
		return HttpResponseRedirect(post_delete_redirect)
	return delete_object(request, Simulation, post_delete_redirect,
			     object_id=sim_id)
delete = login_required(delete)

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

def logout_view(request):
	logout(request)
	return HttpResponseRedirect('/specfem3dglobe/login/')

def events_txt(request, sim_id):
	from cig.seismo.events import CMTSolution

	response = HttpResponse(mimetype='text/plain')

	simulation = get_object_or_404(Simulation, id=sim_id)

	count = simulation.events.source_set.count()
	for event in simulation.events.source_set.all():
		cmtSolution = CMTSolution.createFromDBModel(event)

		# NYI: Specfem requires this to be zero for one
		# source, and positive for others.  Set it to zero for
		# single-source events.  For multi-source kinematic
		# ruptures, for now we assume the user uploaded a
		# proper CMTSOLUTION file.
		if count == 1:
			cmtSolution.timeShift = 0.0
		
		response.write(str(cmtSolution))

	return response

def stations_txt(request, sim_id):
	response = HttpResponse(mimetype='text/plain')

	simulation = get_object_or_404(Simulation, id=sim_id)
	stations = simulation.stations.station_set.all()

	response.write("%d\n" % len(stations))
	for station in stations:
		response.write("%-5s %-3s %8.4f %10.4f %6.1f %6.1f\n" %
			       (station.code, station.network.code,
				station.latitude, station.longitude, station.elevation, station.bur))

	return response

def output_tar_gz(request, sim_id):
	import shutil
	response = HttpResponse(mimetype='application/x-gtar')
	response['Content-Encoding'] = 'gzip'
	response['Content-Disposition'] = 'attachment; filename=simulation-%s.tar.gz' % sim_id
	filename = os.path.join(OUTPUT_ROOT, sim_id + '.tar.gz')
	stream = open(filename, 'rb')
	shutil.copyfileobj(stream,  response)
	stream.close()
	return response

def update_simulation_status(request, sim_id):
	from forms import SimulationStatusManipulator

	manipulator = SimulationStatusManipulator()

	simulation = get_object_or_404(Simulation, id=sim_id)

	if request.method == 'POST':
		response = HttpResponse(mimetype='text/plain')
		new_data = request.POST.copy()
		new_data.update(request.FILES)
		errors = manipulator.get_validation_errors(new_data)
		if errors:
			response.write(repr(errors))
		else:
			manipulator.do_html2python(new_data)
			output = new_data['output']
			if output:
				content = output['content']
				try:
					os.makedirs(OUTPUT_ROOT)
				except OSError: # Directory probably already exists.
					pass
				filename = os.path.join(OUTPUT_ROOT, sim_id + '.tar.gz')
				stream = open(filename, 'wb')
				stream.write(content)
				stream.close()
			simulation.status = new_data['status']
			simulation.save()
			response.write(repr('OK'))
		return response
	
	errors = {}
	new_data = {'status': simulation.status}
	form = forms.FormWrapper(manipulator, new_data, errors)
	return render_to_response('Specfem3DGlobe/simulation_status.html',
				  {'form': form},
				  RequestContext(request, {}))

def registration(request):
	user = request.user
	if user.is_anonymous():
		manipulator = RegistrationAddManipulator()
		template = 'Specfem3DGlobe/register.html'
	else:
		manipulator = RegistrationChangeManipulator(user)
		template = 'Specfem3DGlobe/userinfo_form.html'

	if request.method == 'POST':
		new_data = request.POST.copy()
		errors = manipulator.get_validation_errors(new_data)
		if not errors:
			manipulator.do_html2python(new_data)
			errors = manipulator.save(new_data, request)
			if not errors:
				return HttpResponseRedirect('/specfem3dglobe/')
	else:
		# Populate new_data with a 'flattened' version of the current data.
		new_data = manipulator.flatten_data()
		errors = {}

	# Populate the FormWrapper.
	form = forms.FormWrapper(manipulator, new_data, errors, edit_inline = True)
	
	return render_to_response(template, { 'form': form }, RequestContext(request, {}))
