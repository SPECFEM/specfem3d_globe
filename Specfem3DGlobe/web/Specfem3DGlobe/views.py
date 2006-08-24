# Create your views here.

from django import forms
from django.contrib.auth import logout
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.models import User
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.shortcuts import render_to_response, get_object_or_404
from django.template.context import RequestContext
from datetime import datetime
from Specfem3DGlobe.web.Specfem3DGlobe.models import Mesh, Model, Simulation, UserInfo
from cig.web.forms import TeeManipulator
from cig.web.seismo.events.models import Event


# Create our own version of 'login_required' which redirects to our login page.
login_required = user_passes_test(lambda u: not u.is_anonymous(), "/specfem3dglobe/login")


class SimulationTypeManipulator(forms.Manipulator):
	def __init__(self):
		forms.Manipulator.__init__(self)
		from models import MESH_TYPES, SIMULATION_TYPES
		self.fields = (
			forms.RadioSelectField(field_name='mesh__type', choices=MESH_TYPES),
			forms.RadioSelectField(field_name='simulation_type', choices=SIMULATION_TYPES),
			)

def index(request):
	prev_simulations = Simulation.objects.filter(user=request.user)
	manipulator = SimulationTypeManipulator()
	new_data = {'mesh__type': '1', 'simulation_type': '1' }
	form = forms.FormWrapper(manipulator, new_data, {})

	return render_to_response('Specfem3DGlobe/home.html',
				  {'prev_simulations': prev_simulations,
				   'form': form},
				  RequestContext(request, {}))
index = login_required(index)



	

class SimulationWizardManipulator(TeeManipulator):
	
	def __init__(self, mesh_type):
		TeeManipulator.__init__(self,
					{'mesh__': Mesh.AddManipulator(),
					 'model__': Model.AddManipulator(),
					 '': Simulation.AddManipulator()},
					['user',
					 'mesh',
					 'mesh__user',
					 'mesh__nchunks', # see below
					 'model',
					 'model__user'])
		# Create our own custom fields.
		from models import NCHUNKS_CHOICES
		if mesh_type == '1':
			nchunks_choices = NCHUNKS_CHOICES
		else:
			nchunks_choices = NCHUNKS_CHOICES[:-1]
		regional_field = forms.RadioSelectField(field_name='mesh__nchunks',
							choices=nchunks_choices)
		self.manipulators['mesh__'].fields.append(regional_field)
		self.fields.append(regional_field)
		return
	
	def save(self, new_data):
		self._revert_field_names(new_data)
		mesh = self.manipulators['mesh__'].save(new_data)
		new_data['mesh'] = mesh.id
		model = self.manipulators['model__'].save(new_data)
		new_data['model'] = model.id
		# Save the simulation last.
		self.manipulators[''].save(new_data)
		return


def create_simulation(request):

	mesh_type = request.POST['mesh__type']
	if mesh_type == '1':
		template = 'Specfem3DGlobe/simulation_form_global.html'
		absorbing_conditions = False
	elif mesh_type == '2':
		template = 'Specfem3DGlobe/simulation_form_regional.html'
		absorbing_conditions = True
	else:
		raise RuntimeError()

	manipulator = SimulationWizardManipulator(mesh_type)

	if request.POST.has_key('blank'):
		# Return a new, blank form of the requested type.
		new_data = manipulator.flatten_data()
		new_data['mesh__type'] = mesh_type
		new_data['mesh__nchunks'] = '1'
		new_data['simulation_type'] = request.POST['simulation_type']
		new_data['absorbing_conditions'] = absorbing_conditions
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
		return render_to_response(template, { 'form': form }, RequestContext(request, {}))

	manipulator.do_html2python(new_data)
	
	# Fill-in user data.
	new_data['user'] = request.user.id
	new_data['mesh__user'] = request.user.id
	new_data['model__user'] = request.user.id

	# Fill-in automatic fields.
	new_data['absorbing_conditions'] = absorbing_conditions
	
	manipulator.save(new_data)
	
	return HttpResponseRedirect('/specfem3dglobe/')

create_simulation = login_required(create_simulation)


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


def registration(request):
	user = request.user
	follow = {
		# Suppress these fields.
		'password': False,
		'date_joined': False,
		'last_login': False,
		'is_staff': False,
		'is_active': False,
		'is_superuser': False,
		}
	if user.is_anonymous():
		manipulator = User.AddManipulator(follow)
		template = 'Specfem3DGlobe/register.html'
	else:
		# Create UserInfo if it doesn't exist.
		try:
			userInfo = user.userinfo
		except UserInfo.DoesNotExist:
			userInfo = UserInfo()
			user.userinfo = userInfo
			user.save()
			userInfo.save()
		manipulator = User.ChangeManipulator(user.id, follow)
		template = 'Specfem3DGlobe/userinfo_form.html'

	if request.POST:
		new_data = request.POST.copy()
		errors = manipulator.get_validation_errors(new_data)
		manipulator.do_html2python(new_data)
		new_data['userinfo.0.user'] = user.id # Not sure why this is this necessary.
		new_data['userinfo.0.user_id'] = user.id # Not sure why this is this necessary.
		if not errors:
			manipulator.save(new_data)
			return HttpResponseRedirect('/specfem3dglobe/')
	else:
		# Populate new_data with a 'flattened' version of the current data.
		new_data = manipulator.flatten_data()
		errors = {}

	# Populate the FormWrapper.
	form = forms.FormWrapper(manipulator, new_data, errors, edit_inline = True)
	
	return render_to_response(template, { 'form': form }, RequestContext(request, {}))
