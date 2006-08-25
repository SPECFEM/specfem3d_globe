# Create your views here.

from django import forms
from django.contrib.auth import authenticate, login, logout
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.models import User
from django.core import validators
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Registration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


usernameTaken = "This username is already taken."


def isNotExistingUser(field_data, all_data):
	try:
		User.objects.get(username = field_data)
	except User.DoesNotExist:
		return
	raise validators.ValidationError(usernameTaken)


class RegistrationManipulator(forms.Manipulator):
	
	def __init__(self):
		self.fields = [
			# User
			forms.TextField('first_name',   maxlength=30,  is_required=True),
			forms.TextField('last_name',    maxlength=30,  is_required=True),
			forms.EmailField('email',                      is_required=True),
			# UserInfo
			forms.TextField('institution',  maxlength=100, is_required=True),
			forms.TextField('address1',     maxlength=100),
			forms.TextField('address2',     maxlength=100),
			forms.TextField('address3',     maxlength=100),
			forms.PhoneNumberField('phone'),
		]

	def usernameValidatorList(self):
		return [validators.isAlphaNumeric]


class RegistrationAddManipulator(RegistrationManipulator):
	
	def __init__(self):
		super(RegistrationAddManipulator, self).__init__()
		self.fields.extend([
			forms.TextField('username',     maxlength=30,  is_required=True, validator_list=self.usernameValidatorList()),
			forms.PasswordField('password', maxlength=128, is_required=True),
			])
		
	def save(self, new_data, request):
		user, created = User.objects.get_or_create(
			username = new_data['username'],
			defaults = {'first_name': new_data['first_name'],
				    'last_name':  new_data['last_name'],
				    'email':      new_data['email']})
		if not created:
			# Race: the username was just taken!
			return {'username': [usernameTaken]}
		user.set_password(new_data['password'])
		user.save()
		UserInfo.objects.create(user        = user,
					institution = new_data['institution'],
					address1    = new_data['address1'],
					address2    = new_data['address2'],
					address3    = new_data['address3'],
					phone       = new_data['phone'])
		# Log-in the new user.
		user = authenticate(username=new_data['username'], password=new_data['password'])
		if user is not None:
			login(request, user)
		return {}

	def flatten_data(self):
		return {}
	
	def usernameValidatorList(self):
		validator_list = super(RegistrationAddManipulator, self).usernameValidatorList()
		validator_list.append(isNotExistingUser)
		return validator_list
	
	
class RegistrationChangeManipulator(RegistrationManipulator):

	def __init__(self, user):
		super(RegistrationChangeManipulator, self).__init__()
		self.user = user
	
	def flatten_data(self):
		new_data = {}
		new_data.update(self.user.__dict__)
		try:
			userInfo = self.user.userinfo
		except UserInfo.DoesNotExist:
			pass
		else:
			new_data.update(userInfo.__dict__)
		return new_data
	
	def save(self, new_data, request):
		# Create UserInfo if it doesn't exist.
		user = self.user
		try:
			userInfo = user.userinfo
		except UserInfo.DoesNotExist:
			userInfo = UserInfo()
			user.userinfo = userInfo
		# Save the new user.
		user.first_name      = new_data['first_name']
		user.last_name       = new_data['last_name']
		user.email           = new_data['email']
		userInfo.institution = new_data['institution']
		userInfo.address1    = new_data['address1']
		userInfo.address2    = new_data['address2']
		userInfo.address3    = new_data['address3']
		userInfo.phone       = new_data['phone']
		user.save()
		userInfo.save()
		return {}


def registration(request):
	user = request.user
	if user.is_anonymous():
		manipulator = RegistrationAddManipulator()
		template = 'Specfem3DGlobe/register.html'
	else:
		manipulator = RegistrationChangeManipulator(user)
		template = 'Specfem3DGlobe/userinfo_form.html'

	if request.POST:
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
