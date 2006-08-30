
from cig.web.forms import TeeManipulator
from django import forms
from django.contrib.auth import authenticate, login
from django.contrib.auth.models import User
from django.core import validators
from models import Mesh, Model, Simulation, UserInfo


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


class SimulationTypeManipulator(forms.Manipulator):
	def __init__(self):
		forms.Manipulator.__init__(self)
		from models import NCHUNKS_CHOICES, SIMULATION_TYPES
		self.fields = (
			forms.RadioSelectField(field_name='mesh__nchunks', choices=NCHUNKS_CHOICES),
			forms.RadioSelectField(field_name='simulation_type', choices=SIMULATION_TYPES),
			)

	
class MeshAddManipulator(Mesh.AddManipulator):
	def __init__(self):
		super(MeshAddManipulator, self).__init__()
		# replace generic fields with custom fields
		del self['user']
		del self['nchunks']
		self.fields.extend([
			forms.HiddenField(field_name='user', is_required=True),
			forms.HiddenField(field_name='nchunks', is_required=True),
			])


class OneOr2ChunkMeshAddManipulator(MeshAddManipulator):
	
	def __init__(self):
		super(OneOr2ChunkMeshAddManipulator, self).__init__()
		# replace generic fields with custom fields
		del self['nex_xi']
		del self['nex_eta']
		from models import NEX_C_CHOICES
		self.fields.extend([
			forms.SelectField(field_name='nex_xi_c', choices=NEX_C_CHOICES, is_required=True),
			forms.SelectField(field_name='nex_eta_c', choices=NEX_C_CHOICES, is_required=True),
			])

	def save(self, new_data):
		# Compute derived mesh values.  This must be done
		# (redundantly) on the server, in case JavaScript is
		# unavailable on the client.
		nex_xi_c = int(new_data['nex_xi_c'])
		nex_eta_c = int(new_data['nex_eta_c'])
		nproc_xi = int(new_data['nproc_xi'])
		nproc_eta = int(new_data['nproc_eta'])
		new_data['nex_xi'] = str(16 * nex_xi_c * nproc_xi)
		new_data['nex_eta'] = str(16 * nex_eta_c * nproc_eta)
		return super(OneOr2ChunkMeshAddManipulator, self).save(new_data)


class ThreeOr6ChunkMeshAddManipulator(MeshAddManipulator):
	
	def __init__(self):
		super(ThreeOr6ChunkMeshAddManipulator, self).__init__()
		
		# replace generic fields with custom fields
		del self['nproc_xi']
		del self['nproc_eta']
		del self['nex_xi']
		del self['nex_eta']
		
		from models import NPROC_CHOICES, NEX_C_CHOICES
		
		# Force the user to pick a value, so that the
		# JavaScript will compute.  This isn't necessary (we
		# perform the same computations server-side); it's
		# just that question marks look silly in a completed
		# form.
		nproc_choices = (('', '---------'),) + NPROC_CHOICES
		
		self.fields.extend([
			forms.SelectField(field_name='nproc', choices=nproc_choices, is_required=True),
			forms.SelectField(field_name='nex_c', choices=NEX_C_CHOICES, is_required=True),
			])

	def save(self, new_data):
		# Compute derived mesh values.  This must be done
		# (redundantly) on the server, in case JavaScript is
		# unavailable on the client.
		nex_c = int(new_data['nex_c'])
		nproc = int(new_data['nproc'])
		new_data['nproc_xi'] = str(nproc)
		new_data['nproc_eta'] = str(nproc)
		new_data['nex_xi'] = str(16 * nex_c * nproc)
		new_data['nex_eta'] = str(16 * nex_c * nproc)
		return super(ThreeOr6ChunkMeshAddManipulator, self).save(new_data)


class ModelAddManipulator(Model.AddManipulator):
	def __init__(self):
		super(ModelAddManipulator, self).__init__()
		# replace generic fields with custom fields
		del self['user']
		self.fields.extend([
			forms.HiddenField(field_name='user', is_required=True),
			])


class SimulationAddManipulator(Simulation.AddManipulator):
	
	def __init__(self):
		super(SimulationAddManipulator, self).__init__()
		# replace generic fields with custom fields
		del self['user']
		del self['mesh']
		del self['model']
		del self['simulation_type']
		del self['absorbing_conditions']
		self.fields.extend([
			forms.HiddenField(field_name='user', is_required=True),
			forms.HiddenField(field_name='mesh', is_required=True),
			forms.HiddenField(field_name='model', is_required=True),
			forms.HiddenField(field_name='simulation_type', is_required=True),
			forms.HiddenField(field_name='absorbing_conditions'),
			])
		
	def do_html2python(self, new_data):
		super(SimulationAddManipulator, self).do_html2python(new_data)
		# This field requires special love and care, I guess
		# because it's a HiddenField instead of a
		# CheckboxField.
		new_data['absorbing_conditions'] = {'True': True, 'False': False}[new_data['absorbing_conditions']]


class SimulationWizardManipulator(TeeManipulator):
	
	def __init__(self, nchunks):
		if nchunks < 3:
			MeshAddManip = OneOr2ChunkMeshAddManipulator
		else:
			MeshAddManip = ThreeOr6ChunkMeshAddManipulator
		TeeManipulator.__init__(self,
					{'mesh__': MeshAddManip(),
					 'model__': ModelAddManipulator(),
					 '': SimulationAddManipulator()})
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


class SimulationStatusManipulator(forms.Manipulator):

	def __init__(self):
		from models import STATUS_CHOICES
		self.fields = [
			forms.SelectField(field_name='status', choices=STATUS_CHOICES, is_required=True),
			forms.FileUploadField(field_name='output'),
			]
		return


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

# end of file
