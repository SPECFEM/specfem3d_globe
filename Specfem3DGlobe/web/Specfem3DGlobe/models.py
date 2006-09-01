
from django.db import models
from django.contrib.auth.models import User
from cig.web.models import CurrentUser
from cig.web.seismo.events.models import Event
from cig.web.seismo.stations.models import Station
from Specfem3DGlobe.Daemon import STATUS_CHOICES


MODEL_TYPES = (
	(1, 'isotropic prem'),
	(2, 'transversely isotropic prem'),
	(3, 'iaspei'),
	(4, 'ak135'),
	(5, '3d isotropic'),
	(6, '3d anisotropic'),
	(7, '3d attenuation'),
)

SIMULATION_TYPES = (
	(1, 'forward'),
	(2, 'adjoint'),
	(3, 'both forward and adjoint'),
)


NCHUNKS_CHOICES = (
	(1, 'regional with 1 chunk'),
	(2, 'regional with 2 chunks'),
	(3, 'regional with 3 chunks'),
	(6, 'global (6 chunks)'),
)

oneThruNine = tuple([(i, str(i)) for i in xrange(1, 10)])
NPROC_CHOICES = oneThruNine
NEX_C_CHOICES = oneThruNine


class UserInfo(models.Model):
	user = models.OneToOneField(User)
	institution = models.CharField(maxlength=100, core=True)
	address1 = models.CharField(maxlength=100, null=True, blank=True)
	address2 = models.CharField(maxlength=100, null=True, blank=True)
	address3 = models.CharField(maxlength=100, null=True, blank=True)
	phone = models.PhoneNumberField(null=True, blank=True)
	class Admin:
		pass


class Mesh(models.Model):

	user = models.ForeignKey(User)
	created = models.DateTimeField(auto_now_add=True, editable=False)
        modified = models.DateTimeField(auto_now=True, editable=False)
	
	nchunks = models.IntegerField(core=True, choices=NCHUNKS_CHOICES, default=1)
	nproc_xi = models.IntegerField(core=True, choices=NPROC_CHOICES)
	nproc_eta = models.IntegerField(core=True, choices=NPROC_CHOICES)
	nex_xi = models.IntegerField(core=True)
	nex_eta = models.IntegerField(core=True)
	save_files = models.BooleanField(core=True)

	# this is for regional only (when type == 2), and when global, all these values are fixed
	angular_width_eta = models.FloatField(max_digits=19, decimal_places=10, core=True)
	angular_width_xi = models.FloatField(max_digits=19, decimal_places=10, core=True)
	center_latitude = models.FloatField(max_digits=19, decimal_places=10, core=True)
	center_longitude = models.FloatField(max_digits=19, decimal_places=10, core=True)
	gamma_rotation_azimuth = models.FloatField(max_digits=19, decimal_places=10, core=True)

        class Admin:
            pass


class Model(models.Model):
	
	user = models.ForeignKey(User)
	created = models.DateTimeField(auto_now_add=True, editable=False)
        modified = models.DateTimeField(auto_now=True, editable=False)
	
	type = models.IntegerField(choices=MODEL_TYPES, core=True, default=False)
	oceans = models.BooleanField(core=True, default=False)
	gravity = models.BooleanField(core=True, default=False)
	attenuation = models.BooleanField(core=True, default=False)
	topography = models.BooleanField(core=True, default=False)
	rotation = models.BooleanField(core=True, default=False)
	ellipticity = models.BooleanField(core=True, default=False)

	def get_type_id(self):
		return {
			1: 'isotropic_prem',
			2: 'transversly_isotropic_prem',
			3: 'iaspei',
			4: 'ak135',
			5: '3D_isotropic',
			6: '3D_anisotropic',
			7: '3D_attenuation'}[self.type]

        class Admin:
            pass


class Simulation(models.Model):
	#
	# general information about the simulation
	#
	user = models.ForeignKey(User)
	name = models.CharField(maxlength=100, unique=True, core=True)
	created = models.DateTimeField(auto_now_add=True, editable=False)
        modified = models.DateTimeField(auto_now=True, editable=False)
	started = models.DateTimeField(editable=False, null=True)
	finished = models.DateTimeField(editable=False, null=True)
	
	mesh = models.ForeignKey(Mesh)
	model = models.ForeignKey(Model)
	status = models.CharField(maxlength=100, choices=STATUS_CHOICES, default='SimStatusNew', editable=False)

	#
	# specific information starts here
	#
	record_length = models.FloatField(max_digits=19, decimal_places=10, core=True, default=10.0)
	receivers_can_be_buried = models.BooleanField(core=True)
	print_source_time_function = models.BooleanField(core=True)
	save_forward = models.BooleanField(core=True, default=False)

	movie_surface = models.BooleanField(core=True)
	movie_volume = models.BooleanField(core=True)

	# CMTSOLUTION
	events = models.ManyToManyField(Event,
					limit_choices_to = {'user__exact' : CurrentUser()})
	# STATIONS
	stations = models.ManyToManyField(Station,
					  limit_choices_to = {'user__exact' : CurrentUser()})

	# need to find out what the fields are for...
	# hdur_movie:
	hdur_movie = models.FloatField(max_digits=19, decimal_places=10, core=True, default=0.0, blank=True)
	# absorbing_conditions: set to true for regional, and false for global
	absorbing_conditions = models.BooleanField(core=True)
	# ntstep_between_frames: typical value is 100 time steps
	ntstep_between_frames = models.IntegerField(core=True, default=100)
	# ntstep_between_output_info: typical value is 100 time steps
	ntstep_between_output_info = models.IntegerField(core=True, default=100)
	# ntstep_between_output_seismos : typical value is 5000
	ntstep_between_output_seismos = models.IntegerField(core=True, default=5000)
	# simulation_type:
	simulation_type = models.IntegerField(choices=SIMULATION_TYPES, default=1)

	def get_simulation_type_id(self):
		return {1: 'forward', 2: 'adjoint', 3: 'both'}[self.simulation_type]

	class Admin:
		pass


