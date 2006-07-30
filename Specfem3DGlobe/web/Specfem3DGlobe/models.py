
from django.db import models

MESH_TYPES = (
    (1, 'global'),
    (2, 'regional'),
)

MODEL_TYPES = (
    (1, 'isotropic prem'),
    (2, 'transversely isotropic prem'),
    (3, 'iaspei'),
    (4, 'ak135'),
    (5, '3d isotropic'),
    (6, '3d anisotropic'),
    (7, '3d attenuation'),
)

STATUS_TYPES = (
    (1, 'preparing'),
    (2, 'ready'),
    (3, 'pending'),
    (4, 'running'),
    (5, 'done'),
)
    
class UserInfo(models.Model):
    userid = models.CharField(maxlength=100, core=True)
    lastname = models.CharField(maxlength=100, core=True)
    firstname = models.CharField(maxlength=100, core=True)
    email = models.CharField(maxlength=100, null=True)
    institution = models.CharField(maxlength=100, null=True)
    address1 = models.CharField(maxlength=100, null=True)
    address2 = models.CharField(maxlength=100, null=True)
    address3 = models.CharField(maxlength=100, null=True)
    phone = models.CharField(maxlength=20, null=True)
    class Admin:
        pass

class Mesh(models.Model):
    nchunks = models.IntegerField(core=True)
    nproc_xi = models.IntegerField(core=True)
    nproc_eta = models.IntegerField(core=True)
    nex_xi = models.IntegerField(core=True)
    nex_eta = models.IntegerField(core=True)
    save_files = models.BooleanField(core=True)
    result_file = models.CharField(maxlength=100, core=True, null=True)
    type = models.IntegerField(choices=MESH_TYPES, core=True)

    # what should be the default value for this?
    dry = models.BooleanField(core=True, default=False)

    # this is for regional only (when type == 2), and when global, all these values are null
    angular_width_eta = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)
    angular_width_xi = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)
    center_latitude = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)
    center_longitude = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)
    gamma_rotation_azimuth = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)

class Model(models.Model):
    type = models.IntegerField(choices=MODEL_TYPES, core=True, default=False)
    oceans = models.BooleanField(core=True, default=False)
    gravity = models.BooleanField(core=True, default=False)
    attenuation = models.BooleanField(core=True, default=False)
    topography = models.BooleanField(core=True, default=False)
    rotation = models.BooleanField(core=True, default=False)
    ellipticity = models.BooleanField(core=True, default=False)

class Simulation(models.Model):
    #
    # general information about the simulation
    #
    user = models.ForeignKey(UserInfo, edit_inline=models.TABULAR, num_in_admin=1)
    date = models.DateTimeField('simuation date')
    mesh = models.ForeignKey(Mesh, edit_inline=models.TABULAR, num_in_admin=1)
    model = models.ForeignKey(Model, edit_inline=models.TABULAR, num_in_admin=1)
    status = models.IntegerField(choices=STATUS_TYPES, default=1)

    #
    # specific information starts here
    #
    record_length = models.FloatField(max_digits=19, decimal_places=10, core=True)
    receivers_can_be_buried = models.BooleanField(core=True)
    print_source_time_function = models.BooleanField(core=True)
    save_forward = models.BooleanField(core=True, default=False)

    movie_surface = models.BooleanField(core=True)
    movie_volume = models.BooleanField(core=True)

    # this is the solution file path
    cmt_solution = models.CharField(maxlength=255, core=True, null=True)
    # header file path
    header_file = models.CharField(maxlength=255, core=True, null=True)
    # scratch_seismogram_archive
    scratch_seismogram_archive = models.CharField(maxlength=255, core=True, null=True)
    # seismogram_archive
    seismogram_archive = models.CharField(maxlength=255, core=True, null=True)
    # stations
    stations = models.CharField(maxlength=255, core=True, null=True)
    # output_file
    output_file = models.CharField(maxlength=255, core=True, null=True)

    # need to find out what the fields are for...
    # hdur_movie:
    hdur_movie = models.FloatField(max_digits=19, decimal_places=10, core=True, null=True)
    # absorbing_conditions: set to true for regional, and false for global
    absorbing_conditions = models.BooleanField(core=True)
    # number_of_this_run:
    number_of_this_run = models.IntegerField(core=True, default=1)
    # number_of_runs: choose 1 for a run without restart
    number_of_runs = models.IntegerField(core=True, default=1)
    # ntstep_between_frames: typical value is 100 time steps
    ntstep_between_frames = models.IntegerField(core=True, default=100)
    # simulation_type:
    simulation_type = models.CharField(maxlength=255, core=True, null=True)
    # dry:
    dry = models.BooleanField(core=True, default=False)
    # ntstep_between_output_seismos:
    ntstep_between_output_seismos = models.IntegerField(core=True, default=100)
    # ntstep_between_read_adjsrc:
    ntstep_between_read_adjsrc = models.IntegerField(core=True, default=100)
    # ntstep_between_output_info:
    ntstep_between_output_info = models.IntegerField(core=True, default=100)

    class Admin:
        pass


