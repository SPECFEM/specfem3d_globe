
from django.db import models

class Mesher(models.Model):
    dry = models.BooleanField()
    center_longitude = models.FloatField(max_digits=19, decimal_places=10)
    angular_width_eta = models.FloatField(max_digits=19, decimal_places=10)
    nex_eta = models.IntegerField()
    # skipped output-file <cig.addyndum.properties.OutputFile object at 0x18727b0c>
    center_latitude = models.FloatField(max_digits=19, decimal_places=10)
    angular_width_xi = models.FloatField(max_digits=19, decimal_places=10)
    nproc_eta = models.IntegerField()
    gamma_rotation_azimuth = models.FloatField(max_digits=19, decimal_places=10)
    nchunks = models.IntegerField()
    nproc_xi = models.IntegerField()
    save_files = models.BooleanField()
    nex_xi = models.IntegerField()

class Model(models.Model):
    oceans = models.BooleanField()
    gravity = models.BooleanField()
    attenuation = models.BooleanField()
    topography = models.BooleanField()
    rotation = models.BooleanField()
    ellipticity = models.BooleanField()

class Solver(models.Model):
    receivers_can_be_buried = models.BooleanField()
    # skipped cmt-solution <cig.addyndum.properties.InputFile object at 0x187345cc>
    hdur_movie = models.FloatField(max_digits=19, decimal_places=10)
    ntstep_between_read_adjsrc = models.IntegerField()
    movie_volume = models.BooleanField()
    # skipped header-file <cig.addyndum.properties.OutputFile object at 0x1873462c>
    record_length = models.FloatField(max_digits=19, decimal_places=10)
    save_forward = models.BooleanField()
    absorbing_conditions = models.BooleanField()
    number_of_this_run = models.IntegerField()
    number_of_runs = models.IntegerField()
    # skipped seismogram-archive <cig.addyndum.properties.OutputFile object at 0x18734f6c>
    ntstep_between_frames = models.IntegerField()
    simulation_type = models.CharField(maxlength=255)
    movie_surface = models.BooleanField()
    dry = models.BooleanField()
    # skipped scratch-seismogram-archive <cig.addyndum.properties.ScratchFile object at 0x18734eec>
    ntstep_between_output_seismos = models.IntegerField()
    # skipped stations <cig.addyndum.properties.InputFile object at 0x187345ec>
    print_source_time_function = models.BooleanField()
    # skipped output-file <cig.addyndum.properties.OutputFile object at 0x1873460c>
    ntstep_between_output_info = models.IntegerField()
