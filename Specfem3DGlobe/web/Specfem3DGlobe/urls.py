

from django.conf.urls.defaults import *
from models import Mesh, Model, Simulation, UserInfo


# Mesh

mesh_list_detail_args = {
	'queryset': Mesh.objects.all(),
	'allow_empty': True,
}

mesh_create_update_args = {
	'model': Mesh,
	'post_save_redirect': '/specfem3dglobe/meshes/',
}

mesh_delete_args = {
	'model': Mesh,
	'post_delete_redirect': '/specfem3dglobe/meshes/',
}



# Model

model_list_detail_args = {
	'queryset': Model.objects.all(),
	'allow_empty': True,
}

model_create_update_args = {
	'model': Model,
	'post_save_redirect': '/specfem3dglobe/models/',
}

model_delete_args = {
	'model': Model,
	'post_delete_redirect': '/specfem3dglobe/models/',
}



# Simulation

simulation_list_detail_args = {
	'queryset': Simulation.objects.all(),
	'allow_empty': True,
}

simulation_create_update_args = {
	'model': Simulation,
	'post_save_redirect': '/specfem3dglobe/simulations/',
}

simulation_delete_args = {
	'model': Simulation,
	'post_delete_redirect': '/specfem3dglobe/simulations/',
}


# Registration

password_change_args = {
    'template_name': 'Specfem3DGlobe/password_change_form.html',
    }



# URLs

urlpatterns = patterns('',
	(r'^$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.index'),
	(r'^detail/(?P<sim_id>\d+)/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.detail'),
	(r'^delete/(?P<sim_id>\d+)/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.delete'),
	(r'^wizard/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.create_simulation'),
	(r'^info/(?P<info_str>\w+)/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.info'),

	(r'^login/$', 'django.contrib.auth.views.login'),
	(r'^logout/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.logout_view'),
	(r'^registration/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.registration'),
	(r'^registration/password/$', 'django.contrib.auth.views.password_change', password_change_args),
	(r'^registration/password/done/$', 'django.views.generic.simple.direct_to_template', { 'template': 'Specfem3DGlobe/password_changed.html'}),

	(r'^meshes/$', 'django.views.generic.list_detail.object_list', mesh_list_detail_args),
	(r'^meshes/create/$', 'django.views.generic.create_update.create_object', mesh_create_update_args),
	(r'^meshes/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', mesh_create_update_args),
	(r'^meshes/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', mesh_delete_args),

	(r'^models/$', 'django.views.generic.list_detail.object_list', model_list_detail_args),
	(r'^models/create/$', 'django.views.generic.create_update.create_object', model_create_update_args),
	(r'^models/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', model_create_update_args),
	(r'^models/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', model_delete_args),

	(r'^simulations/$', 'django.views.generic.list_detail.object_list', simulation_list_detail_args),
	(r'^simulations/create/$', 'django.views.generic.create_update.create_object', simulation_create_update_args),
	(r'^simulations/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', simulation_create_update_args),
	(r'^simulations/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', simulation_delete_args),
	(r'^simulations/(?P<sim_id>\d+)/parameters.pml$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.simulation_pml'),
	(r'^simulations/(?P<sim_id>\d+)/stations.txt$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.stations_txt'),
	(r'^simulations/(?P<sim_id>\d+)/events.txt$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.events_txt'),

	(r'^events/', include('cig.web.seismo.events.urls')),
	(r'^stations/', include('cig.web.seismo.stations.urls')),

)
