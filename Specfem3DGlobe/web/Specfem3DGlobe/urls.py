

from django.conf.urls.defaults import *
from models import  Mesh, Model, Simulation

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



# URLs

urlpatterns = patterns('',
    (r'^$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.index'),
    (r'^setparam/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.setparam'),
    (r'^detail/(?P<sim_id>\d+)/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.detail'),
    (r'^delete/(?P<sim_id>\d+)/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.delete'),
    (r'^create_simulation/$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.create_simulation'),

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
    (r'^simulations/(?P<sim_id>\d+).pml$', 'Specfem3DGlobe.web.Specfem3DGlobe.views.simulation_pml'),


)
