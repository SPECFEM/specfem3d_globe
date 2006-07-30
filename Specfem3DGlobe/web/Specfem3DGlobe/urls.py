

from django.conf.urls.defaults import *
from models import  Mesh, Model, Simulation

# Mesh

mesh_list_detail_args = {
    'queryset': Mesh.objects.all(),
    'allow_empty': True,
}

mesh_create_update_args = {
    'model': Mesh,
    'post_save_redirect': '/specfem3dglobe/mesh/',
    }

mesh_delete_args = {
    'model': Mesh,
    'post_delete_redirect': '/specfem3dglobe/mesh/',
    }



# Model

model_list_detail_args = {
    'queryset': Model.objects.all(),
    'allow_empty': True,
}

model_create_update_args = {
    'model': Model,
    'post_save_redirect': '/specfem3dglobe/model/',
    }

model_delete_args = {
    'model': Model,
    'post_delete_redirect': '/specfem3dglobe/model/',
    }



# Simulation

simulation_list_detail_args = {
    'queryset': Simulation.objects.all(),
    'allow_empty': True,
}

simulation_create_update_args = {
    'model': Simulation,
    'post_save_redirect': '/specfem3dglobe/simulation/',
    }

simulation_delete_args = {
    'model': Simulation,
    'post_delete_redirect': '/specfem3dglobe/simulation/',
    }



# URLs

urlpatterns = patterns('',
    (r'^$', 'mysite.Specfem3DGlobe.web.Specfem3DGlobe.views.index'),
    (r'^setparam/$', 'mysite.Specfem3DGlobe.web.Specfem3DGlobe.views.setparam'),
    (r'^detail/(?P<sim_id>\d+)/$', 'mysite.Specfem3DGlobe.web.Specfem3DGlobe.views.detail'),
    (r'^delete/(?P<sim_id>\d+)/$', 'mysite.Specfem3DGlobe.web.Specfem3DGlobe.views.delete'),
    (r'^create_simulation/$', 'mysite.Specfem3DGlobe.web.Specfem3DGlobe.views.create_simulation'),

    (r'^mesh/$', 'django.views.generic.list_detail.object_list', mesh_list_detail_args),
    (r'^mesh/create/$', 'django.views.generic.create_update.create_object', mesh_create_update_args),
    (r'^mesh/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', mesh_create_update_args),
    (r'^mesh/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', mesh_delete_args),


    (r'^model/$', 'django.views.generic.list_detail.object_list', model_list_detail_args),
    (r'^model/create/$', 'django.views.generic.create_update.create_object', model_create_update_args),
    (r'^model/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', model_create_update_args),
    (r'^model/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', model_delete_args),


    (r'^simulation/$', 'django.views.generic.list_detail.object_list', simulation_list_detail_args),
    (r'^simulation/create/$', 'django.views.generic.create_update.create_object', simulation_create_update_args),
    (r'^simulation/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', simulation_create_update_args),
    (r'^simulation/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', simulation_delete_args),


)
