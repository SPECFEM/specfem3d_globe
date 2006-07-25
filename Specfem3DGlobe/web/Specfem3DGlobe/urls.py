

from django.conf.urls.defaults import *
from models import  Mesher, Model, Solver

# Mesher

mesher_list_detail_args = {
    'queryset': Mesher.objects.all(),
    'allow_empty': True,
}

mesher_create_update_args = {
    'model': Mesher,
    'post_save_redirect': '/specfem3dglobe/mesher/',
    }

mesher_delete_args = {
    'model': Mesher,
    'post_delete_redirect': '/specfem3dglobe/mesher/',
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



# Solver

solver_list_detail_args = {
    'queryset': Solver.objects.all(),
    'allow_empty': True,
}

solver_create_update_args = {
    'model': Solver,
    'post_save_redirect': '/specfem3dglobe/solver/',
    }

solver_delete_args = {
    'model': Solver,
    'post_delete_redirect': '/specfem3dglobe/solver/',
    }



# URLs

urlpatterns = patterns('',
    (r'^$', 'django.views.generic.simple.direct_to_template', { 'template': 'Specfem3DGlobe/home.html' }),

    (r'^mesher/$', 'django.views.generic.list_detail.object_list', mesher_list_detail_args),
    (r'^mesher/create/$', 'django.views.generic.create_update.create_object', mesher_create_update_args),
    (r'^mesher/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', mesher_create_update_args),
    (r'^mesher/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', mesher_delete_args),


    (r'^model/$', 'django.views.generic.list_detail.object_list', model_list_detail_args),
    (r'^model/create/$', 'django.views.generic.create_update.create_object', model_create_update_args),
    (r'^model/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', model_create_update_args),
    (r'^model/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', model_delete_args),


    (r'^solver/$', 'django.views.generic.list_detail.object_list', solver_list_detail_args),
    (r'^solver/create/$', 'django.views.generic.create_update.create_object', solver_create_update_args),
    (r'^solver/(?P<object_id>\d+)/$', 'django.views.generic.create_update.update_object', solver_create_update_args),
    (r'^solver/(?P<object_id>\d+)/delete/$', 'django.views.generic.create_update.delete_object', solver_delete_args),


)
