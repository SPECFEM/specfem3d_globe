# Create your views here.

from django.shortcuts import render_to_response, get_object_or_404
from datetime import datetime
from django.http import HttpResponse, HttpResponseRedirect
from Specfem3DGlobe.web.Specfem3DGlobe.models import Mesh, Model, Simulation, UserInfo

def create_test_user():
    _userid = 'test_user'
    _lastname = 'test_user'
    _firstname = 'test_user'
    _email = ''
    _institution = ''
    _address1 = ''
    _address2 = ''
    _address3 = ''
    _phone = ''
    userinfo = UserInfo(
                userid      = _userid,
                lastname    = _lastname,
                firstname   = _firstname,
                email       = _email,
                institution = _institution,
                address1    = _address1,
                address2    = _address2,
                address3    = _address3,
                phone       = _phone
                )
    userinfo.save()

def getres_checkbox(request,key):
    if request.has_key(key) and request.POST[key] == 'on':
        return True
    else: 
        return False

def index(request):
    try:
        user = UserInfo.objects.get(userid='test_user')
    except UserInfo.DoesNotExist:
        create_test_user()
        user = UserInfo.objects.get(userid='test_user')
    prev_simulations = Simulation.objects.filter(user=user)
    
    return render_to_response('Specfem3DGlobe/home.html',
                             {'prev_simulations': prev_simulations})

def setparam(request):
    mesh_type = request.POST['mesh_type']
    if mesh_type == 'regional':
        template = 'Specfem3DGlobe/simulation_form_regional.html'
    else:
        template = 'Specfem3DGlobe/simulation_form_global.html'
    return render_to_response(template, 
                             {'simulation_cmt_solution': request.POST['simulation_cmt_solution'],
                              'simulation_type': request.POST['simulation_type']})

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

def create_simulation(request):
    # 
    # first, do some parameter validity checking!!!
    #
    # Do some checking here!



    #
    # mesh information
    #
    _nchunks = request.POST['mesh_nchunks']
    _nproc_xi = request.POST['mesh_nproc_xi']
    _nproc_eta = request.POST['mesh_nproc_eta']
    _nex_xi = request.POST['mesh_nex_xi']
    _nex_eta = request.POST['mesh_nex_eta']
    _save_files = getres_checkbox(request,'mesh_save_files')
    _type = 1
    if request.POST['mesh_type'] == 'regional':
        _type = 2
    _angular_width_eta = 0.0
    _angular_width_xi = 0.0
    _center_latitude = 0.0
    _center_longitude = 0.0
    _gamma_rotation_azimuth = 0.0
    if request.POST['mesh_type'] == 'regional':
        _angular_width_eta = request.POST['mesh_angular_width_xi']
        _angular_width_xi = request.POST['mesh_angular_width_eta']
        _center_latitude = request.POST['mesh_center_latitude']
        _center_longitude = request.POST['mesh_center_longitude']
        _gamma_rotation_azimuth = request.POST['mesh_gamma_rotation_azimuth']
        
    mesh = Mesh(    
                    nchunks                         = _nchunks, 
                    nproc_xi                        = _nproc_xi,
                    nproc_eta                       = _nproc_eta,
                    nex_xi                          = _nex_xi,
                    nex_eta                         = _nex_eta,
                    save_files                      = _save_files,
                    type                            = _type,
                    angular_width_eta               = _angular_width_eta,
                    angular_width_xi                = _angular_width_xi,
                    center_latitude                 = _center_latitude,
                    center_longitude                = _center_longitude,
                    gamma_rotation_azimuth          = _gamma_rotation_azimuth 
                )
    mesh.save()

    # 
    # model information
    #
    _type = request.POST['model_type']
    _oceans = getres_checkbox(request,'model_oceans')
    _gravity = getres_checkbox(request,'model_gravity')
    _attenuation = getres_checkbox(request,'model_attenuation')
    _topography = getres_checkbox(request,'model_topography')
    _rotation = getres_checkbox(request,'model_rotation')
    _ellipticity = getres_checkbox(request,'model_ellipticity')
    model = Model(  
                    type                            = _type,
                    oceans                          = _oceans,
                    gravity                         = _gravity,
                    attenuation                     = _attenuation,
                    topography                      = _topography,
                    rotation                        = _rotation,
                    ellipticity                     = _ellipticity 
                )
    model.save()

    # 
    # simulation information
    #
    _date = datetime.now()
    _status = 2
    _record_length = request.POST['simulation_record_length']
    _receivers_can_be_buried = getres_checkbox(request,'simulation_receivers_can_be_buried')
    _print_source_time_function = getres_checkbox(request,'simulation_print_source_time_function')
    _save_forward = False
    _movie_surface = getres_checkbox(request,'simulation_movie_surface')
    _movie_volume = getres_checkbox(request,'simulation_movie_volume')
    _cmt_solution = request.POST['simulation_cmt_solution']
    _stations = ''
    _absorbing_conditions = False
    if request.POST['mesh_type'] == 'regional':
        _absorbing_conditions = True
    _ntstep_between_frames = 100
    _simulation_type = request.POST['simulation_type']
    try:
        user = UserInfo.objects.get(userid='test_user')
    except UserInfo.DoesNotExist:
        create_test_user()
        user = UserInfo.objects.get(userid='test_user')
         
    simulation = Simulation(
                    user                            = user,
                    date                            = _date,
                    mesh                            = mesh,
                    model                           = model,
                    status                          = _status,
                    record_length                   = _record_length,
                    receivers_can_be_buried         = _receivers_can_be_buried,
                    print_source_time_function      = _print_source_time_function,
                    save_forward                    = _save_forward,
                    movie_surface                   = _movie_surface,
                    movie_volume                    = _movie_volume,
                    cmt_solution                    = _cmt_solution,
                    stations                        = _stations,
                    absorbing_conditions            = _absorbing_conditions,
                    ntstep_between_frames           = _ntstep_between_frames,
                    simulation_type                 = _simulation_type,
                )
    simulation.save()

    return HttpResponseRedirect('/specfem3dglobe/')

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
