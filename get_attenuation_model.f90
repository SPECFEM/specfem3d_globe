!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_attenuation_model(myrank,iregion_attenuation,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
         tau_mu,tau_sigma,beta,one_minus_sum_beta,factor_scale)

! return attenuation mechanisms Q_mu in PREM using standard linear solids
! the Tau values computed by Jeroen's code are used
! number of relaxation mechanisms: N_SLS = 3
! in the future when more memory is available on computers
! it would be more accurate to use four mechanisms instead of three

  implicit none

  include "constants.h"

  integer iregion_attenuation,myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  double precision, dimension(N_SLS) :: tau_mu,tau_sigma,beta
  double precision one_minus_sum_beta

  integer i

  double precision Q_mu,T_c_source,T_c_source_nondim,f_c_source,w_c_source
  double precision f_0_prem,factor_scale_mu0,factor_scale_mu,factor_scale
  double precision a_val,b_val,big_omega
  double precision scale_t

! check number of SLS is okay
  if(N_SLS /= 3) call exit_MPI(myrank,'wrong number of SLS for attenuation, must be 3')

! clear arrays
  tau_mu(:) = 0.d0
  tau_sigma(:) = 0.d0

  if(MAX_ATTENUATION_PERIOD == 200 .and. MIN_ATTENUATION_PERIOD == 1) then

! period range: 1.000000 -- 200.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   70.710678118654755d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =        31.83098861837910931172d0
  tau_sigma(2) =         2.25079079039276752638d0
  tau_sigma(3) =         0.15915494309189551214d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!! DK DK UGLY this for regional simulations at high frequency
!! DK DK UGLY !--- inner core, not used, but needs to be there even for regional code
!! DK DK UGLY !--- because fictitious mesh is created in inner core
!! DK DK UGLY !--- therefore we just use fictitious values

!--- CMB -> d670, target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670,IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =        31.98831234987152072335d0
    tau_mu(2) =         2.26372027015198540312d0
    tau_mu(3) =         0.16007126864753609685d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =        32.10936720723712056724d0
    tau_mu(2) =         2.27962803170602423819d0
    tau_mu(3) =         0.16117116738629347350d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =        32.15354690558234551645d0
    tau_mu(2) =         2.30398847208382084872d0
    tau_mu(3) =         0.16280124738435630682d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =        31.91896317140586347705d0
    tau_mu(2) =         2.25745362644383806838d0
    tau_mu(3) =         0.15962974749732541935d0
    Q_mu =       600.0000000000d0

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

!! DK DK UGLY this for global simulations at longer periods
!! DK DK UGLY should be merged with regional attenuation parameters at higher frequency above
!! DK DK UGLY in one cleaner general statement

  else if(MAX_ATTENUATION_PERIOD == 1000 .and. MIN_ATTENUATION_PERIOD == 20) then

! period range: 20.000000 -- 1000.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /    7.071067811865475d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =       159.15494309189548971517d0
  tau_sigma(2) =        22.50790790392767704020d0
  tau_sigma(3) =         3.18309886183791013181d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =       161.87151795592419034620d0
    tau_mu(2) =        22.83321303277247693586d0
    tau_mu(3) =         3.24415023587751250034d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =       159.94161523178539141554d0
    tau_mu(2) =        22.59299057282909117816d0
    tau_mu(3) =         3.19929580037721894570d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =       160.82504929484198896716d0
    tau_mu(2) =        22.69650297758172996510d0
    tau_mu(3) =         3.21878389936411624106d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =       162.01122326212433222281d0
    tau_mu(2) =        22.85290281779909804527d0
    tau_mu(3) =         3.24776851553577827758d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =       159.56830719443664179380d0
    tau_mu(2) =        22.55186839448693802979d0
    tau_mu(3) =         3.19148738794593134216d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 1000 .and. MIN_ATTENUATION_PERIOD == 8) then

! period range: 8.000000 -- 1000.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   11.180339887498947d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =       159.15494309189548971517d0
  tau_sigma(2) =        14.23525086834355768417d0
  tau_sigma(3) =         1.27323954473516454122d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =       161.46934711308600185475d0
    tau_mu(2) =        14.51296625431915998661d0
    tau_mu(3) =         1.29975705687209908135d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =       159.96474054645611317937d0
    tau_mu(2) =        14.30733506027169532615d0
    tau_mu(3) =         1.28028958294749606317d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =       160.75072907573880343080d0
    tau_mu(2) =        14.39559888678218868563d0
    tau_mu(3) =         1.28875837639885348906d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =       161.54392370131370171293d0
    tau_mu(2) =        14.52993576942056641599d0
    tau_mu(3) =         1.30132356306765761822d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =       159.59222087574624993067d0
    tau_mu(2) =        14.27243840130370422514d0
    tau_mu(3) =         1.27689204530675759486d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 1000 .and. MIN_ATTENUATION_PERIOD == 5) then

! period range: 5.000000 -- 1000.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   14.142135623730949d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =       159.15494309189548971517d0
  tau_sigma(2) =        11.25395395196383852010d0
  tau_sigma(3) =         0.79577471545947742193d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =       160.78443818062038417338d0
    tau_mu(2) =        11.50450764923400370776d0
    tau_mu(3) =         0.81299073871856486484d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =       159.94156174935702097173d0
    tau_mu(2) =        11.31860135075993589737d0
    tau_mu(3) =         0.80035634323768012344d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =       160.54683603618144616121d0
    tau_mu(2) =        11.39814015853017181712d0
    tau_mu(3) =         0.80585583693146656259d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =       160.76773452790945384550d0
    tau_mu(2) =        11.51994236041913133306d0
    tau_mu(3) =         0.81400623692178109003d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =       159.59481585703019845823d0
    tau_mu(2) =        11.28726813221918057195d0
    tau_mu(3) =         0.79814873748662718000d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 1000 .and. MIN_ATTENUATION_PERIOD == 4) then

! period range: 4.000000 -- 1000.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   15.811388300841896d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =       159.15494309189548971517d0
  tau_sigma(2) =        10.06584242089741820791d0
  tau_sigma(3) =         0.63661977236758215959d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =       160.43705131222463933227d0
    tau_mu(2) =        10.30203571530812034496d0
    tau_mu(3) =         0.65070081804940560488d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =       159.93113725698523808205d0
    tau_mu(2) =        10.12670852224042405965d0
    tau_mu(3) =         0.64036908310685669576d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =       160.44494348323934218570d0
    tau_mu(2) =        10.20166569793829047796d0
    tau_mu(3) =         0.64486747121131060556d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =       160.37362697346378581642d0
    tau_mu(2) =        10.31661392291149503819d0
    tau_mu(3) =         0.65153086897277423528d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =       159.59699380876537588847d0
    tau_mu(2) =        10.09720149757456475470d0
    tau_mu(3) =         0.63856272483804032980d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 400 .and. MIN_ATTENUATION_PERIOD == 20) then

! period range: 20.000000 -- 400.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   11.180339887498947d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =        63.66197723675821862344d0
  tau_sigma(2) =        14.23525086834355768417d0
  tau_sigma(3) =         3.18309886183791013181d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =        64.76417144248684110153d0
    tau_mu(2) =        14.36140908080261269220d0
    tau_mu(3) =         3.24062372548047905596d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =        63.96357795120979972125d0
    tau_mu(2) =        14.26855060142314002292d0
    tau_mu(3) =         3.19833140646331726131d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =        64.31765922582317784872d0
    tau_mu(2) =        14.30875228289421485783d0
    tau_mu(3) =         3.21668529475045783528d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =        64.82653618389130656396d0
    tau_mu(2) =        14.36895778355524200265d0
    tau_mu(3) =         3.24404313208959971249d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =        63.81900994371194002497d0
    tau_mu(2) =        14.25248721242562410794d0
    tau_mu(3) =         3.19098562562620902128d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 400 .and. MIN_ATTENUATION_PERIOD == 8) then

! period range: 8.000000 -- 400.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   17.677669529663682d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =        63.66197723675821862344d0
  tau_sigma(2) =         9.00316316157107365825d0
  tau_sigma(3) =         1.27323954473516409713d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =        64.74582414016121845179d0
    tau_mu(2) =         9.13298209844589337081d0
    tau_mu(3) =         1.29773265247803104572d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =        63.97698091926329055923d0
    tau_mu(2) =         9.03708785215422238934d0
    tau_mu(3) =         1.27973753823006330954d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =        64.32974236408334434145d0
    tau_mu(2) =         9.07838927486644919895d0
    tau_mu(3) =         1.28755597185630898949d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =        64.80118772325930365241d0
    tau_mu(2) =         9.14084899626674207695d0
    tau_mu(3) =         1.29918425004121362853d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =        63.82759066501840550245d0
    tau_mu(2) =         9.02068868832229142640d0
    tau_mu(3) =         1.27660489713826086344d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 400 .and. MIN_ATTENUATION_PERIOD == 5) then

! period range: 5.000000 -- 400.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   22.360679774997898d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =        63.66197723675821862344d0
  tau_sigma(2) =         7.11762543417177706573d0
  tau_sigma(3) =         0.79577471545947742193d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =        64.68991972277899549226d0
    tau_mu(2) =         7.23903487221359842607d0
    tau_mu(3) =         0.81170480859443183697d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =        63.98298909503787257336d0
    tau_mu(2) =         7.14922711793459608742d0
    tau_mu(3) =         0.80000541215927811756d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =        64.32290790448600148466d0
    tau_mu(2) =         7.18783008187737593175d0
    tau_mu(3) =         0.80509173507097531175d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =        64.73532400224027583135d0
    tau_mu(2) =         7.24642788081793653987d0
    tau_mu(3) =         0.81264737181571522484d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =        63.83260919765795904368d0
    tau_mu(2) =         7.13393777784677229903d0
    tau_mu(3) =         0.79796618094443061420d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else if(MAX_ATTENUATION_PERIOD == 400 .and. MIN_ATTENUATION_PERIOD == 4) then

! period range: 4.000000 -- 400.000000 s

! define central period of source in seconds using values from Jeroen's code
  T_c_source = 1000.d0 /   24.999999999999993d0

! tau sigma evenly spaced in log frequency, does not depend on value of Q
  tau_sigma(1) =        63.66197723675821862344d0
  tau_sigma(2) =         6.36619772367582115180d0
  tau_sigma(3) =         0.63661977236758215959d0

! check in which region we are based upon doubling flag

  select case(iregion_attenuation)

!--- inner core, target Q_mu:  84.60

  case(IREGION_ATTENUATION_INNER_CORE)

    tau_mu(1) =        64.65396574106728166953d0
    tau_mu(2) =         6.48254642872048325586d0
    tau_mu(3) =         0.64961295196815083131d0
    Q_mu =        84.60d0

!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.

  case(IREGION_ATTENUATION_CMB_670)

    tau_mu(1) =        63.98539935552845037137d0
    tau_mu(2) =         6.39645572624521818739d0
    tau_mu(3) =         0.64007238467197202780d0
    Q_mu =       312.0000000000d0

!--- d670 -> d220, target Q_mu: 143.

  case(IREGION_ATTENUATION_670_220)

    tau_mu(1) =        64.31665187797257488000d0
    tau_mu(2) =         6.43344580297075729902d0
    tau_mu(3) =         0.64422147261685436259d0
    Q_mu =       143.0000000000d0

!--- d220 -> depth of 80 km, target Q_mu:  80.

  case(IREGION_ATTENUATION_220_80)

    tau_mu(1) =        64.69346058437983515432d0
    tau_mu(2) =         6.48963824090475238648d0
    tau_mu(3) =         0.65038111965209866039d0
    Q_mu =        80.0000000000d0

!--- depth of 80 km -> surface, target Q_mu: 600.

  case(IREGION_ATTENUATION_80_SURFACE)

    tau_mu(1) =        63.83496368779294982687d0
    tau_mu(2) =         6.38181345959871659801d0
    tau_mu(3) =         0.63840836059004624214d0
    Q_mu =       600.0000000000d0

!--- do nothing for fluid outer core (no attenuation there)

  case default

    call exit_MPI(myrank,'wrong attenuation flag in mesh')

  end select

  else
    call exit_MPI(myrank,'incorrect minimum or maximum attenuation period')
  endif

!--- non-dimensionalize the tau values and the period of the source

    scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

    tau_mu(:) = tau_mu(:) / scale_t
    tau_sigma(:) = tau_sigma(:) / scale_t

!--- compute beta
    beta(:) = 1.d0 - tau_mu(:) / tau_sigma(:)

    T_c_source_nondim = T_c_source / scale_t

!--- compute central angular frequency of source (non dimensionalized)
    f_c_source = ONE / T_c_source_nondim
    w_c_source = TWO_PI * f_c_source

!--- non dimensionalize PREM reference of 1 second
    f_0_prem = ONE / ( ONE / scale_t)

!--- quantity by which to scale mu_0 to get mu
    factor_scale_mu0 = ONE + TWO * log(f_c_source / f_0_prem) / (PI * Q_mu)

!--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO
  one_minus_sum_beta = ONE

  do i = 1,N_SLS
    a_val = a_val - w_c_source * w_c_source * tau_mu(i) * &
      (tau_mu(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0)

!--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

!--- total factor by which to scale mu0
  factor_scale = factor_scale_mu * factor_scale_mu0

!--- check that the correction factor is close to one
  if(factor_scale < 0.9 .or. factor_scale > 1.1) &
    call exit_MPI(myrank,'incorrect correction factor in attenuation model')

  end subroutine get_attenuation_model

