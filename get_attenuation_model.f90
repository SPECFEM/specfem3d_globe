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

 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =        32.01772419813112691145d0
     tau_mu(2) =         2.26295014032926911085d0
     tau_mu(3) =         0.16009837572177521015d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =        32.24022626127403157170d0
     tau_mu(2) =         2.27753979117215532568d0
     tau_mu(3) =         0.16124848521135615176d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =        32.56637722770504694836d0
     tau_mu(2) =         2.29935732671327297538d0
     tau_mu(3) =         0.16298642034932270262d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =        31.92809179152154541725d0
     tau_mu(2) =         2.25705689329478609295d0
     tau_mu(3) =         0.15964311991416568759d0
          Q_mu =       600.0000000000d0
 
 !--- do nothing for fluid outer core (no attenuation there)
 
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =       162.14490233390495177446d0
     tau_mu(2) =        22.82419389065037051978d0
     tau_mu(3) =         3.24514126878853792491d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =       159.96238677589397525480d0
     tau_mu(2) =        22.59200103328783626466d0
     tau_mu(3) =         3.19941840854237469216d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =       160.91977214000002049943d0
     tau_mu(2) =        22.69300216006637782584d0
     tau_mu(3) =         3.21918230720303011339d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =       162.31797761220104803215d0
     tau_mu(2) =        22.84287987713201317774d0
     tau_mu(3) =         3.24886786464258259244d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =       159.57440558551263620757d0
     tau_mu(2) =        22.55150522635637955204d0
     tau_mu(3) =         3.19153381437414740418d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =       162.42071935356284484442d0
     tau_mu(2) =        14.49929877579018722145d0
     tau_mu(3) =         1.30045337980719555304d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =       160.03319855314362030185d0
     tau_mu(2) =        14.30549165803078537351d0
     tau_mu(3) =         1.28036835479182142805d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =       161.07811629420274357471d0
     tau_mu(2) =        14.38975485348623095661d0
     tau_mu(3) =         1.28906206911127751980d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =       162.61641330522928683422d0
     tau_mu(2) =        14.51449312134383973216d0
     tau_mu(3) =         1.30227414871675950536d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =       159.61095844748751915176d0
     tau_mu(2) =        14.27166760653473431830d0
     tau_mu(3) =         1.27691138668395165467d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =       162.62983680848708445410d0
     tau_mu(2) =        11.48310665737810865039d0
     tau_mu(3) =         0.81383338355291257038d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =       160.08859189474421214072d0
     tau_mu(2) =        11.31475189095645461634d0
     tau_mu(3) =         0.80049216168943215788d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =       161.20153794984145179114d0
     tau_mu(2) =        11.38755447047194735433d0
     tau_mu(3) =         0.80624514209892939043d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =       162.83211429812064352518d0
     tau_mu(2) =        11.49677225353980247746d0
     tau_mu(3) =         0.81492710793937817026d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =       159.64027081362178250856d0
     tau_mu(2) =        11.28529825778758066690d0
     tau_mu(3) =         0.79821444289159093621d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =       162.72949947544336168903d0
     tau_mu(2) =        10.28008029214220542258d0
     tau_mu(3) =         0.65148722307895390315d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =       160.11586020478787872889d0
     tau_mu(2) =        10.12237643559230093615d0
     tau_mu(3) =         0.64050486320597954659d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =       161.25959323778928933280d0
     tau_mu(2) =        10.19066637823843279875d0
     tau_mu(3) =         0.64523504177037338536d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =       162.93762739129357441925d0
     tau_mu(2) =        10.29287588454070245803d0
     tau_mu(3) =         0.65238830508968748134d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =       159.65390567402135957309d0
     tau_mu(2) =        10.09509723750865362035d0
     tau_mu(3) =         0.63862664710547289992d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =        64.74754934844317233456d0
     tau_mu(2) =        14.36897611143563935343d0
     tau_mu(3) =         3.23900252215628192687d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =        63.95536976467237622046d0
     tau_mu(2) =        14.27086029471458417106d0
     tau_mu(3) =         3.19788010916230369673d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =        64.30294233312994833796d0
     tau_mu(2) =        14.31362245566047697309d0
     tau_mu(3) =         3.21568666661526503248d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =        64.81024764817939853856d0
     tau_mu(2) =        14.37685795447098513478d0
     tau_mu(3) =         3.24234250339832286159d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =        63.81446317410515689517d0
     tau_mu(2) =        14.25369879719033505694d0
     tau_mu(3) =         3.19075364624547574977d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =        64.85795595561907589399d0
     tau_mu(2) =         9.12968077134321731592d0
     tau_mu(3) =         1.29805586567704822620d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =        63.98498337148299697219d0
     tau_mu(2) =         9.03680316300993347056d0
     tau_mu(3) =         1.27976700282590316604d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =        64.36796404728389120464d0
     tau_mu(2) =         9.07719787898379415481d0
     tau_mu(3) =         1.28767253698791295236d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =        64.92718894441718191501d0
     tau_mu(2) =         9.13714835179053253000d0
     tau_mu(3) =         1.29954762980627780422d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =        63.82978503832104166804d0
     tau_mu(2) =         9.02060111267174136174d0
     tau_mu(3) =         1.27661302844902913023d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =        64.91539550103743749787d0
     tau_mu(2) =         7.23443108987387084596d0
     tau_mu(3) =         0.81204067301336624318d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =        63.99919838292981921768d0
     tau_mu(2) =         7.14875688323272839853d0
     tau_mu(3) =         0.80003162995834964377d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =        64.40011413556405273084d0
     tau_mu(2) =         7.18609902209461903766d0
     tau_mu(3) =         0.80520508161117398949d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =        64.98786468038585439899d0
     tau_mu(2) =         7.24140408806958291166d0
     tau_mu(3) =         0.81302347647883421722d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =        63.83708402275915716473d0
     tau_mu(2) =         7.13376182745341047564d0
     tau_mu(3) =         0.79797505252482325844d0
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

 !--- inner core, target Q_mu: 84.60
 
   case(IREGION_ATTENUATION_INNER_CORE)
 
     tau_mu(1) =        64.93473007621230408404d0
     tau_mu(2) =         6.47784347577742991575d0
     tau_mu(3) =         0.64992285325701537602d0
          Q_mu =        84.6000000000d0
 
 !--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.
 
   case(IREGION_ATTENUATION_CMB_670)
 
     tau_mu(1) =        64.00553853880884958016d0
     tau_mu(2) =         6.39587371477192512259d0
     tau_mu(3) =         0.64009931875558101488d0
          Q_mu =       312.0000000000d0
 
 !--- d670 -> d220, target Q_mu: 143.
 
   case(IREGION_ATTENUATION_670_220)
 
     tau_mu(1) =        64.41378421517728725121d0
     tau_mu(2) =         6.43148896649217505228d0
     tau_mu(3) =         0.64432961718299808229d0
          Q_mu =       143.0000000000d0
 
 !--- d220 -> depth of 80 km, target Q_mu:  80.
 
   case(IREGION_ATTENUATION_220_80)
 
     tau_mu(1) =        65.00839889844625929527d0
     tau_mu(2) =         6.48444159417360488362d0
     tau_mu(3) =         0.65073008539433319086d0
          Q_mu =        80.0000000000d0
 
 !--- depth of 80 km -> surface, target Q_mu: 600.
 
   case(IREGION_ATTENUATION_80_SURFACE)
 
     tau_mu(1) =        63.84050388187188929123d0
     tau_mu(2) =         6.38157285455825551423d0
     tau_mu(3) =         0.63841788042370528622d0
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

