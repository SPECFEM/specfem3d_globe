
  integer, parameter :: mxhpar=2
  integer, parameter :: mxkern=200
  integer, parameter :: mxcoef=2000

  character(len=80) refmodel
  character(len=80) kernstri
  character(len=40) desckern(mxkern)
  character(len=80) hsplfile(mxhpar)

  integer ihorpar(mxkern)
  integer ityphpar(mxhpar)
  integer ixlspl(mxcoef,mxhpar)
  integer lmaxhor(mxhpar)
  integer ncoefhor(mxhpar)

  real(kind=4) coef(mxcoef,mxkern)
  real(kind=4) xlaspl(mxcoef,mxhpar)
  real(kind=4) xlospl(mxcoef,mxhpar)
  real(kind=4) xraspl(mxcoef,mxhpar)

  common /tdmodl/ nmodkern,nhorpar,ityphpar, &
    ihorpar,lmaxhor,ncoefhor, &
    xlaspl,xlospl,xraspl,ixlspl,coef, &
    hsplfile,refmodel,kernstri,desckern

