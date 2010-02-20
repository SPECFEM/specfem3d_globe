--------------------------------
readme
--------------------------------

gnuplot
http://www.gnuplot.info/



- contains sources to plot PREM model in directory:
  Gnuplot/display_PREM_model

- contains source file to create gnuplot scripts of regular 2D GLL elements:
  draw_regular_2D_GLL_element_gnuplot.f90
    
    
- basic procedure to execute gnuplot scripts:

  simple example to plot a single seismogram: 
  
  run:

  > gnuplot

  then either use
  
    A) a provided script:
    
      gnuplot> load "PAS.TS.LHZ.gnu"

  or
  
    B) the gnuplot command plot:
    
    gnuplot> plot 'OUTPUT_FILES/PAS.TS.LHZ.sem.ascii' w l
    
    
  Writing the figure to a file may depend on your gnuplot settings.
    
  One example to generate a postscript is shown below:

    gnuplot> plot 'OUTPUT_FILES/PAS.TS.LHZ.sem.ascii' w l
    gnuplot> set term postscript color solid
    
         Terminal type set to 'postscript'
         Options are 'landscape noenhanced color colortext \
            solid dashlength 1.0 linewidth 1.0 defaultplex \
            palfuncparam 2000,0.003 \
            butt "Helvetica" 14'
            
    gnuplot> set output 'my_figure.ps'
    gnuplot> replot
    gnuplot> quit

