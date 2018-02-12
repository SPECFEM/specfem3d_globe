!(doc)TITLE{LIB\_VTK\_IO}
!(doc)SUBTITLE{VTK InPut/OutPut Fortran Library}
!(doc)VERSION{0.2}
!(doc)AUTHOR{Stefano Zaghi}
!(doc)COAUTHORS{Enrico Cavallini and Renato N. Elias}
!(doc)DATE{07-01-2008}

!!\newcommand{\LIBVTKIO}{\MaiuscolettoBS{LIB\_VTK\_IO }}

!(doc)header

!(doc)titlepage

!!\tableofcontents

!!\chapter{Acknowledgements}
!!\label{cap:Acknowledgements}
!!
!!I am very grateful to Renato N. Elias: whitout his support \LIBVTKIO would not born. As a matter of facts \LIBVTKIO is a
!!collection of his tips. Despite the fact that Renato does not write the code he is a \virgo{moral co-author} of the code.
!!
!!I thank Enrico Cavallini for his help in debugging the code. He also develop the MS Windows version
!!of \LIBVTKIO and he is the first co-author that I found.
!!
!!Finally I thank the ParaView mailing list for the great support of its members.
!!
!!\chapter{Introduction}
!!\label{cap:Introduction}
!!\begin{epigraphs}
!! \qitem{\emph{I have not failed. I've just found ten thousand ways that don't work.}}{{\sc Thomas Edison}}
!!\end{epigraphs}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf L}}{IB\_VTK\_IO} is a Fortran library to write and read
!!(actually only to write) data conforming the VTK standard both binary and ascii. Even though there are many
!!wrappers/porting of the VTK source code (C++ code), there is not a Fortran one. This library is not a porting
!!or a wrapper of the VTK code, but it only an exporter/importer of the VTK data format written in pure Fortran
!!language (standard Fortran 95 with some extensions of non standard Fortran 2003) that can be used by Fortran
!!coders (yes, there are still a lot of these brave coders...) without mixing Fortran with C++ language.
!!
!!The library is still in developing and testing, this is first usable release, but there are not all the features
!!of the stable release (the importer is totaly absent and the exporter is not complete). Surely there are a lot of
!!bugs and the progamming style is not the best, but the exporter is usable for the 90\% of the VTK data format.
!!
!!The \LIBVTKIO is an open source project, it is distribuited under the GPL v3 (see appendix \ref{cap:GPL}). Anyone is
!!interest to use, to develop or contribuite to \LIBVTKIO is welcome.
!!
!!\section*{VTK Standard}
!!\label{sec:VTK Standard}
!!VTK, Visualization Toolkit, is an open source software that provides a powerful framework for the computer grafich, for
!!the images processing and for 3D rendering. It is widely used in the world and so it has a very large comunity of users;
!!besides the Kitware\footnote{The Kitware homepage can be found here: \href{http://public.kitware.com}{http://public.kitware.com}.}
!!company provides professional support. The toolkit is written in C++ and a lot of porting/wrappers for Tcl/Tk,
!!Java and Python are provided; unlucky there aren't wrappers for Fortran.
!!
!!Because of its good features the VTK toolkit has been used to develop a large set of open source programs. For my work
!!the most important family of programs is the scientific visualization programs. A lot of high-quality scientific visualization
!!tool are available on the web but for me the best is ParaView: I think that it is one of the best scintific visualization
!!program in the world and it is open source! Paraview is based on VTK.
!!
!!\section*{Paraview}
!!\label{sec:Paraview}
!!ParaView\footnote{The ParaView homepage can be found here: \href{http://www.paraview.org}{http://www.paraview.org}.}
!!is an open source software voted to scientific visualization and able to use the power of parallel architectures. It
!!has an architecture client-server in order to make easy the remote visualization of very large set of data. Because it is based
!!on VTK it inherits all VTK features. ParaView is very useful for Computational Fluid Dynamics visualizations because it provides
!!powerful post-processing tools; it provides a very large set of importers for the most used format like Plot3D and HDF (the list
!!is very large). It is easy to extend ParaView because it supports all the scripting language supported by VTK.
!!
!!\section*{LIB\_VTK\_IO}
!!\label{sec:LIB_VTK_IO}
!!Even though the VTK toolkit is written in C++ and so it is possible to use it in mixed Fortran/c++ code this is not the easiest
!!way. Fortran is still the best language for high performance computing for scientific purpose, like CFD computing. It necessary a
!!tool to deal with VTK standard directly by Fortran code. The library \LIBVTKIO was made to fill this empty: it is a simple
!!Fortran module able to export native Fortran data into VTK data format and to import VTK data into a Fortran code (actually this
!!feature is missing), both in ascii and binary file format.
!!
!!The library provides an automatic way to deal with VTK data format: all the formatting processes is nested into the library and
!!the users comunicate with it by a simple API passing only native Fortran data (native Fortran scalar, vector and matrix).
!!
!!The library \LIBVTKIO is distribuited under the GNU GPL v3 license (see appendix \ref{cap:GPL}). Beyond to the source code there
!!are some precompiled binaries for GNU-Linux (amd x86, amd x86\_64, intel x86, intel x86\_64) and WindowsXP (amd x86, intel x86).
!!
!!Actually the library is still in developing/testing phase (a lot of features are missing); this is not a stable release, but the
!!exporter is quite complete and its API is quite stable. The exporter is usable and I use it for my work.
!!
!!It can be found at: \href{http://stefano.zaghi.googlepages.com/lib\_vtk\_io}{http://stefano.zaghi.googlepages.com/lib\_vtk\_io}.
!!
!!\chapter{News and Changes}
!!\label{cap:NewsChanges}
!!
!!\section*{Version v0.2}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} version v0.2 is the second testing release. From version v0.1 there are
!!only minor changes; this new version does not introduce new features and does not fix bugs: it is a simple code-cleaning. The
!!character variables are now case-insensitive; the names of some variables have been changed. The comments have been translated in
!!English (very poor translation...).
!!
!!\begin{boxred}{List of changes from v0.1}
!!\begin{enumerate1Red}
!!  \item variable {\color{Maroon}formato} is changed in {\color{Maroon}output\_format} and now appears only in VTK\_INI and
!!        in VTK\_INI\_XML.
!!  \item variable {\color{Maroon}nomefile} is changed in {\color{Maroon}filename}.
!!  \item variable {\color{Maroon}titolo} is changed in {\color{Maroon}title}.
!!  \item variable {\color{Maroon}topologia} is changed in {\color{Maroon}mesh\_topology} and now appears only in VTK\_INI and
!!        in VTK\_INI\_XML.
!!  \item variable {\color{Maroon}NCelle} is changed in {\color{Maroon}NC}.
!!  \item variable {\color{Maroon}Nnodi} is changed in {\color{Maroon}NN}.
!!  \item variable {\color{Maroon}tipo} in VTK\_CON and VTK\_CON\_XML is changed in {\color{Maroon}cell\_type}.
!!  \item variable {\color{Maroon}tipo} in VTK\_DAT and VTK\_DAT\_XML is changed in {\color{Maroon}var\_location}.
!!  \item variable {\color{Maroon}azione} in VTK\_DAT\_XML is changed in {\color{Maroon}var\_block\_action}.
!!  \item variable {\color{Maroon}tipo} in VTK\_VAR is changed in {\color{Maroon}vec\_type}.
!!  \item variable {\color{Maroon}nomevar} is changed in {\color{Maroon}varname}.
!!\end{enumerate1Red}
!!\end{boxred}
!!
!!The only relevant news in the v0.2 version is about this guide: now the guide is integrated in the code. The code has particular
!!comments: if the code is processed by the program FortranDOC\footnote{FortranDOC is an open-source Fortran code available at:
!!\href{http://stefano.zaghi.googlepages.com/Fortrandoc}{http://stefano.zaghi.googlepages.com/Fortrandoc}. This code processing a
!!free-format Fortran code generates a corresponding pretty-latex documentation file of the code structure.} a latex source of
!!this guide will be made; compiling the latex file with \virgo{pdflatex} you will obtain this guide in PDF.
!!
!!\section*{Version v0.1}
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} version v0.1 is the first testing release. There are not news and changes.
!!
!!\mainmatter
!!
!!\part{Compile and Install LIB\_VTK\_IO}
!!\label{part:Compile and Install}
!!
!!\chapter{Compile LIB\_VTK\_IO}
!!\label{cap:Compiling Library}
!!\minitoc
!!\vspace*{3mm}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} \LIBVTKIO is open source and so anyone is encouraged to use the source
!!code and to \virgo{patch} it.
!!
!!The code is written in Fortran: the standard adopted is the Fortran 95 standard that is a minor upgrade to the Fortran 90 standard
!!and that is widely supported by the almost all compilers actually available. Unluckily Fortran 95 does not allow the creation of
!!C-binary file (Fortran inserts some bytes before and after each records despite the C standard) that is the standard adopted by
!!VTK. Therefore in order to create binary files that are compatible whit VTK standard the only way is to use a non-standard 95
!!instructions. At today only Fortran 2003 can create C-binary file, but there are not any compilers that completely implement this
!!standard. In the next year (2008) maybe a new minor upgrade of Fortran standard (unofficial named Fortran 2008) will be born
!!and so the support to Fortran 2003/2008 probably will be improved. Luckily we need to use only some features of Fortran 2003
!!that are supported by many compilers.
!!
!!The Fortran 2003 instructions are focused on the opening of the binary file, in particular in the functions
!!\MaiuscolettoBS{VTK\_INI} and \MaiuscolettoBS{VTK\_INI\_XML}. In these functions there are opening instructions like the following:
!!
!!\begin{boxred}{Fortran 2003 instructions}
!!\begin{verbatim}
!!open(unit       = ..., &
!!     file       = ..., &
!!     form       = ..., &
!!     access     = ..., &
!!     action     = ..., &
!!     convert    = 'BIG_ENDIAN', &
!!     recordtype = 'STREAM', &
!!     buffered   = 'YES', &
!!     iostat     = ...)
!!\end{verbatim}
!!\end{boxred}
!!
!!The specifiers \MaiuscolettoBS{convert}, \MaiuscolettoBS{recordtype} and \MaiuscolettoBS{buffered} are non standard for Fortran 95.
!!The \MaiuscolettoBS{buffered} specifier is not necessary and so can be commented or eliminated. The specifiers
!!\MaiuscolettoBS{convert} and \MaiuscolettoBS{recordtype} are instead necessary to write binary file but can be replaced by other
!!specifiers/instructions. In particular an alternative is opening the file with the specifier
!!\MaiuscolettoBS{form = BINARY}\footnote{Remember that also the value \MaiuscolettoBS{BINARY} for form specifier is non standard
!!for Fortran 95.} and using a compiler's option\footnote{Each compilers adopt differents option to achieve conversion of bytes
!!order (if it allows conversion). See the user guide of your compiler. Intel Fortran allows the conversion both by open specifier
!!and by compiling option.} to ensure the \MaiuscolettoBS{BIG\_ENDIAN} encoding. \MaiuscolettoBS{BIG\_ENDIAN} encoding is strictly
!!necessary only for legacy binary file; for XML binary file one can choice also the \MaiuscolettoBS{LITTLE\_ENDIAN} and so the
!!conversion is not necessary.
!!
!!Actually there is also another instruction that is non-standard for Fortran 95: the instruction \MaiuscolettoBS{sizeof}. This
!!instruction is used to comptuing the number of bytes of the saved data in the XML binary files. Maybe there are others
!!alternatives that are Fortran 95 compatible but at the moment I have not the time to implement them.
!!
!!Before you compile \LIBVTKIO ensure that your compiler allows these Fortran 2003 extensions. I use the Intel Fortran
!!Compiler\footnote{\href{http://www.intel.com}{http://www.intel.com}.} that is free for non-commercial use and it has a strong
!!support for Fortran 2003.
!!
!!\section{Compile under GNU/Linux}
!!\label{sec:CompileLinux}
!!
!!\LIBVTKIO can be compiled as a stand-alone library or it can be integrated directly in your code. It is a self-contained module
!!that can be safely included into others Fortran codes. There are no any advices for compile \LIBVTKIO excluding the above non
!!standard instructions.
!!
!!For the GNU/Linux users there is available a makefile already set to compile \LIBVTKIO both as static and dynamic library with
!!Intel Fortran. The makefile has only one option: \MaiuscolettoBS{SHARED}. This variable (default set to \virgo{no}) can assume
!!two values:
!!\begin{enumerate1Blu}
!!\item {\color{RoyalBlue}\MaiuscolettoBS{no}}:  makefile creates a \MaiuscolettoBS{static} library
!!\item {\color{RoyalBlue}\MaiuscolettoBS{yes}}: makefile creates a \MaiuscolettoBS{dynamic} library
!!\end{enumerate1Blu}
!!
!!\section{Compile under MS Windows}
!!\label{sec:CompileWin}
!!
!!For MS Windows users there is not any support at the moment. As soon as I have the time I will make available a MS Visual Studio
!!Project to compile \LIBVTKIO with Intel Visual Fortran for Windows.
!!
!!\clearpage
!!
!!\chapter{Install and Link (Pre)Compiled LIB\_VTK\_IO}
!!\label{cap:Install and Linking}
!!\minitoc
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf T}}{he} \LIBVTKIO is distribuited in two different version (other than source code): the first is a static linking version (extensions are \emph{.a} and \emph{.lib}) and the second is dynamic linking version (extensions are \emph{.so} and \emph{.dll}). The use of these two version is different and it depends on the OS used. The library is been tested only on GNU/Linux (several different distro) and on MS Windows (Windows XP).
!!
!!The library is distribuited with two different archive: \MaiuscolettoBS{LIB\_VTK\_IO-bin-x.x}.tar for GNU/Linux systems and
!!\MaiuscolettoBS{LIB\_VTK\_IO-bin-x.x}.zip for MS Windows systems. Into the archives there is the source code of the library
!!(\MaiuscolettoBS{LIB\_VTK\_IO}.f90), there are both static and dynamic version of the librabry and there is also this guide
!!(\MaiuscolettoBS{LIB\_VTK\_IO\_Guide}.pdf).
!!
!!\section{GNU/Linux}
!!
!!\subsection{Static Library}
!!\label{sec:Linux Static}
!!The static version of the precompiled library (\MaiuscolettoBS{LIB\_VTK\_IO}.a) does not require any kind of installations. It is
!!enough to link against it in the linking phase. It is important to use the interface module \emph{lib\_vtk\_io.mod} distribuited
!!with the library: this is the interface of the subroutines and functions that constitute the library.
!!
!!To use the functions and subroutines of the library it is mandatory to \MaiuscolettoBS{USE} the module. Suppose one has a program
!!(or subprogram) named \emph{test} that use the library; the correct \MaiuscolettoBS{USE} is:
!!
!!\begin{boxred}{The \LIBVTKIO must to be loaded with the USE statement}
!!\begin{verbatim}
!!program test
!!USE LIB_VTK_IO
!!...
!!...
!!...
!!endprogram test
!!\end{verbatim}
!!\end{boxred}
!!
!!With the instruction \verb|USE LIB\_VTK\_IO| the program \emph{test} can use the functions and subroutines of the library. To
!!compile, without link, this code one must give the module interface \emph{lib\_vtk\_io.mod} to the compiler:
!!
!!\begin{boxred}{Static Compiling Phase}
!!\begin{verbatim}
!!ifort -c lib_vtk_io.mod test.f90 -o test.o
!!\end{verbatim}
!!\end{boxred}
!!
!!In this example \emph{ifort} is the Intel Fortran Compiler\footnote{Da aggiungere.} and the \verb|-c| flag compiles preventing
!! linking; the compiler must \virgo{see} the module interface: the file \emph{lib\_vtk\_io.mod} must be placed in a folder visible
!!by the compiler.
!!
!!In the linking phase one simply give the library to the compiler:
!!
!!\begin{boxred}{Static Linking Phase}
!!\begin{verbatim}
!!ifort test.o LIB_VTK_IO.a -o test.out
!!\end{verbatim}
!!\end{boxred}
!!
!!The library must be placed in a folder visible by the compiler.
!!
!!\subsection{Dynamic Library}
!!\label{sec:Linux Dynamic}
!!The dynamic version of the precompiled library must be installed. The operating system must know where is the library so it is
!!necessary to install the library in a folder where the OS search its shared objects. In the most of the GNU/Linux distro the
!!folder \emph{/usr/lib/} is scanned to find shared objects. After you have copied the \MaiuscolettoBS{LIB\_VTK\_IO}.so file in
!!this folder, update the list of the shared objects with the command \verb|ldconfig -v| and the OS is ready to use the library.
!!
!!After you set your OS the compiling and linking phase is identical to the previous (remember to you the module interface at
!!the compiling phase). The only difference is to use the dynamic library at the linking phase:
!!
!!\begin{boxred}{Dynamic Linking Phase}
!!\begin{verbatim}
!!ifort test.o LIB_VTK_IO.so -o test.out
!!\end{verbatim}
!!\end{boxred}
!!
!!\section{MS Windows}
!!
!!Unluckily for MS Windows there is not any support at the moment. As soon as I have the time, I make some instructions on how
!!use \LIBVTKIO with MS Visual Studio and Intel Visual Fortran for MS Windows.
!!
!!\chapter{LIB\_VTK\_IO Programming Style}
!!\label{cap:Programming Style}
!!
!!\lettrine[lines=2,loversize=-.1,lraise=0.2]{{\bf A}}{ll} the \LIBVTKIO functions are \MaiuscolettoBS{4-byte integer functions}:
!!the output of these functions is an integer that is $0$ if the function calling has been done right while it is $> 0$  if some
!!errors occur (the error handling is only at its embryonal phase). Therefore the functions calling must be done in the following
!!way:
!!
!!\begin{boxred}{functions Calling}
!!\begin{verbatim}
!!...
!!integer(4):: E_IO
!!...
!!E_IO = VTK_INI(....
!!...
!!\end{verbatim}
!!\end{boxred}
!!
!!The \LIBVTKIO programming style is based on two main principles: \MaiuscolettoBS{portable kind-precision} of reals and integers
!!variables and \MaiuscolettoBS{dynamic dispatching}. In the appendix \ref{cap:kind precision} and \ref{cap:Dynamic Dispatching}
!!there are more details about these choices. I just remark some consequences of these choices. Using \MaiuscolettoBS{dynamic
!!dispatching} the \LIBVTKIO has a simple API. The user calls a generic procedure (VTK\_INI, VTK\_GEO,...) and the library,
!!depending on the type of the inputs passed, calls the correct internal function (i.e. VTK\_GEO for 8-byte real type if the input
!!passed is 8-byte real type). By this interface only few functions are used whitout the necessity of calling a different function
!!for every different inputs type. \MaiuscolettoBS{Dynamic dispatching} is valid also for the different kind of topology and
!!variables-data-dimensions; the function VTK\_GEO is the same for all topologies, just the inputs passed to the functions change
!!as the topology changes. Also the dimensions of variables-data use the \MaiuscolettoBS{dynamic dispatching}: the function
!!(VTK\_VAR) used to save vectorial data is identical to the one used for scalar data, depending on the dimensions of the data
!!\LIBVTKIO calls the correct internal function. \MaiuscolettoBS{Dynamic dispatching} is based on the internal kind-precision
!!selecting convention: Fortran 90/95 standard has some useful functions to achive the portability of reals and integers precision
!!and \LIBVTKIO uses these functions to define portable kind-precision; because it is important to make portable the code on
!!different architectures I suggest to use this programming style.
!!
!!The data handled by \LIBVTKIO can be classified into two main categories:
!!
!!\begin{enumerate1Red}
!!\item Geometric Data. These are the geometric informations of the mesh and they can be of different kind and different number
!!      depending on the topology choiced. The mesh points coordinates type must be of 4-byte real type or 8-byte real type.
!!\item Variable Data. These are the scalar or vectorial variables appended to the mesh points (both at the cell-nodes and the
!!      cell-centers of the mesh). The type of these data can be of 8-byte real type, 4-byte real type and 4-byte integer type
!!      (for the XML output there are also the 8-byte integer type, 2-byte integer type and 1-byte integer type).
!!\end{enumerate1Red}
!!
!!In the following chapters theare the details of \LIBVTKIO API.
!!
!!\part{LIB\_VTK\_IO API}
!!\label{part:LIBVTKIO API}
module LIB_VTK_IO
!----------------------------------------------------------------------------------------------------------------------------------
!!\LIBVTKIO is a library of functions for Input and Output pure Fortran data (both ascii and binary) in VTK format.
!!
!!The VTK standard can be separated into two main catagories: the \MaiuscolettoBS{VTK Legacy Standard} and the
!!\MaiuscolettoBS{VTK XML Standard}. The latter is more powerful and will has a stronger support from VTk comunity than legacy
!!standard; XML file format would to be preferred despite the legacy one.
!!
!!At the present only a few functions of the final library have been implemented. The InPut functions are totaly absent, but the
!!OutPut functions are almost complete (the \virgo{polydata} functions are the only missing).
!!
!!The functions actually present are:
!!
!!\begin{boxred}{functions for Legacy VTK file format}
!!\begin{enumerate1Red}
!! \item \MaiuscolettoS{VTK\_INI}
!! \item \MaiuscolettoS{VTK\_GEO}
!! \item \MaiuscolettoS{VTK\_CON}
!! \item \MaiuscolettoS{VTK\_DAT}
!! \item \MaiuscolettoS{VTK\_VAR}
!! \item \MaiuscolettoS{VTK\_END}
!!\end{enumerate1Red}
!!\end{boxred}
!!
!!\begin{boxred}{functions for XML VTK file format}
!!\begin{enumerate1Red}
!! \item \MaiuscolettoS{VTK\_INI\_XML}
!! \item \MaiuscolettoS{VTK\_GEO\_XML}
!! \item \MaiuscolettoS{VTK\_CON\_XML}
!! \item \MaiuscolettoS{VTK\_DAT\_XML}
!! \item \MaiuscolettoS{VTK\_VAR\_XML}
!! \item \MaiuscolettoS{VTK\_END\_XML}
!!\end{enumerate1Red}
!!\end{boxred}
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
! functions for VTK LEGACY
public:: VTK_INI
public:: VTK_GEO
public:: VTK_CON
public:: VTK_DAT
public:: VTK_VAR
public:: VTK_END
! functions for VTK XML
public:: VTK_INI_XML
public:: VTK_GEO_XML
public:: VTK_CON_XML
public:: VTK_DAT_XML
public:: VTK_VAR_XML
public:: VTK_END_XML
! portable kind-precision
public:: R16P, FR16P
public:: R8P,  FR8P
public:: R4P,  FR4P
public:: R_P,  FR_P
public:: I8P,  FI8P
public:: I4P,  FI4P
public:: I2P,  FI2P
public:: I1P,  FI1P
public:: I_P,  FI_P
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! overloading of VTK_GEO
interface VTK_GEO
  module procedure VTK_GEO_UNST_R8, & ! real(R8P) UNSTRUCTURED_GRID
                   VTK_GEO_UNST_R4, & ! real(R4P) UNSTRUCTURED_GRID
                   VTK_GEO_STRP_R8, & ! real(R8P) STRUCTURED_POINTS
                   VTK_GEO_STRP_R4, & ! real(R4P) STRUCTURED_POINTS
                   VTK_GEO_STRG_R8, & ! real(R8P) STRUCTURED_GRID
                   VTK_GEO_STRG_R4, & ! real(R4P) STRUCTURED_GRID
                   VTK_GEO_RECT_R8, & ! real(R8P) RECTILINEAR_GRID
                   VTK_GEO_RECT_R4    ! real(R4P) RECTILINEAR_GRID
endinterface
! overloading of VTK_VAR
interface VTK_VAR
  module procedure VTK_VAR_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_TEXT_R8, & ! real(R8P)    vectorial (texture)
                   VTK_VAR_TEXT_R4    ! real(R4P)    vectorial (texture)
endinterface
! overloading of VTK_GEO_XML
interface VTK_GEO_XML
  module procedure VTK_GEO_XML_STRG_R4, & ! real(R4P) StructuredGrid
                   VTK_GEO_XML_STRG_R8, & ! real(R8P) StructuredGrid
                   VTK_GEO_XML_RECT_R8, & ! real(R8P) RectilinearGrid
                   VTK_GEO_XML_RECT_R4, & ! real(R4P) RectilinearGrid
                   VTK_GEO_XML_UNST_R8, & ! real(R8P) UnstructuredGrid
                   VTK_GEO_XML_UNST_R4, & ! real(R4P) UnstructuredGrid
                   VTK_GEO_XML_CLOSEP     ! closing tag "Piece" function
endinterface
! overloading of VTK_VAR_XML
interface VTK_VAR_XML
  module procedure VTK_VAR_XML_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_XML_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_XML_SCAL_I8, & ! integer(I8P) scalar
                   VTK_VAR_XML_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_XML_SCAL_I2, & ! integer(I2P) scalar
                   VTK_VAR_XML_SCAL_I1, & ! integer(I1P) scalar
                   VTK_VAR_XML_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_XML_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_XML_VECT_I8, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I2, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I1    ! integer(I4P) vectorial
endinterface
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!!\LIBVTKIO has a small set of internal variables and parameters some of which have public visibility.
!!
!!The \LIBVTKIO uses a partable kind parameters for real and integer variables. The following are the kind parameters used: these
!!parameters are public and their use is strong encouraged.
!!
!!Real precision definitions:
!!
integer, parameter:: R16P = selected_real_kind(33,4931) ! 33  digits, range $[\pm 10^{-4931}  ,\pm 10^{+4931}   -1]$
integer, parameter:: R8P  = selected_real_kind(15,307)  ! 15  digits, range $[\pm 10^{-307}~~ ,\pm 10^{+307}~~  -1]$
integer, parameter:: R4P  = selected_real_kind(6,37)    ! 6~~~digits, range $[\pm 10^{-37}~~~~,\pm 10^{+37}~~~~ -1]$
integer, parameter:: R_P  = R8P                         ! default real precision
!!Integer precision definitions:
!!
integer, parameter:: I8P  = selected_int_kind(18)       ! range $[-2^{63} ,+2^{63}  -1]$
integer, parameter:: I4P  = selected_int_kind(9)        ! range $[-2^{31} ,+2^{31}  -1]$
integer, parameter:: I2P  = selected_int_kind(4)        ! range $[-2^{15} ,+2^{15}  -1]$
integer, parameter:: I1P  = selected_int_kind(2)        ! range $[-2^{7}~~,+2^{7}~~ -1]$
integer, parameter:: I_P  = I4P                         ! default integer precision
!!
!!Besides the kind parameters there are also the format parameters useful for writing in a well-ascii-format numeric variables.
!!Also these parameters are public.
!!
!! Real output formats:
!!
character(10), parameter:: FR16P = '(E41.33E4)'         ! R16P  output format
character(10), parameter:: FR8P  = '(E23.15E3)'         ! R8P   output format
character(9),  parameter:: FR4P  = '(E14.6E2)'          ! R4P   output format
character(10), parameter:: FR_P  = '(E23.15E3)'         ! R\_P  output format
!! Integer output formats:
!!
character(5), parameter:: FI8P  = '(I21)'               ! I8P  output format
character(5), parameter:: FI4P  = '(I12)'               ! I4P  output format
character(4), parameter:: FI2P  = '(I7)'                ! I2P  output format
character(4), parameter:: FI1P  = '(I5)'                ! I1P  output format
character(5), parameter:: FI_P  = '(I12)'               ! I\_P output format
!!
!!\LIBVTKIO uses a small set of internal variables that are private (not accessible from the outside). The following are
!! private variables:
!!
integer(I4P), parameter:: maxlen       = 500         ! max number of characters os static string
character(1), parameter:: end_rec      = char(10)    ! end-character for binary-record finalize
integer(I4P), parameter:: f_out_ascii  = 0           ! ascii-output-format parameter identifier
integer(I4P), parameter:: f_out_binary = 1           ! binary-output-format parameter identifier
integer(I4P)::            f_out        = f_out_ascii ! current output-format (initialized to ascii format)
character(len=maxlen)::   topology                   ! mesh topology
integer(I4P)::            Unit_VTK                   ! internal logical unit
integer(I4P)::            Unit_VTK_Append            ! internal logical unit for raw binary XML append file
integer(I4P)::            N_Byte                     ! number of byte to be written/read
real(R8P)::               tipo_R8                    ! prototype of R8P real
real(R4P)::               tipo_R4                    ! prototype of R4P real
integer(I8P)::            tipo_I8                    ! prototype of I8P integer
integer(I4P)::            tipo_I4                    ! prototype of I4P integer
integer(I2P)::            tipo_I2                    ! prototype of I2P integer
integer(I1P)::            tipo_I1                    ! prototype of I1P integer
integer(I4P)::            ioffset                    ! offset pointer
integer(I4P)::            indent                     ! indent pointer
!----------------------------------------------------------------------------------------------------------------------------------

!!In the following chapters there is the API reference of all functions of \LIBVTKIO.
contains
  !!\chapter{Auxiliary functions}
  !!\minitoc
  !!\vspace*{8mm}
  !!
  !!\LIBVTKIO uses two auxiliary functions that are not connected with the VTK standard. These functions are private and so they
  !!cannot be called outside the library.
  function GetUnit() result(Free_Unit)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The GetUnit function is used for getting a free logic unit. The users of \LIBVTKIO does not know which is
  !!the logical unit: \LIBVTKIO handels this information without boring the users. The logical unit used is safe-free: if the
  !!program calling \LIBVTKIO has others logical units used \LIBVTKIO will never use these units, but will choice one that is free.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  integer(I4P):: Free_Unit ! free logic unit
  integer(I4P):: n1        ! counter
  integer(I4P):: ios       ! inquiring flag
  logical(4)::   lopen     ! inquiring flag
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of GetUnit function: the units 0, 5, 6, 9 and all non-free units are discarded.
  !!
  !(\doc)codesnippet
  Free_Unit = -1_I4P                                      ! initializing free logic unit
  n1=1_I4P                                                ! initializing counter
  do
    if ((n1 /= 5_I4P) .and. (n1 /= 6_I4P) .and. (n1 /= 9_I4P)) then
      inquire (unit=n1,opened=lopen,iostat=ios)           ! verify logic units
      if (ios == 0_I4P) then
        if (.not. lopen) then
          Free_Unit = n1                                  ! assignment of free logic
          return
        endif
      endif
    endif
    n1=n1+1_I4P                                           ! updating counter
  enddo
  return
  !(doc/)codesnippet
  !!GetUnit function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope to public.
  !--------------------------------------------------------------------------------------------------------------------------------
  end function GetUnit

  function Upper_Case(string)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The Upper\_Case function converts the lower case characters of a string to upper case one. \LIBVTKIO uses this function in
  !!order to achieve case-insensitive: all character variables used within \LIBVTKIO functions are pre-processed by
  !!Uppper\_Case function before these variables are used. So the users can call \LIBVTKIO functions whitout pay attention of the
  !!case of the kwywords passed to the functions: calling the function VTK\_INI with the string \code{E_IO = VTK_INI('Ascii',...)}
  !!or with the string  \code{E_IO = VTK_INI('AscII',...)} is equivalent.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(len=*), intent(IN):: string     ! string to be converted
  character(len=len(string))::   Upper_Case ! converted string
  integer::                      n1         ! characters counter
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  !!The following is the code snippet of Upper\_Case function.
  !!
  !(\doc)codesnippet
  Upper_Case = string
  do n1=1,len(string)
    select case(ichar(string(n1:n1)))
    case(97:122)
      Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! Upper case conversion
    endselect
  enddo
  return
  !(doc/)codesnippet
  !!Upper\_Case function is private and cannot be called outside \LIBVTKIO. If you are interested to use it change its scope
  !!to public.
  !--------------------------------------------------------------------------------------------------------------------------------
  end function Upper_Case

  !!\chapter{VTK LEGACY functions}
  !!\minitoc
  !!\vspace*{8mm}
  !!
  function VTK_INI(output_format,filename,title,mesh_topology) result(E_IO)
  !--------------------------------------------------------------------------------------------------------------------------------
  !!The VTK\_INI function is used for initializing file. This function must be the first to be called.
  !--------------------------------------------------------------------------------------------------------------------------------

  implicit none

  !--------------------------------------------------------------------------------------------------------------------------------
  character(*), intent(IN):: output_format ! output format: ASCII or BINARY
  character(*), intent(IN):: filename      ! name of file
  character(*), intent(IN):: title         ! title
  character(*), intent(IN):: mesh_topology ! mesh topology
  integer(I4P)::             E_IO          ! Input/Output inquiring flag: $0$ if IO is done, $ > 0$ if IO is not done
  !!The VTK\_INI variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}output\_format}] indicates the \virgo{format} of output file. It can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{ascii} (it is case insensitive) $\rightarrow$ creating an ascii output file.
  !!  \item \emph{binary} (it is case insensitive) $\rightarrow$ creating a binary (big\_endian encoding) output file.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}filename}] contains the name (with its path) of the output file.
  !! \item[{\color{RoyalBlue}title}] contains the title of the VTK dataset.
  !! \item[{\color{RoyalBlue}topology}] indicates the topology of the mesh and can assume the following values:
  !! \begin{enumerateABlu}
  !!  \item \emph{STRUCTURED\_POINTS}.
  !!  \item \emph{STRUCTURED\_GRID}.
  !!  \item \emph{UNSTRUCTURED\_GRID}.
  !!  \item \emph{RECTILINEAR\_GRID}.
  !! \end{enumerateABlu}
  !! \item[{\color{RoyalBlue}E\_IO}] contains the inquiring integer flag for error handling.
  !!\end{description}
  !!
  !!The following is an example of VTK\_INI calling:
  !!
  !!\begin{boxred}{VTK\_INI Calling}
  !!\begin{verbatim}
  !!...
  !!E_IO = VTK_INI('Binary','example.vtk','VTK legacy file','UNSTRUCTURED_GRID')
  !!...
  !!\end{verbatim}
  !!\end{boxred}
  !!\noindent Note that the \virgo{.vtk} extension is necessary in the file name.
  !--------------------------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------------------------
  topology = trim(mesh_topology)
  Unit_VTK=GetUnit()
  select case(trim(Upper_Case(output_format)))
  case('ASCII')
    f_out = f_out_ascii
    open(unit     = Unit_VTK, &
         file     = trim(filename), &
         form     = 'FORMATTED', &
         access   = 'SEQUENTIAL', &
         action   = 'WRITE', &
         buffered = 'YES', &
         iostat   = E_IO)
    ! writing header of file
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'# vtk DataFile Version 3.0'
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(title)
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)trim(Upper_Case(output_format))
    write(unit=Unit_VTK,fmt='(A)',iostat=E_IO)'DATASET '//trim(topology)
  case('BINARY')
    f_out = f_out_binary
    open(unit       = Unit_VTK, &
         file       = trim(filename), &
         form       = 'UNFORMATTED', &
         access     = 'SEQUENTIAL', &
         action     = 'WRITE', &
         convert    = 'BIG_ENDIAN', &
         recordtype = 'STREAM', &
         buffered   = 'YES', &
         iostat     = E_IO)
    ! writing header of file
    write(unit=Unit_VTK,iostat=E_IO)'# vtk DataFile Version 3.0'//end_rec
    write(unit=Unit_VTK,iostat=E_IO)trim(title)//end_rec
    write(unit=Unit_VTK,iostat=E_IO)trim(Upper_Case(output_format))//end_rec
    write(unit=Unit_VTK,iostat=E_IO)'DATASET '//trim(topology)//end_rec
  endselect
  return
  !--------------------------------------------------------------------------------------------------------------------------------
  end function VTK_INI

  !!\section{VTK\_GEO}
  !!
  !!VTK\_GEO is an interface to 8 different functions; there are 2 functions for each 4 different topologies actually supported:
  !!one function for mesh coordinates with R8P precision and one for mesh coordinates with R4P precision.
  !!This function must be called after VTK\_INI. It saves the mesh geometry. The inputs that must be passed change depending on
  !!the topologies choiced. Not all VTK topologies have been implemented (\virgo{polydata} topologies are absent). The signatures
  !!for all implemented topologies are now reported.
  !!
  !!\subsection{VTK\_GEO STRUCTURED POINTS}
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=signature,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Points Signature}]
  !! function VTK_GEO(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!The topology \virgo{structured points} is useful for structured grid with uniform discretization steps.
  !!
  !!\begin{boxred}{}
  !!\begin{lstlisting}[style=variables,title=\color{Maroon}\MaiuscolettoBS{VTK\_GEO Structured Points Variables}]
  !!integer(I4P),     intent(IN):: Nx   ! number of nodes in x direction
  !!integer(I4P),     intent(IN):: Ny   ! number of nodes in y direction
  !!integer(I4P),     intent(IN):: Nz   ! number of nodes in z direction
  !!real(R8P or R4P), intent(IN):: X0   ! x coordinate of origin
  !!real(R8P or R4P), intent(IN):: Y0   ! y coordinate of origin
  !!real(R8P or R4P), intent(IN):: Z0   ! z coordinate of origin
  !!real(R8P or R4P), intent(IN):: Dx   ! space step in x
  !!real(R8P or R4P), intent(IN):: Dy   ! space step in y
  !!real(R8P or R4P), intent(IN):: Dz   ! space step in z
  !!integer(I4P)::                 E_IO ! Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done
  !!\end{lstlisting}
  !!\end{boxred}
  !!
  !!Note that the variables \texttt{X0,Y0,Z0,Dx,Dy,Dz} can be passed both as 8-byte real kind and 4-byte real kind; the dynamic
  !!displacement interface will call the correct function. Mixing 8-byte real kind and 4-byte real kind is not allowed: be sure
  !!that all variables are 8-byte real kind or all are 4-byte real kind.
  !!
  !!The VTK\_GEO structured point variables have the following meaning:
  !!
  !!\begin{description}
  !! \item[{\color{RoyalBlue}Nx}] indicates the number of nodes in $X$ direction.
  !! \i
