[![Ready in backlog](https://badge.waffle.io/szaghi/xnplot.png?label=ready&title=Ready)](https://waffle.io/szaghi/xnplot)
[![In Progress](https://badge.waffle.io/szaghi/xnplot.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/xnplot)
[![Open bugs](https://badge.waffle.io/szaghi/xnplot.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/xnplot)

# XnPlot
### <a name="top">Xnavis PostProcessor
XnPlot loads Xnavis mesh and solution files and produces post-processed plot files.

## <a name="toc">Table of Contents

* [Team Members](#team-members)
* [What is XnPlot?](#what)
* [Todos](#todos)
* [Requirements](#requirements)
* [Copyrights](#copyrights)
* [Download XnPlot](#download)
* [Compiling Instructions](#compile)
* [Usage](#usage)
  + [Main Help](#usage-help)
  + [Post-processing only mesh files](#usage-only-mesh)
  + [Post-processing mesh and solutions files](#usage-mesh-sol)
  + [Utilities](#utilities)
* [Version History](#versions)

## <a name="team-members"></a>Team Members
* Stefano Zaghi    <stefano.zaghi@gmail.com>
* Riccardo Broglia

Go to [Top](#top) or [Toc](#toc)
## <a name="what"></a>What is XnPlot?
XnPlot loads Xnavis mesh and solution files and produces post-processed plot files. It is one the many post-processors of Xnavis code.

Go to [Top](#top) or [Toc](#toc)
## <a name="todos"></a>Todos
+ Refactor utilities;
+ any feature request is welcome!

Go to [Top](#top) or [Toc](#toc)
## <a name="download"></a>Download XnPlot
If you use `git` it could be convenient to clone this repository:
```bash
git clone https://github.com/szaghi/XnPlot
```
Other 2 possibilities are:

1. use the GitHub **Download ZIP** button on the right sidebar of this page;
2. download one of the releases in the [release page](https://github.com/szaghi/XnPlot/releases), also listed at the end of this page.

Go to [Top](#top) or [Toc](#toc) 
## <a name="requirements"></a>Requirements
+ Modern Fortran Compiler (standard 2003+);
+ a lot of patience with the author.

XnPlot is developed on a GNU/Linux architecture. For Windows architecture there is no support, however it should be work out-of-the-box.

Go to [Top](#top) or [Toc](#toc)
## <a name="Copyrights"></a>Copyrights
XnPlot is an open source project, it is distributed under the [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html). Anyone is interest to use, to develop or to contribute to XnPlot is welcome.

Go to [Top](#top) or [Toc](#toc)
## <a name="compile"></a>Compiling Instructions
XnPlot has been developed on GNU/Linux architectures. Other OS are not supported (and in general there is no best alternative to GNU/Linux :-).

XnPlot have been successfully compiled with the following compilers:

+ GNU gfortran (version 4.7.0 or higher);
+ Intel Fortran Compiler ifort (version 12.0 or higher)


XnPlot is constituted by several modules. Therefore there are many dependences. The most easy way to compile the code is to start with the provided makefile thus it is necessary that the system has "Make" program (preferably GNU make http://www.gnu.org/software/make/).

The provided makefile has several options. It has one rule that prints all options available and the default settings. Typing in the shell prompt: `code make help` the following output will be printed:

```bash
 Make options of XnPlot code

 Compiler choice: COMPILER=intel => default
  COMPILER=gnu   => GNU gfortran           
  COMPILER=intel => Intel Fortran         

 Compiling options
  DEBUG=yes(no)    => on(off) debug                  (default no)
  F03STD=yes(no)   => on(off) check standard fortran (default no)
  OPTIMIZE=yes(no) => on(off) optimization           (default no)
  OPENMP=yes(no)   => on(off) OpenMP directives      (default no)
  BIGEIN=yes(no)   => on(off) Big Endian input files (default no)

 Preprocessing options
  R16P=yes(no) => on(off) definition of real with "128 bit" (default no)

 Executable directory
  DEXE="your_path" => directory where exe is placed (default ~/bin/)

 External libraries
  TECIO=yes(no) => on(off) Tecplot IO library linking (default )

 Provided Rules
  Defualt rule => ~/bin/XnPlot
  help         => printing this help message
  cleanobj     => cleaning compiled object
  cleanmod     => cleaning .mod files
  cleanmsg     => cleaning make-log massage files
  cleanexe     => cleaning executable files
  clean        => running cleanobj, cleanmod and cleanmsg
  cleanall     => running clean and cleanexe
  tar          => creating a tar archive of the project
```
For example compiling in debug mode with the Intel Fortran compiler you can type:
```bash
make DEBUG=yes COMPILER=intel
```

Go to [Top](#top) or [Toc](#toc)
## <a name="usage"></a>Usage
### <a name="usage-help">Main Help
XnPlot is is a Command Line Tool. To list the available options run it as following:
```bash
./XnPlot -h
```
this help will echoed
```bash
+--> XnPlot, post-processor for Xnavis                                                                                                                                                                                                                                
+--> Parsing Command Line Arguments                                                                                                                                                                                                                                   
|--> The XnPlot Command Line Interface (CLI) has the following options                                                                                                                                                                                                
|-->   XnPlot  -g value [-o value] [-i value] [-s value] [-ngc] [-cell] [-ls] [-nt] [-eq value] [-vordet] [-ascii] [-tec value] [-vtk value] [-proc value] [-os value] [-vb] [--help] [--version]                                                                     
|--> Each Command Line Argument (CLA) has the following meaning:                                                                                                                                                                                                      
|-->   [-g value]                                                                                                                                                                                                                                                     
|-->     Grid file (.grd)                                                                                                                                                                                                                                             
|-->     It is a non optional CLA thus must be passed to CLI                                                                                                                                                                                                          
|-->   [-o value]                                                                                                                                                                                                                                                     
|-->     output file name; default is basename of grd file with the proper extension                                                                                                                                                                                  
|-->     It is a optional CLA which default value is "unset"                                                                                                                                                                                                          
|-->   [-i value]                                                                                                                                                                                                                                                     
|-->     ICC file                                                                                                                                                                                                                                                     
|-->     It is a optional CLA which default value is "unset"                                                                                                                                                                                                          
|-->   [-s value]                                                                                                                                                                                                                                                     
|-->     solution file name; if passed the solution variables are saved                                                                                                                                                                                               
|-->     It is a optional CLA which default value is "unset"                                                                                                                                                                                                          
|-->   [-ngc]                                                                                                                                                                                                                                                         
|-->     mesh without ghosts cells, as geogrd output                                                                                                                                                                                                                  
|-->     It is a optional CLA which default value is ".false."                                                                                                                                                                                                        
|-->   [-cell]                                                                                                                                                                                                                                                        
|-->     variables other than nodes coord. are cell centered                                                                                                                                                                                                          
|-->     It is a optional CLA which default value is ".false."
|-->   [-ls]
|-->     solution with level set
|-->     It is a optional CLA which default value is ".false."
|-->   [-nt]
|-->     no turbulent model used
|-->     It is a optional CLA which default value is ".false."
|-->   [-eq value] with value chosen in: (0,1,2)
|-->     equations turbulent model
|-->     It is a optional CLA which default value is "1"
|-->   [-vordet]
|-->     computing variables for vortices identification
|-->     It is a optional CLA which default value is ".false."
|-->   [-ascii]
|-->     write ascii output files
|-->     It is a optional CLA which default value is ".false."
|-->   [-tec value] with value chosen in: (yes,no)
|-->     write output Tecplot files
|-->     It is a optional CLA which default value is "yes"
|-->   [-vtk value] with value chosen in: (yes,no)
|-->     write output VTK files
|-->     It is a optional CLA which default value is "no"
|-->   [-proc value]
|-->     process number for global block numeration if proc.input is found
|-->     It is a optional CLA which default value is "-1"
|-->   [-os value] with value chosen in: (UIX,WIN)
|-->     type of Operating System
|-->     It is a optional CLA which default value is "UIX"
|-->   [-vb]
|-->     Verbose output
|-->     It is a optional CLA which default value is ".false."
|-->   [--help] or [-h]
|-->     Print this help message
|-->     It is a optional CLA
|-->   [--version] or [-v]
|-->     Print version
|-->     It is a optional CLA
|--> Usage examples:
|-->   -) XnPlot -g xship.grd -o grid
|-->   -) XnPlot -g cc.01.grd -i cc.01 -o mesh.01
|-->   -) XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01
```
### <a name="usage-only-mesh">Post-processing only mesh files
#### GRD files (no ghost cells)
```bash
  XnPlot -g xship.grd -o mesh
```
A file named `mesh.plt` is generated. 
#### Overset/Overott files (with ghost cells)
```bash
  XnPlot -g cc.01.grd -i cc.01 -o mesh.01
```
A file named `mesh.01.plt` is generated. 
### <a name="usage-mesh-sol">Post-processing mesh and solutions files
```bash
  XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01
```
A file named `sol.01.plt` is generated. 
### <a name="utilities">Utilities
To be written.

Go to [Top](#top) or [Toc](#toc)

## <a name="versions"></a>Version History
In the following the changelog of most important releases is reported.
### v0.0.2 
##### Download [ZIP](https://github.com/szaghi/XnPlot/archive/v0.0.2.zip) ball or [TAR](https://github.com/szaghi/XnPlot/archive/v0.0.2.tar.gz) one
XnPlot.f90 module split for improve compiling speed. Fully backward compatible.
### v0.0.1 
##### Download [ZIP](https://github.com/szaghi/XnPlot/archive/v0.0.1.zip) ball or [TAR](https://github.com/szaghi/XnPlot/archive/v0.0.1.tar.gz) one
Stable Release. Fully backward compatible.
