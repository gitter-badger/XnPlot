# XnPlot
### <a name="top">Xnavis PostProcessor

XnPlot loads Xnavis mesh and solution files and produces post-processed plot files.

## <a name="toc">Table of Contents

* [Team Members](#team-members)
* [What is FLAP?](#what)
* [Todos](#todos)
* [Requirements](#requirements)
* [Copyrights](#copyrights)
* [Compiling Instructions](#compile)
* [Usage](#usage)

## <a name="team-members"></a>Team Members
* Stefano Zaghi <stefano.zaghi@gmail.com>

Go to [Top](#top) or [Toc](#toc)
## <a name="what"></a>What is XnPlot?

XnPlot loads Xnavis mesh and solution files and produces post-processed plot files.

Go to [Top](#top) or [Toc](#toc)
## <a name="todos"></a>Todos
+ ...
+ any feature request is welcome!

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
  BIGEIN=yes(no)   => on(off) Big Endian input files (default yes)

 Preprocessing options
  R16P=yes(no) => on(off) definition of real with "128 bit" (default no)

 Executable directory
  DEXE="your_path" => directory where exe is placed (default ~/bin/)

 External libraries
  TECIO=yes(no) => on(off) Tecplot IO library linking (default )

 Provided Rules
  Defualt rule     => ~/bin/XnPatches
  help             => printing this help message
  ~/bin/XnPatches => building OFF code
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
## <a name="usage"></a>Usage
### Basic Usage
XnPlot is is a Command Line Tool. To list the available options run it without arguments:

```bash
./XnPlot
```
this help will echoed
```bash
XnPlot
 Post processing code for Xnavis code
 Usage:
   XnPlot -g file_grd
         [-o file_output
          -i file_icc
          -s file_solution
          -ngc
          -cell
          -ls
          -eq #turbulent_eq_model (0,1,2)                                                                                                                                                                                                                             
          -vordet                                                                                                                                                                                                                                                     
          -ascii                                                                                                                                                                                                                                                      
          -tec yes/no                                                                                                                                                                                                                                                 
          -vtk yes/no                                                                                                                                                                                                                                                 
          -proc #proc                                                                                                                                                                                                                                                 
          -os UIX/WIN]                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                      
 Optional arguments and default values:                                                                                                                                                                                                                               
  -o file_output   => output file name; default is basename of grd file with the proper extension                                                                                                                                                                     
  -i file_icc      => icc file name; if passed the icc variable is saved at cell centers                                                                                                                                                                              
  -s file_solution => solution file name; if passed the solution variables are saved                                                                                                                                                                                  
  -ngc             => mesh without ghosts cells, as geogrd output (default no, grd with ghosts cells)                                                                                                                                                                 
  -cell            => all variables other than nodes coord. are cell centered (default no, node centered)                                                                                                                                                             
  -ls              => solution with level set                                                                                                                                                                                                                         
  -nt              => no turbulent model, laminar flow                                                                                                                                                                                                                
  -eq #0/1/2       => # equations turbulent model (default 1)                                                                                                                                                                                                         
  -vordet          => computing "vordet" variable for vortices identification (default no)                                                                                                                                                                            
  -ascii           => write ascii output file (default no, write binary one)                                                                                                                                                                                          
  -tec yes/no      => write (or not) Tecplot file format (default yes)                                                                                                                                                                                                
  -vtk yes/no      => write (or not) VTK file format (default no)                                                                                                                                                                                                     
  -proc #proc      => if file "proc.input" if found global blocks numeration is used; #proc is the process                                                                                                                                                            
                      number of the current processed file                                                                                                                                                                                                            
  -os UIX/WIN      => type of Operating System write (default *UIX OS type)                                                                                                                                                                                           
                                                                                                                                                                                                                                                                      
 Examples:                                                                                                                                                                                                                                                            
  XnPlot -g xship.grd                       -o mesh.01 (process only grd file)                                                                                                                                                                                        
  XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01  (solution variables are saved)                                                                                                                                                                                 
                                                                                                                                                                                                                                                                      
 Notes:                                                                                                                                                                                                                                                               
   1) the output file name extension is not necessary because it assigned according to the type of output:                                                                                                                                                            
      binary       Tecplot => .plt                                                                                                                                                                                                                                    
      ascii        Tecplot => .dat                                                                                                                                                                                                                                    
      binary/ascii VTK     => .vtm                                                                                                                                                                                                                                    
   2) if a file name "mb.par" is present into the directory where XnPlot is executed the values of ls                                                                                                                                                                 
      turbulence model are loaded from this file thus they can be omitted from the command                                                                                                                                                                            
      line arguments list.                                                                                                                                                                                                                                            
   3) all the variables other than the nodes coordinates are saved at cell center if "-cell" option is used;                                                                                                                                                          
      if blanking is used the blanking mode must be "any corners" or "primary values".
```

### Post-processing only mesh files

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
### Post-processing mesh and solutuions files
```bash
  XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01
```
A file named `sol.01.plt` is generated. 

### Utilities
To be written.

Go to [Top](#top) or [Toc](#toc)
