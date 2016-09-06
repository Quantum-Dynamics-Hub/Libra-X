# Libra-GAMESS_interface

   This file introduces how to execute Libra-GAMESS_interface.

0. Install Libra and GAMESS on your PC or server.
   For installation, the websites below will be helpful:
    Libra:  http://www.acsu.buffalo.edu/~alexeyak/libra/index.html
   GAMESS:  http://www.msg.ameslab.gov/gamess/

1. Create a working directory,say, /home/work (outside this repository). 

2. There, create input files(*.inp).(H2O.inp and 23waters.inp in "run" directory are the simple examples.)
   For more details about how to create that, 
   please see the website http://www.msg.ameslab.gov/gamess/GAMESS_Manual/input.pdf .
   Here, Keep in mind 3 things.
   A. Only semi-empirical methods have been connected to libra so far;
      set GBASIS=MNDO, AM1, PM3, or RM1 in $BASIS section. 
   B. Set RUNTYP=GRADIENT in $CONTROL section.
   C. Use cartesian coordinates and C1 symmetry in $DATA section:

      C1
      C  6.000000 4.377921 -4.769170 -2.758971
      C  6.000000 3.858116 -4.331728 -3.995136
      C  6.000000 2.478331 -4.387937 -4.267327
                           .
                           .
                           .

3. For convinience, copy run_gms.py in "run" directory to the working place.

4. Modify copied run_gms.py for calculation.
   Concretely, set variables for GAMESS, Molecular Dynamics(MD), excited electron dynamics, and debugs.
   See the input manual in "run" directory to know more about the variables.

5. Create "res" and "sd_ham" directories under the working place, where the results will be output.

6. When the precedures above are finished, it is the time to execute this program.
   Here, 2 types of execution can be used.
   A. Only invoke "python run_gms.py" in the working place.
   B. Use queuing system. submit_templ_gms.lsf or submit_templ_gms.slurm in "run" directory are the simple examples for using this.
      Modify the files following your queuing system.   
   
7. After the calculation finished, the results will be set in "res" directory.