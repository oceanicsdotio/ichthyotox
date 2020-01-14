January, 2011
Another thing I plan to do is about boundary bouncing, which is related to the
treatment of INDOMAIN,INWATER,IFOUND,SBOUND,etc. which is totally treated in
`offlag.f90` and `triangle_grid_edge.f90`

1) In mod_var.f90, add variable LAG%np_stop (total number of particles going out of domain from both
OB & SB), while LAG%np_out represents total number of particles going out of domain from OB.
The reason from this distinction is that particles going out of domain from SB (solid boundary)
can be retained in the domain by bounce algorithm.


Keeney changes for cyanotoxin addition (September, 2015)
Search "!--Keeney change"

"offlag.f90" contains...
  program "particle_traj"
    changes:
      file used is "./icthyotox.nc"
  subroutine "set_lag"
    changes:
      add temperature field
      check for existence of state variable file for toxins
  subroutine "lag_update"
  
    
    
1.3) in LAG_UPDATE
    declare variables for start/end of hourly temperature
    change NetCDF file name
    time interpolation of temperature field
    depreciate kinesis, replace with conditional  *_movement subroutines
    update hourly temperature field
    

"mod_var.f90" contains...
  module "mod_lag"
    changes: 
      declaration of LAG_OBJ allocatable
    to do: none
  module "lims"
    changes: none
    to do: none
  module "all_vars"
    changes: none
    to do: none


"data_run.f90" contains...
  subroutine "data_run"
    changes: none
    to do:
      random walk update
    
"alloc_var.f90" contains...
  subroutine "alloc_var"
    changes:
      allocate current and previous temperature arrays
    to do: none
    
"triangle_grid_edges.f90" contains...
  subroutine "triangle_grid_edges"
    changes: none
    to do: none


Huang changes for fish (January, 2011)
look for "!---fish change 1" etc.

1) In ncdio.f90
1.1) specify kb=2, that is using a 2-D field (surface u, v, and S1) (October 2009)
1.2) change "SUBROUTINE NCD_READ(INFILE,UL,VL,WWL,KHL,ELL,S1L,HO)" ouput variables

2) In offlag.f90
2.1) in subroutine "SUBROUTINE SET_LAG", change 6 places
2.2) add subroutine "CALL INTERP_SAL(NDRFT,LAG%HOST,LAG%INDOMAIN,LAG%SBOUND, &
       LAG%XP,LAG%YP,S1NC,LAG%SAL,0)" 
     Note this is only for 2-D case.
     For 3-D case, this subroutine needs to be changed
2.3) in subroutine "SUBROUTINE LAG_UPDATE", change 12 places
     Most important is "CALL TRAJECT(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,LAG%ZPT,LAG%HOST,&
             LAG%INDOMAIN,LAG%SBOUND,UT,U,VT,V,WT,W,H,ET,EL,ST1,S1,LAG%SAL)"
2.4) in "SUBROUTINE TRAJECT(NDRFT,DTI,LAG%XP,LAG%YP,LAG%ZP,LAG%ZPT,LAG%HOST,&
             LAG%INDOMAIN,LAG%SBOUND,UT,U,VT,V,WT,W,H,ET,EL,ST1,S1,LAG%SAL)", change 3 places
2.5) add subroutine "CALL KINESIS"


Haosheng Huang, October 2009:

I did several changes to Martin's code (search "!--Huang change"):

1) Input and output file format
   For every new run, need to change input NETCDF file name in offlag.f90 (4 places)
2) find particle's cell is time consuming in certain cases, I did some change
   (compare offlag.f90 with offlag.f90_1)
3) In triangle_grid_edge.f90 the last few lines, it must be changed for each new case
   (specifying open boundary node)
4) check mod_var.f90, alloc_vars.f90 and ncdio.f90 to see how many variables (T, S etc.) are included

Important!!!
5) Vertical motion of particles is wrong.
Reason: NETCDF input is "ww", which is verticle volocity in z-coordinate, 
        while all calculations use omega (vertical velocity in sigma-coordinate)
        You can compare the dz formula with online version.


Totally there are 12 files for the offline Lagrangian particle tracking code

1) util.f90 is the same for the following three folders:
~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill, ~/FVCOM_source/OFFLINE_LAG_source_Breton, and
~/FVCOM_source/OFFLINE_LAG_source_original

2)  data_run.f90 is the same for 
~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill and ~/FVCOM_source/OFFLINE_LAG_source_Breton.
The difference between the above two folders and ~/FVCOM_source/OFFLINE_LAG_source_original 
is that I change variables "INSTP" and "DTOUT" from integer to real.

3) mod_prec.f90 is the same for the three folders
4) mod_ncd.f90 is the same for the three folders
5) mod_inp.f90 is the same for the three folders

6) mod_var.f90
Between ~/FVCOM_source/OFFLINE_LAG_source_original and ~/FVCOM_source/OFFLINE_LAG_source_Breton
Difference 1 is related to 2) data_run.f90 variables change from integer to real.
Difference 2 is case specific. Since in Breton Sound case no T&S, comment them out to save memory

Between ~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill and ~/FVCOM_source/OFFLINE_LAG_source_Breton
Difference 3: Allow array for 3-D salinity

7) alloc_vars.f90
Between ~/FVCOM_source/OFFLINE_LAG_source_original and ~/FVCOM_source/OFFLINE_LAG_source_Breton
The difference is related to mod_var.f90 Difference 2

Between ~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill and ~/FVCOM_source/OFFLINE_LAG_source_Breton
The difference is related to mod_var.f90 Difference 3

8) makedepends is the same for the three folders
9) makefile is similar for the three folders

10) ncdio.f90
Between ~/FVCOM_source/OFFLINE_LAG_source_original and ~/FVCOM_source/OFFLINE_LAG_source_Breton
Difference 1: find one error TEMP(N,1) should be TEMP(M,1) (for H -- depth)
Difference 2: change from "HT=HO+1" to "HT=HO"
              This is related to the changes I made in offlag.f90
Difference 3: This is related to the fact that Breton case is a 2-D calculation.
              (from "!--Huang change 02" to "!--Huang change 02 end")

Between ~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill and ~/FVCOM_source/OFFLINE_LAG_source_Breton
Keep Difference 1 and Difference 2
In Difference 3,in LONI_oil_spill case I need to read in 3-D u and v, as well as 3-D salinity.
However, I didn't output w in the NETCDF file. So I made some special treatments.

11) triangle_grid_edge.f90
Between ~/FVCOM_source/OFFLINE_LAG_source_original and ~/FVCOM_source/OFFLINE_LAG_source_Breton
Difference: In the last few lines, specifying open boundary node 
            (related to stop calculation near open boundary). It must be manually changed for each new case.

Between ~/FVCOM_source/OFFLINE_LAG_source_LONI_oil_spill and ~/FVCOM_source/OFFLINE_LAG_source_Breton
related to the above Difference

12) offlag.f90
Main changes are in this subroutine.

