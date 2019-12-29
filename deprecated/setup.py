# 
# 
#   
# 
#   # write runfile
#   open(unit=iorun, file="../"//trim(folderprefix)//"/ichthyotox_run.dat", status='replace')
#   write(iorun,"(A)") "INFOFILE = screen"
#   write(iorun,"(A,F4.2)") "DTI = ", 0.02_SP # inner interpolation time step, float"
#   write(iorun,"(A,F4.2)") "INSTP = ", 1.0_SP # time step of physical field data, float
#   write(iorun,"(A,F4.2)") "DTOUT = ", 0.1_SP # output time step, >dti
#   write(iorun,"(A,F4.2)") "DHOR = ", 0.10_SP # horizontal diffusion coefficient
#   write(iorun,"(A,F4.2)") "DTRW = ", 0.02_SP # random walk time step
#   write(iorun,"(A,I4)") "TDRIFT = ", 24*ndays  # number of time steps to iterate, int
#   write(iorun,"(A,I4)") "YEARLAG = ", 2016
#   write(iorun,"(A,I2)") "MONTHLAG = ", 4
#   write(iorun,"(A,I2)") "DAYLAG = ", 1
#   write(iorun,"(A,I2)") "HOURLAG = ", 0
#   write(iorun,"(A,I1)") "IRW = ", 0 # random walk type
#   write(iorun,"(A)") "P_SIGMA = F"
#   write(iorun,"(A)") "OUT_SIGMA = F"
#   write(iorun,"(A)") "F_DEPTH = F" 
#   write(iorun,"(A)") "GEOAREA = box"
#   write(iorun,"(A)") "INPDIR=/"
#   write(iorun,"(A)") "LAGINI=/"
#   write(iorun,"(A)") "OUTDIR=/"
#   
#  