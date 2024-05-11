Program integrating_marcus
use marcus  
implicit none

call Marcus_Integration
call Marcus_plot(rf,rb)
call replace(rf,rb)











end program
