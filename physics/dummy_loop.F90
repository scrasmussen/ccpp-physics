module dummy_loop
  
  implicit none
  
  contains

!> \section arg_table_dummy_loop_run Argument Table
!! \htmlinclude dummy_loop_run.html
!!
  subroutine dummy_loop_run(loop_cnt, errmsg, errflg)
    
    integer,           intent(in   )              :: loop_cnt
        
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
  end subroutine dummy_loop_run
end module dummy_loop