module dummy_loop
  
  implicit none
  
  contains

!> \section arg_table_dummy_loop_init Argument Table
!! \htmlinclude dummy_loop_init.html
!!
  subroutine dummy_loop_init(loop_cnt, loop_max, ncol, errmsg, errflg)

    integer,           intent(in   )              :: loop_cnt
    integer,           intent(in   )              :: loop_max
    integer,           intent(in   )              :: ncol

    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

    call dummy_loop_internal('dummy_loop_init', loop_cnt, loop_max, ncol, -999, -999)

  end subroutine dummy_loop_init

!> \section arg_table_dummy_loop_run Argument Table
!! \htmlinclude dummy_loop_run.html
!!
  subroutine dummy_loop_run(loop_cnt, loop_max, ncol, blk_no, thrd_no, errmsg, errflg)
    
    integer,           intent(in   )              :: loop_cnt
    integer,           intent(in   )              :: loop_max
    integer,           intent(in   )              :: ncol
    integer,           intent(in   )              :: blk_no
    integer,           intent(in   )              :: thrd_no

    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    call dummy_loop_internal('dummy_loop_run', loop_cnt, loop_max, ncol, blk_no, thrd_no)

  end subroutine dummy_loop_run

  subroutine dummy_loop_internal(name, loop_cnt, loop_max, ncol, blk_no, thrd_no)
    character(len=*), intent(in) :: name
    integer, intent(in) :: loop_cnt, loop_max, ncol, blk_no, thrd_no

    !
    write(0,'(3a,i6)') 'YYY ', trim(name), ': loop_cnt = ', loop_cnt
    write(0,'(3a,i6)') 'YYY ', trim(name), ': loop_max = ', loop_max
    write(0,'(3a,i6)') 'YYY ', trim(name), ': ncol     = ', ncol
    write(0,'(3a,i6)') 'YYY ', trim(name), ': blk_no   = ', blk_no
    write(0,'(3a,i6)') 'YYY ', trim(name), ': thrd_no  = ', thrd_no
    !
  end subroutine dummy_loop_internal

end module dummy_loop