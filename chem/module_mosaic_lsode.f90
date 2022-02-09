module module_mosaic_lsode
  
  
  
  implicit none
  private

  public :: MOSAIC_LSODE

contains
  subroutine MOSAIC_LSODE(dtchem)
    use module_data_mosaic_kind

    implicit none


    
    real(r8), intent(in) :: dtchem

    call wrf_error_fatal3("<stdin>",20,&
'*** error - mosaic_lsode has been deactivated ***')    

    return
  end subroutine MOSAIC_LSODE

end module module_mosaic_lsode
