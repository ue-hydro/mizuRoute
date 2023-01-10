interface
    function create_openwq_c() bind(C, name="create_openwq")

        use iso_c_binding
        implicit none
        type(c_ptr) :: create_openwq_c

    end function

    function openwq_decl_c(openWQ) bind(C, name="openwq_decl")
    !function openwq_decl_c(      &
    !    openWQ,                 &
    !    num_hru,                &
    !    nCanopy_2openwq,        &
    !    nSnow_2openwq,          &
    !    nSoil_2openwq,          &
    !    nRunoff_2openwq,        &
    !    nAquifer_2openwq,       &
    !    y_direction,            &
    !    ) bind(C, name="openwq_decl")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_decl_c ! returns a return value of 0 (success) or -1 (failure)
        type(c_ptr), intent(in), value :: openWQ
        !integer(c_int), intent(in), value  :: num_hru
        !integer(c_int), intent(in), value  :: nCanopy_2openwq
        !integer(c_int), intent(in), value  :: nSnow_2openwq
        !integer(c_int), intent(in), value  :: nSoil_2openwq
        !integer(c_int), intent(in), value  :: nAquifer_2openwq
        !integer(c_int), intent(in), value  :: nRunoff_2openwq
        !integer(c_int), intent(in), value  :: y_direction

    end function


end interface