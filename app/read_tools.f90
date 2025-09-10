module read_tools
    use iso_fortran_env, only : i4 => int32, i8 => int64, sp => real32, dp => real64
    implicit none

contains
   
    function read_int(i) result(value) ! recebe um argumento da linha de comando, converte esse argumento para um número inteiro e retorna o valor convertido
        integer(kind=i4), intent(in) :: i
        integer(kind=i4) :: value
        character(len=1000) :: argument 

        call get_command_argument(i, argument)  ! input do terminal - (posicao/quantidade,argumento)
        read(argument, *) value                 ! a função read converte a string armazenada em argument para um número inteiro e armazena em value

    end function read_int

    function aux_name_int(N) result(name)
        ! número inteiro N é convertido para string dentro de aux usando write(aux, '(I0)').
        ! os espaços extras são removidos com trim(adjustl(aux)).
        ! a string final é armazenada em name (que tem tamanho ajustável)
        integer(kind=i4), intent(in) :: N
        character(len=1000) :: aux 
        character(len=:), allocatable :: name

        write(aux, '(I0)') N
        name = trim(adjustl(aux))
    end function aux_name_int
    
    function aux_name_real(N) result(name_real)
        real(kind=dp), intent(in) :: N
        character(len=1000) :: aux 
        character(len=:), allocatable :: name_real

        write(aux, '(F0.3)') N
        name_real = trim(adjustl(aux))
    end function aux_name_real

end module read_tools