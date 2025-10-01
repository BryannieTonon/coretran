module listas_mod
   use iso_fortran_env, only : i4 => int32, i8 => int64, sp => real32, dp => real64
   implicit none

   character(len=*), parameter :: formatador = '(*(g0,x))'

   type :: lista_dinamica
        integer(kind=i4) :: tamanho 
        integer(kind=i4) :: posicao_final
        integer(kind=i4) :: inicio, final ! lista circular
        integer(kind=i4), allocatable :: lista(:)
    contains
        procedure :: inicializar => inicializar_lista
        procedure :: adicionar => adicionar_lista
        procedure :: remover => remover_lista
        procedure :: remover_dupla => remover_dupla
    end type lista_dinamica

    ! exemplo de uso:
    ! type(lista_dinamica) :: lista
    ! call lista%inicializar(10)
    ! lista%lista = (valores nao nulos)

contains

    subroutine inicializar_lista(lista, tamanho)
        class(lista_dinamica), intent(inout) :: lista
        integer(kind=i4), intent(in) :: tamanho

        lista%tamanho = tamanho
        lista%posicao_final = 0
        lista%inicio = 1
        lista%final = 1
        if (.not. allocated(lista%lista)) allocate(lista%lista(tamanho))

        ! assumindo que a lista será totalmente preenchida logo em seguida

    end subroutine inicializar_lista
    
    subroutine adicionar_lista(this, valor)
        class(lista_dinamica), intent(inout) :: this
        integer(kind=i4), intent(in) :: valor

        this%posicao_final = this%posicao_final + 1
        this%lista(this%posicao_final) = valor
    end subroutine

    subroutine remover_lista(this, posicao)
        class(lista_dinamica), intent(inout) :: this
        integer(kind=i4) :: posicao
        
        this%lista(posicao) = this%lista(this%posicao_final)
        this%posicao_final = this%posicao_final - 1

    end subroutine remover_lista

    subroutine remover_dupla(this, posicao1, posicao2)
        class(lista_dinamica), intent(inout) :: this
        integer(kind=i4), intent(in) :: posicao1, posicao2
        integer(kind=i4) :: posicao1_remover, posicao2_remover

        if (posicao1 < posicao2) then
            posicao1_remover = posicao2
            posicao2_remover = posicao1
        else
            posicao1_remover = posicao1
            posicao2_remover = posicao2
        end if

        call this%remover(posicao1_remover) ! posicao1 é maior que posicao2
        call this%remover(posicao2_remover)
    end subroutine

end module