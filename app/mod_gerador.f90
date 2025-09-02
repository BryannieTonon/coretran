module rdn_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
contains

  subroutine seed_with_time()
    ! Sema o gerador com system_clock (para variar a cada execução)
    integer :: n, i
    integer, allocatable :: seed(:)
    integer :: clk, count_rate, count_max

    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clk, count_rate=count_rate, count_max=count_max)
    do i = 1, n
      seed(i) = int(123456 + clk + i*37)  ! maneira simples de criar valores
    end do
    call random_seed(put = seed)
    deallocate(seed)
  end subroutine seed_with_time

  subroutine generate_random_geo(N, L, rmin, rmax, x, y, raio)
    ! Gera arrays aleatórios de coordenadas e raios
    integer, intent(in) :: N
    real(real64), intent(in) :: L, rmin, rmax
    real(real64), allocatable, intent(out) :: x(:), y(:), raio(:)

    allocate(x(N), y(N), raio(N))

    ! Sema o gerador
    call seed_with_time()

    ! Gera valores uniformes (0,1)
    call random_number(x)
    call random_number(y)
    call random_number(raio)

    ! Escala para os intervalos desejados
    x = 0.0_real64 + (L - 0.0_real64) * x
    y = 0.0_real64 + (L - 0.0_real64) * y
    raio = rmin + (rmax - rmin) * raio
  end subroutine generate_random_geo

end module rdn_mod
