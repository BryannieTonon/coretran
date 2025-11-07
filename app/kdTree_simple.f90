program geographic_network_kdtree
use variableKind, only: i32, r64
use m_allocate, only: allocate
use m_deallocate, only: deallocate
use m_random, only: rngUniform
use m_KdTree, only: KdTree, KdTreeSearch
use dArgDynamicArray_class, only: dArgDynamicArray
use rndgen_mod
use read_tools
!use tictoc_mod

implicit none

! COMANDOs PARA COMPILAR 

! gfortran -I./include app/rndgen.f90 app/read_tools.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib
! gfortran -I./include app/rndgen.f90 app/tictoc_mod.f90  app/read_tools.f90 app/kdTree_simple.f90 ./lib/libcoretran.so -o kdTree_simple -Wl,-rpath=./lib  
! ./kdTree_simple

! Mantendo um fixo e o outro variando
    ! for N in {100..700..20}; do     L=500;     ./kdTree_simple $L $N; done

! Variando N e L
    ! for N in {100..700..20}; do for L in {100..1000..100}; do ./kdTree_simple $L $N done done

! Variando N, L e fator_escala
    ! for N in {100..10000..20}; do for L in {100..1000..100}; do ./kdTree_simple $L $N done done
    ! for N in {5500..10000..500}; do     for L in {500..1000..500}; do         for F in {5..20..5}; do             ./kdTree_simple $L $N $F;         done;     done; done


! Parâmetros
integer(i32) :: node_1, node_2, N, fator_escala
integer(i32) :: i, j, k, unidade_arquivo, unidade_tempo, int_L, int_fatorescala, num, id_amostra
integer(i32) :: seed 
integer(i32) :: probabilidade, distribuicao, raio_min, raio_max

!real(r64), parameter :: rmin = 15.0_r64
real(r64) :: rmin
real(r64) :: dx, dy, dist, dist2, dist_min, L, rmax, R, p_max, pj
real(r64), parameter :: gamma = 3.0

real(r64), allocatable :: x(:), y(:), raio(:), a(:)
integer(i32), allocatable :: degree(:)
integer(i32), allocatable :: neighborList(:,:) 

character(len=*), parameter :: formatador = '(*(g0,x))'
character(len=200) :: comando, dist_aux


! Tipos
type(KdTree) :: tree
type(dArgDynamicArray) :: neighbors, da
type(rndgen) :: rnd, rnd_a
type(KdTreeSearch) :: search
!type(tictoc_t)    :: ctimer
  
! INPUT
if (command_argument_count() /= 5) then
    print*, 'Faltam argumentos:'
    print*,"for N in {100..500..10}; do L=15"
    print*, "./kdTree_simple $N $L done"
    stop
end if

L = read_int(1)
N = read_int(2)
fator_escala = read_int(3)
id_amostra = read_int(4)
rmin = read_int(5)

rmax = L / (1.0_dp*fator_escala)

seed = 294727492 + id_amostra
    
! Alocação e sorteio
call rnd%init(seed)
call rnd_a%init(seed)

allocate(x(N), y(N), degree(N), raio(N), a(N))
allocate(neighborList(N,N))
neighborList = 0
degree = 0
raio = 0
int_fatorescala = int(fator_escala)
int_L = int(L)

! Raio de influência
distribuicao = 2

selectcase(distribuicao)

case(1) ! Distribuição Uniforme
    raio = rnd%rnd_array(N,rmin,rmax) ! Uniforme
    dist_aux = "un"
case(2) ! Distribuição Lei de Potência
    do i = 1, N 
        raio(i) = (rnd%rnd() * (rmax**(1.0_dp - gamma) - rmin**(1.0_dp - gamma)) &
        + rmin**(1.0_dp - gamma) )**(1.0_dp / (1.0_dp - gamma))
    end do
    dist_aux = "pw_"//aux_name_real(gamma)//"_"//aux_name_real(rmin)
end select

! Coordenadas dos vértices
x = rnd%rnd_array(N,0.0_r64,L)
y = rnd%rnd_array(N,0.0_r64,L)

! Constrói KD-Tree
tree = KdTree(x, y)


! Casos: (1) simples; (2) depedência com a distância e (3) atratividade
probabilidade = 3
selectcase(probabilidade) 

case(1) ! Simples ################################################################################################################
 
    open(newunit=unidade_arquivo,file='lista-de-pares.dat',action='write',status='unknown')
    !open(newunit=unidade_tempo, file="tempo-CPU-L_"// aux_name_int(int_L)// &
    ! "-N_"// aux_name_int(N)// "-"//aux_name_int(int_fatorescala) // &
    ! "kdtree.dat",position="append",action='write',status='unknown')

    !call ctimer%reset()
    !call ctimer%tic()

    do i = 1, N ! Pega todos os pontos próximos (k=N-1) e filtra pelo raio
        R = raio(i)
        ! vizinhos ordenados pela distância
        neighbors = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), radius = R)
        k = 0
        do j = 1, neighbors%size()
            dx = x(neighbors%i%values(j)) - x(i)
            dy = y(neighbors%i%values(j)) - y(i)
            dist2 = dx*dx + dy*dy
            if ( (dist2 <= R*R) .and. (i /= neighbors%i%values(j))) then
                k = k + 1
                neighborList(i,k) = neighbors%i%values(j)
                node_1 = i
                node_2 = neighbors%i%values(j)
                call swap(node_1, node_2)
                write(unidade_arquivo,formatador) node_1, node_2
                ! call ctimer%toc()
            end if
        end do
        degree(i) = k
    end do
    close(unidade_arquivo)

    ! SAÍDAS

    ! Coordendas e raios para plotar a rede
    ! open(newunit=unidade_arquivo,file='coordenadas-raio-simples-'//trim(dist_aux)// &
    !         '-' //aux_name_int(int_fatorescala)//'.dat',action='write',status='unknown')
    ! do i = 1 , N
    !     write(unidade_arquivo,formatador) degree(i), x(i), y(i), raio(i)
    ! end do
    ! close(unidade_arquivo)

    ! lista de aresta
    !comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat | uniq > lista-de-pares-simples-metricas.dat'
    comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat| uniq > ' // &
            'L_'// aux_name_int(int_L)// &
            '-N_'// aux_name_int(N)// &
            '-' //aux_name_int(int_fatorescala) //'-' //aux_name_int(id_amostra) //'-'//trim(dist_aux)//'-simples.dat'
    call system(trim(comando))

    !write(unidade_tempo, formatador) N, int(L), fator_escala, ctimer%t_tot
    !close(unidade_tempo)

    deallocate(x, y, degree)
    deallocate(neighborList)
    call tree%deallocate()

case(2) ! Dependente da distância ########################################################################################################
        
    open(newunit=unidade_arquivo,file='lista-de-pares.dat',action='write',status='unknown')
    
    !call ctimer%reset()
    !call ctimer%tic()

    do i = 1, N ! Pega todos os pontos próximos (k=N-1) e filtra pelo raio
        R = raio(i)
        ! vizinhos ordenados pela distância
        neighbors = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), radius = R)
        k = 0
        ! Encontra quem é o vizinho mais próximo e calcula a probabilidade máxima
        do j = 1, neighbors%size()
            if (neighbors%i%values(j) /= i) then
                dx = x(neighbors%i%values(j)) - x(i)
                dy = y(neighbors%i%values(j)) - y(i)
                dist2 = dx*dx + dy*dy
                dist_min = sqrt(dist2)
                !Probabilidade máxima
                p_max = 1.0_dp / dist_min
                ! Precisa sair no primeiro vizinho diferente de i que encontrar
                exit 
            end if
        end do

        ! Processo de conexão
        do j = 1, neighbors%size()
            if (neighbors%i%values(j) /= i) then
                dx = x(neighbors%i%values(j)) - x(i)
                dy = y(neighbors%i%values(j)) - y(i)
                dist2 = dx*dx + dy*dy
                dist = sqrt(dist2)
                pj = 1.0_dp / dist 
                if ( rnd%rnd() <= pj/p_max .and. (i /= neighbors%i%values(j))) then
                    k = k + 1
                    neighborList(i,k) = neighbors%i%values(j)
                    node_1 = i
                    node_2 = neighbors%i%values(j)

                    call swap(node_1, node_2)
                    write(unidade_arquivo,formatador) node_1, node_2
                end if
            end if
        end do
        degree(i) = k
    end do
    close(unidade_arquivo)

    ! SAÍDAS

    ! ! Coordendas e raios para plotar a rede
    ! open(newunit=unidade_arquivo,file='coordenadas-raio-inversamente-'//trim(dist_aux) // &
    !         '-' //aux_name_int(int_fatorescala)//'.dat',action='write',status='unknown')
    ! do i = 1 , N
    !     write(unidade_arquivo,formatador) degree(i), x(i), y(i), raio(i)
    ! end do
    ! close(unidade_arquivo)

    ! lista de aresta
    !comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat | uniq > lista-de-pares-simples-metricas.dat'
    comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat| uniq > ' // &
            'L_'// aux_name_int(int_L)// &
            '-N_'// aux_name_int(N)// &
            '-' //aux_name_int(int_fatorescala) //'-' //aux_name_int(id_amostra) //'-'//trim(dist_aux)//'-inversamente.dat'
    call system(trim(comando))

    deallocate(x, y, degree)
    deallocate(neighborList)
    call tree%deallocate()

case(3) ! Atratividade ################################################################################################################
    
    open(newunit=unidade_arquivo,file='lista-de-pares.dat',action='write',status='unknown')

    a = rnd_a%rnd_array(N,0.0_r64,1.0_r64)

    do i = 1, N 
        R = raio(i)
        ! vizinhos ordenados pela distância
        neighbors = search%kNearest(tree, x, y, xQuery = x(i), yQuery = y(i), radius = R)
        k = 0
        do j = 1, neighbors%size()
            pj =  a(neighbors%i%values(j)) 
            if ( rnd%rnd() <= pj .and. (i /= neighbors%i%values(j))) then 
                k = k + 1
                neighborList(i,k) = neighbors%i%values(j)
                node_1 = i
                node_2 = neighbors%i%values(j)

                call swap(node_1, node_2)
                write(unidade_arquivo,formatador) node_1, node_2
            end if
        end do
        degree(i) = k
    end do
    close(unidade_arquivo)

    ! Saídas

    ! ! Coordendas e raios para plotar a rede
    ! open(newunit=unidade_arquivo,file='coordenadas-raio-atratividade-'//trim(dist_aux) // &
    !         '-' //aux_name_int(int_fatorescala)//'.dat',action='write',status='unknown')
    ! do i = 1 , N
    !     write(unidade_arquivo,formatador) degree(i), x(i), y(i), raio(i)
    ! end do
    ! close(unidade_arquivo)

    comando = 'sort -n -k1,1 -k2,2 lista-de-pares.dat| uniq > ' // &
            'L_'// aux_name_int(int_L)// &
            '-N_'// aux_name_int(N)// &
            '-' //aux_name_int(int_fatorescala) //'-' //aux_name_int(id_amostra) //'-'//trim(dist_aux)//'-atratividade.dat'
    call system(trim(comando))
    
    deallocate(x, y, degree)
    deallocate(neighborList)
    call tree%deallocate()

end select

contains

    subroutine swap(a, b)
        integer(i32), intent(inout) :: a, b
        integer(i32) :: temp

        if (a > b) then
            temp = a
            a = b
            b = temp
        end if
    end subroutine swap

end program geographic_network_kdtree
