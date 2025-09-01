! Module for derivative-free optimization of `real`s with a genetic algorithm.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

! Notation follows luke_essentials_2013 p. 31 in most instances.
! Instead of "fitness" or "quality", I am using the term objective function.
! In contrast to luke_essentials_2013, lower objective function values are better.

! Make the simulation fully specified by the chromosome.
! Then no additional `config` type needs to be passed in beyond the genetic algorithm configuration.

module ga

use prec, only: WP
implicit none
private

public :: init_pop, mutate_indiv, cross_two_indivs, select_indiv, evaluate, optimize_ga
public :: constraint_lt, constraint_gt
public :: standard_ga_config

integer, parameter          :: MAX_SAMPLES = 10000
character(len=*), parameter :: GENER_FMT = "(i8)"

type, public :: ga_config
    integer :: n_genes = 0 ! number of genes (default set to zero to catch when not set)
    
    ! martins_engineering_2021 p. 309:
    ! > As a rule of thumb, the population size should be approximately one order of magnitude larger than the number of design
    ! > variables, and this size should be tuned.
    ! I guess `n_pop = 50` is a decent starting point.
    integer :: n_pop = 50 ! number of indivs in population
    
    ! de_jong_analysis_1975 p. 68 (pdf pp. 83):
    ! > a mutation rate of the same order of magnitude as `1/POP_SIZE` seems to be about the best setting
    ! So you could try `p_mutate = min(1.0, p_mutate_factor/real(n_pop, WP))`.
    ! de_jong_analysis_1975 p. 70 (pdf pp. 85): for n_pop = 50, p_mutate = 0.02 is about optimal
    real(WP) :: p_mutate = 0.02_WP
    
    ! de_jong_analysis_1975 pp. 75 (pdf pp. 90): p_cross_indiv = 0.6 is about optimal
    real(WP) :: p_cross_indiv = 0.6_WP
    
    ! `n_select = 2` is most popular according to luke_essentials_2013 p. 45.
    integer :: n_select = 2
    
    ! `p_cross_gene = 0.5_WP` was the first proposed suggestion according to luke_essentials_2013 p. 39.
    real(WP) :: p_cross_gene = 0.5_WP
    
    integer  :: n_gener = 20
    real(WP) :: rel_b   = 0.2_WP ! `b` parameter for Cauchy dist., relative to range of variable determined from `lb` and `ub`
    
    real(WP), allocatable :: lb(:), ub(:) ! lower and upper bounds for each gene
    
    real(WP) :: f_max_all_infeasible   = 0.0_WP
    logical  :: stop_if_all_unfeasible = .true.
    
    logical          :: progress = .true.
    character(len=8) :: f_fmt = "f12.2"
end type ga_config

type, public :: indiv_type
    real(WP), allocatable :: chromo(:)
    
    ! whether the objective function has been set (`.false.` by default)
    logical :: set = .false.
    
    ! objective function value
    real(WP) :: f = huge(1.0_WP)
    
    ! sum of constraint violations
    real(WP) :: sum_g = 0.0_WP
end type indiv_type

type, public :: pop_type
    type(indiv_type), allocatable :: indivs(:)
    
    type(indiv_type) :: best_pop_indiv, best_ever_indiv
end type pop_type

contains

subroutine init_pop(config, rng, pop)
    use purerng, only: rng_type
    use checks, only: assert
    
    type(ga_config), intent(in)    :: config
    type(rng_type), intent(in out) :: rng
    type(pop_type), intent(out)    :: pop
    
    integer :: i_pop, i_gene
    
    call assert(allocated(config%lb), "ga (init_pop): allocated(config%lb) violated")
    call assert(allocated(config%ub), "ga (init_pop): allocated(config%ub) violated")
    
    ! The assertions on `pop` are pointless because `pop` is `intent(out)`.
    
    call assert(.not. allocated(pop%indivs), "ga (init_pop): .not. allocated(pop%indivs) violated")
    allocate(pop%indivs(config%n_pop))
    
    call assert(config%n_genes > 0, "ga (init_pop): config%n_genes > 0 violated")
    
    pop_loop: do i_pop = 1, config%n_pop ! SERIAL
        call assert(.not. allocated(pop%indivs(i_pop)%chromo), "ga (init_pop): .not. allocated(pop%indivs(i_pop)%chromo)")
        allocate(pop%indivs(i_pop)%chromo(config%n_genes))
        
        gene_loop: do i_gene = 1, config%n_genes ! SERIAL
            call assert(config%lb(i_gene) <= config%ub(i_gene), "ga (init_pop): lb > ub?")
            
            call rng%uniform(config%lb(i_gene), config%ub(i_gene), pop%indivs(i_pop)%chromo(i_gene))
            
            call assert(pop%indivs(i_pop)%chromo(i_gene) >= config%lb(i_gene), "ga (init_pop): lower bound violated")
            call assert(pop%indivs(i_pop)%chromo(i_gene) <= config%ub(i_gene), "ga (init_pop): upper bound violated")
        end do gene_loop
    end do pop_loop
    
    pop%best_pop_indiv%f    = huge(1.0_WP)
    pop%best_pop_indiv%set  = .true.
    pop%best_ever_indiv%f   = huge(1.0_WP)
    pop%best_ever_indiv%set = .true.
end subroutine init_pop

pure subroutine mutate_indiv(config, rng, indiv)
    ! Follows luke_essentials_2013 Algorithm 11, p. 23.
    
    use purerng, only: rng_type
    use checks, only: assert
    
    type(ga_config), intent(in)      :: config
    type(rng_type), intent(in out)   :: rng
    type(indiv_type), intent(in out) :: indiv
    
    integer  :: i_gene, i_sample
    real(WP) :: nu, b
    
    call assert(config%n_genes > 0, "ga (mutate_indiv): config%n_genes > 0 violated")
    call assert(allocated(config%lb), "ga (mutate_indiv): allocated(config%lb) violated")
    call assert(allocated(config%ub), "ga (mutate_indiv): allocated(config%ub) violated")
    
    do i_gene = 1, config%n_genes ! SERIAL
        call rng%random_number(nu)
        if (config%p_mutate >= nu) then
            call assert(config%lb(i_gene) <= config%ub(i_gene), "ga (mutation): lb > ub?")
            b = config%rel_b * (config%ub(i_gene) - config%lb(i_gene))
            
            i_sample = 0
            do ! SERIAL
                i_sample = i_sample + 1
                if (i_sample > MAX_SAMPLES) then
                    error stop "ga (mutate_indiv): MAX_SAMPLES exceeded"
                    ! TODO: Alternatively, clip to the bounds.
                end if
                
                call rng%cauchy(0.0_WP, b, nu)
                
                if (((indiv%chromo(i_gene) + nu) >= config%lb(i_gene)) &
                        .and. ((indiv%chromo(i_gene) + nu) <= config%ub(i_gene))) then
                    exit
                end if
            end do
            
            indiv%chromo(i_gene) = indiv%chromo(i_gene) + nu
            indiv%set = .false.
            
            call assert(indiv%chromo(i_gene) >= config%lb(i_gene), "ga (mutate_indiv): lower bound violated")
            call assert(indiv%chromo(i_gene) <= config%ub(i_gene), "ga (mutate_indiv): upper bound violated")
        end if
    end do
end subroutine mutate_indiv

pure subroutine cross_two_indivs(config, rng, indiv_1, indiv_2)
    ! Uniform crossover
    ! Follows luke_essentials_2013 Algorithm 25, p. 39.
    
    use purerng, only: rng_type
    use checks, only: assert, assert_dimension
    
    type(ga_config), intent(in)      :: config
    type(rng_type), intent(in out)   :: rng
    type(indiv_type), intent(in out) :: indiv_1, indiv_2
    
    integer  :: i_gene
    real(WP) :: nu, gene_temp
    
    call assert(config%n_genes > 0, "ga (cross_two_indivs): config%n_genes > 0 violated")
    call assert_dimension(indiv_1%chromo, indiv_2%chromo)
    call assert(allocated(config%lb), "ga (cross_two_indivs): allocated(config%lb) violated")
    call assert(allocated(config%ub), "ga (cross_two_indivs): allocated(config%ub) violated")
    
    do i_gene = 1, config%n_genes ! SERIAL
        call rng%random_number(nu)
        if (config%p_cross_indiv >= nu) then
            call assert(indiv_1%chromo(i_gene) >= config%lb(i_gene), "ga (cross_two_indivs): lower bound violated (1)")
            call assert(indiv_1%chromo(i_gene) <= config%ub(i_gene), "ga (cross_two_indivs): upper bound violated (1)")
            call assert(indiv_2%chromo(i_gene) >= config%lb(i_gene), "ga (cross_two_indivs): lower bound violated (2)")
            call assert(indiv_2%chromo(i_gene) <= config%ub(i_gene), "ga (cross_two_indivs): upper bound violated (2)")
            
            gene_temp              = indiv_1%chromo(i_gene)
            indiv_1%chromo(i_gene) = indiv_2%chromo(i_gene)
            indiv_2%chromo(i_gene) = gene_temp
            
            indiv_1%set = .false.
            indiv_2%set = .false.
        end if
    end do
end subroutine cross_two_indivs

pure subroutine select_indiv(config, rng, pop, indiv)
    ! Tournament selection
    ! Follows luke_essentials_2013 Algorithm 32, p. 45.
    
    use purerng, only: rng_type
    use checks, only: assert, assert_dimension
    
    type(ga_config), intent(in)    :: config
    type(rng_type), intent(in out) :: rng
    type(pop_type), intent(in out) :: pop
    type(indiv_type), intent(out)  :: indiv
    
    integer :: nu, i_pop
    
    call assert(config%n_pop == size(pop%indivs), "ga (select_indiv): config%n_pop == size(pop%indivs) violated")
    
    call rng%int(1, config%n_pop, nu)
    indiv = pop%indivs(nu)
    call assert(indiv%set, "ga (select_indiv): objective function not set (1)")
    do i_pop = 2, config%n_select ! SERIAL
        call rng%int(1, config%n_pop, nu)
        
        call assert(pop%indivs(nu)%set, "ga (select_indiv): objective function not set (2)")
        
        if (pop%indivs(nu)%f < indiv%f) then
            indiv = pop%indivs(nu)
        end if
    end do
end subroutine select_indiv

subroutine evaluate(config, objfun, pop)
    use prec, only: WP
    use checks, only: assert, is_close
    
    type(ga_config), intent(in)    :: config
    type(pop_type), intent(in out) :: pop
    
    integer  :: i_pop
    real(WP) :: f_max
    logical  :: best_pop_indiv_set
    
    interface
        pure subroutine objfun(chromo, f, sum_g)
            use prec, only: WP
            
            ! Passing in all `real`s means that the objective function does not need any of this module's derived types.
            real(WP), intent(in)  :: chromo(:)
            real(WP), intent(out) :: f     ! objective function value
            real(WP), intent(out) :: sum_g ! sum of constraint violations
        end subroutine objfun
    end interface
    
    call assert(config%n_pop == size(pop%indivs), "ga (evaluate): config%n_pop == size(pop%indivs) violated")
    call assert(pop%best_ever_indiv%set, "ga (evaluate): pop%best_ever_indiv%set violated")
    
    do concurrent (i_pop = 1:config%n_pop)
        if (.not. pop%indivs(i_pop)%set) then
            call objfun(pop%indivs(i_pop)%chromo, pop%indivs(i_pop)%f, pop%indivs(i_pop)%sum_g)
        end if
        
        call assert(pop%indivs(i_pop)%sum_g >= 0.0_WP, "ga (evaluate): pop%indivs(i_pop)%sum_g >= 0 violated")
    end do
    
    f_max                = -huge(1.0_WP)
    pop%best_pop_indiv%f = huge(1.0_WP)
    best_pop_indiv_set   = .false.
    do i_pop = 1, config%n_pop ! SERIAL
        if (is_close(pop%indivs(i_pop)%sum_g, 0.0_WP)) then
            f_max      = max(f_max, pop%indivs(i_pop)%f)
            pop%indivs(i_pop)%set = .true.
            
            if (pop%best_pop_indiv%f > pop%indivs(i_pop)%f) then
                pop%best_pop_indiv = pop%indivs(i_pop)
                best_pop_indiv_set = .true.
            end if
        end if
    end do
    
    ! deb_efficient_2000 p. 317: > If no feasible solution exists in a population, $f_max$ is set to zero.
    ! I don't like this as `f` could normally be above 0.
    ! I decided to stop with an error by default in this situation.
    if (.not. best_pop_indiv_set) then
        if (config%stop_if_all_unfeasible) then
            error stop "ga (evaluate): all individuals violate constraints " // &
                            "(can override with config%stop_if_all_unfeasible=.false.)"
        else
            f_max = config%f_max_all_infeasible
        end if
    end if
    
    ! set `f` for indivs that had constraint violations
    ! See deb_efficient_2000 eq. 4.
    do concurrent (i_pop = 1:config%n_pop)
        if (pop%indivs(i_pop)%sum_g > 0.0_WP) then
            ! Infeasible individuals have their fitness recalculated based on the current population.
            ! This is regardless of whether they were `set` before `evaluate` was called.
            ! This avoid issues from the fitness depending on the population.
            pop%indivs(i_pop)%f   = f_max + pop%indivs(i_pop)%sum_g
            pop%indivs(i_pop)%set = .true.
        end if
    end do
    
    ! set best ever individual
    if (pop%best_ever_indiv%f > pop%best_pop_indiv%f) then
        pop%best_ever_indiv = pop%best_pop_indiv
    end if
end subroutine evaluate

subroutine optimize_ga(config, rng, objfun, pop, rc)
    ! Follows luke_essentials_2013 Algorithm 20, p. 37.
    
    use purerng, only: rng_type
    use checks, only: assert
    
    type(ga_config), intent(in)    :: config
    type(rng_type), intent(in out) :: rng
    type(pop_type), intent(in out) :: pop
    integer, intent(out)           :: rc ! TODO: return codes
    
    integer :: i_gener, i_pop
    
    type(pop_type) :: next_pop
    
    interface
        pure subroutine objfun(chromo, f, sum_g)
            use prec, only: WP
            
            ! Passing in all `real`s means that the objective function does not need any of this module's derived types.
            real(WP), intent(in)  :: chromo(:)
            real(WP), intent(out) :: f     ! objective function value
            real(WP), intent(out) :: sum_g ! sum of constraint violations
        end subroutine objfun
    end interface
    
    call assert(allocated(pop%indivs), "ga (optimize_ga): allocated(pop%indivs) violated")
    call assert(config%n_pop == size(pop%indivs), "ga (optimize_ga): config%n_pop == size(pop%indivs) violated")
    call assert(config%n_pop > 0, "ga (optimize_ga): config%n_pop > 0 violated")
    call assert(config%n_genes > 0, "ga (optimize_ga): config%n_genes > 0 violated")
    call assert(config%n_gener > 0, "ga (optimize_ga): config%n_gener > 0 violated")
    call assert(allocated(config%lb), "ga (optimize_ga): allocated(config%lb) violated")
    call assert(allocated(config%ub), "ga (optimize_ga): allocated(config%ub) violated")
    
    allocate(next_pop%indivs(config%n_pop))
    next_pop%best_ever_indiv    = pop%best_ever_indiv
    next_pop%best_pop_indiv%set = .false.
    
    if (config%progress) then
        write(unit=*, fmt="(a)") "   gener    pop best   best ever"
        write(unit=*, fmt=GENER_FMT, advance="no") 0
    end if
    call evaluate(config, objfun, pop)
    if (config%progress) then
        write(unit=*, fmt="(2" // trim(config%f_fmt) // ")") pop%best_pop_indiv%f, pop%best_ever_indiv%f
    end if
    
    rc = 0
    do i_gener = 1, config%n_gener ! SERIAL
        do i_pop = 1, config%n_pop, 2 ! SERIAL
            ! I don't have a parallel RNG right now, so this can't be parallelized at the moment.
            ! Not that parallelizing this would matter much anyway as it's quick.
            
            call select_indiv(config, rng, pop, next_pop%indivs(i_pop))
            call select_indiv(config, rng, pop, next_pop%indivs(i_pop + 1))
            call cross_two_indivs(config, rng, next_pop%indivs(i_pop), next_pop%indivs(i_pop + 1))
            call mutate_indiv(config, rng, next_pop%indivs(i_pop))
            call mutate_indiv(config, rng, next_pop%indivs(i_pop + 1))
        end do
        
        pop%indivs = next_pop%indivs
        
        if (config%progress) then
            write(unit=*, fmt=GENER_FMT, advance="no") i_gener
        end if
        call evaluate(config, objfun, pop)
        if (config%progress) then
            write(unit=*, fmt="(2" // trim(config%f_fmt) // ")") pop%best_pop_indiv%f, pop%best_ever_indiv%f
        end if
    end do
end subroutine optimize_ga

pure subroutine constraint_lt(x, y, sum_g)
    use checks, only: assert
    
    real(WP), intent(in)     :: x, y
    real(WP), intent(in out) :: sum_g
    
    call assert(sum_g >= 0.0_WP, "ga (constraint_lt): sum_g >= 0 violated (1)")
    
    if (.not. (x < y)) then
        sum_g = sum_g + (x - y)
    end if
    
    call assert(sum_g >= 0.0_WP, "ga (constraint_lt): sum_g >= 0 violated (2)")
end subroutine constraint_lt

pure subroutine constraint_gt(x, y, sum_g)
    use checks, only: assert
    
    real(WP), intent(in)     :: x, y
    real(WP), intent(in out) :: sum_g
    
    call assert(sum_g >= 0.0_WP, "ga (constraint_gt): sum_g >= 0 violated (1)")
    
    if (.not. (x > y)) then
        sum_g = sum_g + (y - x)
    end if
    
    call assert(sum_g >= 0.0_WP, "ga (constraint_gt): sum_g >= 0 violated (2)")
end subroutine constraint_gt

pure subroutine standard_ga_config(n_genes, config)
    integer, intent(in)          :: n_genes
    type(ga_config), intent(out) :: config
    
    config%n_genes = n_genes
    
    ! martins_engineering_2021 p. 309:
    ! > As a rule of thumb, the population size should be approximately one order of magnitude larger than the number of design
    ! > variables, and this size should be tuned.
    config%n_pop = 10*n_genes
    
    ! de_jong_analysis_1975 p. 68 (pdf pp. 83):
    ! > a mutation rate of the same order of magnitude as `1/POP_SIZE` seems to be about the best setting
    config%p_mutate = 1.0_WP/real(config%n_pop, WP)
end subroutine standard_ga_config

end module ga
