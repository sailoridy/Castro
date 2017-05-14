
#ifdef CUDA
  attributes(global) &
  subroutine cuda_enforce_consistent_e(lo,hi,state,s_lo,s_hi)

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR => NVAR_d
    use castro_util_module, only: enforce_consistent_e

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    attributes(managed) :: state
    attributes(device) :: lo, hi, s_lo, s_hi

    integer :: my_lo(3), my_hi(3)

    ! Get our spatial index based on the CUDA thread index

    my_lo(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    my_lo(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    my_lo(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    my_hi = my_lo

    if (my_hi(1) .gt. hi(1) .or. my_hi(2) .gt. hi(2) .or. my_hi(3) .gt. hi(3)) return

    call enforce_consistent_e(my_lo, my_hi, state, s_lo, s_hi)

  end subroutine cuda_enforce_consistent_e
#endif



  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi) bind(c, name='ca_enforce_consistent_e')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaDeviceSynchronize, dim3
    use cuda_module, only: threads_and_blocks
#else
    use castro_util_module, only: enforce_consistent_e
#endif

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

#ifdef CUDA

    attributes(managed) :: state

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: s_lo_d(3), s_hi_d(3)

    integer :: cuda_result
    type(dim3) :: numThreads, numBlocks

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    call cuda_enforce_consistent_e<<<numBlocks, numThreads>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

    cuda_result = cudaDeviceSynchronize()

#else

    call enforce_consistent_e(lo, hi, state, s_lo, s_hi)

#endif

  end subroutine ca_enforce_consistent_e
