
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

    attributes(device) :: lo, hi, s_lo, s_hi, state

    integer :: idx(3)

    ! Get our spatial index based on the CUDA thread index

    idx(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    idx(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    idx(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    if (idx(1) .gt. hi(1) .or. idx(2) .gt. hi(2) .or. idx(3) .gt. hi(3)) return

    call enforce_consistent_e(idx, idx, state, s_lo, s_hi)

  end subroutine cuda_enforce_consistent_e
#endif



  subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi,idx) bind(c, name='ca_enforce_consistent_e')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR
    use castro_util_module, only: enforce_consistent_e
#ifdef CUDA
    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, cudaDeviceSynchronize, dim3, cuda_stream_kind
    use cuda_module, only: threads_and_blocks, cuda_streams, max_cuda_streams
#endif

    implicit none

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer, intent(in)     :: idx

#ifdef CUDA

    attributes(device) :: state

    integer, device :: lo_d(3), hi_d(3)
    integer, device :: s_lo_d(3), s_hi_d(3)

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    cuda_result = cudaMemcpyAsync(lo_d, lo, 3, cudaMemcpyHostToDevice)
    cuda_result = cudaMemcpyAsync(hi_d, hi, 3, cudaMemcpyHostToDevice)

    cuda_result = cudaMemcpyAsync(s_lo_d, s_lo, 3, cudaMemcpyHostToDevice)
    cuda_result = cudaMemcpyAsync(s_hi_d, s_hi, 3, cudaMemcpyHostToDevice)

    call threads_and_blocks(lo, hi, numBlocks, numThreads)

    stream = cuda_streams(mod(idx, max_cuda_streams) + 1)

    call cuda_enforce_consistent_e<<<numBlocks, numThreads, 0, stream>>>(lo_d, hi_d, state, s_lo_d, s_hi_d)

#else

    call enforce_consistent_e(lo, hi, state, s_lo, s_hi)

#endif

  end subroutine ca_enforce_consistent_e
