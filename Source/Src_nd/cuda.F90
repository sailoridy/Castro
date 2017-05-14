module cuda_module

  implicit none

contains

  subroutine threads_and_blocks(lo, hi, numBlocks, numThreads)

    use cudafor, only: dim3
    use prob_params_module, only: dim

    implicit none

    integer, intent(in)       :: lo(3), hi(3)
    type(dim3), intent(inout) :: numBlocks, numThreads

    integer :: tile_size(3)

    if (dim .eq. 1) then
       numThreads % x = 256
       numThreads % y = 1
       numThreads % z = 1
    else if (dim .eq. 2) then
       numThreads % x = 16
       numThreads % y = 16
       numThreads % z = 1
    else
       numThreads % x = 8
       numThreads % y = 8
       numThreads % z = 8
    endif

    tile_size = hi - lo + 1

    numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
    numBlocks % y = (tile_size(2) + numThreads % y - 1) / numThreads % y
    numBlocks % z = (tile_size(3) + numThreads % z - 1) / numThreads % z

  end subroutine threads_and_blocks

end module cuda_module
