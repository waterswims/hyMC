#ifndef _ARR_ALLOC_
#define _ARR_ALLOC_

///////////////////////////////////////////////////////////////////////////
/// \brief Allocates a C array in memory.
///
/// \param size_m The size of the array.
/// \return A pointer to the allocated array.
///////////////////////////////////////////////////////////////////////////
template <class T> T* alloc_1darr(int size_m);

///////////////////////////////////////////////////////////////////////////
/// \brief Allocates a two-dimensional C array in memory.
///
/// \param size_m The size of the array in the first dimension.
/// \param size_n The size of the array in the second dimension.
/// \param contig A boolean determining if the 2D array is contiguous in memory.
///               Defaults to true.
/// \return A double pointer too the 2D array.
///////////////////////////////////////////////////////////////////////////
template <class T> T** alloc_2darr(int size_m, int size_n, bool contig=true);

///////////////////////////////////////////////////////////////////////////
/// \brief Allocates a three-dimensional C array in memory.
///
/// \param size_m The size of the array in the first dimension.
/// \param size_n The size of the array in the second dimension.
/// \param size_p The size of the array in the third dimension.
/// \param contig A boolean determining if the 3D array is contiguous in
///               memory. Defaults to true.
/// \return A triple pointer too the 3D array.
///////////////////////////////////////////////////////////////////////////
template <class T> T*** alloc_3darr(int size_m, int size_n, int size_p, bool contig=true);

///////////////////////////////////////////////////////////////////////////
/// \brief Deacllocates a standard C array.
///
/// \param arr The array to be deallocated.
///////////////////////////////////////////////////////////////////////////
template <class T> void dealloc_1darr(T* arr);

///////////////////////////////////////////////////////////////////////////
/// \brief Deacllocates a 2D C array.
///
/// \param size_m The size of the array in the first dimension.
/// \param arr The array to be deallocated.
/// \param contig Was the array allocated contiguously? Defaults to true.
///////////////////////////////////////////////////////////////////////////
template <class T> void dealloc_2darr(int size_m, T** arr, bool contig=true);

///////////////////////////////////////////////////////////////////////////
/// \brief Deacllocates a 3D C array.
///
/// \param size_m The size of the array in the first dimension.
/// \param size_n The size of the array in the second dimension.
/// \param arr The array to be deallocated.
/// \param contig Was the array allocated contiguously? Defaults to true.
///////////////////////////////////////////////////////////////////////////
template <class T> void dealloc_3darr(int size_m, int size_n, T*** arr, bool contig=true);

///////////////////////////////////////////////////////////////////////////
/// \brief Performs a deep copy of a C array, returning a newly allocated
///        array.
///
/// \param size_m The size of the array to be copied.
/// \param arr The array to be copied.
/// \return A pointer to a newly allocated C array where each element
///         matches the input array.
///////////////////////////////////////////////////////////////////////////
template <class T> T* deep_copy_1darr(int size_m, const T* arr);

///////////////////////////////////////////////////////////////////////////
/// \brief Performs a deep copy of a 2D C array, returning a newly
///        allocated array.
///
/// \param size_m The size of the array to be copied in the first
///               dimension.
/// \param size_n The size of the array to be copied in the second
///               dimension.
/// \param arr The array to be copied.
/// \param contig Should the new array be allocated contiguously?
/// \return A pointer to a newly allocated C array where each element
///         matches the input array.
///////////////////////////////////////////////////////////////////////////
template <class T> T** deep_copy_2darr(int size_m, int size_n, const T** arr, bool contig=true);

///////////////////////////////////////////////////////////////////////////
/// \brief Performs a deep copy of a 3D C array, returning a newly
///        allocated array.
///
/// \param size_m The size of the array to be copied in the first
///               dimension.
/// \param size_n The size of the array to be copied in the second
///               dimension.
/// \param size_p The size of the array to be copied in the third
///               dimension.
/// \param arr The array to be copied.
/// \param contig Should the new array be allocated contiguously?
/// \return A pointer to a newly allocated C array where each element
///         matches the input array.
///////////////////////////////////////////////////////////////////////////
template <class T> T*** deep_copy_3darr(int size_m, int size_n, int size_p, const T*** arr, bool contig=true);

#endif
