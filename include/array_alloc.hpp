
#ifndef _ARR_ALLOC_
#define _ARR_ALLOC_

template <class T> T* alloc_1darr(int size_m);

template <class T> T** alloc_2darr(int size_m, int size_n, bool contig=true);

template <class T> T*** alloc_3darr(int size_m, int size_n, int size_p, bool contig=true);

template <class T> void dealloc_1darr(T* arr);

template <class T> void dealloc_2darr(int size_m, T** arr, bool contig=true);

template <class T> void dealloc_3darr(int size_m, int size_n, T*** arr, bool contig=true);

template <class T> T* deep_copy_1darr(int size_m, T* arr);

template <class T> T** deep_copy_2darr(int size_m, int size_n, T** arr, bool contig=true);

template <class T> T*** deep_copy_3darr(int size_m, int size_n, int size_p, T*** arr, bool contig=true);

#endif
