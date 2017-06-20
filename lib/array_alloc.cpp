#include "../include/array_alloc.hpp"
#include <cstdlib>

template <class T>
T* alloc_1darr(int size_m)
{
    return (T*)malloc(sizeof(T)*size_m);
}

template <class T>
T** alloc_2darr(int size_m, int size_n, bool contig)
{
    T** out;
    if(contig)
    {
        T* fullcontig = (T*)malloc(sizeof(T)*size_m*size_n);
        out = (T**)malloc(sizeof(T*)*size_m);
        for(int i=0; i < size_m; i++)
        {
            out[i] = &(fullcontig[i*size_n]);
        }
    }
    else
    {
        out = (T**)malloc(sizeof(T*)*size_m);
        for(int i=0; i < size_m; i++)
        {
            out[i] = alloc_1darr<T>(size_n);
        }
    }
    return out;
}

template <class T>
T*** alloc_3darr(int size_m, int size_n, int size_p, bool contig)
{
    T*** out;
    if(contig)
    {
        T* fullcontig = (T*)malloc(sizeof(T)*size_m*size_n*size_p);
        out = (T***)malloc(sizeof(T*)*size_m);
        for(int j=0; j < size_n; j++)
        {
            T** out1in = (T**)malloc(sizeof(T*)*size_n);
            for(int i=0; i < size_m; i++)
            {
                out1in[i] = &(fullcontig[i*size_n+j*size_n*size_p]);
            }
            out[j] = out1in;
        }
    }
    else
    {
        out = (T***)malloc(sizeof(T**)*size_m);
        for(int i=0; i < size_m; i++)
        {
            out[i] = alloc_2darr<T>(size_n, size_p, contig);
        }
    }
    return out;
}

template <class T>
void dealloc_1darr(T* arr)
{
    free(arr);
}

template <class T>
void dealloc_2darr(int size_m, T** arr, bool contig)
{
    if(!contig)
    {
        for(int i = 0; i < size_m; i++)
        {
            dealloc_1darr(arr[i]);
        }
    }
    else
    {
        free(arr[0]);
    }
    free(arr);
}

template <class T>
void dealloc_3darr(int size_m, int size_n, T*** arr, bool contig)
{
    if(!contig)
    {
        for(int i = 0; i < size_m; i++)
        {
            dealloc_2darr(size_n, arr[i]);
        }
    }
    else
    {
        free(arr[0][0]);
        for(int i = 0; i < size_m; i++)
        {
            free(arr[i]);
        }
    }
    free(arr);
}

template <class T>
T* deep_copy_1darr(int size_m, T* arr)
{
    T* out = alloc_1darr<T>(size_m);
    for(int i = 0; i < size_m; i++)
    {
        out[i] = arr[i];
    }
    return out;
}

template <class T>
T** deep_copy_2darr(int size_m, int size_n, T** arr, bool contig)
{
    T** out = alloc_2darr<T>(size_m, size_n, contig);
    for(int i = 0; i < size_m; i++)
    {
        for(int j = 0; j < size_n; j++)
        {
            out[i][j] = arr[i][j];
        }
    }
    return out;
}

template <class T>
T*** deep_copy_3darr(int size_m, int size_n, int size_p, T*** arr, bool contig)
{
    T*** out = alloc_3darr<T>(size_m, size_n, size_p, contig);
    for(int i = 0; i < size_m; i++)
    {
        for(int j = 0; j < size_n; j++)
        {
            for(int k = 0; k < size_p; k++)
            {
                out[i][j][k] = arr[i][j][k];
            }
        }
    }
    return out;
}

template double* alloc_1darr<double>(int size_m);
template double** alloc_2darr<double>(int size_m, int size_n, bool contig);
template double*** alloc_3darr<double>(int size_m, int size_n, int size_p, bool contig);
template void dealloc_1darr<double>(double* arr);
template void dealloc_2darr<double>(int size_m, double** arr, bool contig);
template void dealloc_3darr<double>(int size_m, int size_n, double*** arr, bool contig);
template double* deep_copy_1darr<double>(int size_m, double* arr);
template double** deep_copy_2darr<double>(int size_m, int size_n, double** arr, bool contig);
template double*** deep_copy_3darr<double>(int size_m, int size_n, int size_p, double*** arr, bool contig);

template int* alloc_1darr<int>(int size_m);
template int** alloc_2darr<int>(int size_m, int size_n, bool contig);
template int*** alloc_3darr<int>(int size_m, int size_n, int size_p, bool contig);
template void dealloc_1darr<int>(int* arr);
template void dealloc_2darr<int>(int size_m, int** arr, bool contig);
template void dealloc_3darr<int>(int size_m, int size_n, int*** arr, bool contig);
template int* deep_copy_1darr<int>(int size_m, int* arr);
template int** deep_copy_2darr<int>(int size_m, int size_n, int** arr, bool contig);
template int*** deep_copy_3darr<int>(int size_m, int size_n, int size_p, int*** arr, bool contig);

template bool* alloc_1darr<bool>(int size_m);
template bool** alloc_2darr<bool>(int size_m, int size_n, bool contig);
template bool*** alloc_3darr<bool>(int size_m, int size_n, int size_p, bool contig);
template void dealloc_1darr<bool>(bool* arr);
template void dealloc_2darr<bool>(int size_m, bool** arr, bool contig);
template void dealloc_3darr<bool>(int size_m, int size_n, bool*** arr, bool contig);
template bool* deep_copy_1darr<bool>(int size_m, bool* arr);
template bool** deep_copy_2darr<bool>(int size_m, int size_n, bool** arr, bool contig);
template bool*** deep_copy_3darr<bool>(int size_m, int size_n, int size_p, bool*** arr, bool contig);
