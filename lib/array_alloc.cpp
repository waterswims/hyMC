#include "../include/array_alloc.hpp"
#include <cstdlib>

template <class T>
T* arr::alloc_1darr(int size_m)
{
    return (T*)malloc(sizeof(T)*size_m);
}

template <class T>
T** arr::alloc_2darr(int size_m, int size_n, bool contig)
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
T*** arr::alloc_3darr(int size_m, int size_n, int size_p, bool contig)
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
void arr::dealloc_1darr(T* arr)
{
    free(arr);
}

template <class T>
void arr::dealloc_2darr(int size_m, T** arr, bool contig)
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
void arr::dealloc_3darr(int size_m, int size_n, T*** arr, bool contig)
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
T* arr::deep_copy_1darr(int size_m, const T* arr)
{
    T* out = alloc_1darr<T>(size_m);
    for(int i = 0; i < size_m; i++)
    {
        out[i] = arr[i];
    }
    return out;
}

template <class T>
T** arr::deep_copy_2darr(int size_m, int size_n, const T** arr, bool contig)
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
T*** arr::deep_copy_3darr(int size_m, int size_n, int size_p, const T*** arr, bool contig)
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

template double* arr::alloc_1darr<double>(int size_m);
template double** arr::alloc_2darr<double>(int size_m, int size_n, bool contig);
template double*** arr::alloc_3darr<double>(int size_m, int size_n, int size_p, bool contig);
template void arr::dealloc_1darr<double>(double* arr);
template void arr::dealloc_2darr<double>(int size_m, double** arr, bool contig);
template void arr::dealloc_3darr<double>(int size_m, int size_n, double*** arr, bool contig);
template double* arr::deep_copy_1darr<double>(int size_m, const double* arr);
template double** arr::deep_copy_2darr<double>(int size_m, int size_n, const double** arr, bool contig);
template double*** arr::deep_copy_3darr<double>(int size_m, int size_n, int size_p, const double*** arr, bool contig);

template int* arr::alloc_1darr<int>(int size_m);
template int** arr::alloc_2darr<int>(int size_m, int size_n, bool contig);
template int*** arr::alloc_3darr<int>(int size_m, int size_n, int size_p, bool contig);
template void arr::dealloc_1darr<int>(int* arr);
template void arr::dealloc_2darr<int>(int size_m, int** arr, bool contig);
template void arr::dealloc_3darr<int>(int size_m, int size_n, int*** arr, bool contig);
template int* arr::deep_copy_1darr<int>(int size_m, const int* arr);
template int** arr::deep_copy_2darr<int>(int size_m, int size_n, const int** arr, bool contig);
template int*** arr::deep_copy_3darr<int>(int size_m, int size_n, int size_p, const int*** arr, bool contig);

template bool* arr::alloc_1darr<bool>(int size_m);
template bool** arr::alloc_2darr<bool>(int size_m, int size_n, bool contig);
template bool*** arr::alloc_3darr<bool>(int size_m, int size_n, int size_p, bool contig);
template void arr::dealloc_1darr<bool>(bool* arr);
template void arr::dealloc_2darr<bool>(int size_m, bool** arr, bool contig);
template void arr::dealloc_3darr<bool>(int size_m, int size_n, bool*** arr, bool contig);
template bool* arr::deep_copy_1darr<bool>(int size_m, const bool* arr);
template bool** arr::deep_copy_2darr<bool>(int size_m, int size_n, const bool** arr, bool contig);
template bool*** arr::deep_copy_3darr<bool>(int size_m, int size_n, int size_p, const bool*** arr, bool contig);
