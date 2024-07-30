#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <cassert>
#include <iostream>
// #ifdef USE_PDQSORT

// #include "pdqsort.h"
// #endif

/*
#define DATA_T int64_t
#define DATA_PRINTF "%ld"
#define RAND_EXPR (rand() % 256 - 128)
*/
#define DATA_T double
#define DATA_PRINTF "%g"
#define RAND_EXPR ((double)rand() / RAND_MAX - 0.5)

// sort algorithm
void swap(DATA_T *xp, DATA_T *yp)
{
    DATA_T temp = *xp;
    *xp = *yp;
    *yp = temp;
}

// Partition function for QuickSort
uint64_t partition(DATA_T *array, uint64_t low, uint64_t high)
{
    DATA_T pivot = array[low];
    uint64_t i = low;
    uint64_t j = high;

    while (i < j)
    {
        while (array[i] <= pivot && i < high)
        {
            i++;
        }
        while (array[j] > pivot && j > low)
        {
            j--;
        }
        if (i < j)
        {
            swap(&array[i], &array[j]);
        }
    }
    swap(&array[low], &array[j]);
    return j;
}

void quickSort(DATA_T *array, uint64_t low, uint64_t high)
{
    if (low < high)
    {
        uint64_t partitionIndex = partition(array, low, high);
        if (partitionIndex > 0)
        {
            quickSort(array, low, partitionIndex - 1); // Only call with partitionIndex > 0
        }
        quickSort(array, partitionIndex + 1, high);
    }
}

void quickSortWrapper(DATA_T *array, uint64_t length)
{
    quickSort(array, 0, length - 1);
}

void merge(DATA_T *array, int left, int mid, int right)
{
    int n1 = mid - left + 1;
    int n2 = right - mid;

    // Dynamically allocate temporary arrays
    DATA_T *leftArray = new DATA_T[n1];
    DATA_T *rightArray = new DATA_T[n2];

    // Copy data to temporary arrays
    for (int i = 0; i < n1; i++)
        leftArray[i] = array[left + i];
    for (int j = 0; j < n2; j++)
        rightArray[j] = array[mid + 1 + j];

    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2)
    {
        if (leftArray[i] <= rightArray[j])
        {
            array[k] = leftArray[i];
            i++;
        }
        else
        {
            array[k] = rightArray[j];
            j++;
        }
        k++;
    }

    // Copy remaining elements of leftArray, if any
    while (i < n1)
    {
        array[k] = leftArray[i];
        i++;
        k++;
    }

    // Copy remaining elements of rightArray, if any
    while (j < n2)
    {
        array[k] = rightArray[j];
        j++;
        k++;
    }

    delete[] leftArray;
    delete[] rightArray;
}

void mergeSort(DATA_T *array, int begin, int end)
{
    if (begin >= end)
        return;

    int mid = begin + (end - begin) / 2;
    mergeSort(array, begin, mid);
    mergeSort(array, mid + 1, end);
    merge(array, begin, mid, end);
}

void mergeSortWrapper(DATA_T *array, uint64_t length)
{
    mergeSort(array, 0, length - 1);
}

void insertionSort(DATA_T *array, uint64_t size)
{

    DATA_T key;
    for (uint64_t i = 1; i < size; i++)
    {
        key = array[i];
        int64_t j = i - 1;
        while (j >= 0 && array[j] > key)
        {
            array[j + 1] = array[j];
            j--;
        }
        array[j + 1] = key;
    }
}

void bubble_sort(DATA_T *array, uint64_t length)
{
    for (uint64_t i = 0; i < length - 1; ++i)
    { // i from 0 to length-2
        for (uint64_t j = i + 1; j < length; ++j)
        { // j from i+1 to length-1
            if (array[i] > array[j])
            {
                // Swap array[i] and array[j]
                DATA_T tmp = array[j];
                array[j] = array[i];
                array[i] = tmp;
            }
        }
    }
}

void selectionSort(DATA_T *array, uint64_t n)
{
    uint64_t i, j, min_idx;
    for (i = 0; i < n - 1; i++)
    {
        min_idx = i;
        for (j = i + 1; j < n; j++)
        {
            if (array[j] < array[min_idx])
                min_idx = j;
        }
        if (min_idx != i)
            swap(&array[min_idx], &array[i]);
    }
}

// Define DATA_T as needed, here using int as an example
// code from lab 9
enum array_ordering
{
    RANDOM,
    SORTED,
    REVERSE_SORTED,
    ALMOST_SORTED
};

void gently_shuffle_array(DATA_T *array, uint64_t length)
{
    uint64_t limit = 10;
    // shuffle array indicies a little, with nearby values
    for (uint64_t i = 0; i < length - limit - 1; i++)
    {
        uint64_t offset = rand() % limit;
        DATA_T tmp = array[i + offset];
        array[i + offset] = array[i];
        array[i] = tmp;
    }
}

DATA_T *create_array(uint64_t length, array_ordering order)
{
    DATA_T *array = (DATA_T *)malloc(length * sizeof(DATA_T));
    if (array == NULL)
    {
        return NULL;
    }
    for (uint64_t i = 0; i < length; i++)
    {
        array[i] = RAND_EXPR;
    }

    // order the array as requested
    switch (order)
    {
    case RANDOM:
        break;
    case SORTED:
        std::sort(array, array + length);
        break;
    case ALMOST_SORTED:
        std::sort(array, array + length);
        gently_shuffle_array(array, length);
        break;
    case REVERSE_SORTED:
        std::sort(array, array + length);
        std::reverse(array, array + length);
        break;
    }
    return array;
}

bool is_sorted(DATA_T *array, uint64_t length)
{
    for (uint64_t i = 0; i < length - 1; i++)
    {
        if (array[i] > array[i + 1])
        {
            return false;
        }
    }
    return true;
}

void time_sort(const char *descr, void (*sort)(DATA_T *, uint64_t), uint64_t length, array_ordering order)
{
    struct timespec start, end;

    DATA_T *array = create_array(length, order);
    if (array == NULL)
    {
        printf("Couldn't allocate.\n");
        return;
    }

    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
    sort(array, length);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
    assert(is_sorted(array, length));

    free(array);
    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("%25s sorted %8lu values in %7.2f ms\n", descr, length, elapsed * 1000);
}
void time_them(uint64_t length)
{
    printf("Timing on array of size: %lu\n", length);

    // Timing different sorting algorithms on random data
    // printf("\n--- Random Data ---\n");
    // time_sort("QuickSort on random data", quickSortWrapper, length, RANDOM);
    // time_sort("MergeSort on random data", mergeSortWrapper, length, RANDOM);
    // time_sort("InsertionSort on random data", insertionSort, length, RANDOM);
    // time_sort("Bubble Sort on random data", bubble_sort, length, RANDOM);
    // time_sort("SelectionSort on random data", selectionSort, length, RANDOM);

    // // Timing different sorting algorithms on sorted data
    // printf("\n--- Sorted Data ---\n");
    // time_sort("QuickSort on sorted data", quickSortWrapper, length, SORTED);
    // time_sort("MergeSort on sorted data", mergeSortWrapper, length, SORTED);
    // time_sort("InsertionSort on sorted data", insertionSort, length, SORTED);
    // time_sort("Bubble Sort on sorted data", bubble_sort, length, SORTED);
    // time_sort("SelectionSort on sorted data", selectionSort, length, SORTED);

    // // Timing different sorting algorithms on reverse sorted data
    // printf("\n--- Reverse Sorted Data ---\n");
    // time_sort("QuickSort on reverse sorted data", quickSortWrapper, length, REVERSE_SORTED);
    // time_sort("MergeSort on reverse sorted data", mergeSortWrapper, length, REVERSE_SORTED);
    // time_sort("InsertionSort on reverse sorted data", insertionSort, length, REVERSE_SORTED);
    // time_sort("Bubble Sort on reverse sorted data", bubble_sort, length, REVERSE_SORTED);
    // time_sort("SelectionSort on reverse sorted data", selectionSort, length, REVERSE_SORTED);

    // // Timing different sorting algorithms on almost sorted data
    // printf("\n--- Almost Sorted Data ---\n");
    // time_sort("QuickSort on almost sorted data", quickSortWrapper, length, ALMOST_SORTED);
    // time_sort("MergeSort on almost sorted data", mergeSortWrapper, length, ALMOST_SORTED);
    // time_sort("InsertionSort on almost sorted data", insertionSort, length, ALMOST_SORTED);
    // time_sort("Bubble Sort on almost sorted data", bubble_sort, length, ALMOST_SORTED);
    // time_sort("SelectionSort on almost sorted data", selectionSort, length, ALMOST_SORTED);
}

void just_sort(void (*sort)(DATA_T *, uint64_t), uint64_t length, array_ordering order)
{
    DATA_T *array = create_array(length, order);
    if (array == NULL)
    {
        printf("Couldn't allocate.\n");
        return;
    }
    sort(array, length);
    assert(is_sorted(array, length));
    free(array);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        printf("Must give array size on command line.\n");
        return 1;
    }
    uint64_t length = atol(argv[1]);
    printf("Array size: %lu kB\n", length * sizeof(DATA_T) / 1024);

    // Warm up to get the CPU out of a low-power state...
    just_sort(quickSortWrapper, length, RANDOM);

    time_them(length);
    return 0;
}

// g++ -Wall -Wpedantic -march=haswell  mini.cpp && ./a.out 10000
// g++ -Wall -Wpedantic -std=c17 -march=haswell -O1 -g mini.cpp &&  valgrind --tool=cachegrind ./a.out 10000