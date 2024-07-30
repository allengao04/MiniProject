#ifndef MINI_H
#define MINI_H


void quickSort(DATA_T arr[], uint64_t low, uint64_t high);
void mergeSort(DATA_T arr[], uint64_t left, uint64_t right);
void insertionSort(DATA_T arr[], uint64_t size);
void bubble_sort(DATA_T arr[], uint64_t length);
void selectionSort(DATA_T arr[], uint64_t n);

// Include these if you haven't defined the functions before they are called
void swap(DATA_T *xp, DATA_T *yp);
uint64_t partition(DATA_T arr[], uint64_t low, uint64_t high);



#endif