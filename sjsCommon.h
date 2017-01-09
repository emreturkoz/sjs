#ifndef _SJSCOMMON_H_
#define _SJSCOMMON_H_

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "omp.h"
#include "time.h"

#define PI 3.14159265359	

/*
This is the common file that includes functions to be used by all the files of the program

*/



template <typename T> 
T ****AllocateDynamic4DArray( int nRows, int nCols, int nSlice, int kSlice)
{
      T ****dynamicArray;

      dynamicArray = new T***[nRows];
      for( int i = 0 ; i < nRows ; i++ ){
		dynamicArray[i] = new T**[nCols];
		for ( int j=0; j<nCols;j++){
			dynamicArray[i][j] = new T*[nSlice];
			for (int k=0; k<nSlice; k++){
				dynamicArray[i][j][k] = new T[kSlice];
				for(int l=0;l<kSlice;l++){
					dynamicArray[i][j][k][l] = 0;
				}
			}
		}
	  }
      return dynamicArray;
}
template <typename T>
T ***AllocateDynamic3DArray(int nRows, int nCols, int nSlice){
      T ***dynamicArray;

      dynamicArray = new T**[nRows];
      for( int i = 0 ; i < nRows ; i++ ){
		dynamicArray[i] = new T*[nCols];
		for ( int j=0; j<nCols;j++){
			dynamicArray[i][j] = new T[nSlice];
			for (int k=0; k<nSlice; k++){
					dynamicArray[i][j][k]= 0;
				}
			}
		}
      return dynamicArray;
}
template <typename T>
T **AllocateDynamic2DArray(int nRows, int nCols){
      T **dynamicArray;

      dynamicArray = new T*[nRows];
      for( int i = 0 ; i < nRows ; i++ ){
		dynamicArray[i] = new T[nCols];
		for ( int j=0; j<nCols;j++){
				dynamicArray[i][j]= 0;
			}
		}
      return dynamicArray;
}
template <typename T>
T *AllocateDynamicVector(int nRows){
      T *dynamicArray;

      dynamicArray = new T[nRows];
      for( int i = 0 ; i < nRows ; i++ ){
			dynamicArray[i]= 0;
		}
      return dynamicArray;
}

#endif // _SJSCOMMON_H_
