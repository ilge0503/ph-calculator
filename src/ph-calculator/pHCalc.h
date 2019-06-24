/*
version : v1.1.7-alpha

MIT License

Copyright (c) 2019 nulLeeKH <i_am@nulleekh.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//Import libraries required for calculation.

#ifndef PH_CALCULATOR_PHCALC_H
#define PH_CALCULATOR_PHCALC_H

#define MAX_NUMBER_LENGTH 255
#define NUMBER_OF_SOLUTE 255
#define NUMBER_OF_DATA 6
#define MAX_DATA_LENGTH 32
#define MAX_SOLUTION_NUMBER 255
#define PH_CALCULATOR_STARTPOINT -3
#define PH_CALCULATOR_ENDPOINT 18
#define PH_CALCULATOR_RESOLUTION 0.001
#define PH_CALCULATOR_INITIALINTERVAL 1
#define RECIPE_FINDER_RESOLUTION 0.0001
//Specify settings used for calculations and operations

#define HFactor(soluteDataBase, index) strtol(soluteDataBase[index][1], NULL, 10)
#define NumIonization(soluteDataBase, index) strtol(soluteDataBase[index][2], NULL, 10)
#define IonizationFactor(soluteDataBase, index, time) -strtod(soluteDataBase[index][2+time], NULL)
//Declare macros approaching the database.

long SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute, char const* name);
double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nAcid, long nBase, long nRest, double const* sAcid, double const* sBase, double const* sRest, long const* iRest, double v);
double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v);
double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v);
double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v);
double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v);
double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, long nRest, double const* sRest, long const* iRest, double pH, double v);
double CalculatePH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nAcid, long nBase, long nRest, double const* sAcid, double const* sBase, double const* sRest, long const* iRest, double v, double nStart, double nEnd, double RESOLUTION, double interval);
//Declare the function prototypes

#endif //PH_CALCULATOR_PHCALC_H