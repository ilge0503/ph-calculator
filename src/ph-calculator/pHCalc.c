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

#include "pHCalc.h"

long SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute, char const* name) {
    long i;
    //Declare the variables needed to operate the function.

    for (i=0; i<nSolute; i++) {
        if (strncmp(&soluteDataBase[i][0][0], &name[0], MAX_DATA_LENGTH) == 0) {
            return i;
        }
    }   //Search the database and return the material's unique number.

    printf("ERR : Invalid solution name\n");
    return -1;
}

double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nAcid, long nBase, long nRest, double const* sAcid, double const* sBase, double const* sRest, long const* iRest, double v) {
    long i;
    double h = 0;
    //Declare the variables needed to operate the function.

    for (i=0; i<nAcid; i++) h+=sAcid[i];
    for (i=0; i<nBase; i++) h-=sBase[i];
    for (i=0; i<nRest; i++) h-=sRest[i]*HFactor(soluteDataBase, iRest[i]);
    h = h/v;
    //Calculate the initial hydrogen ion concentration.

    return h;   //Return the result.
}

double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v) {
    double result;
    //Declare the variables needed to operate the function.

    switch (NumIonization(soluteDataBase, iRest)) {
        case 1 :
            result = CalculateMono(soluteDataBase, sRest, iRest, pH, v);
            break;

        case 2 :
            result = CalculateDi(soluteDataBase, sRest, iRest, pH, v);
            break;

        case 3 :
            result = CalculateTri(soluteDataBase, sRest, iRest, pH, v);
            break;

        default :
            printf("ERR : Error occured while calculating %ld the solution\n", iRest);
            result = 0;
    }   //Classify polyprotic acid and call up the corresponding calculation function.

    return result;  //Return the result.
}

double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v) {
    double cen, x1, k1, result;
    //Declare the variables needed to operate the function.

    cen = sRest/v;
    x1 = pH;
    k1 = pow(10, IonizationFactor(soluteDataBase, iRest, 1));
    result = (cen*k1)/(x1+k1);
    //Calculate the formula for monoprotic acid.

    return result;  //Return the result.
}

double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v) {
    double cen, x1, x2, k1, k2, result;
    //Declare the variables needed to operate the function.

    cen = sRest/v;
    x1 = pH;
    x2 = pow(pH, 2);
    k1 = pow(10, IonizationFactor(soluteDataBase, iRest, 1));
    k2 = pow(10, IonizationFactor(soluteDataBase, iRest, 2));
    result = (cen*((k1*x1)+(2*k1*k2)))/(x2+(k1*x1)+(k1*k2));
    //Calculate the formula for diprotic acid.

    return result;  //Return the result.
}

double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, long iRest, double pH, double v) {
    double cen, x1, x2, x3, k1, k2, k3, result;
    //Declare the variables needed to operate the function.

    cen = sRest/v;
    x1 = pH;
    x2 = pow(pH, 2);
    x3 = pow(pH, 3);
    k1 = pow(10, IonizationFactor(soluteDataBase, iRest, 1));
    k2 = pow(10, IonizationFactor(soluteDataBase, iRest, 2));
    k3 = pow(10, IonizationFactor(soluteDataBase, iRest, 3));
    result = (cen*((k1*x2)+(2*k1*k2*x1)+(3*k1*k2*k3)))/(x3+(k1*x2)+(k1*k2*x1)+(k1*k2*k3));
    //Calculate the formula for triprotic acid.

    return result;  //Return the result.
}

double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, long nRest, double const* sRest, long const* iRest, double pH, double v) {
    long i;
    double result=h+(pow(10.0, -14)/pH);
    //Declare the variables needed to operate the function.

    for (i=0;i<nRest;i++) {
        result += CalculatePolyproticAcid(soluteDataBase, sRest[i], iRest[i], pH, v);
    }   //Calculate the result using the value received.

    return fabs(pH-result); //Return the error between value received and the result.
}

double CalculatePH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nAcid, long nBase, long nRest, double const* sAcid, double const* sBase, double const* sRest, long const* iRest, double v, double nStart, double nEnd, double resolution, double interval) {
    double pH = 0, tmp, lowest, ans = -1, h;
    //Declare the variables needed to operate the function.

    h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, v);

    lowest = CalculateError(soluteDataBase, h, nRest, sRest, iRest, pow(10, -nStart), v);

    while (1) {
        pH += interval;

        tmp = CalculateError(soluteDataBase, h, nRest, sRest, iRest, pow(10, -pH), v);

        if (tmp == -1) return -1;
        else if (tmp<lowest) {
            lowest = tmp;
            ans = pH;
        }

        if (pH >= nEnd) break;
    }   //Locate the point where the error is minimized within the search scope.

    if (interval>resolution) return CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, v, ans-interval, ans+interval, resolution, interval*0.1);
        //If the search interval is greater than the allowable error, narrow down the search interval to re-discover both sides of the point.
    else return ans;
    //If the search interval is less than or equal to an acceptable error, return the result.
}