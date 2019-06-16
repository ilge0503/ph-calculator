/*
version : v1.1.0-alpha

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

#define MAX_NUMBER_LENGTH 255
#define NUMBER_OF_SOLUTE 255
#define NUMBER_OF_DATA 6
#define MAX_DATA_LENGTH 32
#define MAX_SOLUTION_NUMBER 255
#define CALCULATOR_PRECISION 0.001

int SpecifySolute(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], char* name, int nSolute);
double CalcInitialH(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nOther, double* sAcidH, double* sBaseH, double sOther[], int iOther[], double v);
double CalculatePolyproticAcid(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateMono(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateDi(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateTri(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateError(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double pH, double v);
double FindAnswer(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double v, double nStart, double nEnd, double precision, double interval);

int main() {
    int i, j, k;
    int nAcid, nBase, nOther, iOther[MAX_SOLUTION_NUMBER], nSolute;
    double v, h, vol, cncnt, sAcidH[MAX_SOLUTION_NUMBER], sBaseH[MAX_SOLUTION_NUMBER], sOther[MAX_SOLUTION_NUMBER];

    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };

    char temp, soluteList[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH] = { 0, };

    FILE *file;
    file = fopen("../solute.pcd", "r");
    for (i=0; i<NUMBER_OF_SOLUTE; i++) {
        for (j = 0; j < NUMBER_OF_DATA; j++){
            for (k = 0; k < MAX_DATA_LENGTH; k++) {
                temp = fgetc(file);
                if (temp == 0x7C) break;
                else if (temp == 0x20) break;
                else if (temp == EOF) break;
                else soluteList[i][j][k] = temp;
            }
            if (temp == 0x20) break;
            else if (temp == EOF) break;
        }
        if (temp == EOF) break;
    }
    nSolute = i;
    fclose(file);

    printf("Number of strong monoprotic acids : "); scanf("%s", tmp); nAcid=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : invalid nAcid numbers");
        return 1;
    }
    printf("Number of strong bases : "); scanf("%s", tmp); nBase=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : invalid nBase numbers");
        return 1;
    }
    printf("Number of the rests : "); scanf("%s", tmp); nOther=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nOther || nOther>255) {
        printf("ERR : invalid oSolution numbers");
        return 1;
    }
    //Get number of acids and bases

    printf("\n");

    v=0;

    printf("Volume of pure water in liter : "); scanf("%s", tmp); v=atof(tmp); printf("\n");
    //Get total volume(in liter)

    printf("\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strong monoprotic acid in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); v=v+vol;
        printf("Concentration of strong monoprotic acid : "); scanf("%s", tmp); cncnt=atof(tmp); printf("\n"); sAcidH[i]=cncnt*vol;
    }
    for(i=0;i<nBase;i++) {
        printf("Volume of strong base in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); v=v+vol;
        printf("Concentration of strong base : "); scanf("%s", tmp); cncnt=atof(tmp); printf("\n"); sBaseH[i]=cncnt*vol;
    }

    for (i=0;i<nOther;i++) {
        printf("name of other solution : "); scanf("%s", name); iOther[i] = SpecifySolute(soluteList, name, nSolute); printf("\n");
        printf("Volume of another solution in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); v=v+vol;
        printf("Concentration of another solution : "); scanf("%s", tmp); cncnt=atof(tmp); printf("\n"); sOther[i]=vol*cncnt;
    }
    //initial concncnttration&acid constant
    //input ends here

    h = CalcInitialH(soluteList, nAcid, nBase, nOther, sAcidH, sBaseH, sOther, iOther, v);

    printf("%.3lf", FindAnswer(soluteList, h, nOther, sOther, iOther, v, -3, 18, CALCULATOR_PRECISION, 1));

    return 0;
}

int SpecifySolute(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], char* name, int nSolute) {
    int i;
    for (i=0; i<NUMBER_OF_SOLUTE; i++) {
        if (strncmp(&list[i][0][0], &name[0], MAX_DATA_LENGTH) == 0) {
            return i;
        }
    }
    return -1;
}

int HFactor(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index) {
    return atoi(list[index][1]);
}

int NumIonization(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index) {
    return atoi(list[index][2]);
}

double IonizationFactor(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index, int time) {
    return -atof(list[index][2+time]);
}

double CalcInitialH(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nOther, double* sAcidH, double* sBaseH, double sOther[], int iOther[], double v) {
    int i;
    double h = 0;

    for (i=0; i<nAcid; i++) h+=sAcidH[i];
    for (i=0; i<nBase; i++) h-=sBaseH[i];
    for (i=0; i<nOther; i++) h-=sOther[i]*HFactor(list, iOther[i]);
    h = h/v;

    return h;
}

double CalculatePolyproticAcid(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double result;
    switch (NumIonization(list, iOther)) {
        case 1 :
            result = CalculateMono(list, sOther, iOther, pH, v);
            break;
        case 2 :
            result = CalculateDi(list, sOther, iOther, pH, v);
            break;
        case 3 :
            result = CalculateTri(list, sOther, iOther, pH, v);
            break;
        default :
            printf("ERR : error occured while calculating %d other solution", iOther);
            result = 0;
            break;
    }
    return result;
}

double CalculateMono(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cncnt, x1, k1, result;
    cncnt = sOther/v;
    x1 = pH;
    k1 = pow(10, IonizationFactor(list, iOther, 1));
    result = (cncnt*k1)/(x1+k1);
    return result;
}

double CalculateDi(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cncnt, x1, x2, k1, k2, result;
    cncnt = sOther/v;
    x1 = pH;
    x2 = pow(pH, 2);
    k1 = pow(10, IonizationFactor(list, iOther, 1));
    k2 = pow(10, IonizationFactor(list, iOther, 2));
    result = (cncnt*((k1*x1)+(2*k1*k2)))/(x2+(k1*x1)+(k1*k2));
    return result;
}

double CalculateTri(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cncnt, x1, x2, x3, k1, k2, k3, result;
    cncnt = sOther/v;
    x1 = pH;
    x2 = pow(pH, 2);
    x3 = pow(pH, 3);
    k1 = pow(10, IonizationFactor(list, iOther, 1));
    k2 = pow(10, IonizationFactor(list, iOther, 2));
    k3 = pow(10, IonizationFactor(list, iOther, 2));
    result = (cncnt*((k1*x2)+(2*k1*k2*x1)+(3*k1*k2*k3)))/(x3+(k1*x2)+(k1*k2*x1)+(k1*k2*k3));
    return result;
}

double CalculateError(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double pH, double v) {
    int i;
    double result=h+(pow(10.0, -14)/pH);
    for (i=0;i<nOther;i++) {
        result += CalculatePolyproticAcid(list, sOther[i], iOther[i], pH, v);
    }
    return fabs(pH-result);
}

double FindAnswer(char list[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double v, double nStart, double nEnd, double precision, double interval) {
    double pH, tmp, lowest, ans;

    lowest = CalculateError(list, h, nOther, sOther, iOther, pow(10, -pH), v);

    for (pH=nStart+interval;pH<nEnd;pH+=interval) {
        tmp = CalculateError(list, h, nOther, sOther, iOther, pow(10, -pH), v);
        if (tmp == -1) return -1;
        else if (tmp<lowest) {
            lowest = tmp;
            ans = pH;
        }
    }

    if (interval>precision) return FindAnswer(list, h, nOther, sOther, iOther, v, ans-interval, ans+interval, precision, interval*0.1);
    else return ans;
}