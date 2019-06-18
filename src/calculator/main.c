/*
version : v1.1.1-alpha

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

int SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nSolute, char* name);
double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nOther, double* sAcidH, double* sBaseH, double sOther[], int iOther[], double v);
double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v);
double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double pH, double v);
double FindAnswer(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double v, double nStart, double nEnd, double precision, double interval);
int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute);
int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute);

int main() {
    int i, j, k, nSolute;
    char tmp, soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH] = { 0, };

    printf("INFO : Loading solution data\n");

    FILE *file;
    file = fopen("../solute.pcd", "r");
    for (i=0; i<NUMBER_OF_SOLUTE; i++) {
        for (j = 0; j < NUMBER_OF_DATA; j++){
            for (k = 0; k < MAX_DATA_LENGTH; k++) {
                tmp = fgetc(file);
                if (tmp == 0x7C) break;
                else if (tmp == 0x20) break;
                else if (tmp == EOF) break;
                else soluteDataBase[i][j][k] = tmp;
            }
            if (tmp == 0x20) break;
            else if (tmp == EOF) break;
        }
        if (tmp == EOF) break;
    }
    nSolute = i;
    fclose(file);


    printf("INFO : solution data Loaded successfully\n");

    printf("Select calculation type you want [pH-calculator : 0 | graph-generator : 1] : ");scanf("%c", &tmp); printf("\n"); fflush(stdin);
    switch (tmp) {
        case '0' :
            return PhCalculator(soluteDataBase, nSolute);
        case '1' :
            return GraphGenerator(soluteDataBase, nSolute);
        default :
            printf("ERR : invalid calculation type");
            return -1;
    }
}

int SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nSolute, char* name) {
    int i;
    for (i=0; i<nSolute; i++) {
        if (strncmp(&soluteDataBase[i][0][0], &name[0], MAX_DATA_LENGTH) == 0) {
            return i;
        }
    }
    printf("ERR : invalid solution name");
    return -1;
}

int HFactor(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index) {
    return atoi(soluteDataBase[index][1]);
}

int NumIonization(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index) {
    return atoi(soluteDataBase[index][2]);
}

double IonizationFactor(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int index, int time) {
    return -atof(soluteDataBase[index][2+time]);
}

double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nOther, double* sAcidH, double* sBaseH, double sOther[], int iOther[], double v) {
    int i;
    double h = 0;

    for (i=0; i<nAcid; i++) h+=sAcidH[i];
    for (i=0; i<nBase; i++) h-=sBaseH[i];
    for (i=0; i<nOther; i++) h-=sOther[i]*HFactor(soluteDataBase, iOther[i]);
    h = h/v;

    return h;
}

double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double result;
    switch (NumIonization(soluteDataBase, iOther)) {
        case 1 :
            result = CalculateMono(soluteDataBase, sOther, iOther, pH, v);
            break;
        case 2 :
            result = CalculateDi(soluteDataBase, sOther, iOther, pH, v);
            break;
        case 3 :
            result = CalculateTri(soluteDataBase, sOther, iOther, pH, v);
            break;
        default :
            printf("ERR : error occured while calculating %d other solution", iOther);
            result = 0;
            break;
    }
    return result;
}

double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cen, x1, k1, result;
    cen = sOther/v;
    x1 = pH;
    k1 = pow(10, IonizationFactor(soluteDataBase, iOther, 1));
    result = (cen*k1)/(x1+k1);
    return result;
}

double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cen, x1, x2, k1, k2, result;
    cen = sOther/v;
    x1 = pH;
    x2 = pow(pH, 2);
    k1 = pow(10, IonizationFactor(soluteDataBase, iOther, 1));
    k2 = pow(10, IonizationFactor(soluteDataBase, iOther, 2));
    result = (cen*((k1*x1)+(2*k1*k2)))/(x2+(k1*x1)+(k1*k2));
    return result;
}

double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sOther, int iOther, double pH, double v) {
    double cen, x1, x2, x3, k1, k2, k3, result;
    cen = sOther/v;
    x1 = pH;
    x2 = pow(pH, 2);
    x3 = pow(pH, 3);
    k1 = pow(10, IonizationFactor(soluteDataBase, iOther, 1));
    k2 = pow(10, IonizationFactor(soluteDataBase, iOther, 2));
    k3 = pow(10, IonizationFactor(soluteDataBase, iOther, 2));
    result = (cen*((k1*x2)+(2*k1*k2*x1)+(3*k1*k2*k3)))/(x3+(k1*x2)+(k1*k2*x1)+(k1*k2*k3));
    return result;
}

double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double pH, double v) {
    int i;
    double result=h+(pow(10.0, -14)/pH);
    for (i=0;i<nOther;i++) {
        result += CalculatePolyproticAcid(soluteDataBase, sOther[i], iOther[i], pH, v);
    }
    return fabs(pH-result);
}

double FindAnswer(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nOther, double sOther[], int iOther[], double v, double nStart, double nEnd, double precision, double interval) {
    double pH, tmp, lowest, ans;

    lowest = CalculateError(soluteDataBase, h, nOther, sOther, iOther, pow(10, -pH), v);

    for (pH=nStart+interval;pH<nEnd;pH+=interval) {
        tmp = CalculateError(soluteDataBase, h, nOther, sOther, iOther, pow(10, -pH), v);
        if (tmp == -1) return -1;
        else if (tmp<lowest) {
            lowest = tmp;
            ans = pH;
        }
    }

    if (interval>precision) return FindAnswer(soluteDataBase, h, nOther, sOther, iOther, v, ans-interval, ans+interval, precision, interval*0.1);
    else return ans;
}

int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute) {
    int i;
    int nAcid, nBase, nOther, iOther[MAX_SOLUTION_NUMBER];
    double vAll, h, vol, cen, sAcidH[MAX_SOLUTION_NUMBER], sBaseH[MAX_SOLUTION_NUMBER], sOther[MAX_SOLUTION_NUMBER];

    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };

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

    printf("Volume of pure water in liter : "); scanf("%s", tmp); vAll=atof(tmp); printf("\n");
    //Get total volume(in liter)

    printf("\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strong monoprotic acid in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of strong monoprotic acid : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sAcidH[i]=cen*vol;
    }
    for(i=0;i<nBase;i++) {
        printf("Volume of strong base in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of strong base : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sBaseH[i]=cen*vol;
    }

    for (i=0;i<nOther;i++) {
        printf("name of other solution : "); scanf("%s", name); iOther[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n");
        if (iOther[i] == -1) return -1;
        printf("Volume of other solution in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of other solution : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sOther[i]=vol*cen;
    }
    //initial concentration&acid constant
    //input ends here

    h = CalcInitialH(soluteDataBase, nAcid, nBase, nOther, sAcidH, sBaseH, sOther, iOther, vAll);

    printf("%.3lf", FindAnswer(soluteDataBase, h, nOther, sOther, iOther, vAll, -3, 18, CALCULATOR_PRECISION, 1));

    return 0;
}

int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute) {
    int i;
    int nAcid, nBase, nOther, iOther[MAX_SOLUTION_NUMBER];
    double vAll, h, vol, cen, sAcidH[MAX_SOLUTION_NUMBER], sBaseH[MAX_SOLUTION_NUMBER], sOther[MAX_SOLUTION_NUMBER];
    double volAdding, volPerTime, volAdded, cenAdded;

    char result[6];

    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };

    printf("Number of target strong monoprotic acids : "); scanf("%s", tmp); nAcid=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : invalid nAcid numbers");
        return 1;
    }
    printf("Number of target strong bases : "); scanf("%s", tmp); nBase=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : invalid nBase numbers");
        return 1;
    }
    printf("Number of the target rests : "); scanf("%s", tmp); nOther=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nOther || nOther>255) {
        printf("ERR : invalid oSolution numbers");
        return 1;
    }
    //Get number of target acids and bases

    printf("\n");

    printf("Volume of pure water in liter : "); scanf("%s", tmp); vAll=atof(tmp); printf("\n");
    //Get total volume(in liter)

    printf("\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of target strong monoprotic acid in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of target strong monoprotic acid : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sAcidH[i]=cen*vol;
    }
    for(i=0;i<nBase;i++) {
        printf("Volume of target strong base in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of target strong base : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sBaseH[i]=cen*vol;
    }

    for (i=0;i<nOther;i++) {
        printf("name of other target solution : "); scanf("%s", name); iOther[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n");
        if (iOther[i] == -1) return -1;
        printf("Volume of other target solution in liter : "); scanf("%s", tmp); vol=atof(tmp); printf("\n"); vAll=vAll+vol;
        printf("Concentration of other target solution : "); scanf("%s", tmp); cen=atof(tmp); printf("\n"); sOther[i]=vol*cen;
    }
    //initial concentration&acid constant
    //input ends here

    printf("\n");
    printf("\n");

    printf("Type of solution to be added [strong monoprotic acid : 0 | strong base : 1 | the rests : 2] : "); scanf("%s", tmp); printf("\n"); fflush(stdin);
    FILE* file;
    file = fopen("../result.pcd", "w");
    switch (tmp[0]) {
        case '0' :
            nAcid++;
            printf("Volume of adding strong monoprotic acid to be added in liter : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n");
            printf("Volume of adding strong monoprotic acid to be added at a time in liter : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n");
            printf("Concentration of adding strong monoprotic acid : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n");

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sAcidH[nAcid-1]=cenAdded*volPerTime*i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nOther, sAcidH, sBaseH, sOther, iOther, vAll+volPerTime*i);
                gcvt(FindAnswer(soluteDataBase, h, nOther, sOther, iOther, vAll+volPerTime*i, -3, 18, CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }
            break;
        case '1' :
            nBase++;
            printf("Volume of adding strong base to be added in liter : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n");
            printf("Volume of adding strong base to be added at a time in liter : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n");
            printf("Concentration of adding strong base : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n");

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sBaseH[nBase - 1] = cenAdded * volPerTime * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nOther, sAcidH, sBaseH, sOther, iOther, vAll + volPerTime * i);
                gcvt(FindAnswer(soluteDataBase, h, nOther, sOther, iOther, vAll + volPerTime * i, -3, 18, CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }
            break;
        case '2' :
            nOther++;
            printf("name of other adding solution : "); scanf("%s", name); iOther[nOther-1] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n");
            printf("Volume of other adding solution to be added in liter : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n");
            printf("Volume of other adding solution to be added at a time in liter : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n");
            printf("Concentration of other adding solution : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n");

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sOther[nBase - 1] = cenAdded * volPerTime * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nOther, sAcidH, sBaseH, sOther, iOther, vAll + volPerTime * i);
                gcvt(FindAnswer(soluteDataBase, h, nOther, sOther, iOther, vAll + volPerTime * i, -3, 18, CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }
            break;
        default :
            printf("ERR : invalid solution type");
            return 1;
    }

    fclose(file);

    printf("Graph has been generated.");

    return 0;
}