/*
version : v1.1.3-alpha

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
//Import libraries required for calculation and operation.

#define MAX_NUMBER_LENGTH 255
#define NUMBER_OF_SOLUTE 255
#define NUMBER_OF_DATA 6
#define MAX_DATA_LENGTH 32
#define MAX_SOLUTION_NUMBER 255
#define PH_CALCULATOR_STARTPOINT -3
#define PH_CALCULATOR_ENDPOINT 18
#define PH_CALCULATOR_PRECISION 0.0001
#define PH_CALCULATOR_INITIALINTERVAL 1
#define NEUTRALPOINT_FINDER_PRECISION 0.0001
#define NEUTRALPOINT_FINDER_INTERVAL 0.0001
//Specify settings used for calculations and operations

#define HFactor(soluteDataBase, index) atoi(soluteDataBase[index][1])
#define NumIonization(soluteDataBase, index) atoi(soluteDataBase[index][2])
#define IonizationFactor(soluteDataBase, index, time) -atof(soluteDataBase[index][2+time])
//Declare macros approaching the database.

int SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nSolute, char* name);
double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nRest, double* sAcid, double* sBase, double* sRest, int* iRest, double v);
double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v);
double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v);
double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v);
double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v);
double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nRest, double* sRest, int* iRest, double pH, double v);
double CalculatePH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nRest, double* sRest, int* iRest, double v, double nStart, double nEnd, double precision, double interval);
int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute);
int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute);
int RecipeFinder(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute);
//Declare the function prototypes

int main() {
    printf("INFO : Program started\n");

    int i, j, k, nSolute;
    char tmp, soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

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
    //Get the database needed to operate the program.

    printf("INFO : Solution data Loaded successfully\n");

    printf("Select calculation type you want [pH-calculator : 0 | graph-generator : 1 | recipe-finder : 2] : ");scanf("%c", &tmp); printf("\n"); fflush(stdin);
    switch (tmp) {
        case '0' :
            return PhCalculator(soluteDataBase, nSolute);

        case '1' :
            return GraphGenerator(soluteDataBase, nSolute);

        case '2' :
            return RecipeFinder(soluteDataBase, nSolute);

        default :
            printf("ERR : Invalid calculation type\n");
            return -1;
    }   //Specify the function to execute, and call the function.
}

int SpecifySolute(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nSolute, char* name) {
    int i;
    //Declare the variables needed to operate the function.

    for (i=0; i<nSolute; i++) {
        if (strncmp(&soluteDataBase[i][0][0], &name[0], MAX_DATA_LENGTH) == 0) {
            return i;
        }
    }   //Search the database and return the material's unique number.

    printf("ERR : Invalid solution name\n");
    return -1;
}

double CalcInitialH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nAcid, int nBase, int nRest, double* sAcid, double* sBase, double* sRest, int* iRest, double v) {
    int i;
    double h = 0;
    //Declare the variables needed to operate the function.

    for (i=0; i<nAcid; i++) h+=sAcid[i];
    for (i=0; i<nBase; i++) h-=sBase[i];
    for (i=0; i<nRest; i++) h-=sRest[i]*HFactor(soluteDataBase, iRest[i]);
    h = h/v;
    //Calculate the initial hydrogen ion concentration.

    return h;   //Return the result.
}

double CalculatePolyproticAcid(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v) {
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
            printf("ERR : Error occured while calculating %d the solution\n", iRest);
            result = 0;
    }   //Classify polyprotic acid and call up the corresponding calculation function.

    return result;  //Return the result.
}

double CalculateMono(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v) {
    double cen, x1, k1, result;
    //Declare the variables needed to operate the function.

    cen = sRest/v;
    x1 = pH;
    k1 = pow(10, IonizationFactor(soluteDataBase, iRest, 1));
    result = (cen*k1)/(x1+k1);
    //Calculate the formula for monoprotic acid.

    return result;  //Return the result.
}

double CalculateDi(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v) {
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

double CalculateTri(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double sRest, int iRest, double pH, double v) {
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

double CalculateError(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nRest, double* sRest, int* iRest, double pH, double v) {
    int i;
    double result=h+(pow(10.0, -14)/pH);
    //Declare the variables needed to operate the function.

    for (i=0;i<nRest;i++) {
        result += CalculatePolyproticAcid(soluteDataBase, sRest[i], iRest[i], pH, v);
    }   //Calculate the result using the value received.

    return fabs(pH-result); //Return the error between value received and the result.
}

double CalculatePH(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], double h, int nRest, double* sRest, int* iRest, double v, double nStart, double nEnd, double precision, double interval) {
    double pH, tmp, lowest, ans;
    //Declare the variables needed to operate the function.

    lowest = CalculateError(soluteDataBase, h, nRest, sRest, iRest, pow(10, -pH), v);
    for (pH=nStart+interval;pH<nEnd;pH+=interval) {
        tmp = CalculateError(soluteDataBase, h, nRest, sRest, iRest, pow(10, -pH), v);
        if (tmp == -1) return -1;
        else if (tmp<lowest) {
            lowest = tmp;
            ans = pH;
        }
    }   //Locate the point where the error is minimized within the search scope.

    if (interval>precision) return CalculatePH(soluteDataBase, h, nRest, sRest, iRest, v, ans-interval, ans+interval, precision, interval*0.1);
    //If the search interval is greater than the allowable error, narrow down the search interval to re-discover both sides of the point.
    else return ans;
    //If the search interval is less than or equal to an acceptable error, return the result.
}

int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute) {
    printf("INFO : PhCalculator function started\n");

    int i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER];
    double vAll, h, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER];
    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic solutions : "); scanf("%s", tmp); nAcid=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic solutions.
    printf("Number of strongly basic solutions : "); scanf("%s", tmp); nBase=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        return -1;
    }   //Get the number of strongly basic solutions.
    printf("Number of the other solutions : "); scanf("%s", tmp); nRest=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        return -1;
    }   //Get the number of the other solutions.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=atof(tmp); printf("\n"); fflush(stdin);
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic solution (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sAcid[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic solutions.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic solution (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic solutions.

    for (i=0;i<nRest;i++) {
        printf("Name of another solution : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) return -1;
        printf("Volume of the solution (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the solutions.

    printf("\n");

    h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll);
    printf("%.3lf", CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, 1));
    //Calculate and print out the pH value.

    return 0;
}

int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], nSolute) {
    printf("INFO : GraphGenerator function started\n");

    int i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER];
    double vAll, h, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER], volAdding, volPerTime, cenAdded;
    char result[6], tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic titrands : "); scanf("%s", tmp); nAcid=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic titrands.
    printf("Number of strongly basic titrands : "); scanf("%s", tmp); nBase=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        return -1;
    }   //Get the number of strongly basic titrands.
    printf("Number of the other titrands : "); scanf("%s", tmp); nRest=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        return -1;
    }   //Get the number of the other titrands.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=atof(tmp); printf("\n"); fflush(stdin);
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sAcid[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic titrands.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic titrands.
    for (i=0;i<nRest;i++) {
        printf("Name of another titrand : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) return -1;
        printf("Volume of the titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the titrands.

    printf("\n");

    printf("Type of titrant [strong monoprotic acid : 0 | strong base : 1 | the other : 2] : "); scanf("%s", tmp); printf("\n"); fflush(stdin);
    //Get type of titrant.

    FILE* file;
    file = fopen("../result.pcd", "w");
    switch (tmp[0]) {
        case '0' :
            nAcid++;
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n"); fflush(stdin);
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n"); fflush(stdin);
            //Get data of strongly acidic monoprotic titrant.

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sAcid[nAcid-1]=cenAdded*volPerTime*i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+volPerTime*i);
                gcvt(CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll+volPerTime*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        case '1' :
            nBase++;
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n"); fflush(stdin);
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n"); fflush(stdin);
            //Get data of strongly basic titrands.

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sBase[nBase - 1] = cenAdded * volPerTime * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll + volPerTime * i);
                gcvt(CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll + volPerTime * i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        case '2' :
            nRest++;
            printf("Name of the titrant : "); scanf("%s", name); iRest[nRest-1] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
            if (iRest[nRest-1] == -1) return -1;
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volAdding=atof(tmp); printf("\n"); fflush(stdin);
            printf("Volume of the titrant (L) : "); scanf("%s", tmp); volPerTime=atof(tmp); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenAdded=atof(tmp); printf("\n"); fflush(stdin);
            //Get data of the other titrant.

            for (i=0; i<=volAdding/volPerTime; i++) {
                result[0] = "      ";
                sRest[nRest - 1] = cenAdded * volPerTime * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll + volPerTime * i);
                gcvt(CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll + volPerTime * i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, 1), 6, result);
                fputs(result, file);
                fputs("\n", file);
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        default :
            printf("ERR : Invalid solution type\n");
            return -1;
    }
    fclose(file);

    printf("\n\n");

    printf("INFO : Graph has been generated\n");

    return 0;
}

int RecipeFinder(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], int nSolute) {
    printf("INFO : RecipeFinder function started\n");

    int i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER];
    double pH, vAll, h, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER], error, target_pH, cenTitrant;
    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic titrands : "); scanf("%s", tmp); nAcid=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic titrands.
    printf("Number of strongly basic titrands : "); scanf("%s", tmp); nBase=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        return -1;
    }   //Get the number of strongly basic titrands.
    printf("Number of the other titrands : "); scanf("%s", tmp); nRest=atoi(tmp); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        return -1;
    }   //Get the number of the other titrands.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=atof(tmp); printf("\n");
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sAcid[i]=cen*vol;  printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic titrands.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic titrands.
    for (i=0;i<nRest;i++) {
        printf("Name of another titrand : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) return -1;
        printf("Volume of the titrand (L) : "); scanf("%s", tmp); vol=atof(tmp); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=atof(tmp); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the titrands.

    printf("\n");

    printf("Target pH : "); scanf("%s", tmp); target_pH=atof(tmp); printf("\n");
    //Get target pH value.

    printf("\n\n");

    printf("Type of titrant [strong monoprotic acid : 0 | strong base : 1 | the other : 2] : "); scanf("%s", tmp); printf("\n"); fflush(stdin);
    //Get type of titrant.

    i=0;
    switch (tmp[0]) {
        case '0' :
            nAcid++;
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=atof(tmp); printf("\n");
            //Get data of strongly acidic monoprotic titrant.

            do {
                sAcid[nAcid - 1] = cenTitrant * NEUTRALPOINT_FINDER_INTERVAL * vAll * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i);
                pH = CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, PH_CALCULATOR_INITIALINTERVAL);

                error = fabs(pH-target_pH);

                i++;
            } while (error>NEUTRALPOINT_FINDER_PRECISION);
            //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        case '1' :
            nBase++;
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=atof(tmp); printf("\n");
            //Get data of strongly basic titrands.

            do {
                sBase[nBase - 1] = cenTitrant * NEUTRALPOINT_FINDER_INTERVAL * vAll * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i);
                pH = CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, PH_CALCULATOR_INITIALINTERVAL);

                error = fabs(pH-target_pH);

                i++;
            } while (error>NEUTRALPOINT_FINDER_PRECISION);
            //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        case '2' :
            nRest++;
            printf("Name of the titrant : "); scanf("%s", name); iRest[nRest-1] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n");
            if (iRest[i] == -1) return -1;
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=atof(tmp); printf("\n");
            //Get data of the other titrant.

            do {
                sRest[nRest - 1] = cenTitrant * NEUTRALPOINT_FINDER_INTERVAL * vAll * i;
                h = CalcInitialH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i);
                pH = CalculatePH(soluteDataBase, h, nRest, sRest, iRest, vAll + NEUTRALPOINT_FINDER_INTERVAL * vAll * i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_PRECISION, PH_CALCULATOR_INITIALINTERVAL);

                error = fabs(pH - target_pH);

                i++;
            } while (error>NEUTRALPOINT_FINDER_PRECISION);
            //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        default :
            printf("ERR : Invalid solution type\n");
            return -1;
    }

    printf("\n\n");

    printf("volume : %.6lf | pH : %.3lf | error : %.3lf\n", NEUTRALPOINT_FINDER_INTERVAL * vAll * (i-1), pH, error);
    //Print out the founded recipe.

    return 0;
}