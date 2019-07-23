/*
version : v1.1.10-alpha

MIT License

Copyright (c) 2019 Kyung-ha Lee <i_am@nulleekh.com>

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
#include "pHCalc.h"
//Import libraries required for operation.

int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute);
int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute);
int RecipeFinder(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute);
//Declare the function prototypes

int main() {
    printf("INFO : Program started\n");

    long i, j, k, nSolute;
    char tmp = 0x00, soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("INFO : Loading solution data\n");

    FILE *soluteDataBaseFile = fopen("solute.pcd", "r");

    if (soluteDataBaseFile==NULL) {
        printf("ERR : Failed to open DB\n");
        return -1;
    }

    for (i=0; i<NUMBER_OF_SOLUTE; i++) {
        for (j = 0; j < NUMBER_OF_DATA; j++){
            for (k = 0; k < MAX_DATA_LENGTH; k++) {
                tmp = fgetc(soluteDataBaseFile);
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
    fclose(soluteDataBaseFile);
    //Get the database needed to operate the program.

    printf("INFO : Solution data Loaded successfully\n");

    printf("Select calculation type you want [pH-calculator : 0 | graph-generator : 1 | recipe-finder : 2] : "); scanf("%c", &tmp); printf("\n"); fflush(stdin);
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

int PhCalculator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute) {
    printf("INFO : PhCalculator function started\n");

    long i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER];
    double vAll, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER];
    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic solutions : "); scanf("%s", tmp); nAcid=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        printf("INFO : PhCalculator function ended\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic solutions.
    printf("Number of strongly basic solutions : "); scanf("%s", tmp); nBase=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        printf("INFO : PhCalculator function ended\n");
        return -1;
    }   //Get the number of strongly basic solutions.
    printf("Number of the other solutions : "); scanf("%s", tmp); nRest=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        printf("INFO : PhCalculator function ended\n");
        return -1;
    }   //Get the number of the other solutions.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=strtod(tmp, NULL); printf("\n"); fflush(stdin);
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic solution (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sAcid[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic solutions.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic solution (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic solutions.

    for (i=0;i<nRest;i++) {
        printf("Name of another solution : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) {
            printf("INFO : PhCalculator function ended\n");
            return -1;
        }
        printf("Volume of the solution (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the solution (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the solutions.

    printf("\n");

    printf("%.3lf\n\n", CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, 1));
    //Calculate and print out the pH value.

    printf("INFO : PhCalculator function ended\n");
    return 0;
}

int GraphGenerator(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute) {
    printf("INFO : GraphGenerator function started\n");

    long i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER], numAddTitrant;
    double vAll, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER], volPerTime, cenTitrant;
    char result[6], tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic titrands : "); scanf("%s", tmp); nAcid=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        printf("INFO : GraphGenerator function ended\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic titrands.
    printf("Number of strongly basic titrands : "); scanf("%s", tmp); nBase=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        printf("INFO : GraphGenerator function ended\n");
        return -1;
    }   //Get the number of strongly basic titrands.
    printf("Number of the other titrands : "); scanf("%s", tmp); nRest=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        printf("INFO : GraphGenerator function ended\n");
        return -1;
    }   //Get the number of the other titrands.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=strtod(tmp, NULL); printf("\n"); fflush(stdin);
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sAcid[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic titrands.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic titrands.
    for (i=0;i<nRest;i++) {
        printf("Name of another titrand : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) {
            printf("INFO : GraphGenerator function ended\n");
            return -1;
        }
        printf("Volume of the titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the titrands.

    printf("\n");

    printf("Type of titrant [strong monoprotic acid : 0 | strong base : 1 | the other : 2] : "); scanf("%s", tmp); printf("\n"); fflush(stdin);
    //Get type of titrant.

    FILE* resultGraphFile;
    resultGraphFile = fopen("result.pcd", "w");

    if (resultGraphFile==NULL) {
        printf("ERR : Failed to open DB\n");
        printf("INFO : GraphGenerator function ended\n");
        return -1;
    }

    i=0;
    switch (tmp[0]) {
        case '0' :
            nAcid++;

            printf("Number of times to add titrant (times) : "); scanf("%s", tmp); numAddTitrant=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
            printf("Volume of the titrant per time (L) : "); scanf("%s", tmp); volPerTime=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            //Get data of strongly acidic monoprotic titrant.

            while (1) {
                result[0] = 0x00;

                sAcid[nAcid-1]=cenTitrant*volPerTime*i;

                gcvt(CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+volPerTime*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, 1), 6, result);

                fputs(result, resultGraphFile);
                fputs("\n", resultGraphFile);

                i++;

                if (i>numAddTitrant) break;
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        case '1' :
            nBase++;

            printf("Number of times to add titrant (times) : "); scanf("%s", tmp); numAddTitrant=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
            printf("Volume of the titrant per time (L) : "); scanf("%s", tmp); volPerTime=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            //Get data of strongly basic titrands.

            while (1) {
                result[0] = 0x00;

                sBase[nBase - 1] = cenTitrant * volPerTime * i;

                gcvt(CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+volPerTime*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, 1), 6, result);

                fputs(result, resultGraphFile);
                fputs("\n", resultGraphFile);

                i++;

                if (i>numAddTitrant) break;
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        case '2' :
            nRest++;

            printf("Name of the titrant : "); scanf("%s", name); iRest[nRest-1] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
            if (iRest[nRest-1] == -1) {
                printf("INFO : GraphGenerator function ended\n");
                return -1;
            }
            printf("Number of times to add titrant (times) : "); scanf("%s", tmp); numAddTitrant=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
            printf("Volume of the titrant per time (L) : "); scanf("%s", tmp); volPerTime=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n"); fflush(stdin);
            //Get data of the other titrant.

            while (1) {
                result[0] = 0x00;

                sRest[nRest - 1] = cenTitrant * volPerTime * i;

                gcvt(CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+volPerTime*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, 1), 6, result);

                fputs(result, resultGraphFile);
                fputs("\n", resultGraphFile);

                i++;

                if (i>numAddTitrant) break;
            }   //Perform titration experimental simulation as instructed and store results in result.pcd.

            break;

        default :
            printf("ERR : Invalid solution type\n");
            printf("INFO : GraphGenerator function ended\n");
            return -1;
    }
    fclose(resultGraphFile);

    printf("\n\n");

    printf("INFO : Graph has been generated\n\n");

    printf("INFO : GraphGenerator function ended\n");
    return 0;
}

int RecipeFinder(char soluteDataBase[NUMBER_OF_SOLUTE][NUMBER_OF_DATA][MAX_DATA_LENGTH], long nSolute) {
    printf("INFO : RecipeFinder function started\n");

    long i, nAcid, nBase, nRest, iRest[MAX_SOLUTION_NUMBER];
    double vAll, vol, cen, sAcid[MAX_SOLUTION_NUMBER], sBase[MAX_SOLUTION_NUMBER], sRest[MAX_SOLUTION_NUMBER], cenTitrant, titrand_pH, titrant_pH, target_pH, accuracy, error[2], pH[2];
    char tmp[MAX_NUMBER_LENGTH], name[MAX_DATA_LENGTH] = { 0, };
    //Declare the variables needed to operate the function.

    printf("Number of strongly acidic monoprotic titrands : "); scanf("%s", tmp); nAcid=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nAcid || nAcid>255) {
        printf("ERR : Invalid nAcid number\n");
        printf("INFO : RecipeFinder function ended\n");
        return -1;
    }   //Get the number of strongly acidic monoprotic titrands.
    printf("Number of strongly basic titrands : "); scanf("%s", tmp); nBase=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nBase || nBase>255) {
        printf("ERR : Invalid nBase number\n");
        printf("INFO : RecipeFinder function ended\n");
        return -1;
    }   //Get the number of strongly basic titrands.
    printf("Number of the other titrands : "); scanf("%s", tmp); nRest=strtol(tmp, NULL, 10); printf("\n"); fflush(stdin);
    if (0>nRest || nRest>255) {
        printf("ERR : Invalid nRest number\n");
        printf("INFO : RecipeFinder function ended\n");
        return -1;
    }   //Get the number of the other titrands.

    printf("\n\n");

    printf("Volume of pure water (L) : "); scanf("%s", tmp); vAll=strtod(tmp, NULL); printf("\n");
    //Get volume of pure water.

    printf("\n\n");

    for(i=0;i<nAcid;i++) {
        printf("Volume of strongly acidic monoprotic titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sAcid[i]=cen*vol;  printf("\n\n"); fflush(stdin);
    }   //Get data of strongly acidic monoprotic titrands.
    for(i=0;i<nBase;i++) {
        printf("Volume of strongly basic titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sBase[i]=cen*vol; printf("\n\n"); fflush(stdin);
    }   //Get data of strongly basic titrands.
    for (i=0;i<nRest;i++) {
        printf("Name of another titrand : "); scanf("%s", name); iRest[i] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n"); fflush(stdin);
        if (iRest[i] == -1) {
            printf("INFO : RecipeFinder function ended\n");
            return -1;
        }
        printf("Volume of the titrand (L) : "); scanf("%s", tmp); vol=strtod(tmp, NULL); vAll=vAll+vol; printf("\n"); fflush(stdin);
        printf("Concentration of the titrand (mol/L) : "); scanf("%s", tmp); cen=strtod(tmp, NULL); sRest[i]=vol*cen; printf("\n\n"); fflush(stdin);
    }   //Get data of the titrands.

    printf("\n");

    printf("Target pH : "); scanf("%s", tmp); target_pH=strtod(tmp, NULL); printf("\n");
    //Get target pH value.

    printf("\n\n");

    titrand_pH = CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

    printf("Type of titrant [strong monoprotic acid : 0 | strong base : 1 | the other : 2] : "); scanf("%s", tmp); printf("\n"); fflush(stdin);
    //Get type of titrant.

    i=0; error[1]=1024;
    switch (tmp[0]) {
        case '0' :
            nAcid++;
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n");
            //Get data of strongly acidic monoprotic titrant.

            printf("Accuracy of calculation (L) : "); scanf("%s", tmp); accuracy=strtod(tmp, NULL); printf("\n");
            //Get the accuracy required for calculation.

            titrant_pH = CalculatePH(soluteDataBase, 1, 0, 0, &cenTitrant, 0, 0, 0, 1, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

            if (titrand_pH<target_pH) {
                if (titrant_pH<target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            } else if (titrand_pH>target_pH) {
                if (titrant_pH>target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            }

            while (1) {
                error[0] = error[1];
                pH[0] = pH[1];

                sAcid[nAcid - 1] = cenTitrant * accuracy*i;

                pH[1] = CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+accuracy*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

                error[1] = fabs(pH[1] - target_pH);

                i++;

                if (error[1] >= error[0] && (error[0] <= RECIPE_FINDER_RESOLUTION || error[1]-error[0] >= RECIPE_FINDER_RESOLUTION)) break;
            }   //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        case '1' :
            nBase++;
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n");
            //Get data of strongly basic titrant.

            printf("Accuracy of calculation (L) : "); scanf("%s", tmp); accuracy=strtod(tmp, NULL); printf("\n");
            //Get the accuracy required for calculation.

            titrant_pH = CalculatePH(soluteDataBase, 0, 1, 0, 0, &cenTitrant, 0, 0, 1, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

            if (titrand_pH<target_pH) {
                if (titrant_pH<target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            } else if (titrand_pH>target_pH) {
                if (titrant_pH>target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            }

            while (1) {
                error[0] = error[1];
                pH[0] = pH[1];

                sBase[nBase - 1] = cenTitrant * accuracy*i;

                pH[1] = CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+accuracy*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

                error[1] = fabs(pH[1] - target_pH);

                i++;

                if (error[1] >= error[0] && (error[0] <= RECIPE_FINDER_RESOLUTION || error[1]-error[0] >= RECIPE_FINDER_RESOLUTION)) break;
            }   //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        case '2' :
            nRest++;
            printf("Name of the titrant : "); scanf("%s", name); iRest[nRest-1] = SpecifySolute(soluteDataBase, nSolute, name); printf("\n");
            if (iRest[i] == -1) {
                printf("INFO : RecipeFinder function ended\n");
                return -1;
            }
            printf("Concentration of the titrant (mol/L) : "); scanf("%s", tmp); cenTitrant=strtod(tmp, NULL); printf("\n");
            //Get data of the other titrant.

            printf("Accuracy of calculation (L) : "); scanf("%s", tmp); accuracy=strtod(tmp, NULL); printf("\n");
            //Get the accuracy required for calculation.

            titrant_pH = CalculatePH(soluteDataBase, 0, 0, 1, 0, 0, &cenTitrant, &iRest[nRest-1], 1, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

            if (titrand_pH<target_pH) {
                if (titrant_pH<target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            } else if (titrand_pH>target_pH) {
                if (titrant_pH>target_pH) {
                    printf("ERR : Invalid titrant\n");
                    printf("INFO : RecipeFinder function ended\n");
                    return -1;
                }
            }

            while (1) {
                error[0] = error[1];
                pH[0] = pH[1];

                sRest[nRest - 1] = cenTitrant*accuracy*i;

                pH[1] = CalculatePH(soluteDataBase, nAcid, nBase, nRest, sAcid, sBase, sRest, iRest, vAll+accuracy*i, PH_CALCULATOR_STARTPOINT, PH_CALCULATOR_ENDPOINT, PH_CALCULATOR_RESOLUTION, PH_CALCULATOR_INITIALINTERVAL);

                error[1] = fabs(pH[1] - target_pH);

                i++;

                if (error[1] > error[0] && (error[0] <= RECIPE_FINDER_RESOLUTION || error[1]-error[0] > RECIPE_FINDER_RESOLUTION)) break;
            }   //Perform titration experimental simulation execute according to the prescribed rule until reach target pH value.

            break;

        default :
            printf("ERR : Invalid solution type\n");
            printf("INFO : RecipeFinder function ended\n");
            return -1;
    }

    printf("\n\n");

    printf("volume : %.9lf | pH : %.3lf | error : %.3lf\n\n", accuracy*(i-2), pH[0], error[0]);
    //Print out the founded recipe.

    printf("INFO : RecipeFinder function ended\n");
    return 0;
}
