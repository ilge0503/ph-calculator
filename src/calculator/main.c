/*
version : v1.0.1-alpha

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

int main( )
{
    int num, i;
    double x[2] = {10e-8, 10e-8}, h, kw, ka[100], ha[100], cost = -1;
    //Declares the various variables to be used in the calculation

    scanf("%d %lf %lf", &num, &h, &kw);
    //Get the number of solutions participating in the reaction, initial hydrogen ion concentration, and Kw value

    for (i=0; i<num; i++) scanf(" %lf %lf ", &ka[i], &ha[i]);
    //Get ka and ha values of the solutions participating in the reaction

    if (h>=0) {
        while (cost <= -10e-8 || 10e-8 <= cost) {
            x[0] = (h + (kw / x[1]));
            for (i = 0; i < num; i++) x[0] += (ka[i] * ha[i]) / (x[1] + ka[i]);

            cost = x[1] - x[0];
            printf("x = %.15lf | cost = %.15lf\n", x[0], cost);

            x[1] = x[0];
        }
        printf("pH = %.15lf | cost = %.15lf\n", x[0], cost);
    }
        //Calculation of acidity by formula (when the initial hydrogen ion concentration is positive)
    else if (h<0) {
        x[0] = 1; x[1] = 1;
        while (cost <= -10e-8 || 10e-8 <= cost) {
            x[0] = (-h + (kw / x[1]));
            for (i = 0; i < num; i++) x[0] -= (x[1] * ka[i] * ha[i]) / ((kw / x[1]) + ka[i]);

            cost = x[1] - x[0];
            printf("x = %.15lf | cost = %.15lf\n", x[0], cost);

            x[1] = x[0];
        }
        printf("pOH = %.15lf | cost = %.15lf\n", x[0], cost);
    }
        //Calculation of acidity by formula (when the initial hydrogen ion concentration is negative)
    else printf("ERROR : H+ value is not available.");
    //Print error message when calculation is not possible
    return 0;
}