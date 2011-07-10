/*

Minnesota BLUE Results Combiner
Version 1.0
May 10, 2011

J. Pastika and J. Mans

This tar contains version 1.0 of the UMN BLUE results combiner.

This will combine the results of 2 independant measurements each with
their own covariance matrix.  The method is summarized in [1].  The
method is implemented as a root script which accepts the input files
of the seperate measurements as text files. 

[1] Valassi, A. "Combining correlated measurements of several
different physical quantities."  Nuclear instruments & methods in
physics research. Section A, Accelerators, spectrometers, detectors
and associated equipment 500.1-3 (2003):391

 */

#include "TMatrixT.h"
#include "TVectorT.h"
#include <stdio.h>
#include <math.h>
#include <vector>
#include <iostream>

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2] = 0, int ftr = 1, double* ccCovar = NULL);
void readCovMatrix(const char * name, TMatrixT<double>& dv);

void resultCombiner(const char * name_dv1, const char * name_cm1, const int length1, const char * name_dv2, const char * name_cm2, const int length2, int ftr = 1, const char * outputNameStub = "output")
{
    TVectorD dv1(length1), dv2(length2);
    TMatrixT<double> cm1(length1, length1), cm2(length2, length2);
    double binLimits1[length1][2], binLimits2[length2][2], ccCov1[length1], ccCov2[length1];

    readDataVector(name_dv1, dv1, binLimits1, ftr, ccCov1);
    printf("Read data vector 1 (%d)\n",length1);

    readDataVector(name_dv2, dv2, binLimits2, ftr, ccCov2);
    printf("Read data vector 2 (%d)\n",length2);

    readCovMatrix(name_cm1, cm1);
    printf("Read covariance matrix 1\n");

    readCovMatrix(name_cm2, cm2);
    printf("Read covariance matrix 2\n");

    std::vector<double*> binLimits;
    std::vector<std::vector<int > > preU;
    int i1 = 0, i2 = 0;
    while(i1 < length1 || i2 < length2)
    {
        if(i1 < length1 && i2 < length2)
        {
            if((binLimits1[i1][1] + binLimits1[i1][0])/2 > binLimits2[i2][0] && (binLimits1[i1][1] + binLimits1[i1][0])/2 < binLimits2[i2][1])
            {
                binLimits.push_back(binLimits1[i1]);
                std::vector<int> tmp;
                tmp.push_back(i1);
                tmp.push_back(i2);
                preU.push_back(tmp);
                i1++;
                i2++;
            }
            else if((binLimits1[i1][1] + binLimits1[i1][0])/2 <= binLimits2[i2][0])
            {
                binLimits.push_back(binLimits1[i1]);
                std::vector<int> tmp;
                tmp.push_back(i1);
                tmp.push_back(-1);
                preU.push_back(tmp);
                i1++;
            }
            else
            {
                binLimits.push_back(binLimits2[i2]);
                std::vector<int> tmp;
                tmp.push_back(-1);
                tmp.push_back(i2);
                preU.push_back(tmp);
                i2++;
            }
        }
        else if(i1 < length1 && i2 >= length2)
        {
            binLimits.push_back(binLimits1[i1]);
            std::vector<int> tmp;
            tmp.push_back(i1);
            tmp.push_back(-1);
            preU.push_back(tmp);
            i1++;
        }
        else
        {
            binLimits.push_back(binLimits2[i2]);
            std::vector<int> tmp;
            tmp.push_back(-1);
            tmp.push_back(i2);
            preU.push_back(tmp);
            i2++;
        }
    }

    TVectorD dv(length1 + length2);
    for(int i = 0; i < length1 + length2; i++)
    {
        dv[i] = (i < length1) ? dv1[i] : dv2[i - length1];
    }

    TMatrixT<double> cm(length1 + length2, length1 + length2), U(length1 + length2, preU.size());
    for(int i = 0; i < length1; i++)
    {
        for(int j = 0; j < length1; j++)
        {
            cm[i][j] = cm1[i][j];
        }
    }
    for(int i = length1; i < length1 + length2; i++)
    {
        for(int j = length1; j < length1 + length2; j++)
        {
            cm[i][j] = cm2[i - length1][j - length1];
        }
    }

    for(unsigned int i = 0; i < preU.size(); i++)
    {
        if(preU[i][0] >= 0) U[preU[i][0]][i] = 1;
        if(preU[i][1] >= 0) U[preU[i][1] + length1][i] = 1;
        if(ftr > 1 && preU[i][0] >= 0 && preU[i][1] >= 0)  cm[preU[i][0]][preU[i][1] + length1] = cm[preU[i][1] + length1][preU[i][0]] = ccCov1[preU[i][0]]*ccCov2[preU[i][1]];
    }
    
    //    cm.Print();

    TMatrixT<double> Ut(U);
    Ut.T();

    TMatrixT<double> cmInv(cm);
    cmInv.Invert();
    TMatrixT<double> step1 = Ut * cmInv * U;
    TMatrixT<double> step2 = Ut * cmInv;
    TMatrixT<double> lambda = step1.Invert() * step2;
    TVectorD bV = lambda*dv;
    TMatrixT<double> bcm = (Ut * cmInv * U).Invert();

    printf("Done with combination.\n");

    //write output
    FILE *file;
    char bVoutName[128], CMoutName[128];
    sprintf(bVoutName, "%s_data.txt", outputNameStub);

    file = fopen(bVoutName, "w");
    if(file)
    {
        fprintf(file, "#\n#%9s %9s %9s %15s %15s\n", "Bin", "Y_min", "Y_max", "Value", "Uncertainty");
        for(int i = 0; i < bV.GetNoElements(); i++)
        {
            fprintf(file, " %9i %9.2f %9.2f %15e %15e\n", i + 1, binLimits[i][0], binLimits[i][1], bV[i], sqrt(bcm[i][i]));
        }
        fclose(file);
    }

    sprintf(CMoutName, "%s_covMat.txt", outputNameStub);

    file = fopen(CMoutName, "w");
    if(file)
    {
        fprintf(file, "#\n#%9s %9s %15s\n", "Bin i", "Bin j", "Value");
        for(int i = 0; i < bcm.GetNrows(); i++)
        {
            for(int j = 0; j < bcm.GetNcols(); j++)
            {
                fprintf(file, " %9i %9i %15e\n", i + 1, j + 1, bcm[i][j]);
            }
        }
        fclose(file);
    }
    printf("Output complete.\n");
}

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2], int ftr, double* ccCovar)
{
    //read data vectors
    FILE *file;
    char *c, buff[2048], *k;
    int i;
    float d1, d2, d3, f;

    file = fopen(name, "r");
    if(file)
    {
        while(!feof(file))
        {
            c = fgets(buff, 2048, file);
            for(k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';

            if(c != NULL && buff[0] != '#')
            {
                if(sscanf(buff, "%d %f %f %f", &i, &d1, &d2, &d3) == 4)
                {
                    if(binLimits != 0)
                    {
                        binLimits[i - 1][0] = d1;
                        binLimits[i - 1][1] = d2;
                    }
                    dv[i - 1] = d3;
                }
                if(ccCovar != NULL && ftr > 1)
                {
                	ccCovar[i] = 0.0;
                	char delim[128];
                	sprintf(delim, "%%*d");
                	for(int j = 0; j <= ftr; j++)
                	{
                		sprintf(delim, "%s %%*f", delim);
                	}
                	sprintf(delim, "%s %%f", delim);
                	if(sscanf(buff, delim, &f) == 1) 
                	{
                		ccCovar[i - 1] = f;
                	}
                }
            }
        }
        fclose(file);
    }

    return;
}

void readCovMatrix(const char * name, TMatrixT<double>& cm)
{
    //read covariance matrix
    FILE *file;
    char *k, *c, buff[2048];
    int i = 0, j = 0;
    float d;

    file = fopen(name, "r");
    if(file)
    {
        while(!feof(file))
        {
            c = fgets(buff, 2048, file);
            for(k = strchr(buff, ','); k != 0; k = strchr(buff, ',')) *k = ' ';
            if(c != NULL && buff[0] != '#')
            {
                if(sscanf(buff, "%d %d %f", &i, &j, &d) == 3)
                {
                    cm[i - 1][j - 1] = d;
                }
            }
        }
        fclose(file);
    }

    return;
}
