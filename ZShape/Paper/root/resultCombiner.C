//Hello C++

#include "TMatrixT.h"
#include "TVectorT.h"
#include "stdio.h"

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2] = 0);
void readCovMatrix(const char * name, TMatrixT<double>& dv);

void resultCombiner(const char * name_dv1, const char * name_cm1, const char * name_dv2, const char * name_cm2, const int length, const char * outputNameStub = "output")
{
    TVectorD dv1(length), dv2(length);
    TMatrixT<double> cm1(length, length), cm2(length, length);
    double binLimits[length][2];

    readDataVector(name_dv1, dv1, binLimits);
    readDataVector(name_dv2, dv2);

    readCovMatrix(name_cm1, cm1);
    readCovMatrix(name_cm2, cm2);

    TVectorD dv(2 * length);
    for(int i = 0; i < 2 * length; i++)
    {
        dv[i] = (i < length) ? dv1[i] : dv2[i - length];
    }

    TMatrixT<double> cm(2 * length, 2 * length), U(2 * length, length);
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < length; j++)
        {
            cm[i][j] = cm1[i][j];
        }
    }
    for(int i = length; i < 2 * length; i++)
    {
        for(int j = length; j < 2 * length; j++)
        {
            cm[i][j] = cm2[i - length][j - length];
        }
    }
    for(int i = 0; i < length; i++)
    {
        U[i][i] = 1;
        U[i + length][i] = 1;
    }
    TMatrixT<double> Ut(U);
    Ut.T();

    TMatrixT<double> cmInv(cm);
    cmInv.Invert();
    TMatrixT<double> step1 = Ut * cmInv * U;
    TMatrixT<double> step2 = Ut * cmInv;
    TMatrixT<double> lambda = step1.Invert() * step2;
    TVectorD bV = lambda*dv;
    TMatrixT<double> bcm = (Ut * cmInv * U).Invert();

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
}

void readDataVector(const char * name, TVectorD& dv, double binLimits[][2])
{
    //read data vectors
    FILE *file;
    char *c, buff[2048], *k;
    int i;
    float d1, d2, d3;

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
