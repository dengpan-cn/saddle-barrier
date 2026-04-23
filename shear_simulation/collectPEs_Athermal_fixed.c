
#include "AthermRelaxMin_ndim.h"
#include "GeoFamily_ndim.h"
#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"

#define MaxStrain 2E-2
#define StrainStep 1E-7
#define Tmax 1E-10

int main(int argc, char const* argv[]) {
    MSG = "Collecting PEs using Athermal limiting Temperature Min";

    double tic = getTimeStamp();
    char ticStr[32];
    getTimeString(ticStr);
    //=============================================================
    Variable* var = (Variable*)calloc(1, sizeof(Variable));
    readCmdLineVar(var, argc, argv);

    // 0: current; 1: last; 2: conf. before reverse shear.
    Box* box = (Box*)calloc(3, sizeof(Box));
    Particle* particle = (Particle*)calloc(3, sizeof(Particle));
    Update* update = (Update*)calloc(1, sizeof(Update));

    readConf(&box[0], &particle[0], &update[0], var);

    cmdArg* cmd = findVariable(var, "shear");
    if (!cmd || cmd->cmdArgc != 2)
        Abort("--shear yz/xz/xy +/-");
    int stype = -1, sdir = 0;
    stype = (strcmp(cmd->cmdArgv[0], "yz") ? (strcmp(cmd->cmdArgv[0], "xz") ? 2 : 1) : 0);
    sdir = (strcmp(cmd->cmdArgv[1], "+") ? -1 : +1);
    if (stype == -1 || sdir == 0)
        Abort("--shear yz/xz/xy +/-");
    char fname[4096];
    sprintf(fname, "%s/PEs_fixed_%s_%s_%%8lf_%s.bin", var->cwd, (stype == 0 ? "YZ" : (stype == 1 ? "XZ" : "XY")), (sdir == 1 ? "P" : "N"), var->sf);
    FILE* fout = createFileReadWrite(fname);
    sprintf(fname, "%s/cvs_fixed_%s_%s_%%8lf_%s.bin", var->cwd, (stype == 0 ? "YZ" : (stype == 1 ? "XZ" : "XY")), (sdir == 1 ? "P" : "N"), var->sf);
    FILE* fcvs = createFileReadWrite(fname);

    AthermRelaxMin* fire = addAthermRelaxMin(&box[0], &particle[0], &update[0], var);
    GeoFamily* geoFam = addGeoFamily(&box[0], &particle[0], &update[0]);

    md5FingerPrint cmSign[3];
    uptriMat stress[3];
    double ePair[3];
    bool jamState[3];
    double dStrain = StrainStep;
    double strain = 0.0;

    reInitAthermRelaxMin(fire, Tmax);
    minAthermRelax(&box[0], &particle[0], &update[0]);
    assignGeoFamily(&box[0], &particle[0], &update[0]);

    // backup configuration;
    backupSimInfo(&box[0], &particle[0], &box[1], &particle[1]);
    memcpy(&cmSign[1], geoFam->fingerPrint, sizeof(md5FingerPrint));
    uptriMatCpy(stress[1], update->pVirTens);
    ePair[1] = update->ePair;
    jamState[1] = (update->ePair / update->energyUnits < __ZeroEnergy__ ? false : true);
    {
        double data[8];
        data[0] = strain;
        data[1] = update->ePair / update->energyUnits;
        data[2] = update->pVirTens[spaceIdx2voigt(0, 0)] / update->pressureUnits;
        data[3] = update->pVirTens[spaceIdx2voigt(1, 1)] / update->pressureUnits;
        data[4] = update->pVirTens[spaceIdx2voigt(2, 2)] / update->pressureUnits;
        data[5] = update->pVirTens[spaceIdx2voigt(1, 2)] / update->pressureUnits;
        data[6] = update->pVirTens[spaceIdx2voigt(0, 2)] / update->pressureUnits;
        data[7] = update->pVirTens[spaceIdx2voigt(0, 1)] / update->pressureUnits;
        safeFwrite(fcvs, data, sizeof(double), 8);
    }

    addWriteDumpFile(&box[0], &particle[0], update, var);
    void (*simplShear)(Box*, Particle*, Update*, double);
    simplShear = (stype == 0 ? instant_simpShearYz : (stype == 1 ? instant_simpShearXz : instant_simpShearXy));

    while (fabs(strain) < MaxStrain) {
        simplShear(&box[0], &particle[0], update, sdir * dStrain);
        reInitAthermRelaxMin(fire, Tmax);
        minAthermRelax(&box[0], &particle[0], &update[0]);
        assignGeoFamily(&box[0], &particle[0], &update[0]);
        // get current state
        memcpy(&cmSign[0], geoFam->fingerPrint, sizeof(md5FingerPrint));
        uptriMatCpy(stress[0], update->pVirTens);
        ePair[0] = update->ePair;
        jamState[0] = (update->ePair / update->energyUnits < __ZeroEnergy__ ? false : true);
        strain += (sdir * dStrain);

        // discard unjamming steps
        if ((!jamState[1]) || (!jamState[0])) {
            backupSimInfo(&box[0], &particle[0], &box[1], &particle[1]);
            memcpy(&cmSign[1], &cmSign[0], sizeof(md5FingerPrint));
            uptriMatCpy(stress[1], stress[0]);
            ePair[1] = ePair[0];
            jamState[1] = jamState[0];

            {
                double data[8];
                data[0] = strain;
                data[1] = ePair[0];
                data[2] = stress[0][spaceIdx2voigt(0, 0)] / update->pressureUnits;
                data[3] = stress[0][spaceIdx2voigt(1, 1)] / update->pressureUnits;
                data[4] = stress[0][spaceIdx2voigt(2, 2)] / update->pressureUnits;
                data[5] = stress[0][spaceIdx2voigt(1, 2)] / update->pressureUnits;
                data[6] = stress[0][spaceIdx2voigt(0, 2)] / update->pressureUnits;
                data[7] = stress[0][spaceIdx2voigt(0, 1)] / update->pressureUnits;
                safeFwrite(fcvs, data, sizeof(double), 8);
            }

            continue;
        }
        // if the same geometric family (same contact network), do nothing.
        if ((memcmp(&cmSign[0], &cmSign[1], sizeof(md5FingerPrint)) == 0)) {
            backupSimInfo(&box[0], &particle[0], &box[1], &particle[1]);
            memcpy(&cmSign[1], &cmSign[0], sizeof(md5FingerPrint));
            uptriMatCpy(stress[1], stress[0]);
            ePair[1] = ePair[0];
            jamState[1] = jamState[0];

            {
                double data[8];
                data[0] = strain;
                data[1] = ePair[0];
                data[2] = stress[0][spaceIdx2voigt(0, 0)] / update->pressureUnits;
                data[3] = stress[0][spaceIdx2voigt(1, 1)] / update->pressureUnits;
                data[4] = stress[0][spaceIdx2voigt(2, 2)] / update->pressureUnits;
                data[5] = stress[0][spaceIdx2voigt(1, 2)] / update->pressureUnits;
                data[6] = stress[0][spaceIdx2voigt(0, 2)] / update->pressureUnits;
                data[7] = stress[0][spaceIdx2voigt(0, 1)] / update->pressureUnits;
                safeFwrite(fcvs, data, sizeof(double), 8);
            }
            continue;
        }

        // contact network changed: backup
        backupSimInfo(&box[0], &particle[0], &box[2], &particle[2]);
        memcpy(&cmSign[2], &cmSign[0], sizeof(md5FingerPrint));
        uptriMatCpy(stress[2], stress[0]);
        ePair[2] = ePair[0];
        jamState[2] = jamState[0];
        // do reverse shear
        simplShear(&box[0], &particle[0], update, -sdir * dStrain);
        reInitAthermRelaxMin(fire, Tmax);
        minAthermRelax(&box[0], &particle[0], &update[0]);
        assignGeoFamily(&box[0], &particle[0], &update[0]);
        // get state
        memcpy(&cmSign[0], geoFam->fingerPrint, sizeof(md5FingerPrint));
        uptriMatCpy(stress[0], update->pVirTens);
        ePair[0] = update->ePair;
        jamState[0] = (update->ePair / update->energyUnits < __ZeroEnergy__ ? false : true);

        // continuous event: restore the conf. before reverse shear.
        if (memcmp(cmSign[0], cmSign[1], sizeof(md5FingerPrint)) == 0) {
            backupSimInfo(&box[2], &particle[2], &box[1], &particle[1]);
            memcpy(&cmSign[1], &cmSign[2], sizeof(md5FingerPrint));
            uptriMatCpy(stress[1], stress[2]);
            ePair[1] = ePair[2];
            jamState[1] = jamState[2];

            restoreSimInfo(&box[0], &particle[0], &update[0], &box[2], &particle[2]);

            {
                double data[8];
                data[0] = strain;
                data[1] = ePair[1];
                data[2] = stress[1][spaceIdx2voigt(0, 0)] / update->pressureUnits;
                data[3] = stress[1][spaceIdx2voigt(1, 1)] / update->pressureUnits;
                data[4] = stress[1][spaceIdx2voigt(2, 2)] / update->pressureUnits;
                data[5] = stress[1][spaceIdx2voigt(1, 2)] / update->pressureUnits;
                data[6] = stress[1][spaceIdx2voigt(0, 2)] / update->pressureUnits;
                data[7] = stress[1][spaceIdx2voigt(0, 1)] / update->pressureUnits;

                safeFwrite(fcvs, data, sizeof(double), 8);
            }
            continue;
        }

        if (fabs((ePair[1] - ePair[0]) / update->energyUnits) >= 2 * __ZeroEnergy__) {
            double data[8];
            data[0] = strain - (sdir * dStrain);
            data[1] = (ePair[1] - ePair[0]) / update->energyUnits;
            data[2] = (stress[1][spaceIdx2voigt(0, 0)] - stress[0][spaceIdx2voigt(0, 0)]) / update->pressureUnits;
            data[3] = (stress[1][spaceIdx2voigt(1, 1)] - stress[0][spaceIdx2voigt(1, 1)]) / update->pressureUnits;
            data[4] = (stress[1][spaceIdx2voigt(2, 2)] - stress[0][spaceIdx2voigt(2, 2)]) / update->pressureUnits;
            data[5] = (stress[1][spaceIdx2voigt(1, 2)] - stress[0][spaceIdx2voigt(1, 2)]) / update->pressureUnits;
            data[6] = (stress[1][spaceIdx2voigt(0, 2)] - stress[0][spaceIdx2voigt(0, 2)]) / update->pressureUnits;
            data[7] = (stress[1][spaceIdx2voigt(0, 1)] - stress[0][spaceIdx2voigt(0, 1)]) / update->pressureUnits;
            safeFwrite(fout, data, sizeof(double), 8);

            writeDump(&box[1], &particle[1], update);
            writeDump(&box[0], &particle[0], update);
        }

        {
            double data[8];
            data[0] = strain;
            data[1] = ePair[2];
            data[2] = stress[2][spaceIdx2voigt(0, 0)] / update->pressureUnits;
            data[3] = stress[2][spaceIdx2voigt(1, 1)] / update->pressureUnits;
            data[4] = stress[2][spaceIdx2voigt(2, 2)] / update->pressureUnits;
            data[5] = stress[2][spaceIdx2voigt(1, 2)] / update->pressureUnits;
            data[6] = stress[2][spaceIdx2voigt(0, 2)] / update->pressureUnits;
            data[7] = stress[2][spaceIdx2voigt(0, 1)] / update->pressureUnits;

            safeFwrite(fcvs, data, sizeof(double), 8);
        }

        backupSimInfo(&box[2], &particle[2], &box[1], &particle[1]);
        memcpy(&cmSign[1], &cmSign[2], sizeof(md5FingerPrint));
        uptriMatCpy(stress[1], stress[2]);
        ePair[1] = ePair[2];
        jamState[1] = jamState[2];
        restoreSimInfo(&box[0], &particle[0], &update[0], &box[2], &particle[2]);
    }

    delAthermRelaxMin(update);
    delWriteDumpFile(update);
    safeCloseFile(fout);
    safeCloseFile(fcvs);

    //=============================================================
    double toc = getTimeStamp();
    char tocStr[32];
    getTimeString(tocStr);
    safeFprintf(stdout, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    safeFprintf(logFile, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    return EXIT_SUCCESS;
}
