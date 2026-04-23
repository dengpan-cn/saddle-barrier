//
//  AthermalRelaxMin.c
//  TimeDrivenSimLib
//
//  Created by Deng Pan on 2023/9/7.
//

#include "AthermRelaxMin_ndim.h"

#define DELAYSTEP 5
#define DT_GROW 1.1
#define ALPHA_SHRINK 0.99
#define DT_SHRINK 0.5
#define ALPHA0 0.1

void mAR_Integrate(Box *box, Particle *particle, Update *update, AthermRelaxMin *fire) {
    // VV Method: do first half step
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScale(veloc, fire->scaleVeloc, veloc);
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    particle->isForceValid = false;
}
void mAR_DotProduct(Box *box, Particle *particle, Update *update, AthermRelaxMin *fire) {
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    // double dotPrd = 0.0;
    //__float128 sumForce = 0;
    double sumForce = 0;  //?
    // do second half step then calculate dotProduct
    //double totKin = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        //totKin += sNormP2(veloc);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    fire->aveForce = sumForce / particle->nAtom / update->forceUnits;
    
    double Tcur = vdotv / (DIM * particle->nAtom);
    double fact = pow((5E-3 / update->timeUnits) / fire->dt, 2);
    fact = (fact > 1.0 ? 1.0 : fact);
    fire->scaleVeloc = (Tcur < fire->Tlimit * fact ? 1.0 : sqrt(fire->Tlimit * fact / Tcur));
}

AthermRelaxMin *addAthermRelaxMin(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getAthermRelaxMin(update) != NULL) {
        return getAthermRelaxMin(update);
    }
    AthermRelaxMin *fire = (AthermRelaxMin *)calloc(sizeof(AthermRelaxMin), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__AthermRelaxMin__");
    
    if (var == NULL) {
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
    
    cmdArg *cmd = findVariable(var, "atherm");
    if (cmd == NULL) {
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
        
    if(cmd->cmdArgc == 0){
        fire->dtSet = 5E-3;
        fire->Tlimit = 1E-10;
        fire->isInit = false;
        return fire;
    }
    if (cmd->cmdArgc != 1) Abort("--atherm 1E-8");
    
    fire->dtSet = 5E-3;
    fire->Tlimit = atof(cmd->cmdArgv[0]);
    if (fire->Tlimit < 1E-20)
        Abort("upper limits of temperature is too small!");
    
    fire->isInit = true;
    return fire;
}
int reInitAthermRelaxMin(AthermRelaxMin *fire, double T_limit) {
    memset(fire, '\0', sizeof(AthermRelaxMin));
    
    fire->dtSet = 5E-3;
    fire->Tlimit = T_limit;
    if (fire->Tlimit < 1E-20)
        Abort("upper limits of temperature is too small!");
    
    fire->isInit = true;
    return 1;
}
AthermRelaxMin *getAthermRelaxMin(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__AthermRelaxMin__");
    if (whichTool < 0)
        return NULL;
    return (AthermRelaxMin *)update->toolkit.toolkit[whichTool];
}
int minAthermRelax(Box *box, Particle *particle, Update *update) {
    AthermRelaxMin *fire = getAthermRelaxMin(update);
    if (!fire || !fire->isInit)
        Abort("Call addAthermRelaxMin(...) or reInitAthermRelaxMin(...)!");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    update->nebrList.nDelay = 0;
    update->nebrList.skinSet = 0.1;
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    
    double dtmax = 0.5 * update->timeUnits;
    int currStep = 0, last_negative = 0;
    double alpha = ALPHA0;
    fire->dt = fire->dtSet * update->timeUnits;
    double dt0 = fire->dt;
    double ePairLastNeg = 0;
    int rVal = 0;
    
    calcForce(box, particle, update);
    mAR_DotProduct(box, particle, update, fire);
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    fire->scaleVeloc = 1.0;
    ePairLastNeg = update->ePair;
    
    while (!(fire->aveForce <= __ZeroForce__)) {
        // update x
//        if (fire->scaleVeloc < 1.0)
//            printf("%d %g\n", currStep, fire->scaleVeloc);
        mAR_Integrate(box, particle, update, fire);
        
        // update v
        calcForce(box, particle, update);
        mAR_DotProduct(box, particle, update, fire);
        currStep++;
        
        // mixing v and update timestep
        if (fire->vdotf > 0) {
            if (currStep - last_negative > DELAYSTEP) {
                fire->dt *= DT_GROW;
                fire->dt = cpuMin(fire->dt, dtmax);
                alpha *= ALPHA_SHRINK;
            }
            
            double scale1 = 1.0 - alpha;
            double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                doubleVecPtr veloc = particle->veloc[iatom];
                doubleVecPtr force = particle->force[iatom];
                vScale(veloc, scale1, veloc);
                vScaleAdd(veloc, veloc, scale2, force);
            }
        } else {
            int deltaStep = currStep - last_negative;
            if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 0.9;
                dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * update->timeUnits);
            } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 1.1;
                dtmax = cpuMin(dtmax, 0.5 * update->timeUnits);
            }
            
            fire->dt *= DT_SHRINK;
            fire->dt = cpuMax(fire->dt, 0.1 * fire->dtSet * update->timeUnits);
            
            dt0 = fire->dt;
            last_negative = currStep;
            ePairLastNeg = update->ePair;
            alpha = ALPHA0;
            memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
            fire->scaleVeloc = 1.0;
        }
        
        if (update->ePair / update->energyUnits < __ZeroEnergy__) {
            rVal = 1;
            break;
        }
        
        if (currStep >= 1E7) {
            rVal = -1;
            break;
        }
    }
    
    fire->isInit = false;
    return rVal;
}
int delAthermRelaxMin(Update *update) {
    AthermRelaxMin *fire = getAthermRelaxMin(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__AthermRelaxMin__");
    return 0;
}

#if (DIM == 3 || DIM == 2)
//==============
void mAVSR_Integrate(Box *box, Particle *particle, Update *update, AthermVolStressMin *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    double fact = fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0);
    double Qmass = DIM * particle->nAtom * (exp(-fact / fire->stressTol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    double deltaStrain = fire->dt * fire->gVeloc;
    double srate = fabs(deltaStrain) / (fire->dt / update->timeUnits);
    if (srate > maxStrainRate) {
        fire->gVeloc = deltaStrain / fabs(deltaStrain) * maxStrainRate / update->timeUnits;
        deltaStrain = fire->dt * fire->gVeloc;
    }
    fire->deltaStrain += deltaStrain;
    
    box->boxH[spaceIdx2voigt(sdim, gdim)] += deltaStrain * box->boxH[spaceIdx2voigt(gdim, gdim)];
#if (DIM == 3)
    if (sdim == 0 && gdim == 1) {  // xy
        box->boxH[spaceIdx2voigt(0, 2)] += deltaStrain * box->boxH[spaceIdx2voigt(1, 2)];
    }
#endif
    setBoxPara(box);
    
    //=============update particle=========
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        pos[sdim] += deltaStrain * pos[gdim];
        
        vScale(veloc, fire->scaleVeloc, veloc);
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        vScaleAdd(pos, pos, fire->dt, veloc);
    }
    
    update->Edone = update->Pdone = update->Tdone = false;
    particle->isForceValid = false;
}
void mAVSR_DotProduct(Box *box, Particle *particle, Update *update, AthermVolStressMin *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    fire->gForce = (update->pVirTens[spaceIdx2voigt(sdim, gdim)] + fire->stressSet) * box->volume;
    double fact = fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0);
    double Qmass = DIM * particle->nAtom * (exp(-fact / fire->stressTol) + 1.0);
    fire->gVeloc += 0.5 * fire->dt * fire->gForce / Qmass;
    
    double vdotf = 0.0, vdotv = 0.0, fdotf = 0.0;
    double sumForce = 0.0;
    double totKin = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr veloc = particle->veloc[iatom];
        doubleVecPtr force = particle->force[iatom];
        
        vScaleAdd(veloc, veloc, fire->dt * 0.5, force);
        totKin += sNormP2(veloc);
        
        vdotf += sDot(veloc, force);
        vdotv += sDot(veloc, veloc);
        fdotf += sDot(force, force);
        sumForce += sNorm(force);
    }
    fire->vdotf = vdotf;
    fire->vdotv = vdotv;
    fire->fdotf = fdotf;
    
    fire->vdotf += fire->gVeloc * fire->gForce;
    fire->vdotv += fire->gVeloc * fire->gVeloc;
    fire->fdotf += fire->gForce * fire->gForce;
    fire->aveForce = sumForce / update->forceUnits / particle->nAtom;
    
    double Tcur = totKin / (DIM * particle->nAtom);
    double tfact = pow((5E-3 / update->timeUnits) / fire->dt, 2);
    tfact = (tfact > 1.0 ? 1.0 : tfact);
    fire->scaleVeloc = (Tcur < fire->Tlimit * tfact ? 1.0 : sqrt(fire->Tlimit * tfact / Tcur));
}
bool mAVSR_CheckConverge(Box *box, Particle *particle, Update *update, AthermVolStressMin *fire) {
    int sdim = fire->shearDim, gdim = fire->gradDim;
    if (fabs(fire->deltaStrain) > fire->maxDeltaStrain) {
        fire->rtype = -1;
        return true;
    }
    
    if (fire->aveForce > __ZeroForce__)
        return false;
    if (fabs(update->pVirTens[spaceIdx2voigt(sdim, gdim)] / fire->stressSet + 1.0) > fire->stressTol) {
        return false;
    }
    
    fire->rtype = 0;
    return true;
}

AthermVolStressMin *addAthermVolStressMin(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getAthermVolStressMin(update))
        Abort("repetitive addAthermVolStressMin(...)");
    
    AthermVolStressMin *fire = (AthermVolStressMin *)calloc(sizeof(AthermVolStressMin), 1);
    addToolkit(&update->toolkit, (void *)fire, NULL, "__AthermVolStressMin__");
    
    fire->stressTol = __Stol__;
    fire->dtSet = 5E-3;
    
    box->isShapeFixed = false;
    particle->isSizeFixed = true;
    
    cmdArg *cmd = findVariable(var, "csmin");
    if (!cmd)
        Abort("--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit or --csmin atherm");
    if (cmd->cmdArgc == 1) {
        if (strcmp(cmd->cmdArgv[0], "atherm") != 0)
            Abort("--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit or --csmin atherm");
        
        fire->isInit = false;
        return fire;
    }
    
    if (cmd->cmdArgc != 5) {
        Abort("--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit");
    }
    if (strcmp(cmd->cmdArgv[1], "xy") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 1;
    }
#if (DIM == 3)
    else if (strcmp(cmd->cmdArgv[1], "xz") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 2;
    } else if (strcmp(cmd->cmdArgv[1], "yz") == 0) {
        fire->shearDim = 1;
        fire->gradDim = 2;
    }
#endif
    else {
        Abort("--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit");
    }
    fire->stressSet = atof(cmd->cmdArgv[0]);
    fire->maxDeltaStrain = fabs(atof(cmd->cmdArgv[2]));
    fire->deltaStrain = 0;
    
    if (strcmp(cmd->cmdArgv[3], "atherm") != 0)
        Abort("--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit");
    fire->Tlimit = atof(cmd->cmdArgv[4]);
    if (fire->Tlimit < 1E-20)
        Abort("upper limits of temperature is too small!");
    
    fire->isInit = true;
    return fire;
}
int reInitAthermVolStressMin(AthermVolStressMin *fire, double tStress, char *sType, double maxDeltaStrain, double T_limit) {
    memset(fire, '\0', sizeof(AthermVolStressMin));
    
    fire->stressTol = __Stol__;
    fire->dtSet = 5E-3;
    
    if (strcmp(sType, "xy") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 1;
    }
#if (DIM == 3)
    else if (strcmp(sType, "xz") == 0) {
        fire->shearDim = 0;
        fire->gradDim = 2;
    } else if (strcmp(sType, "yz") == 0) {
        fire->shearDim = 1;
        fire->gradDim = 2;
    }
#endif
    else {
        safeFprintf(stderr, "--csmin tStress xy|xz|yz maxDeltaStrain atherm T_limit");
        return -1;
    }
    fire->stressSet = tStress;
    fire->maxDeltaStrain = maxDeltaStrain;
    fire->deltaStrain = 0;
    
    fire->Tlimit = T_limit;
    if (fire->Tlimit < 1E-20) {
        safeFprintf(stderr, "upper limits of temperature is too small!");
        return -1;
    }
    
    fire->isInit = true;
    return 0;
}
AthermVolStressMin *getAthermVolStressMin(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__AthermVolStressMin__");
    if (whichTool < 0)
        return NULL;
    return (AthermVolStressMin *)update->toolkit.toolkit[whichTool];
}
int minAthermVolStressRelax(Box *box, Particle *particle, Update *update) {
    AthermVolStressMin *fire = getAthermVolStressMin(update);
    if (!fire || !fire->isInit) Abort("call addAthermVolStressMin(...) or reInitAthermVolStressMin(...)");
    
    box->isShapeFixed = false;
    particle->isSizeFixed = true;
    update->nebrList.nDelay = 0;
    update->nebrList.skinSet = 0.1;
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    
    double dtmax = 0.5 * update->timeUnits;
    int currStep = 0, last_negative = 0;
    double alpha = ALPHA0;
    fire->dt = fire->dtSet * update->timeUnits;
    double dt0 = fire->dt;
    double stressSet = fire->stressSet;
    fire->stressSet *= update->pressureUnits;  // incorporating units;
    
    calcForce(box, particle, update);
    mAVSR_DotProduct(box, particle, update, fire);
    // remove the side effects of _dotProduct()
    memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
    fire->gVeloc = 0;
    fire->scaleVeloc = 1.0;
    double ePairLastNeg = update->ePair;
    int rVal = 0;
    
    bool isConverge = mAVSR_CheckConverge(box, particle, update, fire);
    while (!isConverge) {
//        if (fire->scaleVeloc < 1.0)
//            printf("%d %g\n", currStep, fire->scaleVeloc);
        
        // update x
        mAVSR_Integrate(box, particle, update, fire);
        
        // update v
        calcForce(box, particle, update);
        mAVSR_DotProduct(box, particle, update, fire);
        isConverge = mAVSR_CheckConverge(box, particle, update, fire);
        currStep++;
        
        // mixing v and update timestep
        if (fire->vdotf > 0) {
            if (currStep - last_negative > DELAYSTEP) {
                fire->dt *= DT_GROW;
                fire->dt = cpuMin(fire->dt, dtmax);
                alpha *= ALPHA_SHRINK;
            }
            
            double scale1 = 1.0 - alpha;
            double scale2 = alpha * sqrt(fire->vdotv / fire->fdotf);
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                doubleVecPtr veloc = particle->veloc[iatom];
                doubleVecPtr force = particle->force[iatom];
                vScale(veloc, scale1, veloc);
                vScaleAdd(veloc, veloc, scale2, force);
            }
            fire->gVeloc = fire->gVeloc * scale1 + fire->gForce * scale2;
        } else {
            int deltaStep = currStep - last_negative;
            if (deltaStep < log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 0.9;
                dtmax = cpuMax(dtmax, 5.0 * fire->dtSet * update->timeUnits);
            } else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + DELAYSTEP) {
                dtmax = dtmax * 1.1;
                dtmax = cpuMin(dtmax, 0.5 * update->timeUnits);
            }
            
            fire->dt *= DT_SHRINK;
            fire->dt = cpuMax(fire->dt, 0.1 * fire->dtSet * update->timeUnits);

            dt0 = fire->dt;
            last_negative = currStep;
            ePairLastNeg = update->ePair;
            alpha = ALPHA0;
            memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
            fire->gVeloc = 0;
        }
        
//        if (currStep % 1000 == 0) {
//            safeFprintf(stdout, "%d %g %g (%g) %.15g\n", currStep, fire->deltaStrain,
//                        update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / update->pressureUnits,
//                        fabs(update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / fire->stressSet + 1.0),
//                        fire->aveForce / update->forceUnits);
//        }
    }
    
//    safeFprintf(stdout, "%d %g %g (%g) %.15g\n", currStep, fire->deltaStrain,
//                update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / update->pressureUnits,
//                fabs(update->pVirTens[spaceIdx2voigt(fire->shearDim, fire->gradDim)] / fire->stressSet + 1.0),
//                fire->aveForce / update->forceUnits);
    
    fire->stressSet = stressSet;
    fire->isInit = false;
    return rVal;
}
int delAthermVolStressMin(Update *update) {
    AthermVolStressMin *fire = getAthermVolStressMin(update);
    if (!fire)
        return -1;
    delToolkit(&update->toolkit, "__AthermVolStressMin__");
    return 0;
}
#endif


#define __Tlimit__  1E-10
void mBBGD_DotProduct(Box *box, Particle *particle, Update *update, minBBGD *gd) {
    //__float128 sumForce = 0;
    double sumForce = 0;  //?
    double dXdotF = 0.0, dFdotF = 0.0, dXdotX = 0.0, normF = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector dF;
        vSub(dF, particle->force[iatom], gd->prevForce[iatom]);
        doubleVector dX;
        vScale(dX, gd->beta, gd->prevForce[iatom]);
        
        dFdotF += sDot(dF, dF);
        dXdotF += sDot(dX, dF);
        dXdotX += sDot(dX, dX);
        
        normF += sNormP2(particle->force[iatom]);
        sumForce += sNorm(particle->force[iatom]);
    }
    gd->dXdotF = dXdotF;
    gd->dFdotF = dFdotF;
    gd->dXdotX = dXdotX;
    gd->normForce = sqrt(normF);
    gd->aveForce = sumForce / particle->nAtom / update->forceUnits;
}
minBBGD *getMinBBGD(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__MinBBGD__");
    if (whichTool < 0)
        return NULL;
    return (minBBGD *)update->toolkit.toolkit[whichTool];
}
minBBGD *addMinBBGD(Box *box, Particle *particle, Update *update, Variable *var) {
    if (getMinBBGD(update) != NULL) {
        return getMinBBGD(update);
    }
    if (!findVariable(var, "BBGD")) {
        Abort("--BBGD: BB iterations of GD, Tmax = 1E-10.")
    }
    minBBGD *gd = (minBBGD *)calloc(sizeof(minBBGD), 1);
    addToolkit(&update->toolkit, (void *)gd, NULL, "__MinBBGD__");
    
    gd->dtSet = 5E-3;
    return gd;
}
int delMinBBGD(Update *update) {
    minBBGD *gd = getMinBBGD(update);
    if (!gd)
        return -1;
    delToolkit(&update->toolkit, "__MinBBGD__");
    return 0;
}
int minBBGDrelax(Box *box, Particle *particle, Update *update) {
    minBBGD *gd = getMinBBGD(update);
    if (!gd)
        Abort("Call addMinBBGD(...)!");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    update->nebrList.nDelay = 0;
    update->nebrList.skinSet = 0.1;
    
    update->Edone = update->Pdone = update->Wdone = update->Tdone = false;
    update->isThermalRun = false;
    update->nebrList.nDelay = 0;
    update->nebrList.nRebuildMax = 0;
    double beta0 = 0.5 * pow(gd->dtSet * update->distanceUnits, 2.0);
    double betaMin = 0.5 * pow(0.1 * gd->dtSet * update->distanceUnits, 2.0);
    double maxBetaNormF = sqrt(__Tlimit__ * DIM * particle->nAtom) * gd->dtSet * update->distanceUnits;
    int currStep = 0;
    
    gd->prevForce = particle->veloc;
    memset(gd->prevForce, '\0', particle->nAtom * sizeof(doubleVector));
    calcForceAllPair(box, particle, update);
    mBBGD_DotProduct(box, particle, update, gd);
    if (gd->aveForce <= __ZeroForce__) {
        return 0;
    }

    // first iteration
    gd->beta = beta0;
    gd->beta = fmax(gd->beta, betaMin);
    if (gd->beta * gd->normForce > maxBetaNormF) {
        gd->beta = maxBetaNormF / gd->normForce;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScaleAdd(particle->pos[iatom], particle->pos[iatom], gd->beta, particle->force[iatom]);
    }
    particle->isForceValid = false;
    // swap force and previous force;
    gd->prevForce = particle->veloc;
    particle->veloc = particle->force;
    particle->force = gd->prevForce;
    gd->prevForce = particle->veloc;
    // calculate force
    calcForceAllPair(box, particle, update);
    mBBGD_DotProduct(box, particle, update, gd);
    currStep++;

    // BB iterations of GD
    while (!(gd->aveForce <= __ZeroForce__)) {
        gd->beta = fmax(fabs(gd->dXdotX/gd->dXdotF), fabs(gd->dXdotX/gd->dXdotF));
        gd->beta = fmax(gd->beta, betaMin);
        if (gd->beta > maxBetaNormF/gd->normForce) {
            gd->beta = maxBetaNormF / gd->normForce;
        }
        // update x
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            vScaleAdd(particle->pos[iatom], particle->pos[iatom], gd->beta, particle->force[iatom]);
        }
        particle->isForceValid = false;
        // swap force and previous force;
        gd->prevForce = particle->veloc;
        particle->veloc = particle->force;
        particle->force = gd->prevForce;
        gd->prevForce = particle->veloc;
        // calculate force
        calcForceAllPair(box, particle, update);
        mBBGD_DotProduct(box, particle, update, gd);
        currStep++;
        
        if (update->ePair / update->energyUnits < __ZeroEnergy__) {
            break;
        }
        
        if (currStep > 2E7) {
            return -1;
        }
    }

    return 0;
}

