//
//  AthermalRelaxMin.h
//  TimeDrivenSimLib
//
//  Created by Deng Pan on 2023/9/7.
//

#ifndef AthermalRelaxMin_h
#define AthermalRelaxMin_h

#include "SimSubFunc_ndim.h"
#include "StructSim_ndim.h"
#include "VectorMath_ndim.h"
#include "FireMin_ndim.h"
#include "GeoFamily_ndim.h"

typedef struct AthermRelaxMin {
    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet, dt;
    double vdotf, vdotv, fdotf;
    double aveForce;

    double Tlimit, scaleVeloc;
    
    bool isInit;
} AthermRelaxMin;

//=================const box shape========
AthermRelaxMin *addAthermRelaxMin(Box *box, Particle *particle, Update *update, Variable *var);
int reInitAthermRelaxMin(AthermRelaxMin *fire, double T_limit);
AthermRelaxMin *getAthermRelaxMin(Update *update);
int minAthermRelax(Box *box, Particle *particle, Update *update);
int delAthermRelaxMin(Update *update);
int minBBGDrelax(Box *box, Particle *particle, Update *update);

typedef struct minBBGD {
    double dtSet;
    
    // criteria: average force amplitude per particle is less than 1E-14
    double dXdotF, dFdotF, dXdotX;
    double aveForce, normForce;
    
    double beta;              // dX = beta*prevForce;
    doubleVector *prevForce;  // pointer of particle->veloc
} minBBGD;

minBBGD *getMinBBGD(Update *update);
minBBGD *addMinBBGD(Box *box, Particle *particle, Update *update, Variable *var);
int delMinBBGD(Update *update);

#if (DIM == 3 || DIM == 2)

typedef struct AthermVolStressMin {
    int shearDim, gradDim;
    double stressSet, stressTol;
    
    double deltaStrain, maxDeltaStrain;
    double gVeloc, gForce;
    
    // criteria: average force amplitude per particle is less than 1E-14
    double dtSet, dt;
    double vdotf, vdotv, fdotf;
    double aveForce;
    
    int rtype;//0: converge, -1: reaching maximum delta strain.
    double Tlimit, scaleVeloc;
    
    bool isInit;
} AthermVolStressMin;

//================const vol and stress=======
AthermVolStressMin *addAthermVolStressMin(Box *box, Particle *particle, Update *update, Variable *var);
int reInitAthermVolStressMin(AthermVolStressMin *fire, double tStress, char *sType, double maxDeltaStrain, double T_limit);
AthermVolStressMin *getAthermVolStressMin(Update *update);
int minAthermVolStressRelax(Box *box, Particle *particle, Update *update);
int delAthermVolStressMin(Update *update);

#endif

#endif /* AthermalRelaxMin_h */
