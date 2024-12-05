#ifndef EVTTIMEFNALBEAM_H
#define EVTTIMEFNALBEAM_H

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

class EvtTimeFNALBeam {
public:
    EvtTimeFNALBeam(); // Constructor
    double TimeOffset();
    int nbatch;

private:
    double fTimeBetweenBuckets;
    double fBucketTimeSigma;
    int fNBucketsPerBatch;
    int fNFilledBucketsPerBatch;
    double fGlobalOffset;
};

#endif // EVTTIMEFNALBEAM_H