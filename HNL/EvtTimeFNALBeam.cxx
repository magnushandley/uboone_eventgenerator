#include "EvtTimeFNALBeam.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"


EvtTimeFNALBeam::EvtTimeFNALBeam() 
    : nbatch(6), 
      fTimeBetweenBuckets(1e9 / 53.103e6),
      fBucketTimeSigma(0.750),
      fNBucketsPerBatch(84),     // NOvA-era 81+3, MINOS-era 81+5
      fNFilledBucketsPerBatch(81), // 81 for both eras 
      fGlobalOffset(0.0){
}

double EvtTimeFNALBeam::TimeOffset() {
  double offset = CLHEP::RandGauss::shoot(0.0,fBucketTimeSigma);

    // pick a bucket within a batch
    // assume all ~ buckets constant in batch until we have another model
  offset +=  fTimeBetweenBuckets * 
              (double)CLHEP::RandFlat::shootInt(fNFilledBucketsPerBatch);

  int ibatch = CLHEP::RandFlat::shootInt(nbatch);
  offset += fTimeBetweenBuckets*(double)fNBucketsPerBatch*(double)ibatch;

  // finally the global offset
  return offset + fGlobalOffset;
}