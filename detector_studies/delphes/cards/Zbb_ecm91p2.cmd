! 1) Settings used in the main program.
Random:setSeed = on
Random:seed = 1                    ! seed
Main:numberOfEvents = 10000       ! number of events
Main:timesAllowErrors = 5          ! how many aborts before run stops

! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showChangedParticleData = off ! list changed particle data
Next:numberCount = 10000           ! print message every n events
Next:numberShowInfo = 1            ! print event information n times
Next:numberShowProcess = 1         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

! 3) Beam parameter settings. Values below agree with default ones.
Beams:idA = 11                     ! first beam, e = 2212, pbar = -2212
Beams:idB = -11                    ! second beam, e = 2212, pbar = -2212

Beams:allowMomentumSpread  = off


! 4) Hard process : Z->all
Beams:eCM = 91.2                  ! CM energy of collision


WeakSingleBoson:ffbar2gmZ = on
23:onMode    = off                 ! switch off Z boson decays
23:onIfAny   = 5                   ! u-quarks 1 2 3 4 5 6

PartonLevel:ISR = on               ! initial-state radiation
PartonLevel:FSR = on               ! final-state radiation
