In order to get the ZFromData code to compile, one must first use the latest and greatest version of the Tag&Probe tools.

To Work with the CMSSW_2_1_X
cvs co -r V01-05-00 PhysicsTools/TagAndProbe


I dunno how important, but these are needed to have a "good" version of isodeposits.
cvs co -r V00-16-07  DataFormats/RecoCandidate
cvs co -r V00-02-05  EgammaAnalysis/EgammaIsolationProducers
cvs co -r V01-02-00  PhysicsTools/IsolationAlgos
cvs co -r V00-02-03  RecoEgamma/EgammaIsolationAlgos

Also, I have an updated IsolationSelector that goes into the PhysicsTools/TagAndProbe. So you need to unpack this tar-ball. (There is 3X development going on in the package and I do not want to mess that up.)

wget http://webusers.physics.umn.edu/~haupt/JasonPhysicsToolsTagAndProve.tar
cd CMSSW_2_1_X/src
tar -xvf JasonPhysicsToolsTagAndProve.tar
----------------------------
Jason