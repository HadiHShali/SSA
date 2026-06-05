DataFile            ARML_1.mom
DataDirectory       ./pre_files
OutputFile          ./fin_files\ARML_1.mom
interpolate         no
TimeUnit              days
estimatepostseismic yes
PhysicalUnit        mm
JSON                yes
ScaleFactor         1.0
RandomiseFirstGuess  yes
periodicsignals     365.25 182.625
estimateoffsets     yes
NoiseModels          FlickerGGM RandomWalkGGM White
GGM_1mphi           6.9e-06
useRMLE             no
Verbose               no
OffsetEpochsFile      OffsetEpochs.txt
