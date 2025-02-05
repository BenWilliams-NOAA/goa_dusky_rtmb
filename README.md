GOA dusky rockfish briding ADMB to RTMB

 - 1. make sure all of the packages in the `utils.r` file are available,
 - 2. the 2024 ADMB data are ported over and used as inputs/comparisons with the RTMB model,
 - 3. first create the data list, then the parameter list, using `map` fixes a parameter
 - 4. build an object using `RTMB::MakeADFun()`, for comparison all parameters can be fixed, then the model doesn't need to be optimized.
