## Test of functional diversity (beta) on permuted data

## Theory:
##    Randomize data, per usual
##      If individual, not a problem; pool and calculate for each ith level
##      If sample based, constrained randomization; randomize i-1 within i+1
##    Find Beta
##      Calculate beta per usual for i-rand at ith level (beta func mult)
##      Calculate beta for s-rand
##        Pool i-1 samples within the ith level
##        Calculate beta funct mult for the ith level
##    Calculate p-value
##      Per usual? - MBM is leaning this way
##        Count of values greater than observed / # randomizations
##      SES seems to be common in functional diversity studies, but this is for
##      ACTUAL FUNCTIONAL DIVERITY VALUES, NOT FUNCTIONAL BETA DIVERSITY
##        Obs - mean(Rand) / SD(Rand)

## Practice:
##    Randomize data; run beta.





