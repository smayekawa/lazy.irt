library(devtools)


Cards="
name type ncat w     v0     v1     v2   v3
  Q1    B    2 1    0       1,,
  Q2    B3   2 1    0       1,,
  Q3    Bn   2 1    0       1,,
  Q4    Bn3  2 1    0       1,,
  Q5    P    3 1    0       1,     2,
  Q6    G    3 1    0       1,     2,
  Q7    Gn   3 1    0       1,     2,
  Q8    P    4 1    0       1,     2,     3
  Q9    G    4 1    0       1,     2,     3
  Q10   Gn   4 1    0       1,     2,     3
 ";
weightA1 <- cards( Cards, header=1 )
use_data(weightA1, overwrite=1)



Cards="
name type ncat  w    v0   v1   v2   v3
Q1    B    2    1    0    1 ,,
Q2    B3   2    1    0    1 ,,
Q3    G    4    1    0    1    2    3
Q4    P    4    1    0    1    2    3
";
weightS1 <- cards( Cards, header=1 )
use_data(weightS1, overwrite=1)

Cards="
name type ncat  w    v0   v1   v2   v3
Q1    B    2    1    0    1 ,,
Q2    B3   2    2    0    1 ,,
Q3    G    4    1    0    1    2    3
Q4    P    4    2    0    1    2    3
";
weightS11 <- cards( Cards, header=1 )
use_data(weightS11, overwrite=1)

Cards="
name type ncat  w    v0   v1   v2   v3
Q1    B    2    1    0    1 ,,
Q2    B3   2    1    0    2 ,,
Q3    G    4    1    0    1    2    5
Q4    P    4    1    0    0    0    1
";
weightS12 <- cards( Cards, header=1 )
use_data(weightS12, overwrite=1)



Cards="
name type ncat  w    v0   v1   v2   v3
Q11    B    2   1    0    1 ,,
Q12    B    2   1    0    1 ,,
Q21    B3   2   1    0    1 ,,
Q22    B3   2   1    0    1 ,,
Q31    G    3   1    0    1    2,
Q32    G    4   1    0    1    2    3
Q41    P    3   1    0    1    2,
Q42    P    4   1    0    1    2    3
";
weightS2 <- cards( Cards, header=1 )
use_data(weightS2, overwrite=1)

Cards="
name type ncat  w    v0   v1   v2   v3
Q11    B    2   1    0    1 ,,
Q12    B    2   2    0    1 ,,
Q21    B3   2   1    0    1 ,,
Q22    B3   2   2    0    1 ,,
Q31    G    3   1    0    1    2,
Q32    G    4   2    0    1    2    3
Q41    P    3   1    0    1    2,
Q42    P    4   2    0    1    2    3";
weightS21 <- cards( Cards, header=1 )
use_data(weightS21, overwrite=1)

Cards="
name type ncat  w    v0   v1   v2   v3
Q11    B    2   1    0    1 ,,
Q12    B    2   1    0    2 ,,
Q21    B3   2   1    0    1 ,,
Q22    B3   2   1    0    3 ,,
Q31    G    3   1    0    1    2,
Q32    G    4   1    0    1    3    5
Q41    P    3   1    0    2    3,
Q42    P    4   1    0    0    0    3
";
weightS22 <- cards( Cards, header=1 )
use_data(weightS22, overwrite=1)
