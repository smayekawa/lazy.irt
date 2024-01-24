Cards="
name type ncat   p1   p2   p3   p4
  Q1    B    2  0.9  0.0,,
  Q2    B3   2  1.0  0.0  0.2,
  Q3    G    4  1.0 -2.0  0.0  2.0
  Q4    P    4  0.8 -2.0  0.0  2.0
 ";
paramS1 <- cards( Cards, header=1 )
use_data(paramS1, overwrite=1)


Cards="
name type ncat   p1   p2   p3   p4
  Q11    B    2  0.9  0.0,,
  Q12    B    2  0.9  1.0,,
  Q21    B3   2  1.0  0.0  0.2,
  Q22    B3   2  1.0  1.0  0.2,
  Q31    G    3  1.0 -1.0  1.0,
  Q32    G    4  1.0 -2.0  0.0  2.0
  Q41    P    3  0.8 -1.0  1.0,
  Q42    P    4  0.8 -2.0  0.0  2.0
 ";
paramS2 <- cards( Cards, header=1 )
use_data(paramS2, overwrite=1)

