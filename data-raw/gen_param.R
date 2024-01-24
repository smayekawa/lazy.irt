library(devtools)


Cards="
name type ncat   p1      p2     p3    p4
  Q1    B    2  1.0     0.0    0.0,
  Q2    B3   2  1.0     0.0    0.2,
  Q3    Bn   2  1.0     0.0    0.0,
  Q4    Bn3  2  1.0     0.0    0.2,
  Q5    P    3  1.0    -1.0    1.0,
  Q6    G    3  1.0664 -1.0288 1.0288,
  Q7    Gn   3  1.0819 -1.0221 1.0221,
  Q8    P    4  0.8    -2.0    0.0    2.0
  Q9    G    4  0.9179 -2.0495 0.000  2.0495
  Q10   Gn   4  0.9396 -2.0245 0.000  2.0245
 ";
paramA1 <- cards( Cards, header=1 )
use_data(paramA1, overwrite=1)


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


Cards="
name type ncat   p1   p2   p3
  Q1    B    2  1.0 -1.0, 0.0
  Q2    B    2  1.0  0.0, 0.0
  Q3    B    2  1.0  1.0  0.0
  Q4    B    2  1.0 -1.0  0.0
  Q5    B    2  1.0  0.0  0.0
  Q6    B    2  1.0  1.0  0.0
  Q7    B    2  1.3 -1.0, 0.0
  Q8    B    2  1.3  0.0, 0.0
  Q9    B    2  1.3  1.0  0.0
  Q10   B    2  1.3 -1.0  0.0
  Q11   B    2  1.3  0.0  0.0
  Q12   B    2  1.3  1.0  0.0
  Q13   B    2  0.7 -1.0, 0.0
  Q14   B    2  0.7  0.0, 0.0
  Q15   B    2  0.7  1.0  0.0
  Q16   B    2  0.7 -1.0  0.0
  Q17   B    2  0.7  0.0  0.0
  Q18   B    2  0.7  1.0  0.0
 ";
paramB1 <- cards( Cards, header=1 )
use_data(paramB1, overwrite=1)


Cards="
name type ncat   p1   p2   p3
  Q1    B3    2  1.0 -1.0  0.2
  Q2    B3    2  1.0  0.0  0.2
  Q3    B3    2  1.0  1.0  0.2
  Q4    B3    2  1.0 -1.0  0.2
  Q5    B3    2  1.0  0.0  0.2
  Q6    B3    2  1.0  1.0  0.2
  Q7    B3    2  1.3 -1.0  0.2
  Q8    B3    2  1.3  0.0  0.2
  Q9    B3    2  1.3  1.0  0.2
  Q10   B3    2  1.3 -1.0  0.2
  Q11   B3    2  1.3  0.0  0.2
  Q12   B3    2  1.3  1.0  0.2
  Q13   B3    2  0.7 -1.0  0.2
  Q14   B3    2  0.7  0.0  0.2
  Q15   B3    2  0.7  1.0  0.2
  Q16   B3    2  0.7 -1.0  0.2
  Q17   B3    2  0.7  0.0  0.2
  Q18   B3    2  0.7  1.0  0.2
 ";
paramB2 <- cards( Cards, header=1 )
use_data(paramB2, overwrite=1)


