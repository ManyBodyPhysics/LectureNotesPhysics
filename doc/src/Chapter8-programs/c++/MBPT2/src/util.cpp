int isASumOfThreeSquares (unsigned n) {
  if (n == 0) {
    return 1;
  }
  while (n % 4 == 0) {
    n /= 4;
  }
  return (n % 8) != 7;
}

unsigned smallestSquareRootAtLeast (unsigned x) {
  unsigned square = 1;
  unsigned delta = 3;
  while(square < x){
    square += delta;
    delta += 2;
  }
  return (delta/2 - 1);
}
